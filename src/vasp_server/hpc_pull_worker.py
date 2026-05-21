import asyncio
import json
import logging
import mimetypes
import os
import subprocess
import tempfile
import threading
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlsplit

import requests

from .failure_guidance import FailureGuidance, derive_failure_guidance, derive_failure_guidance_from_workdir
from .settings import settings
from .vasp_worker import VaspWorker
from .worker_client import ControlPlaneClient


TASK_TYPE_TO_METHOD = {
    "structure_optimization": "run_structure_optimization",
    "scf_calculation": "run_scf_calculation",
    "dos_calculation": "run_dos_calculation",
    "dos_analysis": "run_dos_analysis",
    "md_calculation": "run_md_calculation",
    "md_analysis": "run_md_analysis",
    "band_structure": "run_band_structure_calculation",
    "band_structure_analysis": "run_band_structure_analysis",
    "neb_calculation": "run_neb_calculation",
    "phonon_calculation": "run_phonon_calculation",
    "custom_calculation": "run_custom_calculation",
    "agent_analysis": "run_agent_analysis",
}

TASK_TYPE_TO_ANALYZER = {
    "structure_optimization": "_analyze_results",
    "scf_calculation": "_analyze_scf_results",
    "dos_calculation": "_analyze_dos_results",
    "dos_analysis": "_analyze_dos_results",
    "md_calculation": "_analyze_md_results",
    "md_analysis": "_analyze_md_results",
    "band_structure": "_analyze_band_structure_results",
    "band_structure_analysis": "_analyze_band_structure_results",
    "neb_calculation": "_analyze_neb_results",
    "phonon_calculation": "_analyze_phonon_results",
    "custom_calculation": "_analyze_custom_results",
}

ACTIVE_SCHEDULER_STATES = {"PENDING", "RUNNING", "CONFIGURING", "COMPLETING"}
FAILED_SCHEDULER_STATES = {
    "FAILED",
    "CANCELLED",
    "TIMEOUT",
    "NODE_FAIL",
    "OUT_OF_MEMORY",
    "BOOT_FAIL",
    "PREEMPTED",
}
PUBLIC_ARTIFACT_EXTENSIONS = {".html", ".png", ".jpg", ".jpeg", ".svg", ".webp", ".csv", ".json", ".pdf"}
PUBLIC_ARTIFACT_FILENAMES = {
    "outcar",
    "contcar",
    "poscar",
    "incar",
    "kpoints",
    "oszicar",
    "ibzkpt",
    "result.log",
    "vasp_job.sh",
}
PUBLIC_JSON_PRIVATE_KEYS = {
    "work_directory",
    "source_work_directory",
    "local_path",
}
PUBLIC_JSON_REPORT_KEYS = {
    "analysis_report_html_path",
    "html_analysis_report",
    "md_analysis_report_html_path",
    "dos_analysis_report_html_path",
    "scf_analysis_report_html_path",
    "band_structure_report_html_path",
    "neb_report_html_path",
    "phonon_report_html_path",
    "agent_analysis_report_html_path",
}


class TaskCanceled(Exception):
    pass


logger = logging.getLogger(__name__)


class PullWorker:
    def __init__(
        self,
        control_plane_client: Optional[ControlPlaneClient] = None,
        execution_worker: Optional[VaspWorker] = None,
        worker_id: str = "hpc-worker",
        queue_name: str = "default",
        poll_interval_seconds: Optional[int] = None,
    ) -> None:
        self.worker_id = worker_id
        self.queue_name = queue_name
        self.poll_interval_seconds = poll_interval_seconds or settings.worker_poll_interval_seconds
        self.control_plane_client = control_plane_client or ControlPlaneClient(
            base_url=settings.vasp_server_base_url,
            worker_token=settings.internal_worker_token,
            worker_id=worker_id,
            queue_name=queue_name,
        )
        self.execution_worker = execution_worker or VaspWorker(user_id=worker_id, base_work_dir=settings.vasp_work_root)
        self._control_plane_lock = threading.Lock()

    @staticmethod
    def _build_failure_context(
        *,
        scheduler_job_id: Optional[str],
        work_directory: Optional[str],
        extras: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        context: dict[str, Any] = {}
        if scheduler_job_id:
            context["scheduler_job_id"] = str(scheduler_job_id)
            if str(scheduler_job_id).isdigit():
                context["process_id"] = int(str(scheduler_job_id))
        if work_directory:
            context["work_directory"] = str(work_directory)
        if extras:
            context.update(extras)
        return context

    @staticmethod
    def _append_failure_guidance(error_message: str, guidance: FailureGuidance) -> str:
        parts = [error_message or "任务执行失败"]
        if guidance.reason and guidance.reason not in parts[0]:
            parts.append(f"诊断: {guidance.reason}")
        if guidance.suggested_action and guidance.suggested_action not in parts[0]:
            parts.append(f"建议给智能体的下一步: {guidance.suggested_action}")
        return "\n".join(parts)

    def _derive_failure_guidance_for_report(
        self,
        *,
        error_message: str,
        task_type: str,
        work_directory: Optional[str],
    ) -> FailureGuidance:
        guidance = derive_failure_guidance(error_message, task_type=task_type)
        if guidance.failure_type:
            return guidance

        if work_directory:
            try:
                return derive_failure_guidance_from_workdir(Path(str(work_directory)), task_type=task_type)
            except OSError:
                return FailureGuidance()
        return FailureGuidance()

    def _build_guided_failure_payload(
        self,
        *,
        task_id: str,
        task_type: str,
        error_message: str,
        scheduler_job_id: Optional[str],
        work_directory: Optional[str],
        extras: Optional[dict[str, Any]] = None,
    ) -> dict[str, Any]:
        guidance = self._derive_failure_guidance_for_report(
            error_message=error_message,
            task_type=task_type,
            work_directory=work_directory,
        )
        failure_extras = dict(extras or {})
        if guidance.failure_type:
            failure_extras["failure_type"] = guidance.failure_type
        if guidance.reason:
            failure_extras["failure_reason"] = guidance.reason
        if guidance.suggested_action:
            failure_extras["suggested_action"] = guidance.suggested_action

        return {
            "error_message": self._append_failure_guidance(error_message, guidance),
            "failure_code": guidance.failure_type,
            "failure_context": self._build_failure_context(
                scheduler_job_id=scheduler_job_id,
                work_directory=work_directory,
                extras=failure_extras,
            ),
            "artifact_manifest": self._collect_and_stage_public_artifacts(task_id, {"work_directory": work_directory}),
        }

    def run_loop(self, iterations: Optional[int] = None) -> None:
        completed = 0
        while iterations is None or completed < iterations:
            self.run_once()
            completed += 1
            time.sleep(self.poll_interval_seconds)

    def run_once(self) -> bool:
        if not self._drain_pending_completions():
            return False

        reconciled = self._reconcile_runtime_tasks()
        if reconciled is not None:
            return reconciled

        try:
            task = self._call_control_plane_with_retries(
                "claim",
                self.control_plane_client.claim,
            )
        except Exception as exc:
            logger.warning("worker %s claim failed after retries: %s", self.worker_id, exc)
            return False

        if not task:
            return False

        task_id = task["task_id"]
        lease_token = task["lease_token"]
        task_type = task["task_type"]
        params = dict(task.get("params") or {})
        self._upsert_runtime_state(
            task_id,
            task_type=task_type,
            params=params,
            queue_name=self.queue_name,
            lease_token=lease_token,
            work_directory=str(self._default_work_directory(task_id)),
            started_at=datetime.now(timezone.utc).isoformat(),
        )

        scheduler_job_id: Optional[str] = None
        cancel_requested = False
        heartbeat_state = {
            "scheduler_job_id": None,
            "cancel_requested": False,
        }

        async def progress_callback(progress: int, message: str, pid: Optional[int] = None):
            nonlocal scheduler_job_id, cancel_requested
            if pid is not None:
                scheduler_job_id = str(pid)
                heartbeat_state["scheduler_job_id"] = scheduler_job_id
            self._upsert_runtime_state(
                task_id,
                lease_token=lease_token,
                scheduler_job_id=scheduler_job_id,
                last_progress=progress,
                last_message=message,
            )
            heartbeat = self._call_control_plane_with_retries(
                "heartbeat",
                lambda: self.control_plane_client.heartbeat(task_id, lease_token),
            )
            if heartbeat.get("cancel_requested"):
                cancel_requested = True
                heartbeat_state["cancel_requested"] = True
            if cancel_requested and scheduler_job_id:
                subprocess.run(["scancel", scheduler_job_id], check=False)
                self._ack_cancel(task_id, lease_token)
                raise TaskCanceled(message)
            if cancel_requested and pid is None:
                return
            if cancel_requested:
                if scheduler_job_id:
                    subprocess.run(["scancel", scheduler_job_id], check=False)
                    self._ack_cancel(task_id, lease_token)
                raise TaskCanceled(message)

        try:
            self._call_control_plane_with_retries(
                "mark_running",
                lambda: self.control_plane_client.mark_running(task_id, lease_token),
            )
            heartbeat_stop_event, heartbeat_thread = self._start_background_heartbeat(
                task_id,
                lease_token,
                heartbeat_state,
            )
            result = self._run_task(task_type, task_id, params, progress_callback)
            result.setdefault("work_directory", str(self._default_work_directory(task_id)))
            self._stop_background_heartbeat(heartbeat_stop_event, heartbeat_thread)
            cancel_requested = cancel_requested or bool(heartbeat_state.get("cancel_requested"))
            if cancel_requested:
                self._ack_cancel(task_id, lease_token)
                return True
            succeeded, error_message, failure_extras = self._classify_result(task_type, result)
            if not succeeded:
                work_directory = result.get("work_directory") or str(self._default_work_directory(task_id))
                payload = self._build_guided_failure_payload(
                    task_id=task_id,
                    task_type=task_type,
                    error_message=error_message or "任务执行失败",
                    scheduler_job_id=scheduler_job_id,
                    work_directory=work_directory,
                    extras=failure_extras,
                )
                self._report_failure(task_id, lease_token, payload)
                return True

            payload = {
                "result_data": result,
                "artifact_manifest": self._collect_and_stage_public_artifacts(task_id, result),
            }
            self._report_completion(task_id, lease_token, payload)
            return True
        except TaskCanceled:
            return True
        except Exception as exc:
            work_directory = str(self._default_work_directory(task_id))
            self._report_failure(
                task_id,
                lease_token,
                self._build_guided_failure_payload(
                    task_id=task_id,
                    task_type=task_type,
                    error_message=str(exc),
                    scheduler_job_id=scheduler_job_id,
                    work_directory=work_directory,
                ),
            )
            return True
        finally:
            if "heartbeat_stop_event" in locals():
                self._stop_background_heartbeat(heartbeat_stop_event, locals().get("heartbeat_thread"))

    def _run_task(self, task_type: str, task_id: str, params: dict[str, Any], progress_callback) -> dict[str, Any]:
        method_name = TASK_TYPE_TO_METHOD.get(task_type)
        if not method_name:
            raise ValueError(f"unsupported task type: {task_type}")

        worker_method = getattr(self.execution_worker, method_name)
        return asyncio.run(worker_method(task_id, params, progress_callback))

    @staticmethod
    def _has_nonempty_file(path: Path) -> bool:
        return path.exists() and path.is_file() and path.stat().st_size > 0

    def _validate_success_result(self, task_type: str, result: dict[str, Any]) -> tuple[Optional[str], dict[str, Any]]:
        work_directory = result.get("work_directory")
        if not work_directory:
            return None, {}

        work_dir = Path(str(work_directory))
        if task_type == "scf_calculation":
            missing_outputs = [
                filename
                for filename in ("CHGCAR",)
                if not self._has_nonempty_file(work_dir / filename)
            ]
            if missing_outputs:
                return (
                    "SCF任务未生成有效的 CHGCAR，当前结果不能作为后续 DOS / 能带计算的上游输入",
                    {
                        "validation_stage": "result_validation",
                        "missing_outputs": missing_outputs,
                    },
                )

        return None, {}

    def _classify_result(self, task_type: str, result: dict[str, Any]) -> tuple[bool, Optional[str], dict[str, Any]]:
        succeeded = result.get("success") is True or (
            task_type == "md_calculation" and result.get("convergence")
        )
        if not succeeded:
            return (
                False,
                result.get("error") or result.get("error_message") or "任务执行失败",
                {},
            )

        validation_error, extras = self._validate_success_result(task_type, result)
        if validation_error:
            return False, validation_error, extras
        return True, None, {}

    def _is_success(self, task_type: str, result: dict[str, Any]) -> bool:
        succeeded, _, _ = self._classify_result(task_type, result)
        return succeeded

    def _call_control_plane_with_retries(
        self,
        operation_name: str,
        func,
        attempts: Optional[int] = None,
    ):
        retry_attempts = attempts or settings.worker_control_plane_retry_attempts
        last_error: Optional[Exception] = None
        for attempt in range(1, retry_attempts + 1):
            try:
                with self._control_plane_lock:
                    return func()
            except Exception as exc:
                last_error = exc
                if not self._is_retryable_control_plane_error(exc) or attempt >= retry_attempts:
                    break
                delay_seconds = self._retry_delay_seconds(attempt)
                logger.warning(
                    "worker %s %s failed on attempt %s/%s: %s; retrying in %.1fs",
                    self.worker_id,
                    operation_name,
                    attempt,
                    retry_attempts,
                    exc,
                    delay_seconds,
                )
                time.sleep(delay_seconds)
        assert last_error is not None
        raise last_error

    def _start_background_heartbeat(
        self,
        task_id: str,
        lease_token: str,
        runtime_state: dict[str, Any],
    ) -> tuple[threading.Event, Optional[threading.Thread]]:
        interval_seconds = max(0.0, float(settings.worker_background_heartbeat_interval_seconds))
        stop_event = threading.Event()
        if interval_seconds <= 0:
            return stop_event, None

        def _pump() -> None:
            while not stop_event.wait(interval_seconds):
                try:
                    heartbeat = self._call_control_plane_with_retries(
                        "heartbeat",
                        lambda: self.control_plane_client.heartbeat(task_id, lease_token),
                    )
                    if heartbeat.get("cancel_requested"):
                        runtime_state["cancel_requested"] = True
                except Exception as exc:
                    if not stop_event.is_set():
                        logger.warning(
                            "worker %s background heartbeat failed for task %s: %s",
                            self.worker_id,
                            task_id,
                            exc,
                        )

        thread = threading.Thread(
            target=_pump,
            name=f"control-plane-heartbeat-{task_id}",
            daemon=True,
        )
        thread.start()
        return stop_event, thread

    @staticmethod
    def _stop_background_heartbeat(
        stop_event: threading.Event,
        heartbeat_thread: Optional[threading.Thread],
    ) -> None:
        stop_event.set()
        if heartbeat_thread is not None:
            heartbeat_thread.join(timeout=1.0)

    def _is_retryable_control_plane_error(self, exc: Exception) -> bool:
        if isinstance(exc, requests.exceptions.HTTPError):
            response = getattr(exc, "response", None)
            if response is None:
                return True
            return response.status_code in {408, 429, 500, 502, 503, 504}
        return isinstance(
            exc,
            (
                requests.exceptions.ConnectionError,
                requests.exceptions.Timeout,
                requests.exceptions.RequestException,
                OSError,
            ),
        )

    def _retry_delay_seconds(self, attempt: int) -> float:
        base = settings.worker_control_plane_retry_backoff_seconds
        max_delay = settings.worker_control_plane_retry_max_backoff_seconds
        return min(base * (2 ** (attempt - 1)), max_delay)

    def _ack_cancel(self, task_id: str, lease_token: str) -> bool:
        try:
            self._call_control_plane_with_retries(
                "cancel_ack",
                lambda: self.control_plane_client.cancel_ack(task_id, lease_token),
                attempts=settings.worker_control_plane_final_retry_attempts,
            )
            self._delete_runtime_state(task_id)
            return True
        except Exception as exc:
            logger.warning("worker %s failed to acknowledge cancel for task %s: %s", self.worker_id, task_id, exc)
            return False

    def _report_completion(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> bool:
        pending_path = self._write_pending_completion(task_id, lease_token, payload)
        try:
            self._call_control_plane_with_retries(
                "complete",
                lambda: self.control_plane_client.complete(task_id, lease_token, payload),
                attempts=settings.worker_control_plane_final_retry_attempts,
            )
        except Exception as exc:
            logger.warning(
                "worker %s failed to report completion for task %s after retries; completion payload kept at %s: %s",
                self.worker_id,
                task_id,
                pending_path,
                exc,
            )
            return False
        self._delete_pending_file(pending_path)
        self._delete_runtime_state(task_id)
        return True

    def _report_failure(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> bool:
        try:
            self._call_control_plane_with_retries(
                "fail",
                lambda: self.control_plane_client.fail(task_id, lease_token, payload),
                attempts=settings.worker_control_plane_final_retry_attempts,
            )
            self._delete_runtime_state(task_id)
            return True
        except Exception as exc:
            logger.warning("worker %s failed to report failure for task %s after retries: %s", self.worker_id, task_id, exc)
            return False

    def _pending_completion_dir(self) -> Path:
        preferred_dir = Path(settings.vasp_work_root) / ".control-plane-pending" / self.worker_id
        try:
            preferred_dir.mkdir(parents=True, exist_ok=True)
            return preferred_dir
        except OSError:
            fallback_dir = Path(tempfile.gettempdir()) / "mcp-ic-tool-control-plane-pending" / self.worker_id
            fallback_dir.mkdir(parents=True, exist_ok=True)
            logger.warning(
                "worker %s could not use pending completion dir %s; falling back to %s",
                self.worker_id,
                preferred_dir,
                fallback_dir,
            )
            return fallback_dir

    def _write_pending_completion(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> Path:
        pending_dir = self._pending_completion_dir()
        pending_path = pending_dir / f"{task_id}.json"
        pending_path.write_text(
            json.dumps(
                {
                    "task_id": task_id,
                    "lease_token": lease_token,
                    "payload": payload,
                },
                ensure_ascii=False,
            ),
            encoding="utf-8",
        )
        return pending_path

    def _drain_pending_completions(self) -> bool:
        pending_dir = self._pending_completion_dir()
        if not pending_dir.exists():
            return True

        for pending_path in sorted(pending_dir.glob("*.json")):
            try:
                record = json.loads(pending_path.read_text(encoding="utf-8"))
            except Exception as exc:
                logger.warning("worker %s dropping unreadable pending completion %s: %s", self.worker_id, pending_path, exc)
                self._delete_pending_file(pending_path)
                continue

            try:
                self._call_control_plane_with_retries(
                    "complete",
                    lambda: self.control_plane_client.complete(
                        record["task_id"],
                        record["lease_token"],
                        record["payload"],
                    ),
                    attempts=settings.worker_control_plane_final_retry_attempts,
                )
            except Exception as exc:
                logger.warning(
                    "worker %s could not replay pending completion %s yet: %s",
                    self.worker_id,
                    pending_path,
                    exc,
                )
                return False
            self._delete_pending_file(pending_path)
            self._delete_runtime_state(str(record.get("task_id")))

        return True

    @staticmethod
    def _delete_pending_file(pending_path: Path) -> None:
        try:
            pending_path.unlink(missing_ok=True)
        except Exception:
            logger.warning("failed to delete pending completion file %s", pending_path, exc_info=True)

    def _runtime_state_dir(self) -> Path:
        preferred_dir = Path(settings.vasp_work_root) / ".control-plane-runtime" / self.worker_id
        try:
            preferred_dir.mkdir(parents=True, exist_ok=True)
            return preferred_dir
        except OSError:
            fallback_dir = Path(tempfile.gettempdir()) / "mcp-ic-tool-control-plane-runtime" / self.worker_id
            fallback_dir.mkdir(parents=True, exist_ok=True)
            logger.warning(
                "worker %s could not use runtime state dir %s; falling back to %s",
                self.worker_id,
                preferred_dir,
                fallback_dir,
            )
            return fallback_dir

    def _runtime_state_path(self, task_id: str) -> Path:
        return self._runtime_state_dir() / f"{task_id}.json"

    def _default_work_directory(self, task_id: str) -> Path:
        base_work_dir = Path(getattr(self.execution_worker, "base_work_dir", settings.vasp_work_root))
        return base_work_dir / task_id

    def _load_runtime_state(self, state_path: Path) -> Optional[dict[str, Any]]:
        try:
            return json.loads(state_path.read_text(encoding="utf-8"))
        except FileNotFoundError:
            return None
        except Exception as exc:
            logger.warning("worker %s dropping unreadable runtime state %s: %s", self.worker_id, state_path, exc)
            self._delete_file(state_path)
            return None

    def _upsert_runtime_state(self, task_id: str, **updates: Any) -> Path:
        state_path = self._runtime_state_path(task_id)
        current = {}
        if state_path.exists():
            current = self._load_runtime_state(state_path) or {}
        current.update({k: v for k, v in updates.items() if v is not None})
        current.setdefault("task_id", task_id)
        state_path.write_text(json.dumps(current, ensure_ascii=False), encoding="utf-8")
        return state_path

    def _delete_runtime_state(self, task_id: str) -> None:
        if not task_id:
            return
        self._delete_file(self._runtime_state_path(task_id))

    @staticmethod
    def _delete_file(path: Path) -> None:
        try:
            path.unlink(missing_ok=True)
        except Exception:
            logger.warning("failed to delete file %s", path, exc_info=True)

    def _reconcile_runtime_tasks(self) -> Optional[bool]:
        runtime_dir = self._runtime_state_dir()
        state_paths = sorted(runtime_dir.glob("*.json"))
        if not state_paths:
            return None

        for state_path in state_paths:
            state = self._load_runtime_state(state_path)
            if state is None:
                continue
            reconciled = self._reconcile_runtime_state(state)
            if reconciled is None:
                continue
            return reconciled
        return None

    def _reconcile_runtime_state(self, state: dict[str, Any]) -> Optional[bool]:
        task_id = str(state["task_id"])
        task_type = str(state.get("task_type") or "")
        params = dict(state.get("params") or {})
        work_dir = Path(str(state.get("work_directory") or self._default_work_directory(task_id)))
        scheduler_job_id = str(state.get("scheduler_job_id") or "") or None

        lease_token = state.get("lease_token")
        cancel_requested = False
        task_status = str(state.get("status") or "running")

        if lease_token:
            try:
                heartbeat = self._call_control_plane_with_retries(
                    "heartbeat",
                    lambda: self.control_plane_client.heartbeat(task_id, str(lease_token)),
                )
                cancel_requested = bool(heartbeat.get("cancel_requested"))
                task_status = str(heartbeat.get("status") or task_status)
            except requests.exceptions.HTTPError as exc:
                response = getattr(exc, "response", None)
                if response is None or response.status_code not in {404, 409}:
                    raise
                lease_token = None

        if not lease_token:
            try:
                resumed = self._call_control_plane_with_retries(
                    "resume",
                    lambda: self.control_plane_client.resume(task_id),
                )
            except requests.exceptions.HTTPError as exc:
                response = getattr(exc, "response", None)
                if response is None or response.status_code not in {404, 409}:
                    raise
                logger.info(
                    "worker %s dropping stale runtime state for finalized task %s after resume returned %s",
                    self.worker_id,
                    task_id,
                    response.status_code,
                )
                self._delete_runtime_state(task_id)
                return None
            lease_token = resumed["lease_token"]
            task_status = str(resumed.get("status") or task_status)
            task_type = str(resumed.get("task_type") or task_type)
            params = dict(resumed.get("params") or params)
            cancel_requested = task_status == "cancel_requested"
            self._upsert_runtime_state(
                task_id,
                task_type=task_type,
                params=params,
                lease_token=lease_token,
                work_directory=str(work_dir),
                scheduler_job_id=scheduler_job_id,
                status=task_status,
            )

        if task_status == "leased":
            self._call_control_plane_with_retries(
                "mark_running",
                lambda: self.control_plane_client.mark_running(task_id, str(lease_token)),
            )
            task_status = "running"
            self._upsert_runtime_state(task_id, status=task_status, lease_token=lease_token)

        if cancel_requested:
            if scheduler_job_id:
                subprocess.run(["scancel", scheduler_job_id], check=False)
            self._ack_cancel(task_id, str(lease_token))
            return True

        scheduler_state = self._query_scheduler_state(scheduler_job_id, work_dir=work_dir)
        self._upsert_runtime_state(task_id, scheduler_state=scheduler_state, lease_token=lease_token)

        if scheduler_state in ACTIVE_SCHEDULER_STATES:
            return True

        if scheduler_state == "COMPLETED":
            result = self._analyze_existing_results(
                task_type=task_type,
                work_dir=work_dir,
                scheduler_job_id=scheduler_job_id,
                started_at=state.get("started_at"),
            )
            result.setdefault("work_directory", str(work_dir))
            payload = {
                "result_data": result,
                "artifact_manifest": self._stage_public_artifacts(
                    task_id,
                    self._collect_artifact_manifest(str(work_dir)),
                ),
            }
            succeeded, error_message, failure_extras = self._classify_result(task_type, result)
            if succeeded:
                self._report_completion(task_id, str(lease_token), payload)
            else:
                self._report_failure(
                    task_id,
                    str(lease_token),
                    self._build_guided_failure_payload(
                        task_id=task_id,
                        task_type=task_type,
                        error_message=error_message or "恢复结果分析失败",
                        scheduler_job_id=scheduler_job_id,
                        work_directory=str(work_dir),
                        extras=failure_extras,
                    ),
                )
            return True

        if scheduler_state in FAILED_SCHEDULER_STATES:
            self._report_failure(
                task_id,
                str(lease_token),
                self._build_guided_failure_payload(
                    task_id=task_id,
                    task_type=task_type,
                    error_message=self._build_scheduler_failure_message(scheduler_state, work_dir),
                    scheduler_job_id=scheduler_job_id,
                    work_directory=str(work_dir),
                ),
            )
            return True

        logger.warning(
            "worker %s could not reconcile task %s yet; scheduler state is %s",
            self.worker_id,
            task_id,
            scheduler_state,
        )
        return False

    def _query_scheduler_state(self, scheduler_job_id: Optional[str], work_dir: Optional[Path] = None) -> str:
        if not scheduler_job_id:
            return "UNKNOWN"

        try:
            process = subprocess.run(
                ["squeue", "-j", str(scheduler_job_id), "--noheader", "--format=%T"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            state = self._normalize_scheduler_state(process.stdout)
            if state:
                return state
        except Exception:
            logger.warning("worker %s failed to query squeue for job %s", self.worker_id, scheduler_job_id, exc_info=True)

        try:
            process = subprocess.run(
                ["sacct", "-j", str(scheduler_job_id), "--format=State", "--noheader", "-X"],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                check=False,
            )
            state = self._normalize_scheduler_state(process.stdout)
            if state:
                return state
        except Exception:
            logger.warning("worker %s failed to query sacct for job %s", self.worker_id, scheduler_job_id, exc_info=True)

        return "UNKNOWN"

    @staticmethod
    def _normalize_scheduler_state(raw_output: Any) -> Optional[str]:
        if not raw_output:
            return None
        for line in str(raw_output).splitlines():
            value = line.strip()
            if not value:
                continue
            normalized = value.split()[0].split("+")[0].upper()
            return normalized
        return None

    def _analyze_existing_results(
        self,
        *,
        task_type: str,
        work_dir: Path,
        scheduler_job_id: Optional[str],
        started_at: Optional[str],
    ) -> dict[str, Any]:
        analyzer_name = TASK_TYPE_TO_ANALYZER.get(task_type)
        if not analyzer_name:
            raise ValueError(f"unsupported task type for reconcile: {task_type}")

        analyzer = getattr(self.execution_worker, analyzer_name)
        computation_time = None
        if started_at:
            try:
                started = datetime.fromisoformat(str(started_at))
                computation_time = max((datetime.now(timezone.utc) - started).total_seconds(), 0.0)
            except ValueError:
                computation_time = None

        run_result = {
            "success": True,
            "process_id": int(scheduler_job_id) if scheduler_job_id and str(scheduler_job_id).isdigit() else None,
            "computation_time": computation_time,
            "work_directory": str(work_dir),
        }
        return asyncio.run(analyzer(work_dir, run_result))

    def _build_scheduler_failure_message(self, scheduler_state: str, work_dir: Path) -> str:
        error_msg = f"SLURM 作业以状态 {scheduler_state} 结束"
        guidance = derive_failure_guidance_from_workdir(work_dir)
        if guidance.reason:
            error_msg += f"\n诊断: {guidance.reason}"
        if guidance.suggested_action:
            error_msg += f"\n建议给智能体的下一步: {guidance.suggested_action}"
        error_files = sorted(work_dir.glob("*.err"))
        if error_files:
            try:
                error_content = error_files[0].read_text(encoding="utf-8", errors="ignore").strip()
            except Exception:
                error_content = ""
            if error_content:
                error_msg += f"\n错误日志:\n{error_content}"
        return error_msg

    def _collect_artifact_manifest(self, work_directory: Any) -> list[dict[str, str]]:
        if not work_directory:
            return []

        work_dir = Path(str(work_directory))
        if not work_dir.exists():
            return []

        manifest = []
        nested_output_dirs = {"data", "plots", "figures", "analysis"}

        def append_file(file_path: Path) -> None:
            rel_path = file_path.relative_to(work_dir).as_posix()
            content_type = mimetypes.guess_type(file_path.name)[0]
            item = {
                "filename": rel_path,
                "local_path": str(file_path),
            }
            if content_type:
                item["content_type"] = content_type
            manifest.append(item)

        for item in sorted(work_dir.iterdir()):
            if item.is_file():
                append_file(item)
                continue
            if not item.is_dir():
                continue
            item_name = item.name.lower()
            if (
                item_name.endswith("_output")
                or item_name.endswith("_analysis")
                or item_name in nested_output_dirs
            ):
                for nested_file in sorted(p for p in item.rglob("*") if p.is_file()):
                    append_file(nested_file)
        return manifest

    def _collect_and_stage_public_artifacts(self, task_id: str, result: dict[str, Any]) -> list[dict[str, Any]]:
        work_directory = result.get("work_directory") or str(self._default_work_directory(task_id))
        return self._stage_public_artifacts(
            task_id,
            self._collect_artifact_manifest(work_directory),
        )

    def _stage_public_artifacts(self, task_id: str, artifact_manifest: list[dict[str, Any]]) -> list[dict[str, Any]]:
        staged: list[dict[str, Any]] = []
        for item in artifact_manifest:
            entry = dict(item)
            local_path = entry.get("local_path")
            filename = str(entry.get("filename") or "")
            if not local_path or not filename:
                staged.append(entry)
                continue

            file_path = Path(str(local_path))
            if not file_path.exists():
                staged.append(entry)
                continue

            entry.setdefault("size_bytes", float(file_path.stat().st_size))
            if self._should_upload_public_artifact(file_path, entry.get("content_type")):
                try:
                    uploaded = self.control_plane_client.upload_public_artifact(
                        task_id=task_id,
                        filename=filename,
                        local_path=str(file_path),
                        content_type=entry.get("content_type"),
                    )
                    entry.update(uploaded)
                    entry["storage_backend"] = "local_public"
                    staged.append(entry)
                    continue
                except Exception as exc:
                    logger.warning(
                        "worker %s failed to upload public artifact %s for task %s: %s",
                        self.worker_id,
                        filename,
                        task_id,
                        exc,
                    )

            entry["storage_backend"] = "local"
            entry["storage_key"] = str(file_path)
            entry.pop("object_key", None)
            staged.append(entry)
        self._rewrite_uploaded_public_json_artifacts(task_id, staged)
        return staged

    @classmethod
    def _rewrite_public_json_string(
        cls,
        value: str,
        *,
        urls_by_local_path: dict[str, str],
        urls_by_filename: dict[str, str],
    ) -> Optional[str]:
        text = str(value or "")
        if not text:
            return text

        if text.startswith("/"):
            try:
                resolved = str(Path(text).resolve())
            except OSError:
                resolved = text
            if resolved in urls_by_local_path:
                return urls_by_local_path[resolved]

        basename = os.path.basename(urlsplit(text).path).lower()
        if basename and basename in urls_by_filename:
            return urls_by_filename[basename]

        if text.startswith("/"):
            return None
        return text

    @classmethod
    def _sanitize_public_json_payload(
        cls,
        payload: Any,
        *,
        urls_by_local_path: dict[str, str],
        urls_by_filename: dict[str, str],
    ) -> Any:
        if isinstance(payload, dict):
            sanitized: dict[str, Any] = {}
            discovered_html_url: Optional[str] = None
            for key, value in payload.items():
                if key in PUBLIC_JSON_PRIVATE_KEYS:
                    continue
                if key in PUBLIC_JSON_REPORT_KEYS:
                    rewritten = None
                    if isinstance(value, str):
                        rewritten = cls._rewrite_public_json_string(
                            value,
                            urls_by_local_path=urls_by_local_path,
                            urls_by_filename=urls_by_filename,
                        )
                    if rewritten and not discovered_html_url:
                        discovered_html_url = rewritten
                    continue

                sanitized_value = cls._sanitize_public_json_payload(
                    value,
                    urls_by_local_path=urls_by_local_path,
                    urls_by_filename=urls_by_filename,
                )
                if sanitized_value is None and isinstance(value, str) and value.startswith("/"):
                    continue
                sanitized[key] = sanitized_value

            if discovered_html_url and not sanitized.get("html_report_url"):
                sanitized["html_report_url"] = discovered_html_url
            return sanitized

        if isinstance(payload, list):
            sanitized_items: list[Any] = []
            for item in payload:
                sanitized_item = cls._sanitize_public_json_payload(
                    item,
                    urls_by_local_path=urls_by_local_path,
                    urls_by_filename=urls_by_filename,
                )
                if sanitized_item is None and isinstance(item, str) and item.startswith("/"):
                    continue
                sanitized_items.append(sanitized_item)
            return sanitized_items

        if isinstance(payload, str):
            return cls._rewrite_public_json_string(
                payload,
                urls_by_local_path=urls_by_local_path,
                urls_by_filename=urls_by_filename,
            )
        return payload

    def _rewrite_uploaded_public_json_artifacts(self, task_id: str, staged: list[dict[str, Any]]) -> None:
        urls_by_local_path: dict[str, str] = {}
        urls_by_filename: dict[str, str] = {}
        for entry in staged:
            if entry.get("storage_backend") != "local_public":
                continue
            download_url = str(entry.get("download_url") or "").strip()
            local_path = str(entry.get("local_path") or "").strip()
            filename = str(entry.get("filename") or "").strip()
            if download_url and local_path:
                try:
                    urls_by_local_path[str(Path(local_path).resolve())] = download_url
                except OSError:
                    pass
            if download_url and filename:
                urls_by_filename.setdefault(os.path.basename(filename).lower(), download_url)

        if not urls_by_local_path and not urls_by_filename:
            return

        for entry in staged:
            filename = str(entry.get("filename") or "").strip()
            local_path = str(entry.get("local_path") or "").strip()
            if entry.get("storage_backend") != "local_public" or not filename.lower().endswith(".json") or not local_path:
                continue

            source_path = Path(local_path)
            if not source_path.exists():
                continue

            try:
                original_payload = json.loads(source_path.read_text(encoding="utf-8"))
            except Exception:
                continue

            sanitized_payload = self._sanitize_public_json_payload(
                original_payload,
                urls_by_local_path=urls_by_local_path,
                urls_by_filename=urls_by_filename,
            )
            if sanitized_payload == original_payload:
                continue

            with tempfile.NamedTemporaryFile("w", encoding="utf-8", suffix=".json", delete=False) as fh:
                json.dump(sanitized_payload, fh, ensure_ascii=False, indent=2)
                fh.write("\n")
                temp_path = Path(fh.name)

            try:
                uploaded = self.control_plane_client.upload_public_artifact(
                    task_id=task_id,
                    filename=filename,
                    local_path=str(temp_path),
                    content_type=entry.get("content_type") or "application/json",
                )
                entry.update(uploaded)
                entry["storage_backend"] = "local_public"
            finally:
                temp_path.unlink(missing_ok=True)

    @staticmethod
    def _should_upload_public_artifact(file_path: Path, content_type: Optional[str]) -> bool:
        if file_path.stat().st_size > settings.public_artifact_max_upload_bytes:
            return False

        filename = file_path.name.lower()
        if filename in PUBLIC_ARTIFACT_FILENAMES:
            return True

        suffix = file_path.suffix.lower()
        if suffix not in PUBLIC_ARTIFACT_EXTENSIONS:
            return False
        if content_type and (
            content_type.startswith("text/")
            or content_type.startswith("image/")
            or content_type in {"application/json", "application/pdf"}
        ):
            return True
        return suffix in PUBLIC_ARTIFACT_EXTENSIONS
