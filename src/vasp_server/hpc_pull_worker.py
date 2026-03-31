import asyncio
import json
import logging
import subprocess
import tempfile
import time
from pathlib import Path
from typing import Any, Optional

import requests

from .settings import settings
from .vasp_worker import VaspWorker
from .worker_client import ControlPlaneClient


TASK_TYPE_TO_METHOD = {
    "structure_optimization": "run_structure_optimization",
    "scf_calculation": "run_scf_calculation",
    "dos_calculation": "run_dos_calculation",
    "md_calculation": "run_md_calculation",
    "band_structure": "run_band_structure_calculation",
    "neb_calculation": "run_neb_calculation",
    "phonon_calculation": "run_phonon_calculation",
    "custom_calculation": "run_custom_calculation",
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

    def run_loop(self, iterations: Optional[int] = None) -> None:
        completed = 0
        while iterations is None or completed < iterations:
            self.run_once()
            completed += 1
            time.sleep(self.poll_interval_seconds)

    def run_once(self) -> bool:
        if not self._drain_pending_completions():
            return False

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

        scheduler_job_id: Optional[str] = None
        cancel_requested = False

        async def progress_callback(progress: int, message: str, pid: Optional[int] = None):
            nonlocal scheduler_job_id, cancel_requested
            if pid is not None:
                scheduler_job_id = str(pid)
            heartbeat = self._call_control_plane_with_retries(
                "heartbeat",
                lambda: self.control_plane_client.heartbeat(task_id, lease_token),
            )
            if heartbeat.get("cancel_requested"):
                cancel_requested = True
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
            result = self._run_task(task_type, task_id, params, progress_callback)
            if cancel_requested:
                self._ack_cancel(task_id, lease_token)
                return True
            if not self._is_success(task_type, result):
                payload = {
                    "error_message": result.get("error") or result.get("error_message") or "任务执行失败",
                    "artifact_manifest": self._collect_artifact_manifest(result.get("work_directory")),
                }
                self._report_failure(task_id, lease_token, payload)
                return True

            payload = {
                "result_data": result,
                "artifact_manifest": self._collect_artifact_manifest(result.get("work_directory")),
            }
            self._report_completion(task_id, lease_token, payload)
            return True
        except TaskCanceled:
            return True
        except Exception as exc:
            self._report_failure(
                task_id,
                lease_token,
                {
                    "error_message": str(exc),
                    "artifact_manifest": [],
                },
            )
            return True

    def _run_task(self, task_type: str, task_id: str, params: dict[str, Any], progress_callback) -> dict[str, Any]:
        method_name = TASK_TYPE_TO_METHOD.get(task_type)
        if not method_name:
            raise ValueError(f"unsupported task type: {task_type}")

        worker_method = getattr(self.execution_worker, method_name)
        return asyncio.run(worker_method(task_id, params, progress_callback))

    def _is_success(self, task_type: str, result: dict[str, Any]) -> bool:
        if result.get("success") is True:
            return True
        if task_type == "md_calculation" and result.get("convergence"):
            return True
        return False

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
        return True

    def _report_failure(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> bool:
        try:
            self._call_control_plane_with_retries(
                "fail",
                lambda: self.control_plane_client.fail(task_id, lease_token, payload),
                attempts=settings.worker_control_plane_final_retry_attempts,
            )
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

        return True

    @staticmethod
    def _delete_pending_file(pending_path: Path) -> None:
        try:
            pending_path.unlink(missing_ok=True)
        except Exception:
            logger.warning("failed to delete pending completion file %s", pending_path, exc_info=True)

    def _collect_artifact_manifest(self, work_directory: Any) -> list[dict[str, str]]:
        if not work_directory:
            return []

        work_dir = Path(str(work_directory))
        if not work_dir.exists():
            return []

        manifest = []
        for item in sorted(work_dir.iterdir()):
            if item.is_file():
                manifest.append(
                    {
                        "filename": item.name,
                        "local_path": str(item),
                    }
                )
        return manifest
