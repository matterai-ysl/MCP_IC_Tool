import asyncio
import subprocess
import time
from pathlib import Path
from typing import Any, Optional

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
        task = self.control_plane_client.claim()
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
            heartbeat = self.control_plane_client.heartbeat(task_id, lease_token)
            if heartbeat.get("cancel_requested"):
                cancel_requested = True
            if cancel_requested and scheduler_job_id:
                subprocess.run(["scancel", scheduler_job_id], check=False)
                self.control_plane_client.cancel_ack(task_id, lease_token)
                raise TaskCanceled(message)
            if cancel_requested and pid is None:
                return
            if cancel_requested:
                if scheduler_job_id:
                    subprocess.run(["scancel", scheduler_job_id], check=False)
                    self.control_plane_client.cancel_ack(task_id, lease_token)
                raise TaskCanceled(message)

        try:
            self.control_plane_client.mark_running(task_id, lease_token)
            result = self._run_task(task_type, task_id, params, progress_callback)
            if cancel_requested:
                self.control_plane_client.cancel_ack(task_id, lease_token)
                return True
            if not self._is_success(task_type, result):
                payload = {
                    "error_message": result.get("error") or result.get("error_message") or "任务执行失败",
                    "artifact_manifest": self._collect_artifact_manifest(result.get("work_directory")),
                }
                self.control_plane_client.fail(task_id, lease_token, payload)
                return True

            payload = {
                "result_data": result,
                "artifact_manifest": self._collect_artifact_manifest(result.get("work_directory")),
            }
            self.control_plane_client.complete(task_id, lease_token, payload)
            return True
        except TaskCanceled:
            return True
        except Exception as exc:
            self.control_plane_client.fail(
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
