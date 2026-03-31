import asyncio
from pathlib import Path

import httpx

from src.vasp_server.hpc_pull_worker import PullWorker
from src.vasp_server.task_manager.models import Task


def _request(app, method: str, url: str, **kwargs):
    async def _send():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(transport=transport, base_url="http://testserver") as client:
            return await client.request(method, url, **kwargs)

    return asyncio.run(_send())


class DirectControlPlaneClient:
    def __init__(self, task_manager, worker_id="hpc-a", queue_name="default", cancel_task_id=None, cancel_user_id="test_user"):
        self.task_manager = task_manager
        self.worker_id = worker_id
        self.queue_name = queue_name
        self.cancel_task_id = cancel_task_id
        self.cancel_user_id = cancel_user_id
        self._cancel_requested = False

    def claim(self):
        claimed = self.task_manager.claim_next_task(worker_id=self.worker_id, queue_name=self.queue_name)
        if claimed is None:
            return None
        return {
            "task_id": claimed.id,
            "task_type": claimed.task_type,
            "lease_token": claimed.lease_token,
            "params": claimed.params,
        }

    def mark_running(self, task_id: str, lease_token: str):
        self.task_manager.mark_task_running(task_id, lease_token, self.worker_id)

    def heartbeat(self, task_id: str, lease_token: str):
        if self.cancel_task_id == task_id and not self._cancel_requested:
            self.task_manager.request_cancel(task_id, self.cancel_user_id)
            self._cancel_requested = True
        self.task_manager.heartbeat_task(task_id, lease_token, self.worker_id)
        task = self.task_manager.get_task(task_id, self.cancel_user_id) or self.task_manager.get_task(task_id, "test_user")
        return {"cancel_requested": getattr(task, "status", "") == "cancel_requested"}

    def complete(self, task_id: str, lease_token: str, payload: dict):
        self.task_manager.complete_execution(
            task_id=task_id,
            lease_token=lease_token,
            worker_id=self.worker_id,
            result_data=payload.get("result_data"),
            artifact_manifest=payload.get("artifact_manifest"),
        )

    def fail(self, task_id: str, lease_token: str, payload: dict):
        self.task_manager.fail_execution(
            task_id=task_id,
            lease_token=lease_token,
            worker_id=self.worker_id,
            error_message=payload.get("error_message", "failed"),
            artifact_manifest=payload.get("artifact_manifest"),
        )

    def cancel_ack(self, task_id: str, lease_token: str):
        self.task_manager.ack_cancel(task_id, lease_token, self.worker_id)


class FakeStructureWorker:
    def __init__(self, root: Path):
        self.root = root

    async def run_structure_optimization(self, task_id, params, progress_callback=None):
        work_dir = self.root / task_id
        work_dir.mkdir(parents=True, exist_ok=True)
        (work_dir / "report.html").write_text("<html>ok</html>", encoding="utf-8")
        (work_dir / "result.log").write_text("done", encoding="utf-8")
        (work_dir / "CONTCAR").write_text("optimized-structure", encoding="utf-8")
        if progress_callback:
            await progress_callback(10, "claimed")
            await progress_callback(50, "running", pid=12345)
        return {"success": True, "work_directory": str(work_dir), "final_energy": -10.5}


class FakeSCFWorker:
    def __init__(self):
        self.seen_params = None

    async def run_scf_calculation(self, task_id, params, progress_callback=None):
        self.seen_params = dict(params)
        if progress_callback:
            await progress_callback(20, "running", pid=23456)
        return {"success": True, "work_directory": "/tmp/scf-work", "total_energy": -20.0}


class FakeCancelableWorker:
    async def run_structure_optimization(self, task_id, params, progress_callback=None):
        if progress_callback:
            await progress_callback(10, "claimed")
            await progress_callback(50, "running", pid=34567)
            await progress_callback(60, "still-running", pid=34567)
        return {"success": True, "work_directory": "/tmp/canceled-work"}


def test_end_to_end_structure_optimization_flow(tmp_path):
    from src.vasp_server.vasp_server_api import app, task_manager

    submit_response = _request(
        app,
        "POST",
        "/vasp/structure-optimization",
        json={"user_id": "test_user", "cif_url": "https://structures.example.com/Li2O.cif"},
    )
    assert submit_response.status_code == 200
    task_id = submit_response.json()["task_id"]

    pull_worker = PullWorker(
        control_plane_client=DirectControlPlaneClient(task_manager),
        execution_worker=FakeStructureWorker(tmp_path),
    )
    pull_worker.run_once()

    status_response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})
    payload = status_response.json()

    assert status_response.status_code == 200
    assert payload["status"] == "completed"
    assert payload["html_report_url"].startswith("https://")
    assert payload["artifacts"][0]["download_url"].startswith("https://")


def test_end_to_end_scf_stays_on_upstream_queue_and_cancel_running_task(tmp_path, db_session):
    from src.vasp_server.vasp_server_api import app, task_manager

    upstream_task_id = task_manager.submit_task(
        user_id="test_user",
        task_type="structure_optimization",
        params={"cif_url": "https://structures.example.com/Li2O.cif", "queue_name": "hpc-a"},
    )
    task = db_session.get(Task, upstream_task_id)
    task.status = "completed"
    task.analysis_status = "completed"
    db_session.add(task)
    db_session.commit()
    submit_response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={"user_id": "test_user", "optimized_task_id": upstream_task_id},
    )
    assert submit_response.status_code == 200

    other_queue_worker = PullWorker(
        control_plane_client=DirectControlPlaneClient(task_manager, worker_id="hpc-b", queue_name="hpc-b"),
        execution_worker=FakeSCFWorker(),
    )
    assert other_queue_worker.run_once() is False

    scf_worker = FakeSCFWorker()
    PullWorker(
        control_plane_client=DirectControlPlaneClient(task_manager, worker_id="hpc-a", queue_name="hpc-a"),
        execution_worker=scf_worker,
    ).run_once()

    assert scf_worker.seen_params is not None
    assert scf_worker.seen_params["optimized_task_id"] == upstream_task_id
    assert scf_worker.seen_params["queue_name"] == "hpc-a"
    assert "upstream_artifact_manifest" not in scf_worker.seen_params

    cancel_task_id = task_manager.submit_task(
        user_id="test_user",
        task_type="structure_optimization",
        params={"cif_url": "https://structures.example.com/Li2O-cancel.cif"},
    )
    PullWorker(
        control_plane_client=DirectControlPlaneClient(
            task_manager,
            cancel_task_id=cancel_task_id,
            cancel_user_id="test_user",
        ),
        execution_worker=FakeCancelableWorker(),
    ).run_once()

    canceled_task = task_manager.get_task(cancel_task_id, "test_user")
    assert canceled_task is not None
    assert canceled_task.status == "canceled"
