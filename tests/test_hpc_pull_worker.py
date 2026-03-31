from pathlib import Path

import pytest
import requests


class FakeControlPlaneClient:
    def __init__(self, tasks=None, cancel_on_heartbeat=False, operation_errors=None):
        self.tasks = list(tasks or [])
        self.cancel_on_heartbeat = cancel_on_heartbeat
        self.operation_errors = {name: list(errors) for name, errors in (operation_errors or {}).items()}
        self.claim_calls = 0
        self.mark_running_calls = 0
        self.heartbeat_calls = 0
        self.complete_calls = 0
        self.fail_calls = 0
        self.cancel_ack_calls = 0
        self.heartbeats = []
        self.marked_running = []
        self.completed_tasks = []
        self.failed_tasks = []
        self.canceled_tasks = []
        self.complete_payloads = []
        self.fail_payloads = []

    def _raise_if_needed(self, operation_name: str):
        errors = self.operation_errors.get(operation_name) or []
        if errors:
            raise errors.pop(0)

    def claim(self):
        self.claim_calls += 1
        self._raise_if_needed("claim")
        if not self.tasks:
            return None
        return self.tasks.pop(0)

    def mark_running(self, task_id: str, lease_token: str):
        self.mark_running_calls += 1
        self._raise_if_needed("mark_running")
        self.marked_running.append((task_id, lease_token))

    def heartbeat(self, task_id: str, lease_token: str):
        self.heartbeat_calls += 1
        self._raise_if_needed("heartbeat")
        self.heartbeats.append((task_id, lease_token))
        return {"cancel_requested": self.cancel_on_heartbeat}

    def complete(self, task_id: str, lease_token: str, payload: dict):
        self.complete_calls += 1
        self._raise_if_needed("complete")
        self.completed_tasks.append(task_id)
        self.complete_payloads.append(payload)

    def fail(self, task_id: str, lease_token: str, payload: dict):
        self.fail_calls += 1
        self._raise_if_needed("fail")
        self.failed_tasks.append(task_id)
        self.fail_payloads.append(payload)

    def cancel_ack(self, task_id: str, lease_token: str):
        self.cancel_ack_calls += 1
        self._raise_if_needed("cancel_ack")
        self.canceled_tasks.append(task_id)


class FakeExecutionWorker:
    def __init__(self, result=None, exc: Exception | None = None):
        self.result = result or {"success": True, "work_directory": "/tmp/work"}
        self.exc = exc

    async def run_structure_optimization(self, task_id, params, progress_callback=None):
        if progress_callback:
            await progress_callback(10, "queued")
            await progress_callback(40, "running", pid=12345)
        if self.exc is not None:
            raise self.exc
        return self.result


def test_pull_worker_claims_runs_uploads_and_completes(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("ok", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {"formula": "Li2O"},
                "upstream_artifact_manifest": [],
            }
        ]
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    assert client.completed_tasks == ["task-1"]
    assert client.marked_running == [("task-1", "lease-1")]
    assert client.complete_payloads[0]["artifact_manifest"][0]["filename"] == "result.log"


def test_pull_worker_sends_heartbeats_while_job_is_running():
    from src.vasp_server.hpc_pull_worker import PullWorker

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {"formula": "Li2O"},
                "upstream_artifact_manifest": [],
            }
        ]
    )
    worker = FakeExecutionWorker()

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    assert len(client.heartbeats) >= 2


def test_fail_path_uploads_logs_then_reports_failed(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("boom", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {"formula": "Li2O"},
                "upstream_artifact_manifest": [],
            }
        ]
    )
    worker = FakeExecutionWorker(result={"success": False, "error": "boom", "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    assert client.failed_tasks == ["task-1"]
    assert client.fail_payloads[0]["artifact_manifest"][0]["filename"] == "result.log"


def test_cancel_path_calls_scancel_and_reports_canceled(monkeypatch):
    from src.vasp_server.hpc_pull_worker import PullWorker

    calls = []

    def fake_run(cmd, check=False):
        calls.append(cmd)
        return None

    monkeypatch.setattr("subprocess.run", fake_run)

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {"formula": "Li2O"},
                "upstream_artifact_manifest": [],
            }
        ],
        cancel_on_heartbeat=True,
    )
    worker = FakeExecutionWorker()

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    assert calls == [["scancel", "12345"]]
    assert client.canceled_tasks == ["task-1"]


def test_pull_worker_poll_loop_claims_repeatedly(monkeypatch):
    from src.vasp_server.hpc_pull_worker import PullWorker

    sleeps = []
    monkeypatch.setattr("time.sleep", lambda seconds: sleeps.append(seconds))

    client = FakeControlPlaneClient(tasks=[])
    pull_worker = PullWorker(control_plane_client=client, execution_worker=FakeExecutionWorker(), poll_interval_seconds=2)
    pull_worker.run_loop(iterations=2)

    assert client.claim_calls == 2
    assert sleeps == [2, 2]


def test_pull_worker_retries_claim_failures(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr("time.sleep", lambda seconds: None)
    monkeypatch.setattr(settings, "worker_control_plane_retry_attempts", 4)
    monkeypatch.setattr(settings, "worker_control_plane_retry_backoff_seconds", 0.01)
    monkeypatch.setattr(settings, "worker_control_plane_retry_max_backoff_seconds", 0.01)

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("ok", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {},
            }
        ],
        operation_errors={
            "claim": [
                requests.exceptions.ConnectTimeout("attempt-1"),
                requests.exceptions.ConnectTimeout("attempt-2"),
            ]
        },
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)

    assert pull_worker.run_once() is True
    assert client.claim_calls == 3
    assert client.completed_tasks == ["task-1"]


def test_pull_worker_retries_complete_before_succeeding(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr("time.sleep", lambda seconds: None)
    monkeypatch.setattr(settings, "worker_control_plane_final_retry_attempts", 5)
    monkeypatch.setattr(settings, "worker_control_plane_retry_backoff_seconds", 0.01)
    monkeypatch.setattr(settings, "worker_control_plane_retry_max_backoff_seconds", 0.01)
    monkeypatch.setattr(settings, "vasp_work_root", str(tmp_path))

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("ok", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {},
            }
        ],
        operation_errors={
            "complete": [
                requests.exceptions.ConnectionError("attempt-1"),
                requests.exceptions.ConnectionError("attempt-2"),
            ]
        },
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker, worker_id="worker-a")

    assert pull_worker.run_once() is True
    assert client.complete_calls == 3
    assert client.completed_tasks == ["task-1"]
    assert not list((tmp_path / ".control-plane-pending" / "worker-a").glob("*.json"))


def test_pull_worker_replays_pending_completion_after_network_failure(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr("time.sleep", lambda seconds: None)
    monkeypatch.setattr(settings, "worker_control_plane_final_retry_attempts", 3)
    monkeypatch.setattr(settings, "worker_control_plane_retry_backoff_seconds", 0.01)
    monkeypatch.setattr(settings, "worker_control_plane_retry_max_backoff_seconds", 0.01)
    monkeypatch.setattr(settings, "vasp_work_root", str(tmp_path))

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("ok", encoding="utf-8")
    failing_client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {},
            }
        ],
        operation_errors={
            "complete": [
                requests.exceptions.ConnectionError("attempt-1"),
                requests.exceptions.ConnectionError("attempt-2"),
                requests.exceptions.ConnectionError("attempt-3"),
            ]
        },
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=failing_client, execution_worker=worker, worker_id="worker-a")

    assert pull_worker.run_once() is True
    pending_files = list((tmp_path / ".control-plane-pending" / "worker-a").glob("*.json"))
    assert len(pending_files) == 1
    assert failing_client.completed_tasks == []
    assert failing_client.failed_tasks == []

    recovering_client = FakeControlPlaneClient(tasks=[])
    recovering_worker = PullWorker(
        control_plane_client=recovering_client,
        execution_worker=FakeExecutionWorker(),
        worker_id="worker-a",
    )

    assert recovering_worker.run_once() is False
    assert recovering_client.completed_tasks == ["task-1"]
    assert not list((tmp_path / ".control-plane-pending" / "worker-a").glob("*.json"))
