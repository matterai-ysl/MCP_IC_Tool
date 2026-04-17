import asyncio
import json
from datetime import datetime, timezone
from pathlib import Path

import pytest
import requests


class FakeControlPlaneClient:
    def __init__(self, tasks=None, cancel_on_heartbeat=False, operation_errors=None, resumed_tasks=None):
        self.tasks = list(tasks or [])
        self.resumed_tasks = dict(resumed_tasks or {})
        self.cancel_on_heartbeat = cancel_on_heartbeat
        self.operation_errors = {name: list(errors) for name, errors in (operation_errors or {}).items()}
        self.claim_calls = 0
        self.resume_calls = 0
        self.mark_running_calls = 0
        self.heartbeat_calls = 0
        self.complete_calls = 0
        self.fail_calls = 0
        self.cancel_ack_calls = 0
        self.heartbeats = []
        self.resumed = []
        self.marked_running = []
        self.completed_tasks = []
        self.failed_tasks = []
        self.canceled_tasks = []
        self.uploaded_artifacts = []
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

    def resume(self, task_id: str):
        self.resume_calls += 1
        self._raise_if_needed("resume")
        self.resumed.append(task_id)
        task = self.resumed_tasks.get(task_id)
        if task is None:
            raise ValueError(f"task not found for resume: {task_id}")
        return task

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

    def upload_public_artifact(self, task_id: str, filename: str, local_path: str, content_type: str | None = None):
        payload_text = None
        if str(filename).lower().endswith(".json"):
            payload_text = Path(local_path).read_text(encoding="utf-8")
        self.uploaded_artifacts.append(
            {
                "task_id": task_id,
                "filename": filename,
                "local_path": local_path,
                "content_type": content_type,
                "payload_text": payload_text,
            }
        )
        return {
            "storage_backend": "local_public",
            "storage_key": f"/srv/mcp-ic-tool/public_artifacts/{filename}",
            "object_key": f"tenant/default/tasks/{task_id}/attempts/1/{filename}",
            "download_url": f"https://www.matterai.cn/dft/tenant/default/tasks/{task_id}/attempts/1/{filename}",
            "content_type": content_type,
            "size_bytes": 12,
        }


def _http_error(status_code: int) -> requests.exceptions.HTTPError:
    response = requests.Response()
    response.status_code = status_code
    return requests.exceptions.HTTPError(f"http {status_code}", response=response)


class FakeExecutionWorker:
    def __init__(self, result=None, exc: Exception | None = None, analysis_result=None):
        self.result = result or {"success": True, "work_directory": "/tmp/work"}
        self.exc = exc
        self.analysis_result = analysis_result or {"success": True, "work_directory": "/tmp/work", "final_energy": -10.0}

    async def run_structure_optimization(self, task_id, params, progress_callback=None):
        if progress_callback:
            await progress_callback(10, "queued")
            await progress_callback(40, "running", pid=12345)
        if self.exc is not None:
            raise self.exc
        return self.result

    async def _analyze_results(self, work_dir, run_result):
        result = dict(self.analysis_result)
        result.setdefault("success", True)
        result["work_directory"] = str(work_dir)
        result["computation_time"] = run_result.get("computation_time")
        result["process_id"] = run_result.get("process_id")
        return result

    async def run_dos_analysis(self, task_id, params, progress_callback=None):
        if progress_callback:
            await progress_callback(20, "analyzing")
        return {
            "success": True,
            "work_directory": self.result.get("work_directory", "/tmp/work"),
            "band_gap": 1.23,
            "summary": {"band_gap": 1.23},
        }


def _write_runtime_state(tmp_path: Path, worker_id: str, payload: dict) -> Path:
    runtime_dir = tmp_path / ".control-plane-runtime" / worker_id
    runtime_dir.mkdir(parents=True, exist_ok=True)
    state_path = runtime_dir / f"{payload['task_id']}.json"
    state_path.write_text(json.dumps(payload, ensure_ascii=False), encoding="utf-8")
    return state_path


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


def test_pull_worker_keeps_heartbeating_during_long_agent_analysis(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr(settings, "worker_background_heartbeat_interval_seconds", 0.01)

    class SlowAgentWorker(FakeExecutionWorker):
        async def run_agent_analysis(self, task_id, params, progress_callback=None):
            work_dir = tmp_path / task_id
            work_dir.mkdir(parents=True, exist_ok=True)
            if progress_callback:
                await progress_callback(35, "执行 Agent 分析...")
            await asyncio.sleep(0.05)
            return {"success": True, "work_directory": str(work_dir), "summary": {"steps": 3}}

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "agent-task",
                "task_type": "agent_analysis",
                "lease_token": "lease-1",
                "params": {"source_task_id": "source-task", "question": "请总结结果"},
            }
        ]
    )

    pull_worker = PullWorker(control_plane_client=client, execution_worker=SlowAgentWorker())
    pull_worker.run_once()

    assert client.completed_tasks == ["agent-task"]
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


def test_fail_path_reports_failure_context_with_scheduler_metadata(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-ctx"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("boom", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-ctx",
                "task_type": "structure_optimization",
                "lease_token": "lease-ctx",
                "params": {"formula": "Li2O"},
                "upstream_artifact_manifest": [],
            }
        ]
    )
    worker = FakeExecutionWorker(result={"success": False, "error": "boom", "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    failure_context = client.fail_payloads[0]["failure_context"]
    assert failure_context["scheduler_job_id"] == "12345"
    assert failure_context["work_directory"] == str(work_dir)


def test_scheduler_failure_message_includes_agent_guidance_for_unstable_scf(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-guidance"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text(
        "\n".join(
            [
                "BRMIX: very serious problems",
                "Error EDDDAV: Call to ZHEGV failed. Returncode = 79 1 80",
                "I REFUSE TO CONTINUE WITH THIS SICK JOB",
            ]
        ),
        encoding="utf-8",
    )

    pull_worker = PullWorker(control_plane_client=FakeControlPlaneClient(), execution_worker=FakeExecutionWorker())

    message = pull_worker._build_scheduler_failure_message("FAILED", work_dir)

    assert "先对该结构做一次 SCF" in message
    assert "CHGCAR" in message
    assert "band_structure" in message or "DOS" in message


def test_scheduler_failure_message_includes_agent_guidance_for_lattice_inconsistency(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-lattice"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text(
        "\n".join(
            [
                "HNFORM: k-point generating vectors and reciprocal lattice are incommensurate.",
                "Inconsistent Bravais lattice types found for crystalline and reciprocal lattice.",
                "I REFUSE TO CONTINUE WITH THIS SICK JOB",
            ]
        ),
        encoding="utf-8",
    )

    pull_worker = PullWorker(control_plane_client=FakeControlPlaneClient(), execution_worker=FakeExecutionWorker())

    message = pull_worker._build_scheduler_failure_message("FAILED", work_dir)

    assert "晶格参数" in message
    assert "重新导出" in message or "规范化" in message


def test_collect_artifact_manifest_recurses_and_preserves_relative_paths(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-1"
    nested_dir = work_dir / "BS_output"
    data_dir = work_dir / "data"
    nested_dir.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("outcar", encoding="utf-8")
    (nested_dir / "band_structure_report.html").write_text("<html/>", encoding="utf-8")
    (data_dir / "band_structure.csv").write_text("x,y", encoding="utf-8")

    pull_worker = PullWorker(control_plane_client=FakeControlPlaneClient(), execution_worker=FakeExecutionWorker())

    manifest = pull_worker._collect_artifact_manifest(str(work_dir))
    filenames = {item["filename"] for item in manifest}

    assert "OUTCAR" in filenames
    assert "BS_output/band_structure_report.html" in filenames
    assert "data/band_structure.csv" in filenames


def test_collect_artifact_manifest_includes_dos_analysis_directory(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-dos"
    dos_analysis_dir = work_dir / "dos_analysis"
    dos_analysis_dir.mkdir(parents=True, exist_ok=True)
    (dos_analysis_dir / "pymatgen_dos_analysis_report.html").write_text("<html/>", encoding="utf-8")
    (dos_analysis_dir / "total_dos.png").write_text("png", encoding="utf-8")
    (dos_analysis_dir / "total_dos.csv").write_text("x,y", encoding="utf-8")

    pull_worker = PullWorker(control_plane_client=FakeControlPlaneClient(), execution_worker=FakeExecutionWorker())

    manifest = pull_worker._collect_artifact_manifest(str(work_dir))
    filenames = {item["filename"] for item in manifest}

    assert "dos_analysis/pymatgen_dos_analysis_report.html" in filenames
    assert "dos_analysis/total_dos.png" in filenames
    assert "dos_analysis/total_dos.csv" in filenames


def test_pull_worker_supports_dos_analysis_tasks(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "analysis-task"
    output_dir = work_dir / "dos_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "pymatgen_dos_analysis_report.html").write_text("<html/>", encoding="utf-8")
    (output_dir / "total_dos.png").write_text("png", encoding="utf-8")
    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "analysis-task",
                "task_type": "dos_analysis",
                "lease_token": "lease-1",
                "params": {"source_task_id": "source-task"},
            }
        ]
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(work_dir)})

    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)
    pull_worker.run_once()

    assert client.completed_tasks == ["analysis-task"]
    filenames = {item["filename"] for item in client.complete_payloads[0]["artifact_manifest"]}
    assert "dos_analysis/pymatgen_dos_analysis_report.html" in filenames


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


def test_pull_worker_reconciles_completed_runtime_task_after_restart(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr(settings, "vasp_work_root", str(tmp_path))
    monkeypatch.setattr("time.sleep", lambda seconds: None)

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "result.log").write_text("done", encoding="utf-8")
    _write_runtime_state(
        tmp_path,
        "worker-a",
        {
            "task_id": "task-1",
            "task_type": "structure_optimization",
            "params": {"cif_url": "https://structures.example.com/Li2O.cif"},
            "lease_token": "lease-old",
            "work_directory": str(work_dir),
            "scheduler_job_id": "12345",
            "started_at": datetime.now(timezone.utc).isoformat(),
        },
    )

    client = FakeControlPlaneClient(
        tasks=[],
        operation_errors={"heartbeat": [_http_error(409)]},
        resumed_tasks={
            "task-1": {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-new",
                "status": "running",
                "queue_name": "default",
                "params": {"cif_url": "https://structures.example.com/Li2O.cif"},
            }
        },
    )
    worker = FakeExecutionWorker(analysis_result={"success": True, "final_energy": -20.0})
    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker, worker_id="worker-a")
    monkeypatch.setattr(pull_worker, "_query_scheduler_state", lambda *_args, **_kwargs: "COMPLETED")

    assert pull_worker.run_once() is True
    assert client.claim_calls == 0
    assert client.resume_calls == 1
    assert client.completed_tasks == ["task-1"]
    assert not list((tmp_path / ".control-plane-runtime" / "worker-a").glob("*.json"))


def test_pull_worker_reconciles_running_runtime_task_without_claiming_new_work(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr(settings, "vasp_work_root", str(tmp_path))
    monkeypatch.setattr("time.sleep", lambda seconds: None)

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    _write_runtime_state(
        tmp_path,
        "worker-a",
        {
            "task_id": "task-1",
            "task_type": "structure_optimization",
            "params": {"cif_url": "https://structures.example.com/Li2O.cif"},
            "lease_token": "lease-old",
            "work_directory": str(work_dir),
            "scheduler_job_id": "12345",
            "started_at": datetime.now(timezone.utc).isoformat(),
        },
    )

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-2",
                "task_type": "structure_optimization",
                "lease_token": "lease-2",
                "params": {},
            }
        ],
        operation_errors={"heartbeat": [_http_error(409)]},
        resumed_tasks={
            "task-1": {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-new",
                "status": "running",
                "queue_name": "default",
                "params": {"cif_url": "https://structures.example.com/Li2O.cif"},
            }
        },
    )
    worker = FakeExecutionWorker()
    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker, worker_id="worker-a")
    monkeypatch.setattr(pull_worker, "_query_scheduler_state", lambda *_args, **_kwargs: "RUNNING")

    assert pull_worker.run_once() is True
    assert client.resume_calls == 1
    assert client.claim_calls == 0
    assert client.completed_tasks == []
    assert list((tmp_path / ".control-plane-runtime" / "worker-a").glob("task-1.json"))


def test_pull_worker_drops_finalized_runtime_state_and_claims_new_task(monkeypatch, tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    monkeypatch.setattr(settings, "vasp_work_root", str(tmp_path))
    monkeypatch.setattr("time.sleep", lambda seconds: None)

    stale_work_dir = tmp_path / "task-stale"
    stale_work_dir.mkdir(parents=True, exist_ok=True)
    _write_runtime_state(
        tmp_path,
        "worker-a",
        {
            "task_id": "task-stale",
            "task_type": "agent_analysis",
            "params": {"source_task_id": "source-task", "question": "请总结结果"},
            "lease_token": "lease-stale",
            "work_directory": str(stale_work_dir),
            "scheduler_job_id": "12345",
            "started_at": datetime.now(timezone.utc).isoformat(),
        },
    )

    fresh_work_dir = tmp_path / "task-new"
    fresh_work_dir.mkdir(parents=True, exist_ok=True)
    (fresh_work_dir / "result.log").write_text("ok", encoding="utf-8")

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-new",
                "task_type": "structure_optimization",
                "lease_token": "lease-new",
                "params": {"input_url": "https://structures.example.com/Li2O.cif"},
            }
        ],
        operation_errors={"heartbeat": [_http_error(409)], "resume": [_http_error(409)]},
    )
    worker = FakeExecutionWorker(result={"success": True, "work_directory": str(fresh_work_dir)})
    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker, worker_id="worker-a")

    assert pull_worker.run_once() is True
    assert client.resume_calls == 1
    assert client.claim_calls == 1
    assert client.completed_tasks == ["task-new"]
    assert not list((tmp_path / ".control-plane-runtime" / "worker-a").glob("task-stale.json"))


def test_stage_public_artifacts_uploads_only_lightweight_outputs(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)
    html = work_dir / "optimization_analysis_report.html"
    png = work_dir / "force_convergence.png"
    outcar = work_dir / "OUTCAR"
    result_log = work_dir / "result.log"
    html.write_text("<html/>", encoding="utf-8")
    png.write_bytes(b"png-bytes")
    outcar.write_text("outcar", encoding="utf-8")
    result_log.write_text("log", encoding="utf-8")

    client = FakeControlPlaneClient()
    worker = FakeExecutionWorker()
    pull_worker = PullWorker(control_plane_client=client, execution_worker=worker)

    staged = pull_worker._stage_public_artifacts(
        "task-1",
        [
            {"filename": "optimization_analysis_report.html", "local_path": str(html), "content_type": "text/html"},
            {"filename": "force_convergence.png", "local_path": str(png), "content_type": "image/png"},
            {"filename": "OUTCAR", "local_path": str(outcar), "content_type": "text/plain"},
            {"filename": "result.log", "local_path": str(result_log), "content_type": "text/plain"},
        ],
    )

    uploaded_filenames = {item["filename"] for item in client.uploaded_artifacts}
    assert uploaded_filenames == {
        "optimization_analysis_report.html",
        "force_convergence.png",
        "OUTCAR",
        "result.log",
    }
    staged_by_name = {item["filename"]: item for item in staged}
    assert staged_by_name["optimization_analysis_report.html"]["storage_backend"] == "local_public"
    assert staged_by_name["force_convergence.png"]["storage_backend"] == "local_public"
    assert staged_by_name["OUTCAR"]["storage_backend"] == "local_public"
    assert staged_by_name["result.log"]["storage_backend"] == "local_public"


def test_stage_public_artifacts_rewrites_public_json_without_local_paths(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker

    work_dir = tmp_path / "task-agent"
    analysis_dir = work_dir / "agent_analysis"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    html = analysis_dir / "agent_analysis_report.html"
    summary = work_dir / "analysis_summary.json"
    html.write_text("<html/>", encoding="utf-8")
    summary.write_text(
        json.dumps(
            {
                "success": True,
                "work_directory": str(work_dir),
                "html_report_url": str(html),
                "agent_analysis_report_html_path": str(html),
            },
            ensure_ascii=False,
        ),
        encoding="utf-8",
    )

    client = FakeControlPlaneClient()
    pull_worker = PullWorker(control_plane_client=client, execution_worker=FakeExecutionWorker())

    pull_worker._stage_public_artifacts(
        "task-agent",
        [
            {"filename": "analysis_summary.json", "local_path": str(summary), "content_type": "application/json"},
            {"filename": "agent_analysis/agent_analysis_report.html", "local_path": str(html), "content_type": "text/html"},
        ],
    )

    summary_uploads = [item for item in client.uploaded_artifacts if item["filename"] == "analysis_summary.json"]
    assert len(summary_uploads) == 2
    rewritten_payload = json.loads(summary_uploads[-1]["payload_text"])
    assert rewritten_payload["html_report_url"].startswith("https://www.matterai.cn/dft/")
    assert "agent_analysis_report_html_path" not in rewritten_payload
    assert "work_directory" not in rewritten_payload
    assert "/task-agent/agent_analysis/agent_analysis_report.html" not in summary_uploads[-1]["payload_text"]


def test_pull_worker_uses_default_work_directory_when_result_omits_work_directory(tmp_path):
    from src.vasp_server.hpc_pull_worker import PullWorker
    from src.vasp_server.settings import settings

    settings.vasp_work_root = str(tmp_path)
    task_dir = tmp_path / "hpc-worker" / "task-1"
    task_dir.mkdir(parents=True, exist_ok=True)
    (task_dir / "optimization_analysis_report.html").write_text("<html/>", encoding="utf-8")
    (task_dir / "CONTCAR").write_text("contcar", encoding="utf-8")

    client = FakeControlPlaneClient(
        tasks=[
            {
                "task_id": "task-1",
                "task_type": "structure_optimization",
                "lease_token": "lease-1",
                "params": {},
            }
        ]
    )
    worker = FakeExecutionWorker(result={"success": True})

    class WorkerWithBaseDir(FakeExecutionWorker):
        def __init__(self):
            super().__init__(result={"success": True})
            self.base_work_dir = tmp_path / "hpc-worker"

        async def run_structure_optimization(self, task_id, params, progress_callback=None):
            if progress_callback:
                await progress_callback(40, "running", pid=12345)
            return {"success": True}

    pull_worker = PullWorker(control_plane_client=client, execution_worker=WorkerWithBaseDir())

    assert pull_worker.run_once() is True
    uploaded_filenames = {item["filename"] for item in client.uploaded_artifacts}
    assert "optimization_analysis_report.html" in uploaded_filenames
    assert "CONTCAR" in uploaded_filenames
