import logging
import uuid

import requests

from src.vasp_server.task_manager.models import Artifact, Task
from tests.conftest import make_mock_worker_result, run_worker_with_mock_result


def _create_task(db_session, **overrides):
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="user-10001",
        task_type="structure_optimization",
        status="leased",
        analysis_status="pending",
        params={"cif_url": "https://structures.example.com/Li2O.cif", "queue_name": "default"},
        input_snapshot={"cif_url": "https://structures.example.com/Li2O.cif"},
        queue_name="default",
        tenant_id="default",
        priority=0,
        retry_count=0,
        max_retries=1,
        worker_id="hpc-a",
        lease_token="lease-1",
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_complete_execution_sends_success_notification(task_manager, db_session, monkeypatch):
    task_id = _create_task(db_session)
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.complete_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        result_data={"computation_time": 320.5, "final_energy": -100.5},
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["user_id"] == "user-10001"
    assert payload["notification_type"] == "tool_complete"
    assert payload["channel"] == "email"
    assert payload["language"] == "zh"
    assert payload["context"]["tool_name"] == "preset-vasp-calculation"
    assert payload["context"]["tool_display_name"] == "结构优化计算"
    assert payload["context"]["status"] == "success"
    assert payload["context"]["execution_time_seconds"] == 320.5
    assert "完成" in payload["context"]["result_summary"]


def test_complete_execution_notification_includes_html_report_url(task_manager, db_session, monkeypatch):
    task_id = _create_task(db_session)
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.complete_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        result_data={
            "computation_time": 12.5,
            "final_energy": -10.2,
            "analysis_report_html_path": "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/1/optimization_analysis_report.html",
        },
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["context"]["html_report_url"] == (
        "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/1/optimization_analysis_report.html"
    )
    assert "https://www.matterai.cn/dft/" in payload["context"]["result_summary"]


def test_complete_execution_notification_prefers_public_artifact_html_report_url(task_manager, db_session, monkeypatch):
    task_id = _create_task(db_session)
    artifact = Artifact(
        id=uuid.uuid4().hex,
        task_id=task_id,
        owner_type="analysis",
        owner_id=uuid.uuid4().hex,
        artifact_type="html_report",
        storage_backend="local_public",
        storage_key="tenant/default/tasks/task-x/attempts/2/agent_analysis/agent_analysis_report.html",
        object_key="tenant/default/tasks/task-x/attempts/2/agent_analysis/agent_analysis_report.html",
        mime_type="text/html",
        content_type="text/html",
    )
    db_session.add(artifact)
    db_session.commit()

    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.complete_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        result_data={
            "computation_time": 12.5,
            "final_energy": -10.2,
            "analysis_report_html_path": "/data/home/ysl9527/vasp_calculations/hpc-c2ln6/task-x/agent_analysis/agent_analysis_report.html",
        },
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["context"]["html_report_url"] == (
        "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/2/agent_analysis/agent_analysis_report.html"
    )
    assert "/data/home/ysl9527" not in payload["context"]["result_summary"]


def test_fail_execution_requeue_does_not_send_notification(task_manager, db_session, monkeypatch):
    task_id = _create_task(db_session, max_retries=2)
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="Connection timeout after 30s",
    )

    assert sent_payloads == []
    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "queued"
    assert stored_task.retry_count == 1


def test_fail_execution_terminal_failure_sends_error_notification(task_manager, db_session, monkeypatch):
    task_id = _create_task(db_session, max_retries=1)
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="Fatal VASP parsing error",
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["notification_type"] == "tool_error"
    assert payload["context"]["status"] == "error"
    assert payload["context"]["error_message"] == "Fatal VASP parsing error"

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "failed"


def test_fail_execution_notification_includes_html_report_url_when_available(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        max_retries=1,
        result_summary={
            "html_report_url": "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/1/report.html",
        },
    )
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="Fatal VASP parsing error",
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["context"]["html_report_url"] == (
        "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/1/report.html"
    )


def test_fail_execution_notification_prefers_public_artifact_html_report_url(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        max_retries=1,
        result_summary={
            "html_report_url": "/data/home/ysl9527/vasp_calculations/hpc-c2ln6/task-x/report.html",
        },
    )
    artifact = Artifact(
        id=uuid.uuid4().hex,
        task_id=task_id,
        owner_type="analysis",
        owner_id=uuid.uuid4().hex,
        artifact_type="html_report",
        storage_backend="local_public",
        storage_key="tenant/default/tasks/task-x/attempts/3/report.html",
        object_key="tenant/default/tasks/task-x/attempts/3/report.html",
        mime_type="text/html",
        content_type="text/html",
    )
    db_session.add(artifact)
    db_session.commit()

    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="Fatal VASP parsing error",
    )

    assert len(sent_payloads) == 1
    payload = sent_payloads[0]
    assert payload["context"]["html_report_url"] == (
        "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/3/report.html"
    )


def test_local_worker_success_sends_notification(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        status="queued",
        worker_id=None,
        lease_token=None,
    )
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    run_worker_with_mock_result(
        task_manager,
        task_id,
        make_mock_worker_result("structure_optimization", success=True),
        "structure_optimization",
    )

    assert len(sent_payloads) == 1
    assert sent_payloads[0]["notification_type"] == "tool_complete"
    assert sent_payloads[0]["context"]["status"] == "success"


def test_local_worker_terminal_failure_sends_notification(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        status="queued",
        worker_id=None,
        lease_token=None,
        max_retries=1,
    )
    sent_payloads = []
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: sent_payloads.append(payload),
    )

    run_worker_with_mock_result(
        task_manager,
        task_id,
        make_mock_worker_result("structure_optimization", success=False),
        "structure_optimization",
    )

    assert len(sent_payloads) == 1
    assert sent_payloads[0]["notification_type"] == "tool_error"
    assert sent_payloads[0]["context"]["status"] == "error"


def test_send_notification_logs_and_swallows_failures(monkeypatch, caplog):
    monkeypatch.setenv("NOTIFICATION_SERVICE_BASE_URL", "https://notify.example.com")
    monkeypatch.setenv("NOTIFICATION_SERVICE_API_KEY", "notify-key")

    def raise_request_error(*args, **kwargs):
        raise requests.RequestException("boom")

    monkeypatch.setattr("src.vasp_server.notification_client.requests.post", raise_request_error)

    from src.vasp_server.notification_client import send_notification

    with caplog.at_level(logging.WARNING):
        send_notification(
            {
                "user_id": "user-10001",
                "notification_type": "tool_complete",
                "channel": "email",
                "language": "zh",
                "context": {"tool_name": "preset-vasp-calculation"},
            }
        )

    assert "Notification delivery failed" in caplog.text


def test_send_notification_async_uses_background_thread(monkeypatch):
    started = []

    class FakeThread:
        def __init__(self, target=None, args=(), daemon=None):
            started.append(
                {
                    "target": target,
                    "args": args,
                    "daemon": daemon,
                }
            )

        def start(self):
            started.append("started")

    monkeypatch.setattr("src.vasp_server.notification_client.threading.Thread", FakeThread)

    from src.vasp_server.notification_client import send_notification_async

    send_notification_async(
        {
            "user_id": "user-10001",
            "notification_type": "tool_complete",
            "channel": "email",
            "language": "zh",
            "context": {"tool_name": "preset-vasp-calculation"},
        }
    )

    assert started[0]["daemon"] is True
    assert started[1] == "started"
