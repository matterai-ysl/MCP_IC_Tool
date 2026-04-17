from datetime import datetime, timedelta, timezone
import uuid

from src.vasp_server.task_manager.models import Task, AnalysisRun, ExecutionAttempt


def _create_task(db_session, **overrides):
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="test_user",
        task_type="structure_optimization",
        status="queued",
        analysis_status="pending",
        params={"formula": "Li2O"},
        input_snapshot={"formula": "Li2O"},
        queue_name="default",
        tenant_id="default",
        priority=0,
        retry_count=0,
        max_retries=1,
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_claim_task_assigns_worker_and_lease_token(task_manager, db_session):
    task_id = _create_task(db_session)

    claimed = task_manager.claim_next_task(worker_id="hpc-a", queue_name="default")

    assert claimed is not None
    assert claimed.id == task_id
    assert claimed.worker_id == "hpc-a"
    assert claimed.lease_token
    assert claimed.status == "leased"

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "leased"
    assert stored_task.worker_id == "hpc-a"
    assert stored_task.lease_token == claimed.lease_token
    assert stored_task.lease_expires_at is not None


def test_only_one_worker_can_claim_same_queued_task(task_manager, db_session):
    _create_task(db_session)

    first_claim = task_manager.claim_next_task(worker_id="hpc-a", queue_name="default")
    second_claim = task_manager.claim_next_task(worker_id="hpc-b", queue_name="default")

    assert first_claim is not None
    assert second_claim is None


def test_requeue_expired_leases_returns_task_to_queue(task_manager, db_session):
    task_id = _create_task(
        db_session,
        status="leased",
        worker_id="hpc-a",
        lease_token="expired-token",
        lease_expires_at=datetime.now(timezone.utc) - timedelta(seconds=30),
        heartbeat_at=datetime.now(timezone.utc) - timedelta(seconds=45),
    )

    requeued = task_manager.requeue_expired_leases()

    assert requeued == 1
    refreshed_task = db_session.get(Task, task_id)
    assert refreshed_task is not None
    assert refreshed_task.status == "queued"
    assert refreshed_task.worker_id is None
    assert refreshed_task.lease_token is None


def test_claim_task_automatically_requeues_expired_leases_before_claiming(task_manager, db_session):
    task_id = _create_task(
        db_session,
        status="leased",
        worker_id="hpc-a",
        lease_token="expired-token",
        lease_expires_at=datetime.now(timezone.utc) - timedelta(seconds=30),
        heartbeat_at=datetime.now(timezone.utc) - timedelta(seconds=45),
    )

    claimed = task_manager.claim_next_task(worker_id="hpc-b", queue_name="default")

    assert claimed is not None
    assert claimed.id == task_id
    assert claimed.worker_id == "hpc-b"
    assert claimed.status == "leased"

    refreshed_task = db_session.get(Task, task_id)
    assert refreshed_task is not None
    assert refreshed_task.status == "leased"
    assert refreshed_task.worker_id == "hpc-b"
    assert refreshed_task.lease_token == claimed.lease_token
    assert refreshed_task.lease_expires_at is not None


def test_cancel_requested_task_is_never_newly_claimed(task_manager, db_session):
    task_id = _create_task(db_session, status="cancel_requested")

    claimed = task_manager.claim_next_task(worker_id="hpc-a", queue_name="default")

    assert claimed is None
    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "cancel_requested"


def test_request_cancel_marks_task_and_timestamp(task_manager, db_session):
    task_id = _create_task(db_session, status="leased", worker_id="hpc-a", lease_token="lease-1")

    success = task_manager.request_cancel(task_id=task_id, user_id="test_user")

    assert success is True
    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "cancel_requested"
    assert stored_task.cancel_requested_at is not None


def test_request_cancel_immediately_cancels_unclaimed_queued_task(task_manager, db_session):
    task_id = _create_task(db_session, status="queued")

    success = task_manager.request_cancel(task_id=task_id, user_id="test_user")

    assert success is True
    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "canceled"
    assert stored_task.finalized_at is not None
    assert stored_task.worker_id is None
    assert stored_task.lease_token is None


def test_submit_task_quota_ignores_unleased_cancel_requested_tasks(task_manager, db_session):
    _create_task(
        db_session,
        status="running",
        worker_id="hpc-a",
        lease_token="lease-1",
        heartbeat_at=datetime.now(timezone.utc),
    )
    _create_task(db_session, status="cancel_requested", cancel_requested_at=datetime.now(timezone.utc))
    _create_task(db_session, status="cancel_requested", cancel_requested_at=datetime.now(timezone.utc))

    task_id = task_manager.submit_task(
        user_id="test_user",
        task_type="structure_optimization",
        params={"cif_url": "https://structures.example.com/Li2O-new.cif", "queue_name": "default"},
    )

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "queued"


def test_submit_task_dedup_is_scoped_by_queue_name(task_manager, db_session):
    first_task_id = task_manager.submit_task(
        user_id="test_user",
        task_type="structure_optimization",
        params={"cif_url": "https://structures.example.com/Li2O.cif", "queue_name": "hpc-a"},
    )
    first_task = db_session.get(Task, first_task_id)
    assert first_task is not None
    first_task.status = "completed"
    db_session.add(first_task)
    db_session.commit()

    second_task_id = task_manager.submit_task(
        user_id="test_user",
        task_type="structure_optimization",
        params={"cif_url": "https://structures.example.com/Li2O.cif", "queue_name": "hpc-b"},
    )

    assert second_task_id != first_task_id
    second_task = db_session.get(Task, second_task_id)
    assert second_task is not None
    assert second_task.queue_name == "hpc-b"


def test_complete_execution_triggers_analysis_summary_for_distributed_worker(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        status="leased",
        worker_id="hpc-a",
        lease_token="lease-1",
    )
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: None,
    )

    task_manager.complete_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        result_data={
            "success": True,
            "final_energy": -11.95,
            "computation_time": 12.3,
            "analysis_report_html_path": "https://www.matterai.cn/dft/tenant/default/tasks/task-x/attempts/1/optimization_analysis_report.html",
        },
    )

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "completed"
    assert stored_task.analysis_status == "completed"
    assert stored_task.result_summary is not None
    assert stored_task.result_summary["success"] is True
    assert stored_task.result_summary["final_energy"] == -11.95
    assert stored_task.result_summary["html_report_url"].startswith("https://www.matterai.cn/dft/")

    runs = db_session.query(AnalysisRun).filter(AnalysisRun.task_id == task_id).all()
    assert len(runs) == 1
    assert runs[0].status == "completed"
    assert runs[0].summary is not None


def test_complete_execution_creates_execution_attempt_for_distributed_worker(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        status="leased",
        worker_id="hpc-a",
        lease_token="lease-1",
    )
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: None,
    )

    task_manager.complete_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        result_data={
            "success": True,
            "final_energy": -10.5,
            "slurm_job_id": "12345",
            "work_directory": "/tmp/distributed-success",
        },
    )

    attempts = db_session.query(ExecutionAttempt).filter(ExecutionAttempt.task_id == task_id).all()
    assert len(attempts) == 1
    assert attempts[0].status == "succeeded"
    assert attempts[0].scheduler_job_id == "12345"
    assert attempts[0].work_directory == "/tmp/distributed-success"
    assert attempts[0].worker_id == "hpc-a"
    assert attempts[0].lease_token == "lease-1"

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.current_execution_attempt_id == attempts[0].id


def test_fail_execution_creates_failed_attempt_for_distributed_worker(task_manager, db_session, monkeypatch):
    task_id = _create_task(
        db_session,
        status="leased",
        worker_id="hpc-a",
        lease_token="lease-1",
    )
    monkeypatch.setattr(
        "src.vasp_server.task_manager.manager.send_notification_async",
        lambda payload: None,
    )

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="fatal electronic divergence",
        failure_context={
            "scheduler_job_id": "54321",
            "work_directory": "/tmp/distributed-fail",
        },
    )

    attempts = db_session.query(ExecutionAttempt).filter(ExecutionAttempt.task_id == task_id).all()
    assert len(attempts) == 1
    assert attempts[0].status == "runtime_failed"
    assert attempts[0].scheduler_job_id == "54321"
    assert attempts[0].work_directory == "/tmp/distributed-fail"
    assert attempts[0].failure_detail == "fatal electronic divergence"
    assert attempts[0].worker_id == "hpc-a"
    assert attempts[0].lease_token == "lease-1"

    stored_task = db_session.get(Task, task_id)
    assert stored_task is not None
    assert stored_task.status == "failed"
    assert stored_task.current_execution_attempt_id == attempts[0].id
