from datetime import datetime, timedelta, timezone
import uuid

from src.vasp_server.task_manager.models import Task


def _create_task(db_session, **overrides):
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="test_user",
        task_type="structure_optimization",
        status="leased",
        analysis_status="pending",
        tenant_id="default",
        queue_name="default",
        params={"formula": "Li2O"},
        input_snapshot={"formula": "Li2O"},
        worker_id="hpc-a",
        lease_token="lease-1",
        retry_count=0,
        max_retries=1,
        heartbeat_at=datetime.now(timezone.utc),
        lease_expires_at=datetime.now(timezone.utc) + timedelta(seconds=30),
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_cancel_requested_task_becomes_canceled_after_worker_ack(task_manager, db_session):
    task_id = _create_task(db_session)

    assert task_manager.request_cancel(task_id, "test_user") is True
    task_manager.ack_cancel(task_id, lease_token="lease-1", worker_id="hpc-a")

    task = db_session.get(Task, task_id)
    assert task is not None
    assert task.status == "canceled"
    assert task.finalized_at is not None


def test_cancel_requested_task_is_not_marked_failed(task_manager, db_session):
    task_id = _create_task(db_session, status="cancel_requested", cancel_requested_at=datetime.now(timezone.utc))

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="timeout talking to scheduler",
    )

    task = db_session.get(Task, task_id)
    assert task is not None
    assert task.status == "canceled"


def test_expired_heartbeat_returns_running_task_to_queue(task_manager, db_session):
    task_id = _create_task(
        db_session,
        status="running",
        heartbeat_at=datetime.now(timezone.utc) - timedelta(minutes=10),
        lease_expires_at=datetime.now(timezone.utc) - timedelta(minutes=10),
    )

    recovered = task_manager.mark_orphaned_running_tasks()

    task = db_session.get(Task, task_id)
    assert recovered == 1
    assert task is not None
    assert task.status == "queued"
    assert task.worker_id is None
    assert task.lease_token is None


def test_retryable_infrastructure_failures_increment_retry_count_and_requeue(task_manager, db_session):
    task_id = _create_task(db_session)

    task_manager.fail_execution(
        task_id=task_id,
        lease_token="lease-1",
        worker_id="hpc-a",
        error_message="scheduler timeout while staging job",
    )

    task = db_session.get(Task, task_id)
    assert task is not None
    assert task.status == "queued"
    assert task.retry_count == 1
