from datetime import datetime, timedelta, timezone
import uuid

from src.vasp_server.task_manager.models import Task


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
