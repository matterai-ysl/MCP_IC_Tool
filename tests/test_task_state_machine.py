from src.vasp_server.schemas import TASK_STATES, TaskStatus
from src.vasp_server.task_manager.models import (
    TASK_STATE_TRANSITIONS,
    ExecutionAttempt,
    Task,
    can_transition_task_status,
    resolve_runtime_failure_status,
)


def test_task_status_enum_matches_pull_worker_contract():
    assert TASK_STATES == {
        "queued",
        "leased",
        "running",
        "uploading",
        "analyzing",
        "completed",
        "failed",
        "cancel_requested",
        "canceled",
    }
    assert {status.value for status in TaskStatus} == TASK_STATES


def test_task_lifecycle_for_pull_workers():
    assert can_transition_task_status("queued", "leased")
    assert can_transition_task_status("leased", "running")
    assert can_transition_task_status("running", "uploading")
    assert can_transition_task_status("uploading", "analyzing")
    assert can_transition_task_status("analyzing", "completed")


def test_cancel_requested_transitions_to_canceled():
    assert can_transition_task_status("queued", "cancel_requested")
    assert can_transition_task_status("leased", "cancel_requested")
    assert can_transition_task_status("running", "cancel_requested")
    assert can_transition_task_status("cancel_requested", "canceled")


def test_expired_lease_returns_task_to_queue():
    assert can_transition_task_status("leased", "queued")
    assert "queued" in TASK_STATE_TRANSITIONS["leased"]


def test_runtime_failure_requeues_only_when_retry_budget_remains():
    assert resolve_runtime_failure_status(retry_count=0, max_retries=1) == "queued"
    assert resolve_runtime_failure_status(retry_count=1, max_retries=1) == "failed"


def test_task_model_supports_pull_worker_fields():
    task_columns = Task.__table__.columns.keys()

    for column_name in [
        "tenant_id",
        "queue_name",
        "priority",
        "worker_id",
        "lease_token",
        "lease_expires_at",
        "heartbeat_at",
        "cancel_requested_at",
        "finalized_at",
        "retry_count",
        "max_retries",
    ]:
        assert column_name in task_columns


def test_execution_attempt_model_supports_pull_worker_metadata():
    attempt_columns = ExecutionAttempt.__table__.columns.keys()

    for column_name in [
        "worker_id",
        "lease_token",
        "scheduler_job_id",
        "artifact_manifest",
    ]:
        assert column_name in attempt_columns
