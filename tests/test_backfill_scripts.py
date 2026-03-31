from datetime import datetime, timedelta, timezone
import uuid

from scripts.backfill_artifacts import backfill_local_artifacts
from scripts.requeue_stuck_tasks import requeue_stuck_tasks
from src.vasp_server.storage import ObjectStorageService
from src.vasp_server.task_manager.models import Artifact, Task
from src.vasp_server.vasp_server_api import _build_task_response


def _create_task(db_session, **overrides):
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="test_user",
        task_type="structure_optimization",
        status="completed",
        analysis_status="completed",
        tenant_id="tenant-a",
        queue_name="default",
        params={"formula": "Li2O"},
        input_snapshot={"formula": "Li2O"},
        result_summary={"html_report_url": "/static/report.html"},
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_backfill_artifacts_converts_local_artifact_rows_to_object_storage_rows(db_session):
    task_id = _create_task(db_session)
    artifact = Artifact(
        id=uuid.uuid4().hex,
        task_id=task_id,
        owner_type="analysis",
        artifact_type="html_report",
        storage_backend="local",
        storage_key="/tmp/report.html",
        mime_type="text/html",
    )
    db_session.add(artifact)
    db_session.commit()

    migrated = backfill_local_artifacts(db_session, ObjectStorageService.from_settings())

    stored = db_session.get(Artifact, artifact.id)
    assert migrated == 1
    assert stored is not None
    assert stored.storage_backend == "oss"
    assert stored.object_key.startswith("tenant/")


def test_requeue_script_only_touches_expired_leases(task_manager, db_session):
    expired_task_id = _create_task(
        db_session,
        status="leased",
        analysis_status="pending",
        worker_id="hpc-a",
        lease_token="expired",
        heartbeat_at=datetime.now(timezone.utc) - timedelta(minutes=10),
        lease_expires_at=datetime.now(timezone.utc) - timedelta(minutes=10),
    )
    valid_task_id = _create_task(
        db_session,
        status="leased",
        analysis_status="pending",
        worker_id="hpc-b",
        lease_token="valid",
        heartbeat_at=datetime.now(timezone.utc),
        lease_expires_at=datetime.now(timezone.utc) + timedelta(minutes=10),
    )

    result = requeue_stuck_tasks(task_manager)

    expired_task = db_session.get(Task, expired_task_id)
    valid_task = db_session.get(Task, valid_task_id)
    assert result["expired_leases"] == 1
    assert expired_task is not None and expired_task.status == "queued"
    assert valid_task is not None and valid_task.status == "leased"


def test_cutover_flag_disables_old_local_artifact_urls(monkeypatch, task_manager, db_session):
    from src.vasp_server import vasp_server_api

    task_id = _create_task(db_session)
    task = task_manager.get_task(task_id, "test_user")

    monkeypatch.setattr(vasp_server_api.settings, "legacy_local_artifact_urls_enabled", False)
    response = _build_task_response(task, task_id)

    assert response.html_report_url is None
