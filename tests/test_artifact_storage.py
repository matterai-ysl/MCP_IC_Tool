import uuid

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
        result_summary={"html_report_url": "/tmp/raw/report.html"},
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_register_artifact_uses_object_storage_key(task_manager, db_session):
    task_id = _create_task(db_session)

    artifact = task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="report.html",
        content_type="text/html",
        size_bytes=128.0,
        attempt_no=1,
    )

    assert artifact.storage_backend == "oss"
    assert artifact.storage_key.startswith("tenant/")
    assert artifact.object_key.startswith("tenant/")


def test_object_storage_service_generates_presigned_urls(task_manager):
    upload_url = task_manager.storage_service.create_upload_url(
        object_key="tenant/tenant-a/tasks/task-1/attempts/1/report.html",
        content_type="text/html",
    )
    download_url = task_manager.storage_service.create_download_url(
        object_key="tenant/tenant-a/tasks/task-1/attempts/1/report.html",
    )

    assert upload_url.startswith("https://")
    assert "signature=" in upload_url
    assert "content_type=text%2Fhtml" in upload_url
    assert download_url.startswith("https://")
    assert "signature=" in download_url


def test_register_artifact_persists_object_storage_metadata(task_manager, db_session):
    task_id = _create_task(db_session)

    artifact = task_manager.register_artifact(
        task_id=task_id,
        artifact_type="outcar",
        owner_type="execution",
        owner_id="attempt-1",
        filename="OUTCAR",
        content_type="text/plain",
        size_bytes=1024.0,
        etag="etag-123",
        sha256="sha-256",
        attempt_no=2,
    )

    stored_artifact = db_session.get(Artifact, artifact.id)
    assert stored_artifact is not None
    assert stored_artifact.bucket == task_manager.storage_service.bucket
    assert stored_artifact.object_key.endswith("/OUTCAR")
    assert stored_artifact.etag == "etag-123"
    assert stored_artifact.sha256 == "sha-256"
    assert stored_artifact.content_type == "text/plain"


def test_task_status_prefers_signed_artifact_urls_not_local_paths(task_manager, db_session):
    task_id = _create_task(db_session)
    task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="report.html",
        content_type="text/html",
        size_bytes=128.0,
        attempt_no=1,
    )

    task = task_manager.get_task(task_id, "test_user")
    response = _build_task_response(task, task_id)

    assert response.html_report_url is not None
    assert response.html_report_url.startswith("https://")
    assert "/tmp/raw/report.html" not in response.html_report_url
    assert response.artifacts is not None
    assert response.artifacts[0].download_url is not None


def test_local_public_artifact_uses_dft_public_url(task_manager, db_session, monkeypatch):
    from src.vasp_server.settings import settings

    task_id = _create_task(db_session)
    monkeypatch.setattr(settings, "public_artifact_base_url", "https://www.matterai.cn/dft")

    artifact = Artifact(
        id=uuid.uuid4().hex,
        task_id=task_id,
        owner_type="analysis",
        owner_id="run-1",
        artifact_type="html_report",
        storage_backend="local_public",
        storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/tenant-a/tasks/task-1/attempts/1/report.html",
        object_key="tenant/tenant-a/tasks/task-1/attempts/1/report.html",
        mime_type="text/html",
        content_type="text/html",
    )
    db_session.add(artifact)
    db_session.commit()

    task = task_manager.get_task(task_id, "test_user")
    response = _build_task_response(task, task_id)

    assert response.html_report_url == "https://www.matterai.cn/dft/tenant/tenant-a/tasks/task-1/attempts/1/report.html"


def test_local_only_artifact_has_no_public_download_url(task_manager, db_session):
    task_id = _create_task(db_session)

    artifact = Artifact(
        id=uuid.uuid4().hex,
        task_id=task_id,
        owner_type="execution",
        owner_id="attempt-1",
        artifact_type="outcar",
        storage_backend="local",
        storage_key="/data/home/ysl9527/vasp_calculations/task-1/OUTCAR",
        object_key=None,
        mime_type="text/plain",
        content_type="text/plain",
    )
    db_session.add(artifact)
    db_session.commit()

    task = task_manager.get_task(task_id, "test_user")
    response = _build_task_response(task, task_id)

    assert response.artifacts is not None
    assert response.artifacts[0].download_url is None
