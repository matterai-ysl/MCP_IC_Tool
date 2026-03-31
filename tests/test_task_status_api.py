import asyncio
import json
import uuid

import httpx
from starlette.routing import Mount

from src.vasp_server.task_manager.models import Task


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
        result_summary={"html_report_url": "/static/reports/report.html"},
        result_data={
            "analysis_report_html_path": "/download/file/tenant-a/report.html",
            "work_directory": "/tmp/should-not-leak",
        },
        result_path="/tmp/should-not-leak",
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def _request(app, method: str, url: str, **kwargs):
    async def _send():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(transport=transport, base_url="http://testserver") as client:
            return await client.request(method, url, **kwargs)

    return asyncio.run(_send())


def test_task_status_returns_signed_artifact_urls_not_local_paths(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="report.html",
        content_type="text/html",
        size_bytes=512.0,
        attempt_no=1,
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].startswith("https://")
    assert payload["result_summary"]["html_report_url"].startswith("https://")
    assert payload["result_data"]["analysis_report_html_path"].startswith("https://")
    serialized = json.dumps(payload)
    assert "/static/" not in serialized
    assert "/download/file/" not in serialized
    assert "/tmp/should-not-leak" not in serialized


def test_static_mount_is_not_registered():
    from src.vasp_server.vasp_server_api import app

    assert not any(isinstance(route, Mount) and route.path == "/static" for route in app.routes)


def test_download_file_endpoint_is_deprecated():
    from src.vasp_server.vasp_server_api import app

    response = _request(app, "GET", "/download/file/tenant-a/report.html")

    assert response.status_code == 410


def test_artifact_list_returns_typed_metadata_and_signed_urls_only(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="report.html",
        content_type="text/html",
        size_bytes=512.0,
        attempt_no=1,
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    artifact = response.json()["artifacts"][0]
    assert artifact["artifact_type"] == "html_report"
    assert artifact["content_type"] == "text/html"
    assert artifact["size_bytes"] == 512.0
    assert artifact["download_url"].startswith("https://")
    assert "storage_key" not in artifact
