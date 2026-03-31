import asyncio
import uuid

import httpx

from src.mcp_ic_tool.client import VaspAPIClient
from src.vasp_server.task_manager.models import ExecutionAttempt, Task


def _request(app, method: str, url: str, **kwargs):
    async def _send():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(transport=transport, base_url="http://testserver") as client:
            return await client.request(method, url, **kwargs)

    return asyncio.run(_send())


class LocalVaspAPIClient(VaspAPIClient):
    def __init__(self, app):
        super().__init__(base_url="http://testserver")
        self.app = app

    async def _apost(self, path: str, json: dict):
        transport = httpx.ASGITransport(app=self.app)
        async with httpx.AsyncClient(transport=transport, base_url=self.base_url) as client:
            response = await client.post(path, json=json)
            response.raise_for_status()
            return response.json()

    async def _aget(self, path: str, params: dict):
        transport = httpx.ASGITransport(app=self.app)
        async with httpx.AsyncClient(transport=transport, base_url=self.base_url) as client:
            response = await client.get(path, params=params)
            response.raise_for_status()
            return response.json()


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
        params={"cif_url": "https://structures.example.com/Li2O.cif"},
        input_snapshot={"cif_url": "https://structures.example.com/Li2O.cif"},
        result_summary={},
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


def test_submit_task_only_creates_db_record_and_returns_task_id(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={"user_id": "test_user", "cif_url": "https://structures.example.com/Li2O.cif"},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["status"] == "queued"

    task = db_session.get(Task, payload["task_id"])
    assert task is not None
    assert task.status == "queued"
    attempts = db_session.query(ExecutionAttempt).filter_by(task_id=payload["task_id"]).all()
    assert attempts == []


def test_submit_task_rejects_formula_input(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={"user_id": "test_user", "formula": "Li2O"},
    )

    assert response.status_code == 422


def test_status_endpoint_hydrates_artifact_urls_from_db_metadata(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="report.html",
        content_type="text/html",
        size_bytes=256.0,
        attempt_no=1,
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].startswith("https://")
    assert payload["artifacts"][0]["download_url"].startswith("https://")


def test_submit_with_upstream_task_id_injects_artifact_manifest(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    upstream_task_id = _create_task(
        db_session,
        task_type="structure_optimization",
        status="completed",
        analysis_status="completed",
    )
    task_manager.register_artifact(
        task_id=upstream_task_id,
        artifact_type="contcar",
        owner_type="execution",
        owner_id="attempt-1",
        filename="CONTCAR",
        content_type="text/plain",
        size_bytes=128.0,
        attempt_no=1,
    )

    response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={"user_id": "test_user", "optimized_task_id": upstream_task_id},
    )

    assert response.status_code == 200
    task_id = response.json()["task_id"]
    task = db_session.get(Task, task_id)
    assert task is not None
    manifest = (task.params or {}).get("upstream_artifact_manifest", [])
    assert manifest
    assert manifest[0]["artifact_type"] == "contcar"
    assert manifest[0]["download_url"].startswith("https://")


def test_mcp_client_still_works_with_public_api_shape():
    from src.vasp_server.vasp_server_api import app

    client = LocalVaspAPIClient(app)

    submit_payload = asyncio.run(
        client.submit_structure_optimization(
            {"user_id": "test_user", "cif_url": "https://structures.example.com/Li2O.cif"}
        )
    )
    status_payload = asyncio.run(client.get_task_status(submit_payload["task_id"], "test_user"))

    assert set(submit_payload) >= {"task_id", "status", "message"}
    assert submit_payload["status"] == "queued"
    assert set(status_payload) >= {"task_id", "status", "progress"}
