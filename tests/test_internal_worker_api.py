import asyncio
from datetime import datetime
import uuid

import httpx
import pytest

from src.vasp_server.task_manager.models import Task


def _create_task(db_session, **overrides):
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="test_user",
        task_type="structure_optimization",
        status="queued",
        analysis_status="pending",
        params={
            "formula": "Li2O",
            "upstream_artifact_manifest": [{"artifact_type": "contcar", "object_key": "tenant/task/CONTCAR"}],
        },
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


@pytest.fixture(scope="module")
def app():
    from src.vasp_server.vasp_server_api import app

    return app


def _worker_headers(token: str = "worker-token") -> dict:
    return {"Authorization": f"Bearer {token}"}


def _request(app, method: str, url: str, **kwargs):
    async def _send():
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(transport=transport, base_url="http://testserver") as client:
            return await client.request(method, url, **kwargs)

    return asyncio.run(_send())


def test_worker_claim_endpoint_returns_leased_task(app, db_session):
    task_id = _create_task(db_session)

    resp = _request(
        app,
        "POST",
        "/internal/workers/claim",
        json={"worker_id": "hpc-a", "queue_name": "default"},
        headers=_worker_headers(),
    )

    assert resp.status_code == 200
    payload = resp.json()
    assert payload["task_id"] == task_id
    assert payload["status"] == "leased"
    assert payload["lease_token"]
    assert payload["params"]["formula"] == "Li2O"
    assert payload["upstream_artifact_manifest"] == [
        {"artifact_type": "contcar", "object_key": "tenant/task/CONTCAR"},
    ]


def test_internal_worker_api_rejects_invalid_token(app, db_session):
    _create_task(db_session)

    resp = _request(
        app,
        "POST",
        "/internal/workers/claim",
        json={"worker_id": "hpc-a", "queue_name": "default"},
        headers=_worker_headers("wrong-token"),
    )

    assert resp.status_code == 401


def test_worker_heartbeat_refreshes_lease_expiry(app, db_session):
    _create_task(db_session)
    claim_resp = _request(
        app,
        "POST",
        "/internal/workers/claim",
        json={"worker_id": "hpc-a", "queue_name": "default"},
        headers=_worker_headers(),
    )
    claim_payload = claim_resp.json()
    task_id = claim_payload["task_id"]
    old_expiry = datetime.fromisoformat(claim_payload["lease_expires_at"])

    heartbeat_resp = _request(
        app,
        "POST",
        f"/internal/tasks/{task_id}/heartbeat",
        json={"worker_id": "hpc-a", "lease_token": claim_payload["lease_token"]},
        headers=_worker_headers(),
    )

    assert heartbeat_resp.status_code == 200
    heartbeat_payload = heartbeat_resp.json()
    assert heartbeat_payload["status"] == "leased"
    assert datetime.fromisoformat(heartbeat_payload["lease_expires_at"]) >= old_expiry


def test_complete_and_fail_endpoints_reject_invalid_lease_tokens(app, db_session):
    _create_task(db_session)
    claim_resp = _request(
        app,
        "POST",
        "/internal/workers/claim",
        json={"worker_id": "hpc-a", "queue_name": "default"},
        headers=_worker_headers(),
    )
    claim_payload = claim_resp.json()
    task_id = claim_payload["task_id"]

    complete_resp = _request(
        app,
        "POST",
        f"/internal/tasks/{task_id}/complete",
        json={"worker_id": "hpc-a", "lease_token": "bad-token", "result_data": {"ok": True}},
        headers=_worker_headers(),
    )
    fail_resp = _request(
        app,
        "POST",
        f"/internal/tasks/{task_id}/fail",
        json={"worker_id": "hpc-a", "lease_token": "bad-token", "error_message": "boom"},
        headers=_worker_headers(),
    )

    assert complete_resp.status_code == 409
    assert fail_resp.status_code == 409
