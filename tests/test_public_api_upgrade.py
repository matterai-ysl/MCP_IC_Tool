import asyncio
import uuid
import json
from pathlib import Path
from types import SimpleNamespace
from urllib.parse import urlsplit
from datetime import datetime, timezone

import httpx
import pytest

from src.mcp_ic_tool.client import VaspAPIClient
from src.vasp_server.settings import settings
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
            return self._handle_response(response)

    async def _aget(self, path: str, params: dict):
        transport = httpx.ASGITransport(app=self.app)
        async with httpx.AsyncClient(transport=transport, base_url=self.base_url) as client:
            response = await client.get(path, params=params)
            return self._handle_response(response)


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


def _register_md_artifact(
    db_session,
    task_id: str,
    tmp_path: Path,
    relative_path: str,
    payload: str,
    *,
    content_type: str,
):
    from src.vasp_server.task_manager.models import Artifact

    local_path = tmp_path / task_id / relative_path
    local_path.parent.mkdir(parents=True, exist_ok=True)
    local_path.write_text(payload, encoding="utf-8")
    db_session.add(
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="json" if relative_path.endswith(".json") else "csv",
            storage_backend="local_public",
            storage_key=str(local_path),
            object_key=f"tenant/default/tasks/{task_id}/attempts/1/{relative_path}",
            mime_type=content_type,
            content_type=content_type,
        )
    )


def _register_public_artifact(
    db_session,
    task_id: str,
    filename: str,
    *,
    artifact_type: str,
    content_type: str = "text/plain",
):
    from src.vasp_server.task_manager.models import Artifact

    db_session.add(
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="execution",
            owner_id="attempt-1",
            artifact_type=artifact_type,
            storage_backend="local_public",
            storage_key=f"/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/{task_id}/attempts/1/{filename}",
            object_key=f"tenant/default/tasks/{task_id}/attempts/1/{filename}",
            mime_type=content_type,
            content_type=content_type,
        )
    )
    db_session.commit()


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


def test_submit_scf_persists_custom_incar(db_session):
    from src.vasp_server.vasp_server_api import app

    custom_incar = {
        "ENCUT": 300,
        "LREAL": "Auto",
        "LELF": False,
        "LWAVE": ".FALSE.",
    }

    response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={
            "user_id": "test_user",
            "cif_url": "https://structures.example.com/Li2O.cif",
            "custom_incar": custom_incar,
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["custom_incar"] == custom_incar


def test_local_public_artifact_root_defaults_to_absolute_project_path():
    assert Path(settings.local_public_artifact_root).is_absolute()
    assert Path(settings.local_public_artifact_root).name == "public_artifacts"


def test_submit_band_structure_accepts_input_url_list_and_persists_it(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/band-structure",
        json={
            "user_id": "test_user",
            "input_url": [
                "https://example.com/POSCAR",
                "https://example.com/POTCAR",
                "https://example.com/CHGCAR",
            ],
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["input_url"] == [
        "https://example.com/POSCAR",
        "https://example.com/POTCAR",
        "https://example.com/CHGCAR",
    ]
    assert "cif_url" not in task.params


def test_submit_band_structure_no_longer_accepts_cif_url(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/band-structure",
        json={
            "user_id": "test_user",
            "cif_url": "https://example.com/structure.cif",
        },
    )

    assert response.status_code == 422


def test_analyze_md_multi_can_aggregate_from_uploaded_diffusion_results(db_session, tmp_path):
    from src.vasp_server.vasp_server_api import app

    task_ids = []
    temperatures = [300.0, 400.0, 500.0]
    diffusivities = [1.0e-8, 2.0e-8, 4.0e-8]

    for idx, (temp, diffusivity) in enumerate(zip(temperatures, diffusivities), start=1):
        task_id = _create_task(
            db_session,
            user_id="test_user",
            task_type="md_calculation",
            result_summary={"total_md_steps": 1000},
        )
        task_ids.append(task_id)
        _register_md_artifact(
            db_session,
            task_id,
            tmp_path,
            "MD_output/data/diffusion_results.json",
            json.dumps(
                {
                    "temperature_K": temp,
                    "mobile_species": ["Li"],
                    "diffusion_by_element": {
                        "Li": {
                            "D_m2_per_s": diffusivity,
                            "D_cm2_per_s": diffusivity * 1e4,
                        }
                    },
                    "ionic_conductivity_S_per_m": {"Li": idx * 10.0},
                }
            ),
            content_type="application/json",
        )
        _register_md_artifact(
            db_session,
            task_id,
            tmp_path,
            "MD_output/data/rdf_summary.json",
            json.dumps(
                {
                    "coordination_numbers": {"Li-O": 4.0 + idx},
                    "peak_positions": {"Li-O": [2.0 + idx * 0.1]},
                }
            ),
            content_type="application/json",
        )
        _register_md_artifact(
            db_session,
            task_id,
            tmp_path,
            "MD_output/data/time_series.csv",
            "\n".join(
                [
                    "time_s,energy_eV,temperature_K,pressure_kB,a_A,b_A,c_A,volume_A3,density_kg_m3",
                    f"0.0,{-10.0 - idx},{temp - 5},1.0,4.1,4.1,4.1,{70 + idx},{2300 + idx * 10}",
                    f"1.0,{-10.2 - idx},{temp + 5},1.2,4.1,4.1,4.1,{71 + idx},{2310 + idx * 10}",
                ]
            ),
            content_type="text/csv",
        )
    db_session.commit()

    response = _request(
        app,
        "POST",
        "/vasp/analyze/md-multi",
        json={"user_id": "test_user", "task_ids": task_ids},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["success"] is True
    assert payload["analysis_type"] == "md_multi"
    assert payload["summary"]["task_count"] == 3
    assert payload["summary"]["arrhenius"]["fit"]["Ea_eV"] is not None
    assert payload["summary"]["stability"]["task_count"] == 3
    assert payload["summary"]["rdf"]["task_count"] == 3
    assert payload["html_report_url"].startswith("https://www.matterai.cn/dft/")
    rel_path = urlsplit(payload["html_report_url"]).path.split("/dft/", 1)[1]
    html_text = (Path(settings.local_public_artifact_root).resolve() / rel_path).read_text(encoding="utf-8")
    assert "Arrhenius" in html_text
    assert "稳定性摘要" in html_text
    assert "RDF 摘要" in html_text
    assert "<svg" in html_text


def test_analyze_md_multi_reports_missing_diffusion_results_cleanly(db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        user_id="test_user",
        task_type="md_calculation",
        result_summary={"total_md_steps": 1000},
    )

    response = _request(
        app,
        "POST",
        "/vasp/analyze/md-multi",
        json={"user_id": "test_user", "task_ids": [task_id]},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["success"] is False
    assert "diffusion_results.json" in payload["error_message"]


def test_analyze_dos_task_id_submits_remote_analysis_task(monkeypatch, db_session):
    import src.vasp_server.vasp_server_api as api

    source_task_id = _create_task(
        db_session,
        user_id="test_user",
        task_type="dos_calculation",
        status="completed",
        queue_name="hpc-dos",
        result_path="/data/home/ysl9527/vasp_calculations/hpc-c2ln6/source-dos",
        analysis_status="completed",
    )
    source_task = db_session.get(Task, source_task_id)
    assert source_task is not None

    captured = {}
    analysis_task = SimpleNamespace(
        id="analysis-dos-1",
        user_id="test_user",
        task_type="dos_analysis",
        status="queued",
        progress=0,
        analysis_status="pending",
        result_summary=None,
        params={"source_task_id": source_task_id},
        result_path=None,
        external_job_id=None,
        process_id=None,
        error_message=None,
        result_data=None,
        progress_message=None,
        created_at=datetime.now(timezone.utc),
        updated_at=datetime.now(timezone.utc),
    )

    def fake_submit_task(user_id, task_type, params):
        captured["user_id"] = user_id
        captured["task_type"] = task_type
        captured["params"] = params
        return analysis_task.id

    def fake_get_task(task_id, user_id):
        if task_id == source_task_id:
            return source_task
        if task_id == analysis_task.id:
            return analysis_task
        return None

    monkeypatch.setattr(api.task_manager, "submit_task", fake_submit_task)
    monkeypatch.setattr(api.task_manager, "get_task", fake_get_task)

    response = _request(
        api.app,
        "POST",
        "/vasp/analyze/dos",
        json={"user_id": "test_user", "task_id": source_task_id},
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["success"] is True
    assert payload["analysis_task_id"] == analysis_task.id
    assert payload["status"] == "queued"
    assert captured["task_type"] == "dos_analysis"
    assert captured["params"]["source_task_id"] == source_task_id
    assert captured["params"]["source_task_type"] == "dos_calculation"
    assert captured["params"]["queue_name"] == "hpc-dos"
    assert captured["params"]["source_work_directory"] == "/data/home/ysl9527/vasp_calculations/hpc-c2ln6/source-dos"


def test_agent_analyze_task_id_submits_remote_analysis_task(monkeypatch, db_session):
    import src.vasp_server.vasp_server_api as api

    source_task_id = _create_task(
        db_session,
        user_id="test_user",
        task_type="band_structure",
        status="completed",
        queue_name="hpc-band",
        result_path="/data/home/ysl9527/vasp_calculations/hpc-c2ln6/source-band",
        analysis_status="completed",
    )
    source_task = db_session.get(Task, source_task_id)
    assert source_task is not None

    captured = {}
    analysis_task = SimpleNamespace(
        id="agent-analysis-1",
        user_id="test_user",
        task_type="agent_analysis",
        status="queued",
        progress=0,
        analysis_status="pending",
        result_summary=None,
        params={"source_task_id": source_task_id},
        result_path=None,
        external_job_id=None,
        process_id=None,
        error_message=None,
        result_data=None,
        progress_message=None,
        created_at=datetime.now(timezone.utc),
        updated_at=datetime.now(timezone.utc),
    )

    def fake_submit_task(user_id, task_type, params):
        captured["user_id"] = user_id
        captured["task_type"] = task_type
        captured["params"] = params
        return analysis_task.id

    def fake_get_task(task_id, user_id):
        if task_id == source_task_id:
            return source_task
        if task_id == analysis_task.id:
            return analysis_task
        return None

    monkeypatch.setattr(api.task_manager, "submit_task", fake_submit_task)
    monkeypatch.setattr(api.task_manager, "get_task", fake_get_task)

    response = _request(
        api.app,
        "POST",
        "/vasp/agent/analyze",
        json={
            "user_id": "test_user",
            "task_id": source_task_id,
            "question": "请分析带隙与费米能级",
            "model": "qwen3.5-plus",
        },
    )

    assert response.status_code == 200
    payload = response.json()
    assert payload["success"] is True
    assert payload["analysis_task_id"] == analysis_task.id
    assert payload["status"] == "queued"
    assert captured["task_type"] == "agent_analysis"
    assert captured["params"]["source_task_id"] == source_task_id
    assert captured["params"]["queue_name"] == "hpc-band"
    assert captured["params"]["question"] == "请分析带隙与费米能级"
    assert captured["params"]["model"] == "qwen3.5-plus"


def test_agent_analyze_task_id_attaches_source_artifact_manifest(db_session):
    from src.vasp_server.vasp_server_api import app

    source_task_id = _create_task(
        db_session,
        user_id="test_user",
        task_type="structure_optimization",
        status="completed",
        queue_name="hpc-analysis",
        result_path=None,
        result_data={},
    )
    _register_public_artifact(db_session, source_task_id, "OUTCAR", artifact_type="outcar")

    response = _request(
        app,
        "POST",
        "/vasp/agent/analyze",
        json={
            "user_id": "test_user",
            "task_id": source_task_id,
            "question": "解释这个优化结果",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["analysis_task_id"])
    assert task is not None
    manifest = task.params["source_upstream_artifact_manifest"]
    assert manifest[0]["artifact_type"] == "outcar"
    assert manifest[0]["download_url"].endswith("/OUTCAR")


def test_submit_initial_task_accepts_explicit_queue_name(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/structure-optimization",
        json={
            "user_id": "test_user",
            "cif_url": "https://structures.example.com/Li2O.cif",
            "queue_name": "hpc-a",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.queue_name == "hpc-a"


def test_submit_structure_optimization_defaults_to_lighter_kpoint_density(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/structure-optimization",
        json={
            "user_id": "test_user",
            "cif_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["kpoint_density"] == 15.0


def test_submit_structure_optimization_persists_custom_incar(db_session):
    from src.vasp_server.vasp_server_api import app

    custom_incar = {
        "EDIFFG": -0.05,
        "ENCUT": 350,
        "PREC": "Normal",
        "EDIFF": 1e-4,
    }

    response = _request(
        app,
        "POST",
        "/vasp/structure-optimization",
        json={
            "user_id": "test_user",
            "cif_url": "https://structures.example.com/Li2O.cif",
            "custom_incar": custom_incar,
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["custom_incar"] == custom_incar


def test_submit_dos_defaults_to_lighter_template_settings(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/dos-calculation",
        json={
            "user_id": "test_user",
            "input_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["kpoint_density"] == 20.0
    assert task.params["kpoint_multiplier"] == 1.5
    assert task.params["precision"] == "High"


def test_submit_dos_accepts_input_url_list_and_persists_it(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/dos-calculation",
        json={
            "user_id": "test_user",
            "input_url": [
                "https://example.com/POSCAR",
                "https://example.com/POTCAR",
                "https://example.com/CHGCAR",
                "https://example.com/INCAR",
            ],
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["input_url"] == [
        "https://example.com/POSCAR",
        "https://example.com/POTCAR",
        "https://example.com/CHGCAR",
        "https://example.com/INCAR",
    ]
    assert "cif_url" not in task.params


def test_submit_dos_no_longer_accepts_cif_url(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/dos-calculation",
        json={
            "user_id": "test_user",
            "cif_url": "https://example.com/structure.cif",
        },
    )

    assert response.status_code == 422


def test_submit_band_structure_defaults_to_lighter_template_settings(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/band-structure",
        json={
            "user_id": "test_user",
            "input_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["kpoint_density"] == 20.0
    assert task.params["line_density"] == 20
    assert task.params["precision"] == "High"


def test_submit_band_structure_accepts_input_url_and_persists_it(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/band-structure",
        json={
            "user_id": "test_user",
            "input_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["input_url"] == "https://structures.example.com/Li2O.cif"
    assert "cif_url" not in task.params


def test_submit_neb_defaults_to_lighter_template_settings(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/neb-calculation",
        json={
            "user_id": "test_user",
            "initial_cif_url": "https://structures.example.com/Li2O-initial.cif",
            "final_cif_url": "https://structures.example.com/Li2O-final.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["kpoint_density"] == 15.0
    assert task.params["n_images"] == 5


def test_submit_neb_task_ids_attach_endpoint_artifact_manifests(db_session):
    from src.vasp_server.vasp_server_api import app

    initial_task_id = _create_task(db_session, queue_name="hpc-neb")
    final_task_id = _create_task(db_session, queue_name="hpc-neb")
    _register_public_artifact(db_session, initial_task_id, "CONTCAR", artifact_type="contcar")
    _register_public_artifact(db_session, final_task_id, "CONTCAR", artifact_type="contcar")

    response = _request(
        app,
        "POST",
        "/vasp/neb-calculation",
        json={
            "user_id": "test_user",
            "initial_task_id": initial_task_id,
            "final_task_id": final_task_id,
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["initial_upstream_artifact_manifest"][0]["artifact_type"] == "contcar"
    assert task.params["final_upstream_artifact_manifest"][0]["artifact_type"] == "contcar"
    assert task.params["initial_upstream_artifact_manifest"][0]["download_url"].endswith("/CONTCAR")
    assert task.queue_name == "hpc-neb"


def test_submit_md_accepts_input_url_and_persists_it(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/md-calculation",
        json={
            "user_id": "test_user",
            "input_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["input_url"] == "https://structures.example.com/Li2O.cif"
    assert "cif_url" not in task.params


def test_submit_md_rejects_input_url_list(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/md-calculation",
        json={
            "user_id": "test_user",
            "input_url": [
                "https://example.com/POSCAR",
                "https://example.com/POTCAR",
            ],
        },
    )

    assert response.status_code == 422


def test_submit_md_no_longer_accepts_cif_url(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/md-calculation",
        json={
            "user_id": "test_user",
            "cif_url": "https://example.com/structure.cif",
        },
    )

    assert response.status_code == 422


def test_submit_phonon_defaults_to_lighter_template_settings(db_session):
    from src.vasp_server.vasp_server_api import app

    response = _request(
        app,
        "POST",
        "/vasp/phonon-calculation",
        json={
            "user_id": "test_user",
            "cif_url": "https://structures.example.com/Li2O.cif",
        },
    )

    assert response.status_code == 200
    task = db_session.get(Task, response.json()["task_id"])
    assert task is not None
    assert task.params["kpoint_density"] == 15.0
    assert task.params["displacement"] == 0.015


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
    assert payload["artifacts"] is None


def test_submit_with_upstream_task_id_inherits_upstream_queue_and_only_matching_queue_can_claim(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    upstream_task_id = _create_task(
        db_session,
        task_type="structure_optimization",
        status="completed",
        analysis_status="completed",
        queue_name="hpc-a",
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
    assert task.queue_name == "hpc-a"
    assert "upstream_artifact_manifest" not in (task.params or {})

    other_queue_claim = task_manager.claim_next_task(worker_id="worker-b", queue_name="hpc-b")
    assert other_queue_claim is None

    matching_queue_claim = task_manager.claim_next_task(worker_id="worker-a", queue_name="hpc-a")
    assert matching_queue_claim is not None
    assert matching_queue_claim.id == task_id


def test_submit_with_conflicting_queue_name_and_upstream_task_is_rejected(db_session):
    from src.vasp_server.vasp_server_api import app

    upstream_task_id = _create_task(
        db_session,
        task_type="structure_optimization",
        status="completed",
        analysis_status="completed",
        queue_name="hpc-a",
    )

    response = _request(
        app,
        "POST",
        "/vasp/scf-calculation",
        json={
            "user_id": "test_user",
            "optimized_task_id": upstream_task_id,
            "queue_name": "hpc-b",
        },
    )

    assert response.status_code == 400
    assert "queue_name" in response.json()["detail"]


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


def test_submit_structure_optimization_returns_structured_quota_error(db_session):
    from src.vasp_server.vasp_server_api import app

    for idx in range(3):
        _create_task(
            db_session,
            id=uuid.uuid4().hex,
            user_id="test_user",
            status="running" if idx == 0 else "queued",
            worker_id="hpc-a" if idx == 0 else None,
            lease_token="lease-1" if idx == 0 else None,
        )

    response = _request(
        app,
        "POST",
        "/vasp/structure-optimization",
        json={"user_id": "test_user", "cif_url": "https://structures.example.com/Li2O.cif"},
    )

    assert response.status_code == 429
    payload = response.json()
    assert payload["error"]["code"] == "USER_CONCURRENT_TASK_LIMIT"
    assert payload["error"]["message"] == "已达到并发任务上限"
    assert payload["error"]["retryable"] is False
    assert payload["error"]["details"]["limit"] == 3
    assert payload["error"]["details"]["active_count"] == 3
    assert len(payload["error"]["details"]["active_task_ids"]) == 3


def test_mcp_client_surfaces_structured_error_for_agents(db_session):
    from src.vasp_server.vasp_server_api import app
    from src.mcp_ic_tool.client import VaspAPIError

    client = LocalVaspAPIClient(app)

    for idx in range(3):
        _create_task(
            db_session,
            id=uuid.uuid4().hex,
            user_id="test_user",
            status="running" if idx == 0 else "queued",
            worker_id="hpc-a" if idx == 0 else None,
            lease_token="lease-1" if idx == 0 else None,
        )

    with pytest.raises(VaspAPIError) as exc_info:
        asyncio.run(
            client.submit_structure_optimization(
                {"user_id": "test_user", "cif_url": "https://structures.example.com/Li2O.cif"}
            )
        )

    message = str(exc_info.value)
    assert "USER_CONCURRENT_TASK_LIMIT" in message
    assert "已达到并发任务上限" in message
    assert "active_count=3" in message
