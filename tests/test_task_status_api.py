import asyncio
import json
import uuid

import httpx
from starlette.routing import Mount

from src.vasp_server.task_manager.models import Artifact, Task, ExecutionAttempt


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
    assert "html_report_url" not in payload["result_summary"]
    assert "analysis_report_html_path" not in payload["result_data"]
    serialized = json.dumps(payload)
    assert "/static/" not in serialized
    assert "/download/file/" not in serialized
    assert "/tmp/should-not-leak" not in serialized


def test_task_status_returns_suggested_action_for_scf_like_numeric_failure(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        status="failed",
        analysis_status="pending",
        result_summary=None,
        error_message="结构优化计算失败: BRMIX: very serious problems; Error EDDDAV: Call to ZHEGV failed.",
        result_data=None,
        result_path=None,
    )
    db_session.add(
        ExecutionAttempt(
            id=uuid.uuid4().hex,
            task_id=task_id,
            attempt_no=1,
            executor_type="slurm",
            status="runtime_failed",
            scheduler_job_id="12345",
            work_directory="/tmp/runtime-failure",
            failure_detail="BRMIX: very serious problems\nError EDDDAV: Call to ZHEGV failed.",
        )
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["failure_type"] == "electronic_divergence"
    assert "先对该结构做一次 SCF" in payload["suggested_action"]
    assert "CHGCAR" in payload["suggested_action"]


def test_task_status_returns_suggested_action_for_lattice_inconsistency(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        status="failed",
        analysis_status="pending",
        result_summary=None,
        error_message="结构优化计算失败: HNFORM: k-point generating vectors and reciprocal lattice are incommensurate.",
        result_data=None,
        result_path=None,
    )
    db_session.add(
        ExecutionAttempt(
            id=uuid.uuid4().hex,
            task_id=task_id,
            attempt_no=1,
            executor_type="slurm",
            status="runtime_failed",
            scheduler_job_id="12346",
            work_directory="/tmp/lattice-failure",
            failure_detail=(
                "Inconsistent Bravais lattice types found for crystalline and reciprocal lattice.\n"
                "HNFORM: k-point generating vectors and reciprocal lattice are incommensurate."
            ),
        )
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["failure_type"] == "lattice_inconsistency"
    assert "晶格参数" in payload["suggested_action"]
    assert "重新导出" in payload["suggested_action"] or "规范化" in payload["suggested_action"]


def test_remote_analysis_response_includes_failure_type(db_session):
    from src.vasp_server.vasp_server_api import _build_remote_analysis_response

    task_id = _create_task(
        db_session,
        task_type="md_analysis",
        status="failed",
        analysis_status="pending",
        result_summary=None,
        error_message="结构优化计算失败: BRMIX: very serious problems; Error EDDDAV: Call to ZHEGV failed.",
        result_data=None,
        result_path=None,
    )
    db_session.add(
        ExecutionAttempt(
            id=uuid.uuid4().hex,
            task_id=task_id,
            attempt_no=1,
            executor_type="slurm",
            status="runtime_failed",
            scheduler_job_id="12347",
            work_directory="/tmp/runtime-failure",
            failure_detail="BRMIX: very serious problems\nError EDDDAV: Call to ZHEGV failed.",
        )
    )
    db_session.commit()

    task = db_session.get(Task, task_id)
    response = _build_remote_analysis_response(task, "md")

    assert response.failure_type == "electronic_divergence"
    assert response.suggested_action


def test_task_status_falls_back_to_workdir_logs_for_failure_type(task_manager, db_session, tmp_path):
    from src.vasp_server.vasp_server_api import app

    work_dir = tmp_path / "vasp-failure"
    work_dir.mkdir()
    (work_dir / "result.log").write_text(
        "BRMIX: very serious problems\nError EDDDAV: Call to ZHEGV failed.\n",
        encoding="utf-8",
    )

    task_id = _create_task(
        db_session,
        status="failed",
        analysis_status="pending",
        result_summary=None,
        error_message="结构优化计算失败: VASP计算执行失败: 作业最终状态为 FAILED",
        result_data=None,
        result_path=None,
    )
    db_session.add(
        ExecutionAttempt(
            id=uuid.uuid4().hex,
            task_id=task_id,
            attempt_no=1,
            executor_type="slurm",
            status="runtime_failed",
            scheduler_job_id="12348",
            work_directory=str(work_dir),
            failure_detail="结构优化计算失败: VASP计算执行失败: 作业最终状态为 FAILED",
        )
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["failure_type"] == "electronic_divergence"
    assert "先对该结构做一次 SCF" in payload["suggested_action"]


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
    payload = response.json()
    assert payload["html_report_url"].startswith("https://")
    assert payload["artifacts"] is None


def test_task_status_artifact_list_omits_non_public_local_artifacts(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    db_session.add_all(
        [
            Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="analysis",
                owner_id="run-1",
                artifact_type="html_report",
                storage_backend="local_public",
                storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-1/attempts/1/report.html",
                object_key="tenant/default/tasks/task-1/attempts/1/report.html",
                mime_type="text/html",
                content_type="text/html",
            ),
            Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="execution",
                owner_id="attempt-1",
                artifact_type="potcar",
                storage_backend="local",
                storage_key="/tmp/private/POTCAR",
                object_key=None,
                mime_type="text/plain",
                content_type="text/plain",
            ),
        ]
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["artifacts"] is None
    assert payload["html_report_url"].endswith("/report.html")


def test_task_status_prefers_signed_html_artifact_even_for_legacy_file_rows(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    db_session.add(
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="file",
            storage_backend="oss",
            storage_key="tenant/tenant-a/tasks/task-x/attempts/1/optimization_analysis_report.html",
            object_key="tenant/tenant-a/tasks/task-x/attempts/1/optimization_analysis_report.html",
            mime_type="text/html",
            content_type="text/html",
        )
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].startswith("https://")
    assert "api.matterai.tech" not in payload["html_report_url"]
    assert "html_report_url" not in payload["result_summary"]


def test_task_status_exposes_preview_images_and_data_downloads(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    artifacts = [
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="file",
            storage_backend="oss",
            storage_key="tenant/tenant-a/tasks/task-x/attempts/1/plots/force_convergence.png",
            object_key="tenant/tenant-a/tasks/task-x/attempts/1/plots/force_convergence.png",
            mime_type="image/png",
            content_type="image/png",
        ),
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="file",
            storage_backend="oss",
            storage_key="tenant/tenant-a/tasks/task-x/attempts/1/data/force_convergence.csv",
            object_key="tenant/tenant-a/tasks/task-x/attempts/1/data/force_convergence.csv",
            mime_type="text/csv",
            content_type="text/csv",
        ),
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="file",
            storage_backend="oss",
            storage_key="tenant/tenant-a/tasks/task-x/attempts/1/data/summary.json",
            object_key="tenant/tenant-a/tasks/task-x/attempts/1/data/summary.json",
            mime_type="application/json",
            content_type="application/json",
        ),
    ]
    db_session.add_all(artifacts)
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["preview_images"][0]["filename"] == "plots/force_convergence.png"
    assert payload["preview_images"][0]["download_url"].startswith("https://")
    returned_filenames = {item["filename"] for item in payload["data_downloads"]}
    assert returned_filenames == {"data/force_convergence.csv", "data/summary.json"}
    assert all(item["download_url"].startswith("https://") for item in payload["data_downloads"])


def test_task_status_keeps_report_assets_out_of_generic_artifact_list(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(db_session)
    artifacts = [
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="html_report",
            storage_backend="local_public",
            storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/optimization_analysis_report.html",
            object_key="tenant/default/tasks/task-x/attempts/1/optimization_analysis_report.html",
            mime_type="text/html",
            content_type="text/html",
        ),
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="image",
            storage_backend="local_public",
            storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/force_convergence.png",
            object_key="tenant/default/tasks/task-x/attempts/1/force_convergence.png",
            mime_type="image/png",
            content_type="image/png",
        ),
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="analysis",
            owner_id="run-1",
            artifact_type="csv",
            storage_backend="local_public",
            storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/force_convergence.csv",
            object_key="tenant/default/tasks/task-x/attempts/1/force_convergence.csv",
            mime_type="text/csv",
            content_type="text/csv",
        ),
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="execution",
            owner_id="attempt-1",
            artifact_type="contcar",
            storage_backend="local_public",
            storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/CONTCAR",
            object_key="tenant/default/tasks/task-x/attempts/1/CONTCAR",
            mime_type="text/plain",
            content_type="text/plain",
        ),
    ]
    db_session.add_all(artifacts)
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].endswith("/optimization_analysis_report.html")
    assert {item["filename"] for item in payload["preview_images"]} == {"force_convergence.png"}
    assert {item["filename"] for item in payload["data_downloads"]} == {"force_convergence.csv"}
    assert {item["filename"] for item in payload["artifacts"]} == {"CONTCAR"}


def test_task_status_rewrites_legacy_public_download_urls(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        result_data={
            "optimized_structure_download_url": "https://api.matterai.tech/download/file/hpc-c2ln6/task-1/CONTCAR",
            "band_structure_path": "https://api.matterai.tech/static/hpc-c2ln6/task-1/bands.png",
        },
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["result_data"]["optimized_structure_download_url"].startswith("https://www.matterai.cn/")
    assert "api.matterai.tech" not in payload["result_data"]["optimized_structure_download_url"]
    assert payload["result_data"]["band_structure_path"].startswith("https://www.matterai.cn/")


def test_task_status_hides_same_host_static_html_report_urls_without_public_artifacts(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        result_summary=None,
        result_data={
            "analysis_report_html_path": "https://www.matterai.cn/static/hpc-c2ln6/task-1/optimization_analysis_report.html",
            "work_directory": "/tmp/should-not-leak",
        },
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"] is None
    assert "analysis_report_html_path" not in payload["result_data"]


def test_task_status_exposes_html_report_url_only_once(task_manager, db_session):
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
    assert "html_report_url" not in (payload.get("result_summary") or {})
    assert "analysis_report_html_path" not in (payload.get("result_data") or {})


def test_task_status_hides_agent_analysis_local_report_paths(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        task_type="agent_analysis",
        result_summary={"html_report_url": "/data/home/ysl9527/vasp_calculations/report.html"},
        result_data={
            "html_report_url": "/data/home/ysl9527/vasp_calculations/report.html",
            "agent_analysis_report_html_path": "/data/home/ysl9527/vasp_calculations/report.html",
            "work_directory": "/data/home/ysl9527/vasp_calculations/task-1",
        },
    )
    task_manager.register_artifact(
        task_id=task_id,
        artifact_type="html_report",
        owner_type="analysis",
        owner_id="run-1",
        filename="agent_analysis/agent_analysis_report.html",
        content_type="text/html",
        size_bytes=1024.0,
        attempt_no=1,
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].startswith("https://")
    assert "html_report_url" not in (payload.get("result_summary") or {})
    assert "html_report_url" not in (payload.get("result_data") or {})
    assert "agent_analysis_report_html_path" not in (payload.get("result_data") or {})
    serialized = json.dumps(payload, ensure_ascii=False)
    assert "/data/home/ysl9527/vasp_calculations/" not in serialized


def test_task_status_prefers_public_contcar_artifact_for_optimized_structure_download(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        result_summary=None,
        result_data={
            "optimized_structure_download_url": "https://www.matterai.cn/download/file/hpc-c2ln6/task-1/CONTCAR",
        },
    )
    db_session.add(
        Artifact(
            id=uuid.uuid4().hex,
            task_id=task_id,
            owner_type="execution",
            owner_id="attempt-1",
            artifact_type="contcar",
            storage_backend="local_public",
            storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-1/attempts/1/CONTCAR",
            object_key="tenant/default/tasks/task-1/attempts/1/CONTCAR",
            mime_type="text/plain",
            content_type="text/plain",
        )
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["result_data"]["optimized_structure_download_url"] == (
        "https://www.matterai.cn/dft/tenant/default/tasks/task-1/attempts/1/CONTCAR"
    )


def test_task_status_hides_large_dos_payloads_and_exposes_report_assets(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        task_type="dos_calculation",
        result_summary=None,
        result_data={
            "dos_analysis_report_html_path": "/tmp/private/pymatgen_dos_analysis_report.html",
            "dos_data": {
                "total_dos": {
                    "energy": [0.0, 1.0],
                    "dos_total": [1.0, 2.0],
                }
            },
            "analysis_data": {
                "raw_dos_data": {
                    "TDOS": {
                        "energy": [0.0, 1.0],
                        "dos": [1.0, 2.0],
                    }
                }
            },
        },
    )
    db_session.add_all(
        [
            Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="analysis",
                owner_id="run-1",
                artifact_type="html_report",
                storage_backend="local_public",
                storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/dos_analysis/pymatgen_dos_analysis_report.html",
                object_key="tenant/default/tasks/task-x/attempts/1/dos_analysis/pymatgen_dos_analysis_report.html",
                mime_type="text/html",
                content_type="text/html",
            ),
            Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="analysis",
                owner_id="run-1",
                artifact_type="file",
                storage_backend="local_public",
                storage_key="/srv/mcp-ic-tool/public_artifacts/tenant/default/tasks/task-x/attempts/1/dos_analysis/total_dos.csv",
                object_key="tenant/default/tasks/task-x/attempts/1/dos_analysis/total_dos.csv",
                mime_type="text/csv",
                content_type="text/csv",
            ),
        ]
    )
    db_session.commit()

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert payload["html_report_url"].endswith("/dos_analysis/pymatgen_dos_analysis_report.html")
    assert {item["filename"] for item in payload["data_downloads"]} == {"dos_analysis/total_dos.csv"}
    assert "dos_data" not in (payload.get("result_data") or {})
    assert "analysis_data" not in (payload.get("result_data") or {})


def test_task_status_hides_md_trajectory_payloads(task_manager, db_session):
    from src.vasp_server.vasp_server_api import app

    task_id = _create_task(
        db_session,
        task_type="md_calculation",
        result_summary=None,
        result_data={
            "trajectory_data": {
                "total_steps": 1000,
                "step_intervals": list(range(1, 1001)),
            },
            "final_energy": -10.0,
        },
    )

    response = _request(app, "GET", f"/vasp/task/{task_id}", params={"user_id": "test_user"})

    assert response.status_code == 200
    payload = response.json()
    assert "trajectory_data" not in (payload.get("result_data") or {})
