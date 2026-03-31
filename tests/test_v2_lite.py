"""
v2-lite refactoring tests.

Tests the core requirements:
1. Task creation works
2. ExecutionAttempt is created during worker run
3. Execution success auto-triggers AnalysisRun
4. Analysis success → get_task_status returns summary + html_report_url
5. Analysis failure preserves execution record
6. MD tasks don't fail due to missing 'success' field
7. A task can have multiple artifacts
8. Execution failure sets correct states
"""
import uuid
from unittest.mock import patch

import pytest

from src.vasp_server.task_manager.models import (
    Task, ExecutionAttempt, AnalysisRun, Artifact,
)
from tests.conftest import (
    make_mock_worker_result,
    run_worker_with_mock_result,
)


def _create_task(db_session, task_type="structure_optimization", **overrides):
    """Helper to insert a task directly."""
    task_id = uuid.uuid4().hex
    defaults = dict(
        id=task_id,
        user_id="test_user",
        task_type=task_type,
        status="queued",
        analysis_status="pending",
        params={"formula": "Li2O", "calc_type": "OXC"},
    )
    defaults.update(overrides)
    task = Task(**defaults)
    db_session.add(task)
    db_session.commit()
    return task_id


# ------------------------------------------------------------------ #
#  Test 1: Task creation
# ------------------------------------------------------------------ #
class TestTaskCreation:
    def test_submit_task_creates_task_row(self, task_manager, db_session):
        with patch.object(task_manager, "_run_task_worker"):
            task_id = task_manager.submit_task(
                user_id="test_user",
                task_type="structure_optimization",
                params={"formula": "Li2O", "calc_type": "OXC"},
            )

        task = db_session.get(Task, task_id)
        assert task is not None
        assert str(task.user_id) == "test_user"
        assert str(task.task_type) == "structure_optimization"
        assert str(task.status) == "queued"
        assert str(task.analysis_status) == "pending"
        assert task.input_snapshot is not None

    def test_submit_task_returns_task_id(self, task_manager):
        with patch.object(task_manager, "_run_task_worker"):
            task_id = task_manager.submit_task(
                user_id="u1", task_type="scf_calculation", params={}
            )
        assert task_id and len(task_id) == 32


# ------------------------------------------------------------------ #
#  Test 2: ExecutionAttempt creation
# ------------------------------------------------------------------ #
class TestExecutionAttempt:
    def test_worker_creates_execution_attempt(self, task_manager, db_session):
        task_id = _create_task(db_session, "structure_optimization")
        mock_result = make_mock_worker_result("structure_optimization", success=True)

        run_worker_with_mock_result(task_manager, task_id, mock_result, "structure_optimization")

        attempts = db_session.query(ExecutionAttempt).filter_by(task_id=task_id).all()
        assert len(attempts) == 1
        assert str(attempts[0].status) == "succeeded"
        assert attempts[0].attempt_no == 1


# ------------------------------------------------------------------ #
#  Test 3: Auto-analysis trigger
# ------------------------------------------------------------------ #
class TestAutoAnalysis:
    def test_execution_success_triggers_analysis_run(self, task_manager, db_session):
        task_id = _create_task(db_session, "scf_calculation")
        mock_result = make_mock_worker_result("scf_calculation", success=True)

        run_worker_with_mock_result(task_manager, task_id, mock_result, "scf_calculation")

        runs = db_session.query(AnalysisRun).filter_by(task_id=task_id).all()
        assert len(runs) == 1
        assert str(runs[0].status) == "completed"
        assert runs[0].summary is not None

        task = db_session.get(Task, task_id)
        assert str(task.status) == "completed"
        assert str(task.analysis_status) == "completed"
        assert task.result_summary is not None


# ------------------------------------------------------------------ #
#  Test 4: get_task_status returns summary + html_report_url
# ------------------------------------------------------------------ #
class TestGetTaskStatus:
    def test_completed_task_has_summary_and_report_url(self, task_manager, db_session):
        task_id = _create_task(db_session, "structure_optimization")
        mock_result = make_mock_worker_result("structure_optimization", success=True)

        run_worker_with_mock_result(task_manager, task_id, mock_result, "structure_optimization")

        fetched = task_manager.get_task(task_id, "test_user")
        assert fetched is not None
        assert str(fetched.status) == "completed"
        assert fetched.result_summary is not None
        assert fetched.result_summary.get("html_report_url") == "/static/test/report.html"


# ------------------------------------------------------------------ #
#  Test 5: Analysis failure preserves execution record
# ------------------------------------------------------------------ #
class TestAnalysisFailure:
    def test_analysis_failure_keeps_execution_succeeded(self, task_manager, db_session):
        task_id = _create_task(db_session, "structure_optimization")
        mock_result = make_mock_worker_result("structure_optimization", success=True)

        # Make analysis summary extraction blow up
        with patch.object(
            task_manager, "_extract_analysis_summary",
            side_effect=Exception("analysis exploded"),
        ):
            run_worker_with_mock_result(task_manager, task_id, mock_result, "structure_optimization")

        # Execution should still be succeeded
        attempts = db_session.query(ExecutionAttempt).filter_by(task_id=task_id).all()
        assert len(attempts) == 1
        assert str(attempts[0].status) == "succeeded"

        # Analysis should be failed
        runs = db_session.query(AnalysisRun).filter_by(task_id=task_id).all()
        assert len(runs) == 1
        assert str(runs[0].status) == "failed"

        # Task completed (execution ok), analysis failed
        task = db_session.get(Task, task_id)
        assert str(task.status) == "completed"
        assert str(task.analysis_status) == "failed"


# ------------------------------------------------------------------ #
#  Test 6: MD success field fix
# ------------------------------------------------------------------ #
class TestMDSuccessFix:
    def test_md_single_temp_convergence_counts_as_success(self, task_manager, db_session):
        """MD convergence=True without explicit success field should succeed."""
        task_id = _create_task(db_session, "md_calculation",
                               params={"formula": "Li2O", "calc_type": "OXC", "temperature": 300.0})

        # Result WITHOUT explicit 'success' field
        mock_result = {
            "convergence": True,
            "final_energy": -100.5,
            "average_temperature": 300.0,
            "total_md_steps": 1000,
            "work_directory": "/tmp/test_md",
            "computation_time": 500.0,
            "md_analysis_report_html_path": "/static/test/md_report.html",
        }

        run_worker_with_mock_result(task_manager, task_id, mock_result, "md_calculation")

        task = db_session.get(Task, task_id)
        assert str(task.status) == "completed"  # NOT failed

        attempts = db_session.query(ExecutionAttempt).filter_by(task_id=task_id).all()
        assert str(attempts[0].status) == "succeeded"

    def test_md_single_temp_explicit_success(self, task_manager, db_session):
        """MD with explicit success=True should succeed."""
        task_id = _create_task(db_session, "md_calculation",
                               params={"formula": "Li2O", "calc_type": "OXC", "temperature": 300.0})

        mock_result = {
            "success": True,
            "convergence": True,
            "final_energy": -100.5,
            "average_temperature": 300.0,
            "total_md_steps": 1000,
            "work_directory": "/tmp/test_md",
            "computation_time": 500.0,
            "md_analysis_report_html_path": "/static/test/md_report.html",
        }

        run_worker_with_mock_result(task_manager, task_id, mock_result, "md_calculation")

        task = db_session.get(Task, task_id)
        assert str(task.status) == "completed"

        attempts = db_session.query(ExecutionAttempt).filter_by(task_id=task_id).all()
        assert str(attempts[0].status) == "succeeded"


# ------------------------------------------------------------------ #
#  Test 7: Artifact queries
# ------------------------------------------------------------------ #
class TestArtifacts:
    def test_task_has_html_report_artifact(self, task_manager, db_session):
        """Completed task should have at least html_report artifact."""
        task_id = _create_task(db_session, "structure_optimization")
        mock_result = make_mock_worker_result("structure_optimization", success=True)

        run_worker_with_mock_result(task_manager, task_id, mock_result, "structure_optimization")

        artifacts = db_session.query(Artifact).filter_by(task_id=task_id).all()
        assert len(artifacts) >= 1

        html_arts = [a for a in artifacts if str(a.artifact_type) == "html_report"]
        assert len(html_arts) == 1
        assert html_arts[0].storage_key == "/static/test/report.html"

    def test_get_task_artifacts_method(self, task_manager, db_session):
        """task_manager.get_task_artifacts should return artifacts."""
        task_id = _create_task(db_session, "structure_optimization",
                               status="completed", analysis_status="completed")

        for atype in ["outcar", "contcar", "html_report"]:
            art = Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="execution",
                artifact_type=atype,
                storage_backend="local",
                storage_key=f"/tmp/{atype}",
            )
            db_session.add(art)
        db_session.commit()

        artifacts = task_manager.get_task_artifacts(task_id)
        assert len(artifacts) == 3
        types = {str(a.artifact_type) for a in artifacts}
        assert types == {"outcar", "contcar", "html_report"}


# ------------------------------------------------------------------ #
#  Test 8: Execution failure
# ------------------------------------------------------------------ #
class TestExecutionFailure:
    def test_execution_failure_sets_task_failed(self, task_manager, db_session):
        task_id = _create_task(db_session, "structure_optimization")
        mock_result = make_mock_worker_result("structure_optimization", success=False)

        run_worker_with_mock_result(task_manager, task_id, mock_result, "structure_optimization")

        task = db_session.get(Task, task_id)
        assert str(task.status) == "failed"

        attempts = db_session.query(ExecutionAttempt).filter_by(task_id=task_id).all()
        assert len(attempts) == 1
        assert str(attempts[0].status) == "runtime_failed"

        # No analysis should have been triggered
        runs = db_session.query(AnalysisRun).filter_by(task_id=task_id).all()
        assert len(runs) == 0
