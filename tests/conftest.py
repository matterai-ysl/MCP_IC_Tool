"""
Test fixtures for v2-lite refactoring tests.
Uses in-memory SQLite to avoid needing a real PostgreSQL instance.

IMPORTANT: We must patch the database module before any other project import,
because database.py creates the engine at module level.
"""
import sys
import pytest
from unittest.mock import AsyncMock, MagicMock

from sqlalchemy import create_engine, event
from sqlalchemy.orm import sessionmaker, declarative_base
from sqlalchemy.pool import StaticPool
import types

# ---------------------------------------------------------------------------
# Step 1: Create a test engine + session BEFORE importing project modules
# Use StaticPool so all connections share the same in-memory DB (needed for
# SQLite :memory: when multiple sessions/threads hit the same tables).
# ---------------------------------------------------------------------------
_test_engine = create_engine(
    "sqlite:///:memory:",
    echo=False,
    connect_args={"check_same_thread": False},
    poolclass=StaticPool,
)
_TestSessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=_test_engine)
_TestBase = declarative_base()

# SQLite doesn't enforce FK by default — enable it
@event.listens_for(_test_engine, "connect")
def _set_sqlite_pragma(dbapi_conn, connection_record):
    cursor = dbapi_conn.cursor()
    cursor.execute("PRAGMA foreign_keys=ON")
    cursor.close()


# ---------------------------------------------------------------------------
# Step 2: Monkey-patch the database module so project imports use SQLite
# ---------------------------------------------------------------------------
_fake_db_mod = types.ModuleType("src.vasp_server.task_manager.database")
_fake_db_mod.engine = _test_engine  # type: ignore
_fake_db_mod.SessionLocal = _TestSessionLocal  # type: ignore
_fake_db_mod.Base = _TestBase  # type: ignore
_fake_db_mod.init_db = lambda: _TestBase.metadata.create_all(bind=_test_engine)  # type: ignore
_fake_db_mod.check_and_init_db = lambda: None  # type: ignore

sys.modules["src.vasp_server.task_manager.database"] = _fake_db_mod

# Now safe to import project modules
from src.vasp_server.task_manager.models import Task, ExecutionAttempt, AnalysisRun, Artifact  # noqa: E402
from src.vasp_server.task_manager.manager import TaskManager  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture(autouse=True)
def setup_test_db():
    """Create tables before each test, drop after."""
    _TestBase.metadata.create_all(bind=_test_engine)
    yield
    _TestBase.metadata.drop_all(bind=_test_engine)


@pytest.fixture
def db_session():
    """Provide a fresh DB session for direct assertions."""
    session = _TestSessionLocal()
    try:
        yield session
    finally:
        session.close()


@pytest.fixture
def task_manager():
    """Provide a TaskManager instance."""
    return TaskManager()


# ---------------------------------------------------------------------------
# Mock VaspWorker helper
# ---------------------------------------------------------------------------
class FakeVaspWorker:
    """Replaces VaspWorker for tests — returns pre-set results."""

    def __init__(self, result_map: dict | None = None):
        self.result_map = result_map or {}
        self.base_work_dir = MagicMock()
        # Make Path / operations work
        self.base_work_dir.__truediv__ = lambda s, x: MagicMock(__str__=lambda s: "/tmp/test")

    def _make_method(self, result):
        return AsyncMock(return_value=result)

    def __getattr__(self, name):
        if name.startswith("run_"):
            task_type_key = name  # e.g. run_structure_optimization
            result = self.result_map.get(task_type_key, {"success": True, "work_directory": "/tmp"})
            return AsyncMock(return_value=result)
        raise AttributeError(name)


def run_worker_with_mock_result(task_manager, task_id, mock_result, task_type="structure_optimization"):
    """
    Run _run_task_worker with a mocked VaspWorker that returns mock_result.
    Patches the import inside manager._run_task_worker.
    """
    import threading

    method_name = {
        "structure_optimization": "run_structure_optimization",
        "scf_calculation": "run_scf_calculation",
        "dos_calculation": "run_dos_calculation",
        "md_calculation": "run_md_calculation",
    }[task_type]

    fake_worker = FakeVaspWorker({method_name: mock_result})

    # The manager imports VaspWorker inside _run_task_worker via `from ..vasp_worker import VaspWorker`
    # We patch it at the module level where it gets looked up
    import unittest.mock as um

    original_import = __builtins__.__import__ if hasattr(__builtins__, '__import__') else __import__

    def patched_import(name, *args, **kwargs):
        # When manager tries to import vasp_worker, return a fake module
        if 'vasp_worker' in name:
            mod = types.ModuleType(name)
            mod.VaspWorker = lambda user_id, base_work_dir: fake_worker  # type: ignore
            return mod
        return original_import(name, *args, **kwargs)

    cancel_event = threading.Event()

    with um.patch("src.vasp_server.task_manager.manager.VaspWorker", create=True,
                   new=lambda user_id, base_work_dir: fake_worker):
        # Also need to patch the lazy import inside the method
        import src.vasp_server.task_manager.manager as mgr_mod
        # Inject VaspWorker into the vasp_worker module namespace
        try:
            from src.vasp_server import vasp_worker as vw_mod
            orig_vw = getattr(vw_mod, 'VaspWorker', None)
            vw_mod.VaspWorker = lambda user_id, base_work_dir: fake_worker
        except ImportError:
            # vasp_worker may have heavy deps not available in test
            vw_mod = types.ModuleType("src.vasp_server.vasp_worker")
            vw_mod.VaspWorker = lambda user_id, base_work_dir: fake_worker  # type: ignore
            sys.modules["src.vasp_server.vasp_worker"] = vw_mod
            orig_vw = None

        task_manager._run_task_worker(task_id, cancel_event)

        # Restore
        if orig_vw is not None:
            vw_mod.VaspWorker = orig_vw


# ---------------------------------------------------------------------------
# Result builders
# ---------------------------------------------------------------------------
def make_mock_worker_result(task_type: str, success: bool = True) -> dict:
    base = {
        "success": success,
        "work_directory": "/tmp/test_work_dir",
        "computation_time": 123.4,
    }
    if task_type == "structure_optimization":
        base.update({
            "convergence": success,
            "force_convergence": success,
            "final_max_force": 0.01 if success else 1.0,
            "energy_convergence": success,
            "final_energy": -100.5 if success else None,
            "analysis_report_html_path": "/static/test/report.html" if success else None,
        })
    elif task_type == "scf_calculation":
        base.update({
            "convergence": success,
            "total_energy": -100.5 if success else None,
            "fermi_energy": -2.3 if success else None,
            "band_gap": 1.5 if success else None,
            "scf_analysis_report_html_path": "/static/test/scf_report.html" if success else None,
        })
    elif task_type == "dos_calculation":
        base.update({
            "convergence": success,
            "total_energy": -100.5 if success else None,
            "fermi_energy": -2.3 if success else None,
            "band_gap": 1.5 if success else None,
            "dos_analysis_report_html_path": "/static/test/dos_report.html" if success else None,
        })
    elif task_type == "md_calculation":
        base.update({
            "convergence": success,
            "final_energy": -100.5 if success else None,
            "average_temperature": 300.0,
            "total_md_steps": 1000,
            "md_analysis_report_html_path": "/static/test/md_report.html" if success else None,
        })
    if not success:
        base["error"] = "Test error"
    return base


def make_mock_md_multi_temp_result(completed: int = 2, failed: int = 1) -> dict:
    return {
        "success": completed > 0,
        "is_multi_temperature": True,
        "total_subtasks": completed + failed,
        "completed_subtasks": completed,
        "failed_subtasks": failed,
        "convergence": completed > 0,
        "work_directory": "/tmp/test_work_dir",
        "computation_time": 500.0,
        "md_analysis_report_html_path": "/static/test/multi_md_report.html",
        "subtask_results": [],
    }
