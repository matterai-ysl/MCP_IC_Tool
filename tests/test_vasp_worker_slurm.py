import asyncio

import pytest


class _FakeProcess:
    def __init__(self, returncode=1, stdout=b"", stderr=b"submit failed"):
        self.returncode = returncode
        self._stdout = stdout
        self._stderr = stderr

    async def communicate(self):
        return self._stdout, self._stderr


def test_run_vasp_calculation_defaults_to_absolute_vasp_path(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker
    from src.vasp_server.settings import settings

    monkeypatch.setenv("VASP_PATH", "/opt/vasp/bin/vasp_std")
    monkeypatch.delenv("SLURM_MODULE_LOAD", raising=False)
    monkeypatch.delenv("ONEAPI_ENV_SCRIPT", raising=False)
    monkeypatch.delenv("VASP_SLURM_RUN_LINE", raising=False)
    monkeypatch.setattr(settings, "slurm_module_load", None)
    monkeypatch.setattr(settings, "oneapi_env_script", None)
    monkeypatch.setattr(settings, "vasp_slurm_run_line", None)

    async def fake_create_subprocess_shell(*args, **kwargs):
        return _FakeProcess()

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.asyncio.create_subprocess_shell",
        fake_create_subprocess_shell,
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(Exception, match="SLURM作业提交失败"):
        asyncio.run(worker._run_vasp_calculation(work_dir))

    script = (work_dir / "vasp_job.sh").read_text(encoding="utf-8")
    assert "module load" not in script
    assert "source /data/app/intel/oneapi-2023.2/setvars.sh" not in script
    assert "mpirun -np $SLURM_NPROCS /opt/vasp/bin/vasp_std>result.log 2>&1" in script


def test_run_vasp_calculation_can_export_intel_pmi_and_use_srun_mpi(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker
    from src.vasp_server.settings import settings

    monkeypatch.setenv("VASP_PATH", "/opt/vasp/bin/vasp_std")
    monkeypatch.setenv("I_MPI_PMI_LIBRARY", "/usr/lib64/libpmi2.so")
    monkeypatch.setenv("SLURM_SRUN_MPI", "pmi2")
    monkeypatch.delenv("VASP_SLURM_RUN_LINE", raising=False)
    monkeypatch.setattr(settings, "vasp_slurm_run_line", None)
    monkeypatch.setattr(settings, "slurm_srun_mpi", None, raising=False)
    monkeypatch.setattr(settings, "i_mpi_pmi_library", None, raising=False)

    async def fake_create_subprocess_shell(*args, **kwargs):
        return _FakeProcess()

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.asyncio.create_subprocess_shell",
        fake_create_subprocess_shell,
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "task-pmi"
    work_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(Exception, match="SLURM作业提交失败"):
        asyncio.run(worker._run_vasp_calculation(work_dir))

    script = (work_dir / "vasp_job.sh").read_text(encoding="utf-8")
    assert "export I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so" in script
    assert "srun --mpi=pmi2 /opt/vasp/bin/vasp_std >result.log 2>&1" in script


def test_band_structure_seed_scf_does_not_require_missing_convergence_field(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    task_id = "band-seed-task"
    run_calls = []

    async def fake_prepare_band_structure_files(work_dir, params, progress_callback=None):
        return {"_band_input_mode": "seed_scf"}

    async def fake_prepare_seed_files(work_dir, params):
        (work_dir / "CHGCAR").write_text("seed-chgcar", encoding="utf-8")
        (work_dir / "INCAR").write_text("SYSTEM = BAND_SCF\n", encoding="utf-8")

    async def fake_run_vasp_calculation(work_dir, progress_callback=None):
        run_calls.append(str(work_dir))
        return {
            "success": True,
            "slurm_job_id": f"job-{len(run_calls)}",
            "work_directory": str(work_dir),
        }

    async def fake_generate_band_structure_inputs(work_dir, params, bs_files):
        return None

    async def fake_analyze_band_structure_results(work_dir, vasp_result):
        return {
            "success": True,
            "band_gap": 1.23,
            "work_directory": str(work_dir),
        }

    monkeypatch.setattr(worker, "_prepare_band_structure_files", fake_prepare_band_structure_files)
    monkeypatch.setattr(worker, "_prepare_scf_then_band_structure_files", fake_prepare_seed_files)
    monkeypatch.setattr(worker, "_run_vasp_calculation", fake_run_vasp_calculation)
    monkeypatch.setattr(worker, "_generate_band_structure_inputs", fake_generate_band_structure_inputs)
    monkeypatch.setattr(worker, "_analyze_band_structure_results", fake_analyze_band_structure_results)

    result = asyncio.run(
        worker.run_band_structure_calculation(
            task_id,
            {"input_url": "https://example.com/structure.cif"},
        )
    )

    assert result["success"] is True
    assert run_calls == [
        str(worker.base_work_dir / task_id),
        str(worker.base_work_dir / task_id),
    ]


def test_check_convergence_does_not_match_voluntary_context_switches(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text(
        "General timing and accounting informations for this job:\n"
        "Voluntary context switches: 42\n",
        encoding="utf-8",
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))

    assert worker._check_convergence(outcar_path) is False


def test_analyze_results_uses_real_task_id_for_optimization_report(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "task-42"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("OUTCAR", encoding="utf-8")
    (work_dir / "CONTCAR").write_text("CONTCAR", encoding="utf-8")
    (work_dir / "optimization_analysis_report.html").write_text("<html/>", encoding="utf-8")

    captured = {}

    class FakeAnalyzer:
        def __init__(self, input_path, task_id=None):
            captured["analyzer_task_id"] = task_id

        def analyze(self):
            return {
                "convergence_analysis": {
                    "force_convergence": {"converged": True, "final_max_force": 0.01},
                    "energy_convergence": {"converged": True, "final_energy": -12.34},
                }
            }

    def fake_report(input_path, task_id=None, output_dir=None):
        captured["report_task_id"] = task_id
        return str(work_dir / "optimization_analysis_report.html")

    monkeypatch.setattr("src.vasp_server.optimization_analyzer.OUTCARAnalyzer", FakeAnalyzer)
    monkeypatch.setattr("src.vasp_server.optimization_analyzer.generate_optimization_report", fake_report)

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    monkeypatch.setattr(worker, "_check_convergence", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(worker, "_extract_energy", lambda *_args, **_kwargs: -12.34)
    monkeypatch.setattr(worker, "_extract_forces", lambda *_args, **_kwargs: [])

    result = asyncio.run(
        worker._analyze_results(
            "task-real-id",
            work_dir,
            {"computation_time": 1.23, "process_id": "slurm-1"},
        )
    )

    assert captured["analyzer_task_id"] == "task-real-id"
    assert captured["report_task_id"] == "task-real-id"
    assert result["analysis_report_html_path"]


def test_analyze_results_marks_missing_ionic_force_history_as_failure(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "task-bad-opt"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("OUTCAR", encoding="utf-8")
    (work_dir / "CONTCAR").write_text("CONTCAR", encoding="utf-8")

    class FakeAnalyzer:
        def __init__(self, input_path, task_id=None):
            pass

        def analyze(self):
            return {
                "convergence_analysis": {
                    "force_convergence": {
                        "available": False,
                        "converged": None,
                        "final_max_force": None,
                        "force_history": [],
                    },
                    "energy_convergence": {
                        "available": True,
                        "converged": False,
                        "final_energy": -11.95,
                    },
                },
                "optimization_process": {
                    "total_steps": 0,
                    "data_available": False,
                },
            }

    monkeypatch.setattr("src.vasp_server.optimization_analyzer.OUTCARAnalyzer", FakeAnalyzer)
    monkeypatch.setattr(
        "src.vasp_server.optimization_analyzer.generate_optimization_report",
        lambda *_args, **_kwargs: str(work_dir / "optimization_analysis_report.html"),
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    monkeypatch.setattr(worker, "_check_convergence", lambda *_args, **_kwargs: False)
    monkeypatch.setattr(worker, "_extract_energy", lambda *_args, **_kwargs: -11.95)
    monkeypatch.setattr(worker, "_extract_forces", lambda *_args, **_kwargs: [])

    result = asyncio.run(
        worker._analyze_results(
            "task-bad-opt",
            work_dir,
            {"computation_time": 1.0, "process_id": "slurm-2", "stdout": "MPI startup(): PMI server not found"},
        )
    )

    assert result["success"] is False
    assert "离子步" in result["error_message"]


def test_analyze_results_allows_zero_step_but_tail_converged_optimization(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "task-zero-step-opt"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("OUTCAR", encoding="utf-8")
    (work_dir / "CONTCAR").write_text("CONTCAR", encoding="utf-8")
    (work_dir / "optimization_analysis_report.html").write_text("<html/>", encoding="utf-8")

    class FakeAnalyzer:
        def __init__(self, input_path, task_id=None):
            pass

        def analyze(self):
            return {
                "convergence_analysis": {
                    "tail_check": {"matched": True},
                    "force_convergence": {
                        "available": False,
                        "converged": None,
                        "final_max_force": None,
                        "force_history": [],
                    },
                    "energy_convergence": {
                        "available": False,
                        "converged": None,
                        "final_energy": -11.95,
                        "energy_history": [-11.95],
                    },
                },
                "optimization_process": {
                    "total_steps": 0,
                    "data_available": False,
                },
            }

    monkeypatch.setattr("src.vasp_server.optimization_analyzer.OUTCARAnalyzer", FakeAnalyzer)
    monkeypatch.setattr(
        "src.vasp_server.optimization_analyzer.generate_optimization_report",
        lambda *_args, **_kwargs: str(work_dir / "optimization_analysis_report.html"),
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    monkeypatch.setattr(worker, "_check_convergence", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(worker, "_extract_energy", lambda *_args, **_kwargs: -11.95)
    monkeypatch.setattr(worker, "_extract_forces", lambda *_args, **_kwargs: [])

    result = asyncio.run(
        worker._analyze_results(
            "task-zero-step-opt",
            work_dir,
            {"computation_time": 1.0, "process_id": "slurm-3", "stdout": ""},
        )
    )

    assert result["success"] is True
    assert "error_message" not in result
    assert result["final_energy"] == -11.95
