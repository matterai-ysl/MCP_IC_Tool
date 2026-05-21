import asyncio
from pathlib import Path

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
    monkeypatch.delenv("SLURM_GPUS", raising=False)
    monkeypatch.delenv("ONEAPI_ENV_SCRIPT", raising=False)
    monkeypatch.delenv("VASP_SLURM_RUN_LINE", raising=False)
    monkeypatch.setattr(settings, "slurm_module_load", None)
    monkeypatch.setattr(settings, "oneapi_env_script", None)
    monkeypatch.setattr(settings, "vasp_slurm_run_line", None)
    monkeypatch.setattr(settings, "slurm_gpus", None, raising=False)

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
    assert "vasp_exit_code=$?" in script
    assert 'exit "$vasp_exit_code"' in script
    assert "#SBATCH --gpus=" not in script


def test_run_vasp_calculation_can_request_slurm_gpus(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker
    from src.vasp_server.settings import settings

    monkeypatch.setenv("SLURM_GPUS", "1")
    monkeypatch.setattr(settings, "slurm_gpus", None, raising=False)

    async def fake_create_subprocess_shell(*args, **kwargs):
        return _FakeProcess()

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.asyncio.create_subprocess_shell",
        fake_create_subprocess_shell,
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "task-gpu"
    work_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(Exception, match="SLURM作业提交失败"):
        asyncio.run(worker._run_vasp_calculation(work_dir))

    script = (work_dir / "vasp_job.sh").read_text(encoding="utf-8")
    assert "#SBATCH --gpus=1" in script


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


def test_run_vasp_calculation_accepts_neb_image_outcars(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    commands = []

    async def fake_create_subprocess_shell(command, *args, **kwargs):
        commands.append(command)
        if command.startswith("sbatch"):
            return _FakeProcess(returncode=0, stdout=b"Submitted batch job 12345\n", stderr=b"")
        if command.startswith("squeue"):
            return _FakeProcess(returncode=0, stdout=b"", stderr=b"")
        if command.startswith("sacct"):
            return _FakeProcess(returncode=0, stdout=b"COMPLETED\n", stderr=b"")
        return _FakeProcess(returncode=1, stdout=b"", stderr=b"unexpected command")

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.asyncio.create_subprocess_shell",
        fake_create_subprocess_shell,
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "neb-task"
    (work_dir / "01").mkdir(parents=True)
    (work_dir / "02").mkdir(parents=True)
    (work_dir / "01" / "OUTCAR").write_text("image 1 outcar", encoding="utf-8")
    (work_dir / "02" / "OUTCAR").write_text("image 2 outcar", encoding="utf-8")
    (work_dir / "result.log").write_text("NEB finished", encoding="utf-8")

    result = asyncio.run(worker._run_vasp_calculation(work_dir, output_check="neb"))

    assert result["success"] is True
    assert result["slurm_job_id"] == "12345"
    assert str(work_dir / "01" / "OUTCAR") in result["output_files"]
    assert any(command.startswith("sacct") for command in commands)


def test_run_vasp_calculation_tolerates_gpu_deallocate_after_complete_outputs(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    async def fake_create_subprocess_shell(command, *args, **kwargs):
        if command.startswith("sbatch"):
            return _FakeProcess(returncode=0, stdout=b"Submitted batch job 12345\n", stderr=b"")
        if command.startswith("squeue"):
            return _FakeProcess(returncode=0, stdout=b"", stderr=b"")
        if command.startswith("sacct"):
            return _FakeProcess(returncode=0, stdout=b"FAILED\n", stderr=b"")
        return _FakeProcess(returncode=1, stdout=b"", stderr=b"unexpected command")

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.asyncio.create_subprocess_shell",
        fake_create_subprocess_shell,
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "gpu-deallocate-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("free  energy   TOTEN  =       -33.04916306 eV\n", encoding="utf-8")
    (work_dir / "CHGCAR").write_text("charge density", encoding="utf-8")
    (work_dir / "AECCAR0").write_text("aeccar0", encoding="utf-8")
    (work_dir / "AECCAR2").write_text("aeccar2", encoding="utf-8")
    (work_dir / "result.log").write_text(
        "   1 F= -.33049163E+02 E0= -.33049163E+02  d E =-.123760E-13\n"
        "forrtl: severe (173): A pointer passed to DEALLOCATE points to an object that cannot be deallocated\n",
        encoding="utf-8",
    )

    result = asyncio.run(worker._run_vasp_calculation(work_dir))

    assert result["success"] is True
    assert result["slurm_job_id"] == "12345"


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


def test_band_structure_seed_scf_tolerates_delayed_chgcar_visibility(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    monkeypatch.setattr(worker, "_required_output_wait_timeout_seconds", 0.1, raising=False)
    monkeypatch.setattr(worker, "_required_output_wait_poll_seconds", 0.01, raising=False)

    task_id = "band-delayed-chgcar"
    run_calls = []

    async def fake_prepare_band_structure_files(work_dir, params, progress_callback=None):
        return {"_band_input_mode": "seed_scf"}

    async def fake_prepare_seed_files(work_dir, params):
        (work_dir / "INCAR").write_text("SYSTEM = BAND_SCF\n", encoding="utf-8")

    async def fake_run_vasp_calculation(work_dir, progress_callback=None):
        run_calls.append(str(work_dir))
        if len(run_calls) == 1:
            async def _materialize_chgcar():
                await asyncio.sleep(0.02)
                (work_dir / "CHGCAR").write_text("seed-chgcar", encoding="utf-8")

            asyncio.create_task(_materialize_chgcar())
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
            "band_gap": 0.42,
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


def test_analyze_band_structure_results_uses_low_memory_vasprun_parse(monkeypatch, tmp_path):
    import sys
    import types
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "band-analysis"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("OUTCAR", encoding="utf-8")
    (work_dir / "vasprun.xml").write_text("<xml/>", encoding="utf-8")
    (work_dir / "POSCAR").write_text("POSCAR", encoding="utf-8")

    captured = {}

    class FakeBandStructure:
        def get_band_gap(self):
            return {"energy": 0.8, "direct": False}

        def get_vbm(self):
            return {"energy": 0.0}

        def get_cbm(self):
            return {"energy": 0.8}

    class FakeVasprun:
        def __init__(self, path, parse_projected_eigen=False):
            captured["path"] = path
            captured["parse_projected_eigen"] = parse_projected_eigen

        def get_band_structure(self, line_mode=True):
            captured["line_mode"] = line_mode
            return FakeBandStructure()

    pymatgen_module = types.ModuleType("pymatgen")
    io_module = types.ModuleType("pymatgen.io")
    vasp_module = types.ModuleType("pymatgen.io.vasp")
    electronic_module = types.ModuleType("pymatgen.electronic_structure")
    plotter_module = types.ModuleType("pymatgen.electronic_structure.plotter")
    plotter_module.BSPlotter = object
    vasp_module.Vasprun = FakeVasprun
    monkeypatch.setitem(sys.modules, "pymatgen", pymatgen_module)
    monkeypatch.setitem(sys.modules, "pymatgen.io", io_module)
    monkeypatch.setitem(sys.modules, "pymatgen.io.vasp", vasp_module)
    monkeypatch.setitem(sys.modules, "pymatgen.electronic_structure", electronic_module)
    monkeypatch.setitem(sys.modules, "pymatgen.electronic_structure.plotter", plotter_module)

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    monkeypatch.setattr(worker, "_check_convergence", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(worker, "_extract_energy", lambda *_args, **_kwargs: -22.5)
    monkeypatch.setattr(worker, "_extract_fermi_energy", lambda *_args, **_kwargs: 1.2)
    monkeypatch.setattr(
        "src.vasp_server.analyzers.band_structure.BandStructureAnalyzer",
        lambda *_args, **_kwargs: types.SimpleNamespace(analyze=lambda: {"ok": True}),
    )
    monkeypatch.setattr(
        "src.vasp_server.analyzers.band_structure.generate_band_structure_report",
        lambda *_args, **_kwargs: None,
    )

    result = asyncio.run(
        worker._analyze_band_structure_results(
            work_dir,
            {"success": True, "computation_time": 1.0, "process_id": 123},
        )
    )

    assert captured["parse_projected_eigen"] is False
    assert captured["line_mode"] is True
    assert result["success"] is True
    assert result["band_gap"] == 0.8


def test_generate_band_structure_kpoints_falls_back_when_pymatgen_returns_no_path(monkeypatch, tmp_path):
    import sys
    import types
    from src.vasp_server.vasp_worker import VaspWorker

    class FakeStructure:
        @staticmethod
        def from_file(_path):
            return object()

    class FakeHighSymmKpath:
        def __init__(self, _structure):
            self.kpath = None

    monkeypatch.setitem(sys.modules, "pymatgen.core", types.SimpleNamespace(Structure=FakeStructure))
    monkeypatch.setitem(
        sys.modules,
        "pymatgen.symmetry.bandstructure",
        types.SimpleNamespace(HighSymmKpath=FakeHighSymmKpath),
    )

    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "POSCAR").write_text("fallback-poscar", encoding="utf-8")

    asyncio.run(worker._generate_band_structure_kpoints(work_dir, {"line_density": 7}))

    kpoints_text = (work_dir / "KPOINTS").read_text(encoding="utf-8")
    assert "Line-mode" in kpoints_text
    assert "! G" in kpoints_text
    assert "! X" in kpoints_text


def test_structure_optimization_accepts_contcar_url_as_structure_source(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    task_id = "opt-from-contcar-url"

    async def fail_if_old_cif_path(*args, **kwargs):
        raise AssertionError("legacy cif-only path should not be used")

    async def fake_materialize(work_dir, source_url, params, dest_name="POSCAR"):
        poscar_path = work_dir / dest_name
        poscar_path.write_text("poscar", encoding="utf-8")
        return str(poscar_path)

    async def fake_generate_vasp_inputs(work_dir, params):
        (work_dir / "POTCAR").write_text("potcar", encoding="utf-8")

    async def fake_run_vasp_calculation(work_dir, progress_callback=None):
        return {"success": True, "work_directory": str(work_dir)}

    async def fake_analyze_results(task_id_arg, work_dir, result):
        return {"success": True, "optimized_structure_download_url": str(work_dir / "CONTCAR")}

    monkeypatch.setattr(worker, "_get_cif_file", fail_if_old_cif_path)
    monkeypatch.setattr(worker, "_materialize_url_structure_source_to_poscar", fake_materialize)
    monkeypatch.setattr(worker, "_generate_vasp_inputs", fake_generate_vasp_inputs)
    monkeypatch.setattr(worker, "_run_vasp_calculation", fake_run_vasp_calculation)
    monkeypatch.setattr(worker, "_analyze_results", fake_analyze_results)

    result = asyncio.run(
        worker.run_structure_optimization(
            task_id,
            {"cif_url": "https://example.com/tenant/default/tasks/task-1/attempts/1/CONTCAR"},
        )
    )

    assert result["success"] is True
    assert (worker.base_work_dir / task_id / "POSCAR").read_text(encoding="utf-8") == "poscar"


def test_get_structure_for_scf_accepts_contcar_url_as_structure_source(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    class _Response:
        content = b"contcar-structure"

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(),
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "scf-contcar-url"
    work_dir.mkdir(parents=True, exist_ok=True)

    poscar_path = asyncio.run(
        worker._get_structure_for_scf(
            work_dir,
            {"cif_url": "https://example.com/tenant/default/tasks/task-1/attempts/1/CONTCAR"},
        )
    )

    assert Path(poscar_path).name == "POSCAR"
    assert Path(poscar_path).read_text(encoding="utf-8") == "contcar-structure"


def test_prepare_dos_files_rejects_empty_chgcar_from_scf_task(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    shared_root = tmp_path / "vasp-root"
    upstream_worker_dir = shared_root / "hpc-c2ln6"
    worker = VaspWorker(user_id="hpc-c2ln6-2", base_work_dir=str(shared_root))

    scf_task_id = "scf-empty-chgcar"
    scf_work_dir = upstream_worker_dir / scf_task_id
    scf_work_dir.mkdir(parents=True, exist_ok=True)
    for filename, content in {
        "POSCAR": "poscar",
        "POTCAR": "potcar",
        "CHG": "chg",
        "CHGCAR": "",
        "WAVECAR": "wavecar",
    }.items():
        (scf_work_dir / filename).write_text(content, encoding="utf-8")

    work_dir = worker.base_work_dir / "dos-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(Exception, match="未生成有效的 CHGCAR"):
        asyncio.run(worker._prepare_dos_files(work_dir, {"scf_task_id": scf_task_id}))


def test_modify_incar_for_dos_falls_back_to_gaussian_for_sparse_mesh(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    base_work_dir = tmp_path / "worker-root"
    user_id = "test-user"
    scf_task_id = "scf-task"
    scf_work_dir = base_work_dir / user_id / scf_task_id
    scf_work_dir.mkdir(parents=True, exist_ok=True)
    (scf_work_dir / "INCAR").write_text(
        "\n".join(
            [
                "SYSTEM = SCF",
                "PREC = Accurate",
                "ICHARG = 2",
                "NSW = 0",
                "IBRION = -1",
                "ISMEAR = 0",
                "SIGMA = 0.05",
                "LWAVE = .TRUE.",
                "LCHARG = .TRUE.",
            ]
        )
        + "\n",
        encoding="utf-8",
    )

    worker = VaspWorker(user_id=user_id, base_work_dir=str(base_work_dir))
    work_dir = worker.base_work_dir / "dos-sparse-kmesh"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "KPOINTS").write_text(
        "Automatically generated K-points\n0\nGamma\n1 1 1\n0 0 0\n",
        encoding="utf-8",
    )

    params = {"precision": "High", "scf_task_id": scf_task_id}
    asyncio.run(worker._modify_incar_for_dos(work_dir, params))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "ISMEAR = 0" in incar_text
    assert "ISMEAR = -5" not in incar_text
    assert "SIGMA = 0.05" in incar_text
    assert params["_dos_smearing_mode"] == "gaussian"


def test_analyze_scf_results_rejects_empty_chgcar(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "scf-empty-chgcar"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text("", encoding="utf-8")
    (work_dir / "POSCAR").write_text("poscar", encoding="utf-8")
    (work_dir / "CHGCAR").write_text("", encoding="utf-8")

    result = asyncio.run(
        worker._analyze_scf_results(
            work_dir,
            {"success": True, "work_directory": str(work_dir)},
        )
    )

    assert result["success"] is False
    assert result["missing_outputs"] == ["CHGCAR"]
    assert "未生成有效的 CHGCAR" in result["error"]


def test_analyze_scf_results_infers_convergence_from_oszicar(monkeypatch, tmp_path):
    from src.vasp_server import scf_analyzer
    from src.vasp_server.vasp_worker import VaspWorker

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "scf-gpu-tail"
    work_dir.mkdir(parents=True, exist_ok=True)
    (work_dir / "OUTCAR").write_text(
        "free  energy   TOTEN  =       -33.04916306 eV\n"
        "E-fermi :   6.4390\n",
        encoding="utf-8",
    )
    (work_dir / "INCAR").write_text("EDIFF = 1E-4\n", encoding="utf-8")
    (work_dir / "OSZICAR").write_text(
        "RMM:  10    -0.330491578606E+02   -0.12473E-03   -0.30313E-04    27\n"
        "RMM:  11    -0.330491630599E+02   -0.51993E-05   -0.55633E-05    21\n"
        "   1 F= -.33049163E+02 E0= -.33049163E+02  d E =-.123760E-13\n",
        encoding="utf-8",
    )
    (work_dir / "POSCAR").write_text("poscar", encoding="utf-8")
    (work_dir / "CHGCAR").write_text("charge density", encoding="utf-8")
    (work_dir / "AECCAR0").write_text("aeccar0", encoding="utf-8")
    (work_dir / "AECCAR2").write_text("aeccar2", encoding="utf-8")

    class FakeAnalyzer:
        def __init__(self, input_path, task_id=None):
            pass

        def analyze(self):
            return {"calculation_settings": {"EDIFF": "1E-4"}}

    def fake_report(input_path, task_id=None):
        report_path = work_dir / "scf_analysis_report.html"
        report_path.write_text("<html/>", encoding="utf-8")
        return str(report_path)

    monkeypatch.setattr(worker, "_run_bader_analysis", lambda *_args, **_kwargs: True)
    monkeypatch.setattr(scf_analyzer, "SCFAnalyzer", FakeAnalyzer)
    monkeypatch.setattr(scf_analyzer, "generate_scf_report", fake_report)

    result = asyncio.run(
        worker._analyze_scf_results(
            work_dir,
            {"success": True, "work_directory": str(work_dir)},
        )
    )

    assert result["success"] is True
    assert result["convergence"] is True
    assert result["electronic_steps"] == 2


def test_scf_analyzer_uses_oszicar_for_gpu_tail_convergence(tmp_path):
    from src.vasp_server.analyzers.scf import SCFAnalyzer

    work_dir = tmp_path / "scf-report"
    work_dir.mkdir()
    (work_dir / "OUTCAR").write_text(
        "free  energy   TOTEN  =       -33.04916306 eV\n",
        encoding="utf-8",
    )
    (work_dir / "OSZICAR").write_text(
        "RMM:  10    -0.330491578606E+02   -0.12473E-03   -0.30313E-04    27\n"
        "RMM:  11    -0.330491630599E+02   -0.51993E-05   -0.55633E-05    21\n"
        "   1 F= -.33049163E+02 E0= -.33049163E+02  d E =-.123760E-13\n",
        encoding="utf-8",
    )

    data = SCFAnalyzer(str(work_dir), task_id="scf-tail").analyze()

    convergence = data["electronic_convergence"]["analysis"]
    assert convergence["converged"] is True
    assert convergence["total_electronic_steps"] == 2
    assert data["final_results"]["converged"] is True


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
