import asyncio
from pathlib import Path

from src.vasp_server.input_resolver import UpstreamInputResolver
from src.vasp_server.vasp_worker import VaspWorker


def _write_file(path: Path, content: str) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    return path


def test_scf_uses_upstream_artifact_manifest_instead_of_local_task_dir(tmp_path):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "scf-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    contcar = _write_file(tmp_path / "upstream" / "CONTCAR", "resolved-from-manifest")

    poscar_path = asyncio.run(
        worker._get_structure_for_scf(
            work_dir,
            {
                "optimized_task_id": "optimized-task",
                "upstream_artifact_manifest": [
                    {"artifact_type": "contcar", "local_path": str(contcar), "filename": "CONTCAR"},
                ],
            },
        )
    )

    assert poscar_path is not None
    assert Path(poscar_path).read_text(encoding="utf-8") == "resolved-from-manifest"


def test_dos_resolves_charge_and_wavefunction_from_manifest(tmp_path):
    resolver = UpstreamInputResolver()
    work_dir = tmp_path / "dos-task"
    artifacts = [
        {"artifact_type": "poscar", "local_path": str(_write_file(tmp_path / "src" / "POSCAR", "poscar"))},
        {"artifact_type": "potcar", "local_path": str(_write_file(tmp_path / "src" / "POTCAR", "potcar"))},
        {"artifact_type": "chgcar", "local_path": str(_write_file(tmp_path / "src" / "CHGCAR", "chgcar"))},
        {"artifact_type": "wavecar", "local_path": str(_write_file(tmp_path / "src" / "WAVECAR", "wavecar"))},
        {"artifact_type": "chg", "local_path": str(_write_file(tmp_path / "src" / "CHG", "chg"))},
    ]

    copied = resolver.materialize_inputs("dos_calculation", artifacts, work_dir)

    assert set(copied) >= {"POSCAR", "POTCAR", "CHGCAR", "WAVECAR"}
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "poscar"
    assert (work_dir / "CHGCAR").read_text(encoding="utf-8") == "chgcar"
    assert (work_dir / "WAVECAR").read_text(encoding="utf-8") == "wavecar"


def test_formula_fallback_for_scf_still_works(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "scf-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    cif_path = _write_file(work_dir / "structure.cif", "cif-data")
    poscar_path = work_dir / "POSCAR"

    async def fake_get_cif_file(_work_dir, _params, _progress_callback=None):
        return str(cif_path)

    async def fake_convert(_cif_path, _work_dir, _params):
        poscar_path.write_text("generated-from-cif", encoding="utf-8")
        return str(poscar_path)

    monkeypatch.setattr(worker, "_get_cif_file", fake_get_cif_file)
    monkeypatch.setattr(worker, "_convert_cif_to_poscar", fake_convert)

    resolved = asyncio.run(
        worker._get_structure_for_scf(
            work_dir,
            {"formula": "Li2O"},
        )
    )

    assert resolved == str(poscar_path)
    assert poscar_path.read_text(encoding="utf-8") == "generated-from-cif"
