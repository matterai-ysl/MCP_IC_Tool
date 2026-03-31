import asyncio
from pathlib import Path
import pytest

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


def test_band_structure_uses_upstream_artifact_manifest_instead_of_local_task_dir(tmp_path):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    source_dir = tmp_path / "upstream-band"
    poscar = _write_file(source_dir / "POSCAR", "band-poscar")
    potcar = _write_file(source_dir / "POTCAR", "band-potcar")
    chgcar = _write_file(source_dir / "CHGCAR", "band-chgcar")

    copied = asyncio.run(
        worker._prepare_band_structure_files(
            work_dir,
            {
                "scf_task_id": "scf-task",
                "upstream_artifact_manifest": [
                    {"artifact_type": "poscar", "local_path": str(poscar), "filename": "POSCAR"},
                    {"artifact_type": "potcar", "local_path": str(potcar), "filename": "POTCAR"},
                    {"artifact_type": "chgcar", "local_path": str(chgcar), "filename": "CHGCAR"},
                ],
            },
        )
    )

    assert set(copied) >= {"POSCAR", "POTCAR", "CHGCAR"}
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "band-poscar"
    assert (work_dir / "POTCAR").read_text(encoding="utf-8") == "band-potcar"
    assert (work_dir / "CHGCAR").read_text(encoding="utf-8") == "band-chgcar"


def test_scf_can_download_upstream_artifact_from_download_url(tmp_path, monkeypatch):
    resolver = UpstreamInputResolver()
    work_dir = tmp_path / "scf-task"

    class _Response:
        content = b"downloaded-contcar"

        def raise_for_status(self):
            return None

    def fake_get(url, timeout=60):
        assert url == "https://artifacts.example.com/contcar"
        return _Response()

    monkeypatch.setattr("src.vasp_server.input_resolver.requests.get", fake_get)

    copied = resolver.resolve_for_scf(
        [
            {
                "artifact_type": "contcar",
                "download_url": "https://artifacts.example.com/contcar",
                "filename": "CONTCAR",
            }
        ],
        work_dir,
    )

    assert copied["POSCAR"].endswith("POSCAR")
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "downloaded-contcar"


def test_formula_input_is_no_longer_supported_for_scf(tmp_path):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "scf-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    with pytest.raises(Exception, match="cif_url 或 optimized_task_id"):
        asyncio.run(
            worker._get_structure_for_scf(
                work_dir,
                {"formula": "Li2O"},
            )
        )
