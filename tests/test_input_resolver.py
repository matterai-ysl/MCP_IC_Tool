import asyncio
import io
from pathlib import Path
import tarfile
import zipfile
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


def test_dos_rejects_empty_chgcar_from_manifest(tmp_path):
    resolver = UpstreamInputResolver()
    work_dir = tmp_path / "dos-empty-chgcar"
    artifacts = [
        {"artifact_type": "poscar", "local_path": str(_write_file(tmp_path / "src" / "POSCAR", "poscar"))},
        {"artifact_type": "potcar", "local_path": str(_write_file(tmp_path / "src" / "POTCAR", "potcar"))},
        {"artifact_type": "chgcar", "local_path": str(_write_file(tmp_path / "src" / "CHGCAR", ""))},
    ]

    with pytest.raises(ValueError, match="未生成有效的 CHGCAR"):
        resolver.materialize_inputs("dos_calculation", artifacts, work_dir)


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


def test_source_analysis_materializes_manifest_when_local_task_dir_is_gone(tmp_path):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    source_outcar = _write_file(tmp_path / "public" / "OUTCAR", "outcar-data")

    materialized = worker._resolve_source_work_directory(
        {
            "source_task_id": "old-source-task",
            "source_upstream_artifact_manifest": [
                {
                    "artifact_type": "outcar",
                    "filename": "OUTCAR",
                    "local_path": str(source_outcar),
                }
            ],
        },
        materialize_root=worker.base_work_dir / "analysis-task" / "_source",
    )

    assert materialized.name == "_source"
    assert (materialized / "OUTCAR").read_text(encoding="utf-8") == "outcar-data"


def test_band_structure_can_use_uploaded_bundle_with_chgcar(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-bundle-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    archive_buffer = io.BytesIO()
    with zipfile.ZipFile(archive_buffer, "w") as zf:
        zf.writestr("POSCAR", "bundle-poscar")
        zf.writestr("POTCAR", "bundle-potcar")
        zf.writestr("CHGCAR", "bundle-chgcar")
        zf.writestr("WAVECAR", "bundle-wavecar")
        zf.writestr("INCAR", "SYSTEM = SCF\nPREC = Accurate\nICHARG = 2\n")

    class _Response:
        def __init__(self, content: bytes):
            self.content = content

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(archive_buffer.getvalue()),
    )

    copied = asyncio.run(
        worker._prepare_band_structure_files(
            work_dir,
            {
                "input_url": "https://example.com/band-input.zip",
            },
        )
    )

    assert copied["_band_input_mode"] == "bundle_with_chgcar"
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "bundle-poscar"
    assert (work_dir / "CHGCAR").read_text(encoding="utf-8") == "bundle-chgcar"
    assert (work_dir / "WAVECAR").read_text(encoding="utf-8") == "bundle-wavecar"
    assert Path(copied["_band_seed_incar_path"]).read_text(encoding="utf-8").startswith("SYSTEM = SCF")


def test_band_structure_marks_uploaded_structure_bundle_for_seed_scf(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-seed-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    archive_buffer = io.BytesIO()
    with zipfile.ZipFile(archive_buffer, "w") as zf:
        zf.writestr(
            "POSCAR",
            "\n".join(
                [
                    "Li2O",
                    "1.0",
                    "4.0 0.0 0.0",
                    "0.0 4.0 0.0",
                    "0.0 0.0 4.0",
                    "Li O",
                    "2 1",
                    "Direct",
                    "0.0 0.0 0.0",
                    "0.5 0.5 0.5",
                    "0.25 0.25 0.25",
                ]
            ) + "\n",
        )

    class _Response:
        def __init__(self, content: bytes):
            self.content = content

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(archive_buffer.getvalue()),
    )

    copied = asyncio.run(
        worker._prepare_band_structure_files(
            work_dir,
            {
                "input_url": "https://example.com/structure-only.zip",
            },
        )
    )

    assert copied["_band_input_mode"] == "seed_scf"
    assert (work_dir / "POSCAR").exists()


def test_band_structure_can_use_multiple_input_urls_with_chgcar(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-multi-url-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    payloads = {
        "https://example.com/POSCAR": b"multi-poscar",
        "https://example.com/POTCAR": b"multi-potcar",
        "https://example.com/CHGCAR": b"multi-chgcar",
        "https://example.com/INCAR": b"SYSTEM = SCF\nPREC = Accurate\nICHARG = 2\n",
    }

    class _Response:
        def __init__(self, content: bytes):
            self.content = content

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(payloads[url]),
    )

    copied = asyncio.run(
        worker._prepare_band_structure_files(
            work_dir,
            {
                "input_url": [
                    "https://example.com/POSCAR",
                    "https://example.com/POTCAR",
                    "https://example.com/CHGCAR",
                    "https://example.com/INCAR",
                ],
            },
        )
    )

    assert copied["_band_input_mode"] == "bundle_with_chgcar"
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "multi-poscar"
    assert (work_dir / "CHGCAR").read_text(encoding="utf-8") == "multi-chgcar"
    assert Path(copied["_band_seed_incar_path"]).read_text(encoding="utf-8").startswith("SYSTEM = SCF")


def test_band_structure_can_use_multiple_input_urls_with_structure_only(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "band-multi-structure-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    cif_content = (
        "data_Li2O\n"
        "_cell_length_a 4.0\n"
        "_cell_length_b 4.0\n"
        "_cell_length_c 4.0\n"
        "_cell_angle_alpha 90\n"
        "_cell_angle_beta 90\n"
        "_cell_angle_gamma 90\n"
        "_symmetry_space_group_name_H-M 'P 1'\n"
        "loop_\n"
        "_atom_site_label\n"
        "_atom_site_type_symbol\n"
        "_atom_site_fract_x\n"
        "_atom_site_fract_y\n"
        "_atom_site_fract_z\n"
        "Li1 Li 0 0 0\n"
        "Li2 Li 0.5 0.5 0.5\n"
        "O1 O 0.25 0.25 0.25\n"
    )
    payloads = {
        "https://example.com/structure.cif": cif_content.encode("utf-8"),
        "https://example.com/POTCAR": b"multi-potcar",
    }

    class _Response:
        def __init__(self, content: bytes):
            self.content = content

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(payloads[url]),
    )

    async def fake_convert(cif_path, work_dir, params):
        poscar_path = Path(work_dir) / "POSCAR"
        poscar_path.write_text("converted-poscar", encoding="utf-8")
        return str(poscar_path)

    monkeypatch.setattr(worker, "_convert_cif_to_poscar", fake_convert)

    copied = asyncio.run(
        worker._prepare_band_structure_files(
            work_dir,
            {
                "input_url": [
                    "https://example.com/structure.cif",
                    "https://example.com/POTCAR",
                ],
            },
        )
    )

    assert copied["_band_input_mode"] == "seed_scf"
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "converted-poscar"
    assert (work_dir / "POTCAR").read_text(encoding="utf-8") == "multi-potcar"


def test_dos_can_use_multiple_input_urls_with_chgcar(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "dos-multi-url-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    payloads = {
        "https://example.com/POSCAR": b"dos-poscar",
        "https://example.com/POTCAR": b"dos-potcar",
        "https://example.com/CHGCAR": b"dos-chgcar",
        "https://example.com/INCAR": b"SYSTEM = SCF\nPREC = Accurate\nICHARG = 2\n",
    }

    class _Response:
        def __init__(self, content: bytes):
            self.content = content

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(payloads[url]),
    )

    copied = asyncio.run(
        worker._prepare_dos_files(
            work_dir,
            {
                "input_url": [
                    "https://example.com/POSCAR",
                    "https://example.com/POTCAR",
                    "https://example.com/CHGCAR",
                    "https://example.com/INCAR",
                ],
            },
        )
    )

    assert copied["_dos_input_mode"] == "bundle_with_chgcar"
    assert (work_dir / "POSCAR").read_text(encoding="utf-8") == "dos-poscar"
    assert (work_dir / "CHGCAR").read_text(encoding="utf-8") == "dos-chgcar"
    assert Path(copied["_dos_seed_incar_path"]).read_text(encoding="utf-8").startswith("SYSTEM = SCF")


def test_md_can_use_input_url_structure_file(tmp_path, monkeypatch):
    worker = VaspWorker(user_id="u1", base_work_dir=str(tmp_path / "workspace"))
    work_dir = worker.base_work_dir / "md-input-url-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    class _Response:
        content = (
            "data_Li2O\n"
            "_cell_length_a 4.0\n"
            "_cell_length_b 4.0\n"
            "_cell_length_c 4.0\n"
            "_cell_angle_alpha 90\n"
            "_cell_angle_beta 90\n"
            "_cell_angle_gamma 90\n"
            "_symmetry_space_group_name_H-M 'P 1'\n"
        ).encode("utf-8")

        def raise_for_status(self):
            return None

    monkeypatch.setattr(
        "src.vasp_server.vasp_worker.requests.get",
        lambda url, timeout=60: _Response(),
    )

    async def fake_convert(cif_path, work_dir, params):
        poscar_path = Path(work_dir) / "POSCAR"
        poscar_path.write_text(
            "\n".join(
                [
                    "Li2O",
                    "1.0",
                    "4.0 0.0 0.0",
                    "0.0 4.0 0.0",
                    "0.0 0.0 4.0",
                    "Li O",
                    "2 1",
                    "Direct",
                    "0.0 0.0 0.0",
                    "0.5 0.5 0.5",
                    "0.25 0.25 0.25",
                ]
            )
            + "\n",
            encoding="utf-8",
        )
        return str(poscar_path)

    monkeypatch.setattr(worker, "_convert_cif_to_poscar", fake_convert)
    async def fake_prepare_single_point_md_files(work_dir, params):
        return None

    monkeypatch.setattr(worker, "_prepare_single_point_md_files", fake_prepare_single_point_md_files)

    copied = asyncio.run(
        worker._prepare_md_files(
            work_dir,
            {
                "input_url": "https://example.com/structure.cif",
            },
        )
    )

    assert copied["POSCAR"].endswith("POSCAR")
    assert (work_dir / "POSCAR").exists()


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
