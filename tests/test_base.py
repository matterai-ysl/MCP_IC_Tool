import asyncio
from pathlib import Path


def _write_li_mn_o_poscar(target_dir: Path):
    (target_dir / "POSCAR").write_text(
        """Li Mn O
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 6.0
Li Mn O
8 4 12
Direct
0 0 0
""",
        encoding="utf-8",
    )


def test_cif_to_poscar_uses_pymatgen_without_external_cif2cell(monkeypatch, tmp_path):
    from src.vasp_server.base import cif_to_poscar

    cif_path = tmp_path / "structure.cif"
    cif_path.write_text(
        """data_Li2O
_symmetry_space_group_name_H-M   'P 1'
_cell_length_a   3.0
_cell_length_b   3.0
_cell_length_c   3.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_Int_Tables_number 1
loop_
 _symmetry_equiv_pos_as_xyz
  x,y,z
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 Li1 Li 0.0 0.0 0.0
 Li2 Li 0.5 0.5 0.5
 O1 O 0.25 0.25 0.25
""",
        encoding="utf-8",
    )

    def fail_if_subprocess_used(*args, **kwargs):
        raise AssertionError("cif_to_poscar should not invoke external cif2cell")

    monkeypatch.setattr("subprocess.run", fail_if_subprocess_used)

    poscar_path = Path(cif_to_poscar(str(cif_path), str(tmp_path)))

    assert poscar_path.exists()
    poscar_text = poscar_path.read_text(encoding="utf-8")
    assert "Li O" in poscar_text
    assert "2 1" in poscar_text


def test_cif_to_poscar_supports_cif_without_atom_site_label(tmp_path):
    from src.vasp_server.base import cif_to_poscar

    cif_path = tmp_path / "h2.cif"
    cif_path.write_text(
        """data_H2_iso
_cell_length_a 15.00000
_cell_length_b 15.00000
_cell_length_c 15.00000
_cell_angle_alpha 90.00000
_cell_angle_beta 90.00000
_cell_angle_gamma 90.00000
_symmetry_space_group_name_H-M 'P 1'
loop_
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 H 0.500000 0.500000 0.500000
 H 0.500000 0.500000 0.549333
""",
        encoding="utf-8",
    )

    poscar_path = Path(cif_to_poscar(str(cif_path), str(tmp_path)))

    assert poscar_path.exists()
    poscar_text = poscar_path.read_text(encoding="utf-8")
    assert "H" in poscar_text
    assert "2" in poscar_text


def test_cif_to_poscar_cleans_tiny_lattice_components(tmp_path):
    from src.vasp_server.base import cif_to_poscar

    cif_path = tmp_path / "hex.cif"
    cif_path.write_text(
        """data_hex
_symmetry_space_group_name_H-M 'P 1'
_cell_length_a 3.0
_cell_length_b 3.0
_cell_length_c 10.0
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 120
loop_
 _symmetry_equiv_pos_as_xyz
  x,y,z
loop_
 _atom_site_label
 _atom_site_type_symbol
 _atom_site_fract_x
 _atom_site_fract_y
 _atom_site_fract_z
 Si1 Si 0.0 0.0 0.0
""",
        encoding="utf-8",
    )

    poscar_path = Path(cif_to_poscar(str(cif_path), str(tmp_path)))
    lattice_lines = poscar_path.read_text(encoding="utf-8").splitlines()[2:5]

    first_vector = [float(value) for value in lattice_lines[0].split()]
    second_vector = [float(value) for value in lattice_lines[1].split()]

    assert first_vector[2] == 0.0
    assert second_vector[2] == 0.0


def test_sanitize_structure_for_vasp_cleans_export_noise_above_machine_epsilon():
    from pymatgen.core import Lattice, Structure

    from src.vasp_server.base import _sanitize_structure_for_vasp

    structure = Structure(
        Lattice(
            [
                [6.9190822385541466, 0.0, -0.0001414495007766],
                [0.0000139245219281, 5.1353799097168578, 3.7953479789363485],
                [0.0, 0.0, 8.0732772500000003],
            ]
        ),
        ["Eu", "P"],
        [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]],
    )

    sanitized = _sanitize_structure_for_vasp(structure)
    first_vector = sanitized.lattice.matrix[0]
    second_vector = sanitized.lattice.matrix[1]

    assert first_vector[2] == 0.0
    assert second_vector[0] == 0.0


def test_generate_kpoints_respects_explicit_target_product(tmp_path):
    from src.vasp_server.base import generate_kpoints

    _write_li_mn_o_poscar(tmp_path)

    assert generate_kpoints(str(tmp_path), target_product=15.0) is True

    kpoints_text = (tmp_path / "KPOINTS").read_text(encoding="utf-8")
    assert "  5   5   2" in kpoints_text


def test_generate_incar_uses_lighter_structure_optimization_defaults_and_auto_magmom(monkeypatch, tmp_path):
    from src.vasp_server.base import generate_incar

    _write_li_mn_o_poscar(tmp_path)

    monkeypatch.setattr(
        "src.vasp_server.base.load_u_values",
        lambda _path: {"Li": 0.0, "Mn": 4.0, "O": 0.0},
    )

    assert generate_incar(str(tmp_path)) is True

    incar_text = (tmp_path / "INCAR").read_text(encoding="utf-8")
    assert "PREC = Normal" in incar_text
    assert "ENCUT = 450" in incar_text
    assert "EDIFF = 1E-4" in incar_text
    assert "EDIFFG = -0.02" in incar_text
    assert "ISIF = 2" in incar_text
    assert "MAGMOM = 8*0.0 4*3.0 12*0.0" in incar_text
    assert "LDAUU = 0.0 4.0 0.0" in incar_text


def test_generate_vasp_inputs_passes_structure_optimization_kpoint_density(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    called = {}

    def fake_generate_kpoints(work_dir, target_product=None):
        called["work_dir"] = work_dir
        called["target_product"] = target_product
        (Path(work_dir) / "KPOINTS").write_text("", encoding="utf-8")
        return True

    def fake_generate_potcar(work_dir):
        (Path(work_dir) / "POTCAR").write_text("", encoding="utf-8")
        return True

    def fake_generate_incar(work_dir):
        (Path(work_dir) / "INCAR").write_text("SYSTEM = OPT\n", encoding="utf-8")
        return True

    monkeypatch.setattr("src.vasp_server.base.generate_kpoints", fake_generate_kpoints)
    monkeypatch.setattr("src.vasp_server.base.generate_potcar", fake_generate_potcar)
    monkeypatch.setattr("src.vasp_server.base.generate_incar", fake_generate_incar)

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    work_dir = worker.base_work_dir / "task-1"
    work_dir.mkdir(parents=True, exist_ok=True)

    asyncio.run(worker._generate_vasp_inputs(work_dir, {"kpoint_density": 15.0}))

    assert called["work_dir"] == str(work_dir)
    assert called["target_product"] == 15.0


def test_generate_single_point_dos_incar_uses_lighter_defaults(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "dos-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    _write_li_mn_o_poscar(work_dir)

    monkeypatch.setattr(
        "src.vasp_server.base.load_u_values",
        lambda _path: {"Li": 0.0, "Mn": 4.0, "O": 0.0},
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))

    asyncio.run(worker._generate_single_point_dos_incar(work_dir, {"precision": "High"}))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "PREC = High" in incar_text
    assert "EDIFF = 1E-5" in incar_text
    assert "NELMIN = 2" in incar_text
    assert "NELM = 120" in incar_text
    assert "NEDOS = 1500" in incar_text


def test_generate_single_point_dos_incar_falls_back_to_gaussian_for_sparse_mesh(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "dos-sparse-kmesh"
    work_dir.mkdir(parents=True, exist_ok=True)
    _write_li_mn_o_poscar(work_dir)
    (work_dir / "KPOINTS").write_text(
        "Automatically generated K-points\n0\nGamma\n3 1 1\n0 0 0\n",
        encoding="utf-8",
    )

    monkeypatch.setattr(
        "src.vasp_server.base.load_u_values",
        lambda _path: {"Li": 0.0, "Mn": 4.0, "O": 0.0},
    )

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))
    params = {"precision": "High"}

    asyncio.run(worker._generate_single_point_dos_incar(work_dir, params))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "ISMEAR = 0" in incar_text
    assert "ISMEAR = -5" not in incar_text
    assert "SIGMA = 0.05" in incar_text
    assert params["_dos_smearing_mode"] == "gaussian"


def test_modify_incar_for_band_structure_respects_lighter_precision(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker
    from src.vasp_server.base import generate_incar

    work_dir = tmp_path / "band-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    _write_li_mn_o_poscar(work_dir)

    monkeypatch.setattr(
        "src.vasp_server.base.load_u_values",
        lambda _path: {"Li": 0.0, "Mn": 4.0, "O": 0.0},
    )

    generate_incar(str(work_dir))
    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))

    asyncio.run(worker._modify_incar_for_band_structure(work_dir, {"precision": "High"}))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "SYSTEM = BAND" in incar_text
    assert "PREC = High" in incar_text


def test_modify_incar_for_band_structure_preserves_seed_precision_when_reusing_chgcar(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    base_work_dir = tmp_path / "worker-root"
    user_id = "test-user"
    worker = VaspWorker(user_id=user_id, base_work_dir=str(base_work_dir))
    work_dir = worker.base_work_dir / "band-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    seed_incar_path = work_dir / "seed_INCAR"
    seed_incar_path.write_text(
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

    asyncio.run(
        worker._modify_incar_for_band_structure(
            work_dir,
            {
                "precision": "High",
                "_band_input_mode": "bundle_with_chgcar",
                "_band_seed_incar_path": str(seed_incar_path),
            },
        )
    )

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "SYSTEM = BAND" in incar_text
    assert "PREC = Accurate" in incar_text
    assert "PREC = High" not in incar_text
    assert "ICHARG = 11" in incar_text


def test_modify_incar_for_dos_preserves_scf_precision_when_reusing_chgcar(tmp_path):
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
    work_dir = worker.base_work_dir / "dos-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    asyncio.run(
        worker._modify_incar_for_dos(
            work_dir,
            {"precision": "High", "scf_task_id": scf_task_id},
        )
    )

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert "SYSTEM = DOS" in incar_text
    assert "PREC = Accurate" in incar_text
    assert "PREC = High" not in incar_text
    assert "ICHARG = 11" in incar_text


def test_prepare_band_structure_files_resolves_upstream_from_sibling_worker_directory(tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    shared_root = tmp_path / "vasp-root"
    upstream_worker_dir = shared_root / "hpc-c2ln6"
    worker = VaspWorker(user_id="hpc-c2ln6-2", base_work_dir=str(shared_root))

    scf_task_id = "scf-task"
    scf_work_dir = upstream_worker_dir / scf_task_id
    scf_work_dir.mkdir(parents=True, exist_ok=True)
    for filename, content in {
        "POSCAR": "poscar",
        "POTCAR": "potcar",
        "CHG": "chg",
        "CHGCAR": "chgcar",
        "WAVECAR": "wavecar",
        "INCAR": "SYSTEM = SCF\nPREC = Accurate\n",
    }.items():
        (scf_work_dir / filename).write_text(content, encoding="utf-8")

    work_dir = worker.base_work_dir / "band-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    copied_files = asyncio.run(
        worker._prepare_band_structure_files(work_dir, {"scf_task_id": scf_task_id})
    )

    assert (work_dir / "POSCAR").exists()
    assert (work_dir / "CHGCAR").exists()
    assert copied_files["_band_input_mode"] == "scf_task"
    assert Path(copied_files["_band_seed_incar_path"]).exists()


def test_generate_neb_inputs_uses_lighter_defaults(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "neb-task"
    work_dir.mkdir(parents=True, exist_ok=True)
    for idx in range(7):
        image_dir = work_dir / f"{idx:02d}"
        image_dir.mkdir(parents=True, exist_ok=True)
        _write_li_mn_o_poscar(image_dir)

    called = {}

    def fake_generate_potcar(target):
        (Path(target) / "POTCAR").write_text("POTCAR", encoding="utf-8")
        return True

    def fake_generate_kpoints(target, target_product=None):
        called["target_product"] = target_product
        (Path(target) / "KPOINTS").write_text("KPOINTS", encoding="utf-8")
        return True

    monkeypatch.setattr("src.vasp_server.base.generate_potcar", fake_generate_potcar)
    monkeypatch.setattr("src.vasp_server.base.generate_kpoints", fake_generate_kpoints)

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))

    asyncio.run(worker._generate_neb_inputs(work_dir, {"n_images": 5, "kpoint_density": 15.0}))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert called["target_product"] == 15.0
    assert "NSW = 300" in incar_text
    assert "ENCUT = 450" in incar_text


def test_generate_phonon_inputs_uses_lighter_defaults(monkeypatch, tmp_path):
    from src.vasp_server.vasp_worker import VaspWorker

    work_dir = tmp_path / "phonon-task"
    work_dir.mkdir(parents=True, exist_ok=True)

    called = {}

    def fake_generate_kpoints(target, target_product=None):
        called["target_product"] = target_product
        (Path(target) / "KPOINTS").write_text("KPOINTS", encoding="utf-8")
        return True

    monkeypatch.setattr("src.vasp_server.base.generate_kpoints", fake_generate_kpoints)

    worker = VaspWorker(user_id="test-user", base_work_dir=str(tmp_path))

    asyncio.run(worker._generate_phonon_inputs(work_dir, {"kpoint_density": 15.0}))

    incar_text = (work_dir / "INCAR").read_text(encoding="utf-8")
    assert called["target_product"] == 15.0
    assert "EDIFF = 1E-6" in incar_text
    assert "ENCUT = 450" in incar_text
    assert "PREC = Normal" in incar_text
    assert "LREAL = Auto" in incar_text
