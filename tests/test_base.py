from pathlib import Path


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
