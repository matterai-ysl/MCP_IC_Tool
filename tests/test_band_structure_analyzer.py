import types


def test_band_structure_analyzer_uses_low_memory_vasprun_parse(monkeypatch, tmp_path):
    from src.vasp_server.analyzers.band_structure import BandStructureAnalyzer

    (tmp_path / "OUTCAR").write_text("OUTCAR", encoding="utf-8")
    (tmp_path / "vasprun.xml").write_text("<xml/>", encoding="utf-8")

    captured = {}

    class FakeBandStructure:
        nb_bands = 2
        bands = {"up": []}
        kpoints = [type("KPoint", (), {"label": "G", "frac_coords": (0.0, 0.0, 0.0)})()]
        branches = []
        distance = []
        structure = type("Structure", (), {"sites": []})()

        def get_band_gap(self):
            return {"energy": 1.5, "direct": True, "transition": "G-G"}

        def get_direct_band_gap_dict(self):
            return {"up": {"value": 1.5, "kpoint_index": 0, "band_indices": []}}

        def get_vbm(self):
            return {
                "energy": 0.0,
                "kpoint": None,
                "kpoint_index": [],
                "band_index": {},
                "projections": {},
            }

        def get_cbm(self):
            return {
                "energy": 1.5,
                "kpoint": None,
                "kpoint_index": [],
                "band_index": {},
                "projections": {},
            }

        def is_metal(self):
            return False

    class FakeVasprun:
        def __init__(self, path, parse_projected_eigen=False):
            captured["path"] = path
            captured["parse_projected_eigen"] = parse_projected_eigen
            self.efermi = 0.7
            self.final_energy = -12.3

        def get_band_structure(self, line_mode=True):
            captured["line_mode"] = line_mode
            return FakeBandStructure()

    import sys

    pymatgen_module = types.ModuleType("pymatgen")
    io_module = types.ModuleType("pymatgen.io")
    vasp_module = types.ModuleType("pymatgen.io.vasp")
    vasp_module.Vasprun = FakeVasprun
    monkeypatch.setitem(sys.modules, "pymatgen", pymatgen_module)
    monkeypatch.setitem(sys.modules, "pymatgen.io", io_module)
    monkeypatch.setitem(sys.modules, "pymatgen.io.vasp", vasp_module)

    analyzer = BandStructureAnalyzer(str(tmp_path), task_id="band-low-mem")
    analyzer._analyze_band_structure()

    assert captured["parse_projected_eigen"] is False
    assert captured["line_mode"] is True
    assert analyzer.data["final_results"]["band_gap"] == 1.5
