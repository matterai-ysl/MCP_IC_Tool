import numpy as np
from pymatgen.core import Lattice, Structure


def test_md_rdf_csv_link_accepts_numpy_arrays():
    from src.vasp_server.analyzers.md import MDHTMLReportGenerator

    generator = MDHTMLReportGenerator({})

    result = generator._rdf_csv_link(
        np.array([1.0, 2.0, 3.0]),
        {"Li-O": {"g_r": [0.1, 0.2, 0.3]}},
    )

    assert "CSV" in result


def test_md_rdf_features_falls_back_to_numpy_trapezoid(monkeypatch, tmp_path):
    from src.vasp_server.analyzers import md as md_module

    analyzer = md_module.VASP_MDAnalyzer(str(tmp_path))
    structure = Structure(
        Lattice.cubic(10.0),
        ["Li", "O"],
        [[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]],
    )
    r_centers = np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5], dtype=float)
    g = np.array([0.4, 1.2, 1.9, 1.1, 0.8, 0.6, 0.9], dtype=float)

    monkeypatch.setattr(
        md_module.np,
        "convolve",
        lambda *_args, **_kwargs: np.array([0.6, 1.2, 1.7, 1.3, 0.9, 0.6, 0.8], dtype=float),
    )
    monkeypatch.delattr(md_module.np, "trapz", raising=False)

    calls: dict[str, np.ndarray] = {}

    def fake_trapezoid(y, x=None):
        calls["y"] = np.array(y, dtype=float)
        calls["x"] = np.array(x, dtype=float)
        return 3.25

    monkeypatch.setattr(md_module.np, "trapezoid", fake_trapezoid, raising=False)

    peaks, mins, cn = analyzer._analyze_rdf_features(r_centers, g, structure, "Li", "O")

    assert peaks
    assert mins
    assert "x" in calls
    assert cn == 3.25
