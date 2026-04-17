from src.vasp_server.analyzers.band_structure import (
    BandStructureAnalyzer,
    BandStructureHTMLGenerator,
    _display_branch_name,
    _display_transition,
    estimate_effective_mass_from_parabola,
    export_band_structure_assets,
    summarize_projection_contributions,
)
from src.vasp_server.analyzers.dos import PyMatGenDOSHTMLGenerator
from src.vasp_server.analyzers.md import MDHTMLReportGenerator, VASP_MDAnalyzer
from src.vasp_server.analyzers.neb import NEBHTMLGenerator
from src.vasp_server.analyzers.optimization import OptimizationHTMLGenerator, OUTCARAnalyzer
from src.vasp_server.analyzers.phonon import PhononHTMLGenerator
from src.vasp_server.analyzers.scf import SCFHTMLGenerator


def test_band_structure_report_embeds_png_and_raw_data_download_links():
    generator = BandStructureHTMLGenerator(
        {
            "task_info": {"task_id": "band-1", "composition": "Li2O"},
            "file_info": {},
            "band_structure": {
                "band_gap": 1.2,
                "direct_gap": 1.35,
                "direct_gap_kpoint_label": "Γ",
                "is_direct": True,
                "is_metal": False,
                "nb_bands": 2,
                "transition": "Γ-Γ",
                "nkpoints": 16,
                "n_branches": 3,
                "spin_channels": ["up"],
                "vbm_kpoint_label": "Γ",
                "cbm_kpoint_label": "Γ",
                "vbm_band_degeneracy": 2,
                "cbm_band_degeneracy": 1,
                "vbm_kpoint_degeneracy": 1,
                "cbm_kpoint_degeneracy": 1,
                "fermi_offset_from_midgap": 0.05,
                "vbm_band_widths": [{"band": "up:12", "width": 3.1}],
                "cbm_band_widths": [{"band": "up:13", "width": 4.6}],
                "effective_masses": {
                    "electron": {"m_over_me": 0.42, "branch": "Γ-X", "fit_points": 5},
                    "hole": {"m_over_me": 1.37, "branch": "Γ-X", "fit_points": 5},
                },
                "vbm_character": {
                    "top_species_orbitals": [
                        {"label": "O-p", "weight": 0.62},
                        {"label": "Li-s", "weight": 0.18},
                    ],
                    "top_elements": [
                        {"label": "O", "weight": 0.68},
                        {"label": "Li", "weight": 0.22},
                    ],
                },
                "cbm_character": {
                    "top_species_orbitals": [
                        {"label": "Li-s", "weight": 0.57},
                        {"label": "O-p", "weight": 0.21},
                    ],
                    "top_elements": [
                        {"label": "Li", "weight": 0.61},
                        {"label": "O", "weight": 0.24},
                    ],
                },
            },
            "final_results": {"vbm": -0.1, "cbm": 1.1, "fermi_energy": 0.5, "total_energy": -10.0},
            "visualizations": {"band_structure_plot": "ZmFrZS1wbmc="},
            "raw_band_data": {
                "distance": [0.0, 0.5],
                "bands": {"up": [[-1.0, -0.5], [0.5, 1.0]]},
            },
        }
    )

    html = generator._generate_body_sections()

    assert 'download="band_structure.png"' in html
    assert 'download="band_structure_data.csv"' in html
    assert "带隙跃迁" in html
    assert "路径 k 点数" in html
    assert "位于 Γ 点" in html
    assert "直接跃迁带隙" in html
    assert "费米能级偏离带隙中心" in html
    assert "带边特征" in html
    assert "VBM 主导轨道" in html
    assert "CBM 主导轨道" in html
    assert "有效质量估计" in html
    assert "电子有效质量" in html
    assert "0.42 m" in html
    assert "科研解读" in html
    assert "O-p (62.0%)" in html
    assert "Li-s (57.0%)" in html


def test_band_structure_report_includes_how_to_read_chart_guidance():
    generator = BandStructureHTMLGenerator(
        {
            "task_info": {"task_id": "band-guide", "composition": "Li2O"},
            "file_info": {},
            "band_structure": {
                "band_gap": 1.2,
                "direct_gap": 1.35,
                "is_direct": True,
                "is_metal": False,
                "nb_bands": 2,
                "transition": "Γ-Γ",
            },
            "final_results": {"vbm": -0.1, "cbm": 1.1, "fermi_energy": 0.0, "total_energy": -10.0},
            "visualizations": {"band_structure_plot": "ZmFrZS1wbmc="},
            "raw_band_data": {
                "distance": [0.0, 0.5],
                "bands": {"up": [[-1.0, -0.5], [0.5, 1.0]]},
            },
        }
    )

    html = generator._generate_body_sections()

    assert "怎么看这张图" in html
    assert "先看带隙两侧最靠近 0 eV 的带边" in html
    assert "如果价带顶和导带底出现在同一个 k 点" in html


def test_band_structure_export_writes_preview_png_and_csv(tmp_path):
    export_band_structure_assets(
        {
            "visualizations": {"band_structure_plot": "ZmFrZS1wbmc="},
            "raw_band_data": {
                "distance": [0.0, 0.5],
                "bands": {"up": [[-1.0, -0.5], [0.5, 1.0]]},
            },
        },
        tmp_path,
    )

    assert (tmp_path / "band_structure.png").exists()
    assert (tmp_path / "band_structure_data.csv").exists()


def test_summarize_projection_contributions_aggregates_species_and_orbitals():
    projections = {
        "up": {
            "p": [0.6, 0.1, 0.2],
            "d": [0.1, 0.4, 0.0],
            "s": [0.0, 0.1, 0.2],
        }
    }
    result = summarize_projection_contributions(projections, ["O", "Mn", "O"])

    assert result["top_species_orbitals"][0]["label"] == "O-p"
    assert round(result["top_species_orbitals"][0]["weight"], 3) == 0.471
    assert result["top_elements"][0]["label"] == "O"
    assert round(result["top_elements"][0]["weight"], 3) == 0.647


def test_estimate_effective_mass_from_parabola_handles_electron_and_hole_edges():
    k = [-0.1, 0.0, 0.1]
    electron_energies = [0.0380998212, 0.0, 0.0380998212]
    hole_energies = [-0.0380998212, 0.0, -0.0380998212]

    electron = estimate_effective_mass_from_parabola(k, electron_energies, carrier="electron")
    hole = estimate_effective_mass_from_parabola(k, hole_energies, carrier="hole")

    assert electron is not None
    assert hole is not None
    assert abs(electron["m_over_me"] - 1.0) < 0.05
    assert abs(hole["m_over_me"] - 1.0) < 0.05


def test_display_branch_name_falls_back_when_pymatgen_branch_name_is_none_none():
    class FakeKpoint:
        def __init__(self, label, frac_coords):
            self.label = label
            self.frac_coords = frac_coords

    class FakeBS:
        kpoints = [
            FakeKpoint(None, [0.0, 0.0, 0.0]),
            FakeKpoint(None, [0.5, 0.0, 0.0]),
        ]

    branch = {"name": "None-None", "start_index": 0, "end_index": 1}
    display = _display_branch_name(FakeBS(), branch)

    assert display == "(0.000, 0.000, 0.000) -> (0.500, 0.000, 0.000)"


def test_display_transition_prefers_single_kpoint_text_for_direct_gap():
    text = _display_transition(
        {
            "is_direct": True,
            "transition": "(0.000,0.000,0.000)-(0.000,0.000,0.000)",
            "direct_gap_kpoint_label": "(0.000, 0.000, 0.000)",
        }
    )

    assert text == "位于 (0.000, 0.000, 0.000) 点"


def test_band_structure_report_uses_friendly_projection_placeholder_when_missing():
    generator = BandStructureHTMLGenerator(
        {
            "task_info": {"task_id": "band-2"},
            "file_info": {},
            "band_structure": {
                "band_gap": 2.1,
                "direct_gap": 2.1,
                "is_direct": True,
                "is_metal": False,
                "nb_bands": 10,
                "transition": "(0.000,0.000,0.000)-(0.000,0.000,0.000)",
                "direct_gap_kpoint_label": "(0.000, 0.000, 0.000)",
                "effective_masses": {
                    "electron": {"m_over_me": 1.2, "branch": "None-None", "fit_points": 3, "band": "up:7"},
                    "hole": {"m_over_me": 0.8, "branch": "None-None", "fit_points": 3, "band": "up:6"},
                },
                "vbm_character": {"top_species_orbitals": [], "top_elements": []},
                "cbm_character": {"top_species_orbitals": [], "top_elements": []},
            },
            "final_results": {"vbm": 0.0, "cbm": 2.1, "fermi_energy": 1.0, "total_energy": -1.0},
            "visualizations": {},
            "raw_band_data": {},
        }
    )

    html = generator._generate_body_sections()

    assert "未提供投影数据" in html
    assert "未标注路径（带边附近）" in html
    assert "None-None" not in html


def test_band_structure_analyzer_generate_plots_accepts_axes_like_object(monkeypatch, tmp_path):
    (tmp_path / "OUTCAR").write_text("", encoding="utf-8")
    analyzer = BandStructureAnalyzer(str(tmp_path), task_id="band-plot")
    analyzer._bs = object()

    class FakeFigure:
        def savefig(self, buf, format=None, dpi=None, bbox_inches=None):
            buf.write(b"fake-png")

    class FakeAxes:
        def __init__(self):
            self.figure = FakeFigure()

    class FakePlotter:
        def __init__(self, _bs):
            pass

        def get_plot(self, ylim=None):
            return FakeAxes()

    import sys
    import types

    matplotlib_module = types.ModuleType("matplotlib")
    pyplot_module = types.ModuleType("matplotlib.pyplot")
    matplotlib_module.use = lambda _backend: None
    pyplot_module.close = lambda _fig: None
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib_module)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot_module)

    plotter_module = types.ModuleType("pymatgen.electronic_structure.plotter")
    plotter_module.BSPlotter = FakePlotter
    monkeypatch.setitem(sys.modules, "pymatgen.electronic_structure.plotter", plotter_module)

    analyzer._generate_plots()

    assert analyzer.data["visualizations"]["band_structure_plot"] == "ZmFrZS1wbmc="


def test_band_structure_analyzer_draws_explicit_red_fermi_level(monkeypatch, tmp_path):
    (tmp_path / "OUTCAR").write_text("", encoding="utf-8")
    analyzer = BandStructureAnalyzer(str(tmp_path), task_id="band-fermi-line")
    analyzer._bs = object()

    captured = {}

    class FakeFigure:
        def savefig(self, buf, format=None, dpi=None, bbox_inches=None):
            buf.write(b"fake-png")

    class FakeAxes:
        def __init__(self):
            self.figure = FakeFigure()

        def axhline(self, y, color=None, linestyle=None, linewidth=None, alpha=None, zorder=None):
            captured["line"] = {
                "y": y,
                "color": color,
                "linestyle": linestyle,
                "linewidth": linewidth,
                "alpha": alpha,
                "zorder": zorder,
            }

    class FakePlotter:
        def __init__(self, _bs):
            pass

        def get_plot(self, ylim=None):
            return FakeAxes()

    import sys
    import types

    matplotlib_module = types.ModuleType("matplotlib")
    pyplot_module = types.ModuleType("matplotlib.pyplot")
    matplotlib_module.use = lambda _backend: None
    pyplot_module.close = lambda _fig: None
    monkeypatch.setitem(sys.modules, "matplotlib", matplotlib_module)
    monkeypatch.setitem(sys.modules, "matplotlib.pyplot", pyplot_module)

    plotter_module = types.ModuleType("pymatgen.electronic_structure.plotter")
    plotter_module.BSPlotter = FakePlotter
    monkeypatch.setitem(sys.modules, "pymatgen.electronic_structure.plotter", plotter_module)

    analyzer._generate_plots()

    assert captured["line"]["y"] == 0
    assert captured["line"]["color"] == "red"
    assert captured["line"]["linestyle"] == "--"


def test_neb_report_embeds_png_download_link_alongside_csv():
    generator = NEBHTMLGenerator(
        {
            "task_info": {"task_id": "neb-1", "composition": "Li2O"},
            "file_info": {},
            "neb": {"relative_energies": [0.0, 0.2, 0.1]},
            "visualizations": {"neb_profile": "ZmFrZS1wbmc="},
        }
    )

    html = generator._generate_body_sections()

    assert 'download="neb_profile.png"' in html
    assert 'download="neb_energies.csv"' in html


def test_dos_report_includes_plain_language_interpretation():
    generator = PyMatGenDOSHTMLGenerator(
        {
            "summary": {
                "formula": "MgO",
                "space_group": "Fm-3m",
                "band_gap": 4.8,
                "material_type": "insulator",
                "is_magnetic": False,
                "magnetic_type": "non-magnetic",
            },
            "structure": {
                "formula": "MgO",
                "reduced_formula": "MgO",
                "num_sites": 2,
                "density": 3.5,
                "space_group": "Fm-3m",
                "point_group": "m-3m",
                "elements": ["Mg", "O"],
                "lattice_parameters": {
                    "a": 4.2,
                    "b": 4.2,
                    "c": 4.2,
                    "alpha": 90.0,
                    "beta": 90.0,
                    "gamma": 90.0,
                    "volume": 74.0,
                },
            },
            "dos_analysis": {
                "fermi_energy": 0.0,
                "band_gap": 4.8,
                "cbm_energy": 4.7,
                "vbm_energy": -0.1,
                "material_type": "insulator",
                "is_metal": False,
                "is_spin_polarized": False,
                "orbital_analysis": {
                    "O": {"p": 0.7, "s": 0.1, "d": 0.0, "f": 0.0},
                    "Mg": {"s": 0.5, "p": 0.1, "d": 0.0, "f": 0.0},
                },
            },
            "band_structure": {
                "has_dos_gap_analysis": True,
                "fundamental_gap": 4.8,
                "fundamental_type": "直接带隙",
                "direct_gap": 4.8,
                "is_metal": False,
                "num_bands": 32,
                "fermi_level": 0.0,
            },
            "chemical_properties": {},
            "magnetic_properties": {"is_magnetic": False},
            "visualizations": {},
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "怎么理解这份 DOS 结果" in html
    assert "不等于器件最终性能" in html


def test_dos_visualization_section_includes_per_chart_reading_guidance():
    generator = PyMatGenDOSHTMLGenerator(
        {
            "summary": {"formula": "MgO"},
            "structure": {"formula": "MgO", "elements": ["Mg", "O"]},
            "dos_analysis": {
                "is_spin_polarized": False,
                "orbital_analysis": {
                    "O": {"p": 0.7, "s": 0.1, "d": 0.0, "f": 0.0},
                    "Mg": {"s": 0.5, "p": 0.1, "d": 0.0, "f": 0.0},
                },
            },
            "visualizations": {
                "total_dos_plot": "ZmFrZS1wbmc=",
                "element_dos_plot": "ZmFrZS1wbmc=",
                "spd_dos_plot": "ZmFrZS1wbmc=",
                "band_structure_plot": "ZmFrZS1wbmc=",
            },
            "raw_dos_data": {
                "TDOS": {
                    "energy": [-1.0, 0.0, 1.0],
                    "dos_up": [0.2, 0.0, 0.3],
                },
                "O": {
                    "energy": [-1.0, 0.0, 1.0],
                    "dos_up": [0.1, 0.0, 0.2],
                },
                "Mg": {
                    "energy": [-1.0, 0.0, 1.0],
                    "dos_up": [0.1, 0.0, 0.1],
                },
            },
        }
    )

    html = generator._generate_body_sections()

    assert "怎么看这张图" in html
    assert "先看 0 eV（费米能级）附近" in html
    assert "看哪种元素在带边附近贡献更高" in html
    assert "看 s、p、d、f 哪类轨道主导带边" in html
    assert "如果这里显示的是能带信息" in html


def test_scf_report_includes_plain_language_interpretation():
    generator = SCFHTMLGenerator(
        {
            "electronic_convergence": {
                "analysis": {
                    "energy_convergence": {
                        "converged": True,
                        "ediff_threshold": 1e-5,
                        "final_energy_change": 5e-6,
                        "final_energy": -10.0,
                    },
                    "total_electronic_steps": 12,
                },
                "electronic_steps": [],
                "scf_energies": [-10.0, -10.1, -10.10001],
            },
            "electronic_structure": {
                "fermi_energy": 0.4,
                "band_gap": 1.8,
                "total_electrons": 24,
                "irreducible_kpoints": 12,
            },
            "forces_and_stress": {"forces": [], "stress_tensor": None},
            "magnetic_properties": {"is_magnetic": False},
            "bader_analysis": {},
            "elf_analysis": {},
            "final_results": {
                "converged": True,
                "total_electronic_steps": 12,
                "final_energy": -10.10001,
                "fermi_energy": 0.4,
                "band_gap": 1.8,
                "total_magnetization": 0.0,
            },
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "这一步主要是在把电子状态算稳" in html
    assert "可以作为后续 DOS、能带" in html


def test_scf_sections_include_how_to_read_chart_guidance():
    generator = SCFHTMLGenerator(
        {
            "electronic_convergence": {
                "analysis": {
                    "energy_convergence": {
                        "converged": True,
                        "ediff_threshold": 1e-5,
                        "final_energy_change": 5e-6,
                        "final_energy": -10.0,
                    },
                    "total_electronic_steps": 12,
                },
                "electronic_steps": [{"step": 1, "free_energy": -10.0, "energy_change": -0.1}],
                "scf_energies": [-10.0, -10.1, -10.10001],
            },
            "electronic_structure": {
                "fermi_energy": 0.4,
                "band_gap": 1.8,
                "total_electrons": 24,
                "irreducible_kpoints": 12,
            },
            "forces_and_stress": {
                "forces": [{"magnitude": 0.12, "fx": 0.1, "fy": 0.0, "fz": 0.0}],
                "stress_tensor": None,
                "max_force": 0.12,
                "rms_force": 0.12,
            },
            "magnetic_properties": {
                "is_spin_polarized": True,
                "total_magnetization": 1.2,
                "atom_magnetizations": [{"atom_index": 0, "total_magnitude": 1.2}],
            },
            "bader_analysis": {
                "available": True,
                "summary": {
                    "total_atoms": 1,
                    "total_charge_transfer": 0.2,
                    "avg_charge_transfer": 0.2,
                    "max_charge_transfer": 0.2,
                    "min_charge_transfer": 0.2,
                    "avg_bader_charge": 6.8,
                    "element_summary": {
                        "O": {
                            "count": 1,
                            "avg_bader_charge": 6.8,
                            "avg_charge_transfer": -0.2,
                            "total_charge_transfer": -0.2,
                        }
                    },
                },
                "atom_charges": [
                    {
                        "atom_index": 1,
                        "element": "O",
                        "bader_charge": 6.8,
                        "charge_transfer": -0.2,
                        "position": [0.0, 0.0, 0.0],
                        "zval": 6.0,
                    }
                ],
            },
            "elf_analysis": {},
            "final_results": {
                "converged": True,
                "total_electronic_steps": 12,
                "final_energy": -10.10001,
                "fermi_energy": 0.4,
                "band_gap": 1.8,
                "total_magnetization": 1.2,
            },
        }
    )

    convergence_html = generator._generate_electronic_convergence_section()
    forces_html = generator._generate_forces_stress_section()
    magnetic_html = generator._generate_magnetic_section()
    bader_html = generator._generate_bader_section()

    assert "看柱子是否在后半段逐渐贴近 0" in convergence_html
    assert "最后几步几乎重合" in convergence_html
    assert "少数原子受力特别大" in forces_html
    assert "磁性主要由这些原子贡献" in magnetic_html
    assert "Bader 电荷是否自然分组" in bader_html
    assert "失去电子还是获得电子" in bader_html


def test_md_report_includes_plain_language_interpretation():
    generator = MDHTMLReportGenerator(
        {
            "task_info": {"task_id": "md-1", "output_dir": "."},
            "diffusion": {
                "diffusion_by_element": {"Li": {"D_m2_per_s": 2.5e-10}},
                "ionic_conductivity_S_per_m": {"Li": 1.2},
                "arrhenius_single": {"T_list_K": [900.0]},
            },
            "rdf": {
                "r_A": [0.1, 0.2, 0.3],
                "g_pairs": {
                    "Li-O": {
                        "coordination_number": 4.2,
                        "peaks": [{"r": 2.01}],
                        "method": "rdf",
                        "g_r": [0.0, 0.5, 1.0],
                    }
                },
                "evolution": {"first_peak_r_series_A": [2.0, 2.01, 2.02]},
            },
            "stability": {
                "lattice": {"a": [5.1, 5.1], "b": [5.1, 5.1], "c": [5.1, 5.1], "vol": [120.0, 120.2]},
                "times_s": [0.0, 1e-15],
                "temperature_K": [890.0, 905.0],
                "energy_eV": [-100.0, -99.9],
                "pressure_kB": [1.0, 2.0],
                "density_kg_per_m3": [2500.0, 2498.0],
            },
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "扩散系数越大" in html
    assert "RDF" in html


def test_md_analyzer_embeds_existing_output_images_when_generating_html(tmp_path):
    work_dir = tmp_path / "md-run"
    output_dir = work_dir / "MD_output"
    output_dir.mkdir(parents=True)
    (output_dir / "msd.png").write_bytes(b"fake-png")

    analyzer = VASP_MDAnalyzer(str(work_dir), task_id="md-inline", output_dir=str(output_dir))
    analyzer.analysis_data.update(
            {
                "diffusion": {
                    "diffusion_by_element": {"Li": {"D_m2_per_s": 1.0e-10}},
                    "ionic_conductivity_S_per_m": {"Li": 1.0},
                    "arrhenius_single": {"T_list_K": [600.0]},
                    "msd_by_element": {"Li": {"time_s": [0.0], "msd": [0.0]}},
                },
                "rdf": {"r_A": [0.1], "g_pairs": {}, "evolution": {}},
                "stability": {"times_s": []},
            }
        )

    html_path = analyzer._generate_html_report()
    html = (output_dir / "md_analysis_report.html").read_text(encoding="utf-8")

    assert html_path.endswith("md_analysis_report.html")
    assert "data:image/png;base64," in html
    assert 'download="msd.png"' in html


def test_md_report_uses_single_temperature_arrhenius_explanation_instead_of_empty_placeholder():
    generator = MDHTMLReportGenerator(
        {
            "task_info": {"task_id": "md-arrh", "output_dir": "."},
            "diffusion": {
                "diffusion_by_element": {"Li": {"D_m2_per_s": 2.5e-10}},
                "ionic_conductivity_S_per_m": {"Li": 1.2},
                "arrhenius_single": {"T_list_K": [900.0]},
                "msd_by_element": {"Li": {"time_s": [0.0], "msd": [0.0]}},
            },
            "rdf": {"r_A": [0.1], "g_pairs": {}, "evolution": {}},
            "stability": {"times_s": []},
        }
    )

    html = generator._section_diffusion({})

    assert "Arrhenius（多温聚合时显示）" not in html
    assert "当前是单温度 MD" not in html
    assert "MSD（按元素）" in html


def test_md_report_skips_empty_rdf_pair_cards_when_pair_plot_missing():
    generator = MDHTMLReportGenerator(
        {
            "task_info": {"task_id": "md-rdf", "output_dir": "."},
            "rdf": {
                "r_A": [0.1, 0.2],
                "g_pairs": {
                    "Li-O": {
                        "coordination_number": 4.2,
                        "peaks": [{"r": 2.01}],
                        "method": "rdf",
                        "g_r": [0.0, 0.5],
                    }
                },
                "evolution": {},
            },
        }
    )

    html = generator._generate_rdf_pair_images({}, generator.data["rdf"]["g_pairs"])

    assert "当前没有可展示的元素对 RDF 图像" in html
    assert '<div class="note">无图像</div>' not in html


def test_md_sections_include_how_to_read_chart_guidance():
    generator = MDHTMLReportGenerator(
        {
            "task_info": {"task_id": "md-guide", "output_dir": "."},
            "diffusion": {
                "diffusion_by_element": {"Li": {"D_m2_per_s": 2.5e-10}},
                "ionic_conductivity_S_per_m": {"Li": 1.2},
                "arrhenius_single": {"T_list_K": [700.0, 900.0], "D_list_m2_per_s": [1e-10, 2e-10]},
                "msd_by_element": {"Li": {"time_s": [0.0, 1.0], "msd": [0.0, 0.5]}},
            },
            "rdf": {
                "r_A": [0.1, 0.2],
                "g_pairs": {
                    "Li-O": {
                        "coordination_number": 4.2,
                        "peaks": [{"r": 2.01}],
                        "method": "rdf",
                        "g_r": [0.0, 0.5],
                    }
                },
                "evolution": {},
            },
            "stability": {
                "times_s": [0.0, 1.0],
                "energy_eV": [-10.0, -10.1],
                "temperature_K": [700.0, 702.0],
                "pressure_kB": [1.0, 1.2],
                "lattice": {"a": [4.1, 4.1], "b": [4.1, 4.1], "c": [4.1, 4.1], "vol": [68.9, 68.9]},
                "density_kg_per_m3": [2100.0, 2099.0],
            },
        }
    )

    diffusion_html = generator._section_diffusion({"msd": "ZmFrZS1wbmc=", "arrhenius": "ZmFrZS1wbmc="})
    rdf_html = generator._section_rdf({"rdf_all": "ZmFrZS1wbmc=", "rdf_first_peak_evolution": "ZmFrZS1wbmc=", "rdf_Li-O": "ZmFrZS1wbmc="})
    stability_html = generator._section_stability({"stability": "ZmFrZS1wbmc=", "lattice_density": "ZmFrZS1wbmc="})

    assert "MSD 曲线抬升越快" in diffusion_html
    assert "Arrhenius 图里温度升高后扩散系数如果整体上升" in diffusion_html
    assert "RDF 的峰越尖锐、位置越稳定" in rdf_html
    assert "如果能量、温度、压力只是在小范围波动" in stability_html


def test_phonon_report_includes_plain_language_interpretation():
    generator = PhononHTMLGenerator(
        {
            "phonon": {
                "n_modes": 12,
                "n_imaginary": 0,
                "dynamically_stable": True,
                "max_imaginary_freq_cm1": 0.0,
                "max_real_freq_cm1": 500.0,
                "min_real_freq_cm1": 120.0,
                "frequencies_cm1": [120.0, 200.0, 500.0],
            },
            "visualizations": {},
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "虚频" in html
    assert "动力学稳定" in html


def test_phonon_report_includes_how_to_read_chart_guidance():
    generator = PhononHTMLGenerator(
        {
            "phonon": {
                "n_modes": 12,
                "n_imaginary": 1,
                "dynamically_stable": False,
                "max_imaginary_freq_cm1": -35.0,
                "max_real_freq_cm1": 500.0,
                "min_real_freq_cm1": 120.0,
                "frequencies_cm1": [-35.0, 120.0, 500.0],
            },
            "visualizations": {"phonon_hist": "ZmFrZS1wbmc="},
        }
    )

    html = generator._generate_body_sections()

    assert "怎么看这张图" in html
    assert "先看有没有落在 0 以下的频率" in html
    assert "数量越多、绝对值越大" in html


def test_neb_report_includes_plain_language_interpretation():
    generator = NEBHTMLGenerator(
        {
            "neb": {
                "n_images": 5,
                "relative_energies": [0.0, 0.1, 0.25, 0.12, -0.05],
                "forward_barrier_eV": 0.25,
                "backward_barrier_eV": 0.30,
                "reaction_energy_eV": -0.05,
                "ts_image_index": 2,
            },
            "visualizations": {},
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "能垒越低" in html
    assert "更容易发生" in html


def test_neb_report_includes_how_to_read_chart_guidance():
    generator = NEBHTMLGenerator(
        {
            "neb": {
                "n_images": 5,
                "relative_energies": [0.0, 0.1, 0.25, 0.12, -0.05],
                "forward_barrier_eV": 0.25,
                "backward_barrier_eV": 0.30,
                "reaction_energy_eV": -0.05,
                "ts_image_index": 2,
            },
            "visualizations": {"neb_profile": "ZmFrZS1wbmc="},
        }
    )

    html = generator._generate_body_sections()

    assert "怎么看这张图" in html
    assert "横轴是反应路径上的图像序号" in html
    assert "曲线最高点通常就是最难跨过去的能量山峰" in html


def test_optimization_export_writes_preview_pngs_and_csvs(tmp_path):
    from src.vasp_server.analyzers.optimization import export_optimization_assets

    export_optimization_assets(
        {
            "scf_energies": [-10.0, -10.2, -10.25],
            "convergence_analysis": {
                "force_convergence": {"force_history": [0.2, 0.1, 0.03]},
                "energy_convergence": {"energy_history": [-10.0, -10.2, -10.25]},
            },
            "optimization_process": {
                "energy_profile": [
                    {"step": 1, "energy": -10.0},
                    {"step": 2, "energy": -10.2},
                    {"step": 3, "energy": -10.25},
                ],
                "force_profile": [
                    {"step": 1, "max_force": 0.2, "rms_force": 0.1},
                    {"step": 2, "max_force": 0.1, "rms_force": 0.05},
                    {"step": 3, "max_force": 0.03, "rms_force": 0.01},
                ],
            },
        },
        tmp_path,
    )

    expected = {
        "force_convergence.csv",
        "scf_energy_convergence.csv",
        "energy_evolution.csv",
        "force_evolution.csv",
        "force_convergence.png",
        "scf_energy_convergence.png",
        "energy_evolution.png",
        "force_evolution.png",
    }
    assert expected.issubset({path.name for path in tmp_path.iterdir()})


def test_optimization_report_marks_missing_energy_history_as_unavailable():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-1", "composition": "Li2MnO3"},
            "file_info": {},
            "scf_energies": [],
            "convergence_analysis": {
                "overall_convergence": True,
                "force_convergence": {
                    "converged": True,
                    "threshold": 0.02,
                    "final_max_force": 0.01,
                    "force_history": [0.2, 0.05, 0.01],
                },
                "energy_convergence": {
                    "available": False,
                    "reason": "missing_energy_history",
                },
            },
            "optimization_process": {
                "total_steps": 3,
                "max_steps": 500,
                "completion_ratio": 1.0,
                "is_converged": True,
                "energy_profile": [
                    {"step": 1, "energy": -10.0},
                    {"step": 2, "energy": -10.2},
                    {"step": 3, "energy": -10.25},
                ],
                "force_profile": [
                    {"step": 1, "max_force": 0.2, "rms_force": 0.1},
                    {"step": 2, "max_force": 0.05, "rms_force": 0.03},
                    {"step": 3, "max_force": 0.01, "rms_force": 0.005},
                ],
            },
            "final_results": {
                "final_energy": -10.25,
                "final_step": 3,
                "max_residual_force": 0.01,
                "rms_residual_force": 0.005,
            },
        }
    )

    html = generator._generate_body_sections()

    assert "能量波动分析" in html
    assert "⚪ 数据不足" in html
    assert "最终能量:</strong> --" in html
    assert "平均能量变化:</strong> --" in html
    assert "0.000000 eV" not in html
    assert "0.00e+00 eV" not in html
    assert "最终最大力:</strong> 0.0000 eV/Å" not in html


def test_optimization_convergence_ignores_missing_energy_history_when_force_converged(tmp_path):
    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text("", encoding="utf-8")

    analyzer = OUTCARAnalyzer(str(outcar_path), task_id="opt-2")
    analyzer.data["calculation_settings"] = {"EDIFFG": -0.02}
    analyzer.data["ionic_steps"] = [
        {"step": 1, "max_force": 0.01, "rms_force": 0.005, "energy": -10.0},
    ]
    analyzer.data["scf_energies"] = []

    analyzer._analyze_convergence()

    convergence = analyzer.data["convergence_analysis"]
    assert convergence["force_convergence"]["converged"] is True
    assert convergence["energy_convergence"]["available"] is False
    assert convergence["energy_convergence"]["reason"] == "missing_energy_history"
    assert convergence["overall_convergence"] is True


def test_optimization_convergence_does_not_treat_voluntary_context_switches_as_tail_success(tmp_path):
    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text(
        "General timing and accounting informations for this job:\n"
        "Voluntary context switches: 12345\n",
        encoding="utf-8",
    )

    analyzer = OUTCARAnalyzer(str(outcar_path), task_id="opt-3")
    analyzer.data["ionic_steps"] = []
    analyzer._analyze_convergence()

    convergence = analyzer.data["convergence_analysis"]
    assert convergence["tail_check"]["matched"] is False
    assert convergence["overall_convergence"] is False


def test_optimization_convergence_still_analyzes_energy_without_ionic_steps(tmp_path):
    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text("", encoding="utf-8")

    analyzer = OUTCARAnalyzer(str(outcar_path), task_id="opt-4")
    analyzer.data["calculation_settings"] = {"EDIFFG": -0.02}
    analyzer.data["ionic_steps"] = []
    analyzer.data["scf_energies"] = [-11.0, -11.1, -11.10001]

    analyzer._analyze_convergence()

    convergence = analyzer.data["convergence_analysis"]
    assert convergence["force_convergence"]["available"] is False
    assert convergence["energy_convergence"]["available"] is True
    assert convergence["energy_convergence"]["final_energy"] == -11.10001


def test_optimization_energy_convergence_uses_ediff_threshold(tmp_path):
    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text("", encoding="utf-8")

    analyzer = OUTCARAnalyzer(str(outcar_path), task_id="opt-5")
    analyzer.data["calculation_settings"] = {"EDIFFG": -0.02, "EDIFF": 1e-2}
    analyzer.data["ionic_steps"] = []
    analyzer.data["scf_energies"] = [-11.0, -11.005, -11.009]

    analyzer._analyze_convergence()

    convergence = analyzer.data["convergence_analysis"]
    energy_conv = convergence["energy_convergence"]
    assert energy_conv["available"] is True
    assert energy_conv["threshold"] == 1e-2
    assert energy_conv["converged"] is True


def test_optimization_report_shows_energy_threshold():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-6", "composition": "MgO"},
            "file_info": {},
            "scf_energies": [-11.0, -11.005, -11.009],
            "convergence_analysis": {
                "overall_convergence": True,
                "force_convergence": {
                    "available": False,
                    "converged": None,
                    "threshold": 0.02,
                    "final_max_force": None,
                    "force_history": [],
                },
                "energy_convergence": {
                    "available": True,
                    "converged": True,
                    "threshold": 1e-2,
                    "final_energy": -11.009,
                    "energy_history": [-11.0, -11.005, -11.009],
                    "energy_changes": [0.005, 0.004],
                    "avg_recent_change": 0.0045,
                },
            },
            "optimization_process": {
                "total_steps": 0,
                "max_steps": 500,
                "completion_ratio": 0.0,
                "is_converged": True,
                "energy_profile": [],
                "force_profile": [],
            },
            "final_results": {
                "final_energy": -11.009,
                "final_step": 0,
                "max_residual_force": None,
                "rms_residual_force": None,
            },
        }
    )

    html = generator._generate_body_sections()

    assert "能量波动分析" in html
    assert "收敛阈值:</strong> 1.00e-02 eV" in html


def test_optimization_report_uses_energy_stability_wording():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-7", "composition": "MgO"},
            "file_info": {},
            "scf_energies": [-11.0, -11.001, -11.004],
            "convergence_analysis": {
                "overall_convergence": False,
                "force_convergence": {
                    "available": False,
                    "converged": None,
                    "threshold": 0.02,
                    "final_max_force": None,
                    "force_history": [],
                },
                "energy_convergence": {
                    "available": True,
                    "converged": False,
                    "threshold": 1e-4,
                    "final_energy": -11.004,
                    "energy_history": [-11.0, -11.001, -11.004],
                    "energy_changes": [0.001, 0.003],
                    "avg_recent_change": 0.002,
                },
            },
            "optimization_process": {
                "total_steps": 0,
                "max_steps": 500,
                "completion_ratio": 0.0,
                "is_converged": False,
                "energy_profile": [],
                "force_profile": [],
            },
            "final_results": {
                "final_energy": -11.004,
                "final_step": 0,
                "max_residual_force": None,
                "rms_residual_force": None,
            },
        }
    )

    html = generator._generate_body_sections()

    assert "❌ 波动较大" in html
    assert "未收敛" not in html or "收敛状态" in html


def test_optimization_report_includes_plain_language_interpretation():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-plain", "composition": "MgO"},
            "file_info": {},
            "scf_energies": [-11.0, -11.001, -11.004],
            "convergence_analysis": {
                "overall_convergence": True,
                "force_convergence": {
                    "available": True,
                    "converged": True,
                    "threshold": 0.02,
                    "final_max_force": 0.01,
                    "force_history": [0.03, 0.01],
                },
                "energy_convergence": {
                    "available": True,
                    "converged": False,
                    "threshold": 1e-4,
                    "final_energy": -11.004,
                    "energy_history": [-11.0, -11.001, -11.004],
                    "energy_changes": [0.001, 0.003],
                    "avg_recent_change": 0.002,
                },
            },
            "optimization_process": {
                "total_steps": 2,
                "max_steps": 500,
                "completion_ratio": 0.004,
                "is_converged": True,
                "data_available": True,
                "energy_profile": [{"step": 1, "energy": -11.004}],
                "force_profile": [{"step": 1, "max_force": 0.01, "rms_force": 0.005}],
            },
            "final_results": {
                "final_energy": -11.004,
                "final_step": 2,
                "max_residual_force": 0.01,
                "rms_residual_force": 0.005,
                "force_statistics": {"mean_force_magnitude": 0.004, "std_force_magnitude": 0.001},
            },
        }
    )

    html = generator._generate_body_sections()

    assert "通俗解读" in html
    assert "这份报告最重要的是看结构是否已经基本定住" in html


def test_optimization_process_handles_string_nsw_without_crashing(tmp_path):
    outcar_path = tmp_path / "OUTCAR"
    outcar_path.write_text("", encoding="utf-8")

    analyzer = OUTCARAnalyzer(str(outcar_path), task_id="opt-8")
    analyzer.data["calculation_settings"] = {"NSW": "500", "EDIFFG": "-0.02", "EDIFF": "1E-4"}
    analyzer.data["ionic_steps"] = []
    analyzer.data["scf_energies"] = [-11.0, -11.00005, -11.00008]

    analyzer._analyze_convergence()
    analyzer._analyze_optimization_process()

    process = analyzer.data["optimization_process"]
    assert process["max_steps"] == 500
    assert process["completion_ratio"] == 0.0


def test_optimization_report_avoids_process_monitoring_wording():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-9", "composition": "MgO"},
            "file_info": {},
            "scf_energies": [-11.0, -11.001, -11.004],
            "convergence_analysis": {
                "overall_convergence": False,
                "force_convergence": {
                    "available": True,
                    "converged": True,
                    "threshold": 0.02,
                    "final_max_force": 0.01,
                    "force_history": [0.03, 0.01],
                },
                "energy_convergence": {
                    "available": True,
                    "converged": False,
                    "threshold": 1e-4,
                    "final_energy": -11.004,
                    "energy_history": [-11.0, -11.001, -11.004],
                    "energy_changes": [0.001, 0.003],
                    "avg_recent_change": 0.002,
                },
            },
            "optimization_process": {
                "total_steps": 1,
                "max_steps": 500,
                "completion_ratio": 0.002,
                "is_converged": False,
                "data_available": True,
                "energy_profile": [{"step": 1, "energy": -11.004}],
                "force_profile": [{"step": 1, "max_force": 0.01, "rms_force": 0.005}],
            },
            "final_results": {
                "final_energy": -11.004,
                "final_step": 1,
                "max_residual_force": 0.01,
                "rms_residual_force": 0.005,
            },
        }
    )

    html = generator._generate_body_sections()

    assert "优化过程监控" not in html
    assert "完成状态" not in html
    assert "步骤进度" not in html
    assert "进行中" not in html
    assert "收敛状态" not in html
    assert "任务状态" in html
    assert "已完成" in html
    assert "实际离子步数" in html
    assert "迭代曲线参考" in html


def test_optimization_sections_include_how_to_read_chart_guidance():
    generator = OptimizationHTMLGenerator(
        {
            "task_info": {"task_id": "opt-guide", "composition": "Li2O"},
            "file_info": {},
            "scf_energies": [-10.0, -10.1, -10.11],
            "convergence_analysis": {
                "force_convergence": {
                    "converged": True,
                    "threshold": 0.02,
                    "final_max_force": 0.01,
                    "force_history": [0.2, 0.08, 0.01],
                },
                "energy_convergence": {
                    "converged": True,
                    "threshold": 1e-4,
                    "final_energy": -10.11,
                    "avg_recent_change": 8e-5,
                    "energy_history": [-10.0, -10.1, -10.11],
                },
            },
            "optimization_process": {
                "data_available": True,
                "energy_profile": [{"step": 1, "energy": -10.0}, {"step": 2, "energy": -10.08}],
                "force_profile": [{"step": 1, "max_force": 0.2, "rms_force": 0.1}, {"step": 2, "max_force": 0.08, "rms_force": 0.04}],
            },
        }
    )

    convergence_html = generator._generate_convergence_section()
    process_html = generator._generate_optimization_process_section()

    assert "红线是否逐步压到阈值线下方" in convergence_html
    assert "后几步是否只剩小幅波动" in convergence_html
    assert "整体能量是否逐步下降并趋于平缓" in process_html
    assert "最大力和 RMS 力是否一起往下走" in process_html
