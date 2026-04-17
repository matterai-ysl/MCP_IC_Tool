"""
BandStructureAnalyzer — 从VASP能带结构计算输出中提取数据并生成报告。

依赖 pymatgen 解析 vasprun.xml 获取能带结构数据，
使用 matplotlib 绘制能带图，生成HTML分析报告。
"""

import base64
import csv
import logging
import math
import re
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, Iterable, Optional

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator

logger = logging.getLogger(__name__)


def _format_kpoint_label(kpoint) -> Optional[str]:
    if not kpoint:
        return None
    label = getattr(kpoint, "label", None)
    if label:
        return str(label)
    frac = getattr(kpoint, "frac_coords", None)
    if frac is None:
        return None
    try:
        return "(" + ", ".join(f"{float(v):.3f}" for v in frac) + ")"
    except Exception:
        return None


def _format_kpoint_coords(kpoint) -> Optional[str]:
    if not kpoint:
        return None
    frac = getattr(kpoint, "frac_coords", None)
    if frac is None:
        return None
    try:
        return "(" + ", ".join(f"{float(v):.3f}" for v in frac) + ")"
    except Exception:
        return None


def _display_branch_name(bs: Any, branch: Dict[str, Any]) -> str:
    raw_name = str(branch.get("name") or "").strip()
    if raw_name and "None" not in raw_name:
        return raw_name

    start_kpoint = None
    end_kpoint = None
    try:
        start_kpoint = bs.kpoints[int(branch["start_index"])]
        end_kpoint = bs.kpoints[int(branch["end_index"])]
    except Exception:
        return raw_name or "N/A"

    start_label = _format_kpoint_label(start_kpoint) or _format_kpoint_coords(start_kpoint) or "N/A"
    end_label = _format_kpoint_label(end_kpoint) or _format_kpoint_coords(end_kpoint) or "N/A"
    return f"{start_label} -> {end_label}"


def _normalize_spin_key(spin: Any) -> str:
    return getattr(spin, "name", str(spin)).lower().split(".")[-1]


def _normalize_orbital_key(orbital: Any) -> str:
    if hasattr(orbital, "name"):
        return str(getattr(orbital, "name")).lower()
    label = str(orbital)
    if "." in label:
        label = label.split(".")[-1]
    return label.lower()


def _site_species(structure: Any) -> list[str]:
    species: list[str] = []
    for site in getattr(structure, "sites", []):
        label = getattr(site, "species_string", None)
        if not label:
            specie = getattr(site, "specie", None)
            label = getattr(specie, "symbol", None) or str(specie or "X")
        species.append(str(label))
    return species


def summarize_projection_contributions(
    projections: Dict[Any, Dict[Any, Iterable[float]]],
    site_species: list[str],
    top_n: int = 3,
) -> Dict[str, list[Dict[str, float | str]]]:
    """Aggregate edge projections into element and element-orbital contributions."""
    species_orbital_totals: Dict[str, float] = {}
    element_totals: Dict[str, float] = {}
    total_weight = 0.0

    for orbital_map in (projections or {}).values():
        for orbital, weights in (orbital_map or {}).items():
            orbital_label = _normalize_orbital_key(orbital)
            for idx, raw_weight in enumerate(weights or []):
                try:
                    weight = float(raw_weight)
                except (TypeError, ValueError):
                    continue
                if idx >= len(site_species):
                    continue
                species = site_species[idx]
                total_weight += weight
                element_totals[species] = element_totals.get(species, 0.0) + weight
                key = f"{species}-{orbital_label}"
                species_orbital_totals[key] = species_orbital_totals.get(key, 0.0) + weight

    if total_weight <= 0:
        return {"top_species_orbitals": [], "top_elements": []}

    def _top_items(source: Dict[str, float]) -> list[Dict[str, float | str]]:
        return [
            {"label": label, "weight": value / total_weight}
            for label, value in sorted(source.items(), key=lambda item: item[1], reverse=True)[:top_n]
        ]

    return {
        "top_species_orbitals": _top_items(species_orbital_totals),
        "top_elements": _top_items(element_totals),
    }


def _band_indices_for_display(band_index_map: Dict[Any, Iterable[int]]) -> Dict[str, list[int]]:
    return {
        _normalize_spin_key(spin): [int(idx) + 1 for idx in indices]
        for spin, indices in (band_index_map or {}).items()
        if indices
    }


def _band_width_summaries(bs: Any, band_index_map: Dict[Any, Iterable[int]]) -> list[Dict[str, float | str]]:
    items: list[Dict[str, float | str]] = []
    for spin, indices in (band_index_map or {}).items():
        spin_name = _normalize_spin_key(spin)
        energies = bs.bands.get(spin)
        if energies is None:
            continue
        for idx in indices:
            band = energies[idx]
            width = float(max(band) - min(band))
            items.append({"band": f"{spin_name}:{int(idx) + 1}", "width": width})
    return items


def _dispersion_label(width: float | None) -> str:
    if width is None:
        return "N/A"
    if width < 1.0:
        return "较平坦"
    if width < 3.0:
        return "中等色散"
    return "较强色散"


def _format_percentage_items(items: list[Dict[str, float | str]], empty: str = "N/A") -> str:
    if not items:
        return empty
    formatted = []
    for item in items:
        label = str(item.get("label", "N/A"))
        weight = item.get("weight")
        if isinstance(weight, (int, float)):
            formatted.append(f"{label} ({weight * 100:.1f}%)")
        else:
            formatted.append(label)
    return "，".join(formatted)


def _format_band_widths(items: list[Dict[str, float | str]], empty: str = "N/A") -> str:
    if not items:
        return empty
    return "，".join(
        f"{item['band']} ({float(item['width']):.2f} eV)"
        for item in items
        if item.get("width") is not None
    ) or empty


def _display_transition(bs_data: Dict[str, Any]) -> str:
    transition = str(bs_data.get("transition") or "").strip()
    direct_label = bs_data.get("direct_gap_kpoint_label") or bs_data.get("direct_gap_kpoint_coords")

    if bs_data.get("is_direct") and direct_label:
        return f"位于 {direct_label} 点"
    if transition:
        return transition
    return "N/A"


def _projection_summary_or_placeholder(items: list[Dict[str, float | str]], empty: str = "未提供投影数据") -> str:
    return _format_percentage_items(items, empty=empty)


def _friendly_branch_display(branch: Any) -> str:
    text = str(branch or "").strip()
    if not text:
        return "当前路径"
    if "None" in text:
        return "未标注路径（带边附近）"
    return text


def estimate_effective_mass_from_parabola(
    k_values: Iterable[float],
    energies: Iterable[float],
    carrier: str,
) -> Optional[Dict[str, float]]:
    """Estimate a 1D effective mass from a local parabolic fit in eV vs 1/Å."""
    try:
        import numpy as np
    except ImportError:
        return None

    k = np.asarray(list(k_values), dtype=float)
    e = np.asarray(list(energies), dtype=float)
    if len(k) < 3 or len(e) < 3 or len(k) != len(e):
        return None
    if np.unique(k).size < 3:
        return None

    coeffs = np.polyfit(k, e, 2)
    curvature = float(coeffs[0])
    if math.isclose(curvature, 0.0, abs_tol=1e-12):
        return None

    carrier_lower = carrier.lower()
    if carrier_lower == "electron":
        if curvature <= 0:
            return None
        signed_mass = 3.80998212 / curvature
    elif carrier_lower == "hole":
        if curvature >= 0:
            return None
        signed_mass = 3.80998212 / curvature
    else:
        raise ValueError(f"unsupported carrier type: {carrier}")

    return {
        "m_over_me": abs(signed_mass),
        "curvature_eva2": curvature,
        "fit_points": float(len(k)),
    }


def _estimate_edge_effective_mass(
    bs: Any,
    edge_info: Dict[str, Any],
    carrier: str,
    max_window: int = 2,
) -> Optional[Dict[str, Any]]:
    kpoint_indices = list(edge_info.get("kpoint_index") or [])
    band_index_map = edge_info.get("band_index") or {}
    if not kpoint_indices or not band_index_map:
        return None

    best_result: Optional[Dict[str, Any]] = None
    for branch in getattr(bs, "branches", []) or []:
        start = int(branch["start_index"])
        end = int(branch["end_index"])
        branch_indices = [idx for idx in kpoint_indices if start <= int(idx) <= end]
        if not branch_indices:
            continue
        branch_name = _display_branch_name(bs, branch)
        for edge_kidx in branch_indices:
            for spin, band_indices in band_index_map.items():
                energies = bs.bands.get(spin)
                if energies is None:
                    continue
                for band_idx in band_indices:
                    local_k = []
                    local_e = []
                    for neighbor in range(max(start, edge_kidx - max_window), min(end, edge_kidx + max_window) + 1):
                        local_k.append(float(bs.distance[neighbor]))
                        local_e.append(float(energies[band_idx][neighbor]))
                    result = estimate_effective_mass_from_parabola(local_k, local_e, carrier=carrier)
                    if not result:
                        continue
                    candidate = {
                        "m_over_me": result["m_over_me"],
                        "curvature_eva2": result["curvature_eva2"],
                        "fit_points": int(result["fit_points"]),
                        "branch": branch_name,
                        "band": f"{_normalize_spin_key(spin)}:{int(band_idx) + 1}",
                        "kpoint_label": _format_kpoint_label(bs.kpoints[edge_kidx]) or _format_kpoint_coords(bs.kpoints[edge_kidx]),
                    }
                    if best_result is None or candidate["fit_points"] > best_result["fit_points"]:
                        best_result = candidate
    return best_result


def _build_research_notes(bs_data: Dict[str, Any], results: Dict[str, Any]) -> list[str]:
    notes: list[str] = []
    transition_text = _display_transition(bs_data)
    if bs_data.get("is_metal"):
        notes.append("该体系存在穿越费米能级的能带，重点应转向费米面附近的带色散、自旋分裂与可能的导电通道。")
    elif bs_data.get("is_direct"):
        notes.append(
            f"基态表现为直接带隙，最低直接跃迁{transition_text}；在不考虑跃迁矩阵元时，通常更有利于直接光学跃迁。"
        )
    else:
        notes.append(
            f"基态表现为间接带隙，基本带隙路径为 {transition_text}；光学激发往往需要额外的声子参与。"
        )

    direct_gap = bs_data.get("direct_gap")
    band_gap = bs_data.get("band_gap")
    if isinstance(direct_gap, (int, float)) and isinstance(band_gap, (int, float)) and direct_gap > band_gap + 1e-3:
        notes.append(
            f"最低直接带隙约为 {direct_gap:.4f} eV，高于基本带隙 {band_gap:.4f} eV，说明最低能激发与最低光学跃迁并不完全重合。"
        )

    fermi_offset = bs_data.get("fermi_offset_from_midgap")
    if isinstance(fermi_offset, (int, float)):
        if abs(fermi_offset) <= 0.1:
            notes.append("费米能级接近带隙中心，电子结构更接近本征半导体。")
        elif fermi_offset > 0:
            notes.append(f"费米能级相对带隙中心上移 {fermi_offset:.4f} eV，更接近导带侧，可继续结合缺陷/掺杂判断 n 型趋势。")
        else:
            notes.append(f"费米能级相对带隙中心下移 {abs(fermi_offset):.4f} eV，更接近价带侧，可继续结合缺陷/掺杂判断 p 型趋势。")

    vbm_character = ((bs_data.get("vbm_character") or {}).get("top_species_orbitals") or [])
    cbm_character = ((bs_data.get("cbm_character") or {}).get("top_species_orbitals") or [])
    if vbm_character:
        top = vbm_character[0]
        notes.append(
            f"价带顶主要由 {top['label']} 贡献（约 {float(top['weight']) * 100:.1f}%），这通常决定空穴端的轨道性质与局域化倾向。"
        )
    if cbm_character:
        top = cbm_character[0]
        notes.append(
            f"导带底主要由 {top['label']} 贡献（约 {float(top['weight']) * 100:.1f}%），这通常决定电子端的轨道来源与导电通道特征。"
        )
    if not vbm_character and not cbm_character:
        notes.append("当前结果未包含可用的投影本征值，暂时无法可靠判断 VBM/CBM 的主导元素与轨道组成。")

    vbm_widths = bs_data.get("vbm_band_widths") or []
    cbm_widths = bs_data.get("cbm_band_widths") or []
    if vbm_widths:
        notes.append(
            f"沿当前高对称路径，价带边代表带宽约 {float(vbm_widths[0]['width']):.2f} eV，表现为{_dispersion_label(float(vbm_widths[0]['width']))}。"
        )
    if cbm_widths:
        notes.append(
            f"沿当前高对称路径，导带边代表带宽约 {float(cbm_widths[0]['width']):.2f} eV，表现为{_dispersion_label(float(cbm_widths[0]['width']))}。"
        )

    spin_channels = bs_data.get("spin_channels") or []
    if len(spin_channels) > 1:
        notes.append("体系包含多个自旋通道，建议结合自旋分辨 DOS 或投影带进一步判断自旋极化与交换劈裂。")

    masses = bs_data.get("effective_masses") or {}
    electron_mass = masses.get("electron")
    hole_mass = masses.get("hole")
    if electron_mass and isinstance(electron_mass.get("m_over_me"), (int, float)):
        notes.append(
            f"沿 {_friendly_branch_display(electron_mass.get('branch'))}，电子有效质量约为 {float(electron_mass['m_over_me']):.2f} mₑ。"
        )
    if hole_mass and isinstance(hole_mass.get("m_over_me"), (int, float)):
        notes.append(
            f"沿 {_friendly_branch_display(hole_mass.get('branch'))}，空穴有效质量约为 {float(hole_mass['m_over_me']):.2f} mₑ。"
        )

    return notes


class BandStructureAnalyzer(BaseAnalyzer):
    """从 vasprun.xml / EIGENVAL / OUTCAR 中提取能带结构信息。"""

    def __init__(self, input_path: str, task_id: Optional[str] = None):
        super().__init__(input_path, task_id=task_id or "band")

    def _init_data(self) -> Dict[str, Any]:
        return {
            'file_info': {},
            'calculation_settings': {},
            'final_results': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {},
            'band_structure': {},
            'visualizations': {},
        }

    def _do_analyze(self):
        """执行能带结构分析。"""
        self._analyze_band_structure()
        self._generate_plots()

    def _analyze_band_structure(self):
        """使用 pymatgen 解析能带结构。"""
        vasprun_path = self.work_dir / "vasprun.xml"

        if not vasprun_path.exists():
            logger.warning("vasprun.xml 不存在，跳过能带结构解析")
            self.data['band_structure'] = {'error': 'vasprun.xml 不存在'}
            return

        try:
            from pymatgen.io.vasp import Vasprun

            vasprun = Vasprun(str(vasprun_path), parse_projected_eigen=True)
            bs = vasprun.get_band_structure(line_mode=True)

            gap_info = bs.get_band_gap()
            band_gap = gap_info['energy']
            is_direct = gap_info['direct']
            transition = gap_info.get('transition', '')

            vbm_info = bs.get_vbm()
            cbm_info = bs.get_cbm()
            direct_gap = None
            direct_gap_kpoint_label = None
            direct_gap_kpoint_coords = None
            direct_gap_band_indices = {}
            if not bs.is_metal():
                direct_gap_dict = bs.get_direct_band_gap_dict()
                best_spin, best_info = min(direct_gap_dict.items(), key=lambda item: item[1]["value"])
                direct_gap = float(best_info["value"])
                direct_gap_kpoint = bs.kpoints[int(best_info["kpoint_index"])]
                direct_gap_kpoint_label = _format_kpoint_label(direct_gap_kpoint)
                direct_gap_kpoint_coords = _format_kpoint_coords(direct_gap_kpoint)
                direct_gap_band_indices = {
                    _normalize_spin_key(best_spin): [int(idx) + 1 for idx in best_info.get("band_indices", [])]
                }

            site_species = _site_species(getattr(bs, "structure", None))
            vbm_character = summarize_projection_contributions(vbm_info.get("projections") or {}, site_species)
            cbm_character = summarize_projection_contributions(cbm_info.get("projections") or {}, site_species)
            vbm_band_indices = _band_indices_for_display(vbm_info.get("band_index") or {})
            cbm_band_indices = _band_indices_for_display(cbm_info.get("band_index") or {})
            vbm_band_widths = _band_width_summaries(bs, vbm_info.get("band_index") or {})
            cbm_band_widths = _band_width_summaries(bs, cbm_info.get("band_index") or {})
            effective_masses = {
                "electron": _estimate_edge_effective_mass(bs, cbm_info, carrier="electron"),
                "hole": _estimate_edge_effective_mass(bs, vbm_info, carrier="hole"),
            }
            fermi_energy = float(vasprun.efermi)
            fermi_offset_from_midgap = None
            if vbm_info and cbm_info and vbm_info.get("energy") is not None and cbm_info.get("energy") is not None:
                midgap = 0.5 * (float(vbm_info["energy"]) + float(cbm_info["energy"]))
                fermi_offset_from_midgap = fermi_energy - midgap

            self.data['band_structure'] = {
                'band_gap': band_gap,
                'direct_gap': direct_gap,
                'direct_gap_kpoint_label': direct_gap_kpoint_label,
                'direct_gap_kpoint_coords': direct_gap_kpoint_coords,
                'direct_gap_band_indices': direct_gap_band_indices,
                'is_direct': is_direct,
                'transition': transition,
                'vbm': vbm_info['energy'] if vbm_info else None,
                'cbm': cbm_info['energy'] if cbm_info else None,
                'is_metal': bs.is_metal(),
                'nb_bands': bs.nb_bands,
                'nkpoints': len(getattr(bs, 'kpoints', []) or []),
                'n_branches': len(getattr(bs, 'branches', []) or []),
                'spin_channels': [
                    getattr(spin, 'name', str(spin)).lower().split('.')[-1]
                    for spin in bs.bands.keys()
                ],
                'vbm_kpoint_label': _format_kpoint_label(vbm_info.get('kpoint') if vbm_info else None),
                'vbm_kpoint_coords': _format_kpoint_coords(vbm_info.get('kpoint') if vbm_info else None),
                'cbm_kpoint_label': _format_kpoint_label(cbm_info.get('kpoint') if cbm_info else None),
                'cbm_kpoint_coords': _format_kpoint_coords(cbm_info.get('kpoint') if cbm_info else None),
                'vbm_band_indices': vbm_band_indices,
                'cbm_band_indices': cbm_band_indices,
                'vbm_band_degeneracy': sum(len(indices) for indices in vbm_band_indices.values()),
                'cbm_band_degeneracy': sum(len(indices) for indices in cbm_band_indices.values()),
                'vbm_kpoint_degeneracy': len(vbm_info.get('kpoint_index') or []),
                'cbm_kpoint_degeneracy': len(cbm_info.get('kpoint_index') or []),
                'fermi_offset_from_midgap': fermi_offset_from_midgap,
                'vbm_band_widths': vbm_band_widths,
                'cbm_band_widths': cbm_band_widths,
                'effective_masses': effective_masses,
                'vbm_character': vbm_character,
                'cbm_character': cbm_character,
            }

            self.data['final_results']['band_gap'] = band_gap
            self.data['final_results']['is_direct'] = is_direct
            self.data['final_results']['vbm'] = vbm_info['energy'] if vbm_info else None
            self.data['final_results']['cbm'] = cbm_info['energy'] if cbm_info else None
            self.data['final_results']['is_metal'] = bs.is_metal()

            # 提取费米能级和总能量
            self.data['final_results']['fermi_energy'] = fermi_energy
            self.data['final_results']['total_energy'] = vasprun.final_energy
            self.data['raw_band_data'] = {
                'distance': [float(v) for v in getattr(bs, 'distance', [])],
                'bands': {
                    getattr(spin, 'name', str(spin)).lower().split('.')[-1]: energies.tolist()
                    for spin, energies in bs.bands.items()
                },
            }

            # 保存 band structure 对象供绘图使用
            self._bs = bs
            self._vasprun = vasprun
            self.data['band_structure']['research_notes'] = _build_research_notes(
                self.data['band_structure'],
                self.data['final_results'],
            )

            logger.info("能带结构解析完成: gap=%.4f eV, direct=%s, metal=%s",
                        band_gap, is_direct, bs.is_metal())

        except ImportError:
            logger.error("pymatgen 未安装")
            self.data['band_structure'] = {'error': 'pymatgen 未安装'}
        except Exception as e:
            logger.error("能带结构解析失败: %s", e)
            self.data['band_structure'] = {'error': str(e)}
            self._parse_eigenval_fallback()

    def _parse_eigenval_fallback(self):
        """EIGENVAL 回退解析（当 vasprun.xml 解析失败时）。"""
        eigenval_path = self.work_dir / "EIGENVAL"
        if not eigenval_path.exists():
            return

        try:
            content = self._get_outcar_content()
            if content:
                # 提取费米能级
                m = re.search(r'E-fermi\s*:\s*([-\d.]+)', content)
                if m:
                    self.data['final_results']['fermi_energy'] = float(m.group(1))

                # 提取带隙
                m = re.search(r'band gap\s*=\s*([\d.]+)', content, re.IGNORECASE)
                if m:
                    self.data['final_results']['band_gap'] = float(m.group(1))
        except Exception as e:
            logger.warning("EIGENVAL 回退解析失败: %s", e)

    def _generate_plots(self):
        """生成能带结构图。"""
        if not hasattr(self, '_bs'):
            return

        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt
            from pymatgen.electronic_structure.plotter import BSPlotter

            plotter = BSPlotter(self._bs)

            # 能带结构图
            plot_obj = plotter.get_plot(ylim=(-5, 5))
            fig = plot_obj if hasattr(plot_obj, 'savefig') else getattr(plot_obj, 'figure', None)
            if fig is None:
                raise TypeError(f"unexpected plot object type: {type(plot_obj)!r}")
            ax = plot_obj if hasattr(plot_obj, 'axhline') else (fig.axes[0] if getattr(fig, 'axes', None) else None)
            if ax is not None:
                # BSPlotter aligns energies to E_F = 0; draw an explicit guide line so the
                # rendered plot matches the report wording and is easy to spot visually.
                ax.axhline(0, color='red', linestyle='--', linewidth=1.2, alpha=0.9, zorder=10)
            buf = BytesIO()
            fig.savefig(buf, format='png', dpi=150, bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            self.data['visualizations']['band_structure_plot'] = base64.b64encode(buf.read()).decode()

            logger.info("能带结构图已生成")

        except Exception as e:
            logger.warning("能带结构图生成失败: %s", e)


class BandStructureHTMLGenerator(BaseHTMLGenerator):
    """生成能带结构分析 HTML 报告。"""

    report_title = "能带结构分析报告"
    report_emoji = ""
    report_heading = "能带结构分析报告"
    report_attribution = "由VASP API能带结构分析模块生成"
    header_gradient = "linear-gradient(135deg, #4facfe 0%, #00f2fe 100%)"
    use_chartjs = False

    def _generate_body_sections(self) -> str:
        sections = []
        sections.append(self._section_summary())
        sections.append(self._section_band_structure())
        sections.append(self._section_effective_masses())
        sections.append(self._section_band_edge_character())
        sections.append(self._section_edge_projection_bars())
        sections.append(self._section_research_notes())
        sections.append(self._section_calculation_settings())
        return "\n".join(sections)

    def _generate_javascript(self) -> str:
        return ""

    def _section_summary(self) -> str:
        bs_data = self.data.get('band_structure', {})
        results = self.data.get('final_results', {})

        band_gap = bs_data.get('band_gap', 'N/A')
        is_direct = bs_data.get('is_direct')
        is_metal = bs_data.get('is_metal', False)
        transition = _display_transition(bs_data)
        vbm = results.get('vbm', 'N/A')
        cbm = results.get('cbm', 'N/A')
        fermi = results.get('fermi_energy', 'N/A')
        total_e = results.get('total_energy', 'N/A')
        nkpoints = bs_data.get('nkpoints', 'N/A')
        n_branches = bs_data.get('n_branches', 'N/A')
        spin_channels = ", ".join(bs_data.get('spin_channels') or []) or 'N/A'
        vbm_kpoint = bs_data.get('vbm_kpoint_label', 'N/A')
        cbm_kpoint = bs_data.get('cbm_kpoint_label', 'N/A')
        direct_gap = bs_data.get('direct_gap', 'N/A')
        direct_gap_kpoint = bs_data.get('direct_gap_kpoint_label', 'N/A')
        fermi_offset = bs_data.get('fermi_offset_from_midgap', 'N/A')

        gap_type = "直接带隙" if is_direct else "间接带隙"
        material_type = "金属" if is_metal else f"半导体/绝缘体 ({gap_type})"

        def fmt(v, unit=""):
            if v is None or v == 'N/A':
                return 'N/A'
            return f"{v:.4f} {unit}" if isinstance(v, (int, float)) else str(v)

        return f"""
        <div class="card">
            <h2>计算摘要</h2>
            <table>
                <tr><td>材料类型</td><td><strong>{material_type}</strong></td></tr>
                <tr><td>带隙</td><td><strong>{fmt(band_gap, 'eV')}</strong></td></tr>
                <tr><td>带隙类型</td><td>{gap_type if not is_metal else 'N/A (金属)'}</td></tr>
                <tr><td>价带顶 (VBM)</td><td>{fmt(vbm, 'eV')}</td></tr>
                <tr><td>导带底 (CBM)</td><td>{fmt(cbm, 'eV')}</td></tr>
                <tr><td>费米能级</td><td>{fmt(fermi, 'eV')}</td></tr>
                <tr><td>总能量</td><td>{fmt(total_e, 'eV')}</td></tr>
                <tr><td>能带数</td><td>{bs_data.get('nb_bands', 'N/A')}</td></tr>
                <tr><td>带隙跃迁</td><td>{transition}</td></tr>
                <tr><td>直接跃迁带隙</td><td>{fmt(direct_gap, 'eV')}</td></tr>
                <tr><td>最低直接跃迁 k 点</td><td>{direct_gap_kpoint or 'N/A'}</td></tr>
                <tr><td>VBM 所在 k 点</td><td>{vbm_kpoint or 'N/A'}</td></tr>
                <tr><td>CBM 所在 k 点</td><td>{cbm_kpoint or 'N/A'}</td></tr>
                <tr><td>路径 k 点数</td><td>{nkpoints}</td></tr>
                <tr><td>路径段数</td><td>{n_branches}</td></tr>
                <tr><td>自旋通道</td><td>{spin_channels}</td></tr>
                <tr><td>费米能级偏离带隙中心</td><td>{fmt(fermi_offset, 'eV')}</td></tr>
            </table>
        </div>
        """

    def _section_band_structure(self) -> str:
        vis = self.data.get('visualizations', {})
        bs_plot = vis.get('band_structure_plot')

        if not bs_plot:
            return """
            <div class="card">
                <h2>能带结构图</h2>
                <p>能带结构图生成失败或数据不可用。</p>
            </div>
            """

        png_dl = self._png_download_link('band_structure.png', bs_plot)
        csv_dl = self._band_csv_link()

        return f"""
        <div class="card">
            <h2>能带结构图</h2>
            <div style="text-align: center;">
                <img src="data:image/png;base64,{bs_plot}"
                     alt="Band Structure" style="max-width: 100%; height: auto;" />
                {png_dl}
                {csv_dl}
            </div>
            <p style="text-align: center; color: #666; margin-top: 10px;">
                能带结构 E(k) 色散关系。红色虚线为费米能级。
            </p>
            <div class="info-box" style="margin-top: 16px;">
                <strong>怎么看这张图：</strong>
                先看带隙两侧最靠近 0 eV 的带边，它们就是价带顶和导带底。
                如果价带顶和导带底出现在同一个 k 点，通常就是直接带隙；如果不在同一个 k 点，通常就是间接带隙。
                另外，曲线越分散，往往说明电子或空穴更容易移动；曲线越平，通常意味着载流子更“重”、更局域。
            </div>
        </div>
        """

    def _band_csv_link(self) -> str:
        raw = self.data.get('raw_band_data', {})
        distances = raw.get('distance') or []
        bands = raw.get('bands') or {}
        if not distances or not bands:
            return ''

        headers = ['k_distance']
        series = []
        for spin_name, spin_bands in bands.items():
            for idx, band in enumerate(spin_bands, start=1):
                headers.append(f'{spin_name}_band_{idx}')
                series.append(band)

        rows = []
        for i, distance in enumerate(distances):
            rows.append([distance] + [band[i] if i < len(band) else '' for band in series])

        return self._csv_download_link(
            'band_structure_data.csv',
            headers,
            rows,
            label='下载能带数据 (CSV)',
        )

    def _section_calculation_settings(self) -> str:
        settings = self.data.get('calculation_settings', {})
        if not settings:
            return ""

        rows = "".join(
            f"<tr><td>{k}</td><td>{v}</td></tr>"
            for k, v in sorted(settings.items())
        )
        return f"""
        <div class="card">
            <h2>计算设置</h2>
            <table>{rows}</table>
        </div>
        """

    def _section_effective_masses(self) -> str:
        masses = (self.data.get('band_structure') or {}).get('effective_masses') or {}

        def _fmt_mass(label: str, item: Dict[str, Any] | None) -> str:
            if not item or not isinstance(item.get("m_over_me"), (int, float)):
                return f"<tr><td>{label}</td><td>未能从当前 line-mode 路径稳定拟合</td></tr>"
            branch_label = _friendly_branch_display(item.get('branch'))
            return (
                f"<tr><td>{label}</td><td>"
                f"{float(item['m_over_me']):.2f} m<sub>e</sub>"
                f"（路径 {branch_label}，带 {item.get('band') or 'N/A'}，"
                f"拟合点数 {item.get('fit_points') or 'N/A'}）"
                f"</td></tr>"
            )

        return f"""
        <div class="card">
            <h2>有效质量估计</h2>
            <table>
                {_fmt_mass('电子有效质量', masses.get('electron'))}
                {_fmt_mass('空穴有效质量', masses.get('hole'))}
            </table>
            <p style="color: #666; margin-top: 10px;">
                注：这里给出的是沿当前高对称路径、带边附近的一维局部抛物线拟合结果，适合做快速科研判断，不等同于完整各向异性输运有效质量张量。
            </p>
        </div>
        """

    def _section_band_edge_character(self) -> str:
        bs_data = self.data.get('band_structure', {})

        def _fmt_indices(index_map: Dict[str, list[int]] | None) -> str:
            if not index_map:
                return 'N/A'
            return "；".join(
                f"{spin}: {', '.join(str(v) for v in values)}"
                for spin, values in index_map.items()
            )

        return f"""
        <div class="card">
            <h2>带边特征</h2>
            <table>
                <tr><td>VBM 能带编号</td><td>{_fmt_indices(bs_data.get('vbm_band_indices'))}</td></tr>
                <tr><td>CBM 能带编号</td><td>{_fmt_indices(bs_data.get('cbm_band_indices'))}</td></tr>
                <tr><td>VBM 带简并度</td><td>{bs_data.get('vbm_band_degeneracy', 'N/A')}</td></tr>
                <tr><td>CBM 带简并度</td><td>{bs_data.get('cbm_band_degeneracy', 'N/A')}</td></tr>
                <tr><td>VBM k 点简并度</td><td>{bs_data.get('vbm_kpoint_degeneracy', 'N/A')}</td></tr>
                <tr><td>CBM k 点简并度</td><td>{bs_data.get('cbm_kpoint_degeneracy', 'N/A')}</td></tr>
                <tr><td>VBM 分数坐标</td><td>{bs_data.get('vbm_kpoint_coords') or 'N/A'}</td></tr>
                <tr><td>CBM 分数坐标</td><td>{bs_data.get('cbm_kpoint_coords') or 'N/A'}</td></tr>
                <tr><td>VBM 主导轨道</td><td>{_projection_summary_or_placeholder((bs_data.get('vbm_character') or {}).get('top_species_orbitals') or [])}</td></tr>
                <tr><td>CBM 主导轨道</td><td>{_projection_summary_or_placeholder((bs_data.get('cbm_character') or {}).get('top_species_orbitals') or [])}</td></tr>
                <tr><td>VBM 主导元素</td><td>{_projection_summary_or_placeholder((bs_data.get('vbm_character') or {}).get('top_elements') or [])}</td></tr>
                <tr><td>CBM 主导元素</td><td>{_projection_summary_or_placeholder((bs_data.get('cbm_character') or {}).get('top_elements') or [])}</td></tr>
                <tr><td>VBM 边带宽</td><td>{_format_band_widths(bs_data.get('vbm_band_widths') or [])}</td></tr>
                <tr><td>CBM 边带宽</td><td>{_format_band_widths(bs_data.get('cbm_band_widths') or [])}</td></tr>
            </table>
        </div>
        """

    def _section_edge_projection_bars(self) -> str:
        bs_data = self.data.get('band_structure') or {}

        def _bar_list(title: str, items: list[Dict[str, float | str]]) -> str:
            if not items:
                return f"<div><h3>{title}</h3><p>未提供投影数据</p></div>"
            rows = []
            for item in items:
                label = str(item.get("label", "N/A"))
                weight = float(item.get("weight", 0.0))
                rows.append(
                    f"""
                    <div style="margin-bottom: 10px;">
                        <div style="display:flex;justify-content:space-between;font-size:0.95rem;">
                            <span>{label}</span><span>{weight * 100:.1f}%</span>
                        </div>
                        <div style="background:#e8eef8;border-radius:999px;height:10px;overflow:hidden;">
                            <div style="width:{max(0.0, min(weight * 100, 100.0)):.1f}%;background:linear-gradient(90deg,#4facfe,#00c6ff);height:100%;"></div>
                        </div>
                    </div>
                    """
                )
            return f"<div><h3>{title}</h3>{''.join(rows)}</div>"

        vbm_char = bs_data.get("vbm_character") or {}
        cbm_char = bs_data.get("cbm_character") or {}
        return f"""
        <div class="card">
            <h2>带边轨道组成图</h2>
            <div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(260px,1fr));gap:18px;">
                {_bar_list('VBM 主导轨道', vbm_char.get('top_species_orbitals') or [])}
                {_bar_list('CBM 主导轨道', cbm_char.get('top_species_orbitals') or [])}
                {_bar_list('VBM 主导元素', vbm_char.get('top_elements') or [])}
                {_bar_list('CBM 主导元素', cbm_char.get('top_elements') or [])}
            </div>
        </div>
        """

    def _section_research_notes(self) -> str:
        bs_data = self.data.get('band_structure') or {}
        notes = bs_data.get('research_notes') or _build_research_notes(
            bs_data,
            self.data.get('final_results') or {},
        )
        if not notes:
            return ""
        items = "".join(f"<li>{note}</li>" for note in notes)
        return f"""
        <div class="card">
            <h2>科研解读</h2>
            <ul>
                {items}
            </ul>
        </div>
        """


def generate_band_structure_report(work_dir: str, task_id: str = "band") -> Optional[str]:
    """生成能带结构分析 HTML 报告并返回文件路径。"""
    try:
        analyzer = BandStructureAnalyzer(work_dir, task_id=task_id)
        data = analyzer.analyze()

        output_dir = Path(work_dir) / "BS_output"
        output_dir.mkdir(exist_ok=True)
        export_band_structure_assets(data, output_dir)
        output_path = str(output_dir / "band_structure_report.html")

        generator = BandStructureHTMLGenerator(data)
        generator.generate_html_report(output_path)

        logger.info("能带结构 HTML 报告已保存: %s", output_path)
        return output_path

    except Exception as e:
        logger.error("生成能带结构报告失败: %s", e)
        return None


def export_band_structure_assets(data: Dict[str, Any], output_dir: Path) -> None:
    """将能带分析中的轻量图像/数据导出为独立文件，便于前端直连展示。"""
    output_dir.mkdir(parents=True, exist_ok=True)

    plot_b64 = ((data.get("visualizations") or {}).get("band_structure_plot"))
    if plot_b64:
        png_path = output_dir / "band_structure.png"
        png_path.write_bytes(base64.b64decode(plot_b64))

    raw = data.get("raw_band_data") or {}
    distances = raw.get("distance") or []
    bands = raw.get("bands") or {}
    if not distances or not bands:
        return

    headers = ["k_distance"]
    series = []
    for spin_name, spin_bands in bands.items():
        for idx, band in enumerate(spin_bands, start=1):
            headers.append(f"{spin_name}_band_{idx}")
            series.append(list(band))

    csv_path = output_dir / "band_structure_data.csv"
    with csv_path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for i, distance in enumerate(distances):
            row = [distance]
            for band in series:
                row.append(band[i] if i < len(band) else "")
            writer.writerow(row)
