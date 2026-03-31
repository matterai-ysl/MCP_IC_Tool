"""
BandStructureAnalyzer — 从VASP能带结构计算输出中提取数据并生成报告。

依赖 pymatgen 解析 vasprun.xml 获取能带结构数据，
使用 matplotlib 绘制能带图，生成HTML分析报告。
"""

import base64
import logging
import re
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, Optional

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator

logger = logging.getLogger(__name__)


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

            self.data['band_structure'] = {
                'band_gap': band_gap,
                'is_direct': is_direct,
                'transition': transition,
                'vbm': vbm_info['energy'] if vbm_info else None,
                'cbm': cbm_info['energy'] if cbm_info else None,
                'is_metal': bs.is_metal(),
                'nb_bands': bs.nb_bands,
            }

            self.data['final_results']['band_gap'] = band_gap
            self.data['final_results']['is_direct'] = is_direct
            self.data['final_results']['vbm'] = vbm_info['energy'] if vbm_info else None
            self.data['final_results']['cbm'] = cbm_info['energy'] if cbm_info else None
            self.data['final_results']['is_metal'] = bs.is_metal()

            # 提取费米能级和总能量
            self.data['final_results']['fermi_energy'] = vasprun.efermi
            self.data['final_results']['total_energy'] = vasprun.final_energy

            # 保存 band structure 对象供绘图使用
            self._bs = bs
            self._vasprun = vasprun

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
            fig = plotter.get_plot(ylim=(-5, 5))
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
        vbm = results.get('vbm', 'N/A')
        cbm = results.get('cbm', 'N/A')
        fermi = results.get('fermi_energy', 'N/A')
        total_e = results.get('total_energy', 'N/A')

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

        return f"""
        <div class="card">
            <h2>能带结构图</h2>
            <div style="text-align: center;">
                <img src="data:image/png;base64,{bs_plot}"
                     alt="Band Structure" style="max-width: 100%; height: auto;" />
            </div>
            <p style="text-align: center; color: #666; margin-top: 10px;">
                能带结构 E(k) 色散关系。红色虚线为费米能级。
            </p>
        </div>
        """

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


def generate_band_structure_report(work_dir: str, task_id: str = "band") -> Optional[str]:
    """生成能带结构分析 HTML 报告并返回文件路径。"""
    try:
        analyzer = BandStructureAnalyzer(work_dir, task_id=task_id)
        data = analyzer.analyze()

        output_dir = Path(work_dir) / "BS_output"
        output_dir.mkdir(exist_ok=True)
        output_path = str(output_dir / "band_structure_report.html")

        generator = BandStructureHTMLGenerator(data)
        generator.generate_html_report(output_path)

        logger.info("能带结构 HTML 报告已保存: %s", output_path)
        return output_path

    except Exception as e:
        logger.error("生成能带结构报告失败: %s", e)
        return None
