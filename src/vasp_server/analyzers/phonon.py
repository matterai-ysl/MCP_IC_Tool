"""
PhononAnalyzer — 从 VASP IBRION=6 声子计算 OUTCAR 中解析 Gamma 点声子频率。

使用 IBRION=6（有限位移法+对称性），VASP 直接在 OUTCAR 中输出声子本征频率。
支持识别虚频（动力学不稳定）并生成频率分布直方图。
"""

import base64
import logging
import re
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator

logger = logging.getLogger(__name__)


def _parse_phonon_frequencies(outcar_content: str) -> List[Tuple[float, bool]]:
    """
    从 OUTCAR 内容中提取声子频率列表。

    VASP 输出格式示例：
        2 f  =   14.520 THz     91.244 2PiTHz    484.267 cm-1     60.028 meV
        3 f/i=    0.003 THz      0.017 2PiTHz      0.090 cm-1      0.011 meV

    返回: [(freq_cm1, is_imaginary), ...]
    """
    # 匹配 "N f  = ..." 或 "N f/i= ..."
    pattern = re.compile(
        r'^\s*\d+\s+f(/i)?\s*=\s*[\d.]+\s+THz\s+[\d.]+\s+2PiTHz\s+([\d.]+)\s+cm-1',
        re.MULTILINE,
    )
    freqs = []
    for m in pattern.finditer(outcar_content):
        is_imaginary = m.group(1) == '/i'
        freq_cm1 = float(m.group(2))
        freqs.append((freq_cm1, is_imaginary))
    return freqs


class PhononAnalyzer(BaseAnalyzer):
    """解析 VASP 声子计算（IBRION=5/6）OUTCAR，提取 Gamma 点声子频率。"""

    def __init__(self, input_path: str, task_id: Optional[str] = None):
        super().__init__(input_path, task_id=task_id or "phonon")

    def _init_data(self) -> Dict[str, Any]:
        return {
            'file_info': {},
            'calculation_settings': {},
            'final_results': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {},
            'phonon': {},
            'visualizations': {},
        }

    def _do_analyze(self):
        self._extract_frequencies()
        self._generate_frequency_plot()

    def _extract_frequencies(self):
        """从 OUTCAR 中解析声子频率。"""
        content = self._get_outcar_content()
        if not content:
            self.data['phonon'] = {'error': 'OUTCAR 不存在或为空'}
            return

        freqs = _parse_phonon_frequencies(content)
        if not freqs:
            self.data['phonon'] = {'error': 'OUTCAR 中未找到声子频率（确认 IBRION=5 或 6）'}
            return

        real_freqs = [f for f, imag in freqs if not imag]
        imag_freqs = [f for f, imag in freqs if imag]
        all_cm1 = [f for f, _ in freqs]
        signs = [-f if imag else f for f, imag in freqs]  # 虚频取负值，方便绘图

        self.data['phonon'] = {
            'n_modes': len(freqs),
            'n_imaginary': len(imag_freqs),
            'dynamically_stable': len(imag_freqs) == 0,
            'max_imaginary_freq_cm1': max(imag_freqs) if imag_freqs else 0.0,
            'max_real_freq_cm1': max(real_freqs) if real_freqs else 0.0,
            'min_real_freq_cm1': min(real_freqs) if real_freqs else 0.0,
            'frequencies_cm1': signs,          # 虚频标为负值
            'frequencies_raw': all_cm1,        # 绝对值
        }
        self.data['final_results'].update({
            'n_modes': len(freqs),
            'n_imaginary': len(imag_freqs),
            'dynamically_stable': len(imag_freqs) == 0,
            'max_imaginary_freq_cm1': max(imag_freqs) if imag_freqs else 0.0,
            'max_real_freq_cm1': max(real_freqs) if real_freqs else 0.0,
        })

        stability = "稳定" if not imag_freqs else f"不稳定（{len(imag_freqs)} 个虚频）"
        logger.info("声子分析完成: %d 个模式, %s", len(freqs), stability)

    def _generate_frequency_plot(self):
        """生成声子频率分布直方图。"""
        phonon = self.data.get('phonon', {})
        if 'error' in phonon or 'frequencies_cm1' not in phonon:
            return
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            freqs = phonon['frequencies_cm1']
            n_imag = phonon['n_imaginary']

            fig, ax = plt.subplots(figsize=(9, 5))

            real_f = [f for f in freqs if f >= 0]
            imag_f = [abs(f) for f in freqs if f < 0]

            if real_f:
                ax.hist(real_f, bins=30, color='steelblue', alpha=0.8, label='实频')
            if imag_f:
                ax.hist(imag_f, bins=10, color='tomato', alpha=0.8, label=f'虚频 ({n_imag})')

            ax.set_xlabel('Frequency (cm⁻¹)')
            ax.set_ylabel('Count')
            ax.set_title('Phonon Frequency Distribution (Γ point)')
            ax.legend()
            ax.grid(True, alpha=0.3)

            stability_label = 'Dynamically Stable' if n_imag == 0 else f'UNSTABLE: {n_imag} imaginary modes'
            color = 'green' if n_imag == 0 else 'red'
            ax.text(0.98, 0.95, stability_label, transform=ax.transAxes,
                    ha='right', va='top', color=color, fontsize=11, fontweight='bold')

            buf = BytesIO()
            plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            self.data['visualizations']['phonon_hist'] = base64.b64encode(buf.read()).decode()
        except Exception as e:
            logger.warning("生成声子频率图失败: %s", e)


class PhononHTMLGenerator(BaseHTMLGenerator):
    report_title = "声子计算分析报告"
    report_emoji = "🎵"
    report_heading = "声子频率分析（Γ 点）"
    report_attribution = "由 VASP 声子分析模块生成"
    header_gradient = "linear-gradient(135deg, #43e97b 0%, #38f9d7 100%)"
    use_chartjs = False

    def _generate_body_sections(self) -> str:
        phonon = self.data.get('phonon', {})
        if 'error' in phonon:
            return (f'<div class="section"><p style="color:red">'
                    f'声子分析错误: {phonon["error"]}</p></div>')

        stable = phonon.get('dynamically_stable', True)
        stability_badge = (
            '<span class="convergence-status converged">动力学稳定</span>'
            if stable else
            f'<span class="convergence-status not-converged">'
            f'动力学不稳定（{phonon.get("n_imaginary", 0)} 个虚频）</span>'
        )

        html = f"""
        <div class="section">
          <h2>声子计算结果摘要</h2>
          <p>动力学稳定性: {stability_badge}</p>
          <table class="data-table">
            <tr><th>参数</th><th>值</th></tr>
            <tr><td>声子模式总数</td><td>{phonon.get("n_modes", "N/A")}</td></tr>
            <tr><td>虚频数量</td><td>{phonon.get("n_imaginary", "N/A")}</td></tr>
            <tr><td>最大虚频</td><td>{phonon.get("max_imaginary_freq_cm1", 0):.2f} cm⁻¹</td></tr>
            <tr><td>最高实频</td><td>{phonon.get("max_real_freq_cm1", 0):.2f} cm⁻¹</td></tr>
            <tr><td>最低实频</td><td>{phonon.get("min_real_freq_cm1", 0):.2f} cm⁻¹</td></tr>
          </table>
        </div>
        """

        hist_img = self.data.get('visualizations', {}).get('phonon_hist')
        if hist_img:
            html += f"""
        <div class="section">
          <h2>声子频率分布图</h2>
          <img src="data:image/png;base64,{hist_img}" style="max-width:100%;">
        </div>
            """

        freqs = phonon.get('frequencies_cm1', [])
        if freqs:
            html += self._csv_download_link(
                'phonon_frequencies.csv',
                ['mode_index', 'frequency_cm1', 'type'],
                [[str(i + 1), f'{abs(f):.4f}', 'imaginary' if f < 0 else 'real']
                 for i, f in enumerate(freqs)],
                label='下载频率数据 (CSV)',
            )

        return html

    def _generate_javascript(self) -> str:
        return ""


def generate_phonon_report(work_dir: str, task_id: str = "phonon",
                           output_dir: Optional[str] = None) -> str:
    """生成声子分析 HTML 报告，返回 HTML 文件路径。"""
    import os

    analyzer = PhononAnalyzer(work_dir, task_id=task_id)
    data = analyzer.analyze()

    out_dir = output_dir or work_dir
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    generator = PhononHTMLGenerator(data)
    html_path = os.path.join(out_dir, f"phonon_report_{task_id}.html")
    generator.generate_html_report(html_path)
    return html_path
