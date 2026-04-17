"""
NEBAnalyzer — 从 VASP NEB 计算各图像目录的 OUTCAR 中提取能量曲线与势垒信息。
"""

import base64
import logging
import re
from io import BytesIO
from pathlib import Path
from typing import Any, Dict, List, Optional

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator

logger = logging.getLogger(__name__)


def _parse_outcar_final_energy(outcar_path: Path) -> Optional[float]:
    """从 OUTCAR 中解析最终 TOTEN 能量（eV）。"""
    if not outcar_path.exists():
        return None
    try:
        energy = None
        with open(outcar_path, 'r', errors='ignore') as f:
            for line in f:
                if 'TOTEN' in line:
                    m = re.search(r'TOTEN\s*=\s*([-\d.]+)', line)
                    if m:
                        energy = float(m.group(1))
        return energy
    except Exception:
        return None


class NEBAnalyzer(BaseAnalyzer):
    """解析 VASP NEB 各图像 OUTCAR，提取能量曲线与势垒高度。"""

    def __init__(self, input_path: str, task_id: Optional[str] = None):
        super().__init__(input_path, task_id=task_id or "neb")

    def _init_data(self) -> Dict[str, Any]:
        return {
            'file_info': {},
            'calculation_settings': {},
            'final_results': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {},
            'neb': {},
            'visualizations': {},
        }

    def _do_analyze(self):
        self._extract_image_energies()
        self._generate_profile_plot()

    def _extract_image_energies(self):
        """遍历数字命名子目录（00, 01, …）读取各图像能量。"""
        image_dirs = sorted(
            [d for d in self.work_dir.iterdir() if d.is_dir() and d.name.isdigit()],
            key=lambda d: int(d.name),
        )
        if not image_dirs:
            self.data['neb'] = {'error': '未找到 NEB 图像目录（00, 01, ...）'}
            return

        energies: List[Optional[float]] = []
        for img_dir in image_dirs:
            e = _parse_outcar_final_energy(img_dir / "OUTCAR")
            energies.append(e)
            logger.info("NEB 图像 %s: E = %s eV", img_dir.name, e)

        valid = [(i, e) for i, e in enumerate(energies) if e is not None]
        if len(valid) < 2:
            self.data['neb'] = {'error': '有效能量数据不足', 'raw_energies': energies}
            return

        e0 = valid[0][1]
        rel_energies = [e - e0 for _, e in valid]
        max_rel = max(rel_energies)
        ts_idx = rel_energies.index(max_rel)

        self.data['neb'] = {
            'n_images': len(image_dirs),
            'absolute_energies': [e for _, e in valid],
            'relative_energies': rel_energies,
            'forward_barrier_eV': max_rel,
            'backward_barrier_eV': max_rel - rel_energies[-1],
            'reaction_energy_eV': rel_energies[-1],
            'ts_image_index': ts_idx,
        }
        self.data['final_results'].update({
            'forward_barrier_eV': max_rel,
            'backward_barrier_eV': max_rel - rel_energies[-1],
            'reaction_energy_eV': rel_energies[-1],
            'n_images': len(image_dirs),
        })

    def _generate_profile_plot(self):
        """生成 NEB 能量曲线 PNG，存为 base64。"""
        neb = self.data.get('neb', {})
        if 'error' in neb or 'relative_energies' not in neb:
            return
        try:
            import matplotlib
            matplotlib.use('Agg')
            import matplotlib.pyplot as plt

            rel_e = neb['relative_energies']
            x = list(range(len(rel_e)))
            ts_idx = neb['ts_image_index']

            fig, ax = plt.subplots(figsize=(8, 5))
            ax.plot(x, rel_e, 'o-', color='steelblue', linewidth=2, markersize=8)
            ax.plot(ts_idx, rel_e[ts_idx], 'r*', markersize=16,
                    label=f'TS: {rel_e[ts_idx]:.3f} eV')
            ax.axhline(0, color='gray', linestyle='--', linewidth=0.8)
            ax.set_xlabel('Image Index')
            ax.set_ylabel('Relative Energy (eV)')
            ax.set_title('NEB Energy Profile')
            ax.legend()
            ax.grid(True, alpha=0.3)

            buf = BytesIO()
            plt.savefig(buf, format='png', dpi=100, bbox_inches='tight')
            plt.close(fig)
            buf.seek(0)
            self.data['visualizations']['neb_profile'] = base64.b64encode(buf.read()).decode()
        except Exception as e:
            logger.warning("生成 NEB 能量曲线图失败: %s", e)


class NEBHTMLGenerator(BaseHTMLGenerator):
    report_title = "NEB 过渡态分析报告"
    report_emoji = "⚡"
    report_heading = "NEB 过渡态能量曲线"
    report_attribution = "由 VASP NEB 分析模块生成"
    header_gradient = "linear-gradient(135deg, #f093fb 0%, #f5576c 100%)"
    use_chartjs = False

    def _generate_body_sections(self) -> str:
        neb = self.data.get('neb', {})
        if 'error' in neb:
            return (f'<div class="section"><p style="color:red">'
                    f'NEB 分析错误: {neb["error"]}</p></div>')

        def fmt(v):
            return f'{v:.4f} eV' if isinstance(v, float) else 'N/A'

        fwd = neb.get('forward_barrier_eV')
        bwd = neb.get('backward_barrier_eV')
        rxn = neb.get('reaction_energy_eV')

        html = f"""
        <div class="section">
          <h2>NEB 计算结果摘要</h2>
          <table class="data-table">
            <tr><th>参数</th><th>值</th></tr>
            <tr><td>图像数（含端点）</td><td>{neb.get("n_images", "N/A")}</td></tr>
            <tr><td>正向势垒（活化能）</td><td>{fmt(fwd)}</td></tr>
            <tr><td>逆向势垒</td><td>{fmt(bwd)}</td></tr>
            <tr><td>反应能（ΔE）</td><td>{fmt(rxn)}</td></tr>
            <tr><td>过渡态图像序号</td><td>{neb.get("ts_image_index", "N/A")}</td></tr>
          </table>
        </div>
        """

        if isinstance(fwd, float):
            if fwd < 0.3:
                barrier_text = f"正向能垒约 {fwd:.3f} eV，通常说明这条迁移/反应路径相对容易发生。"
            elif fwd < 0.8:
                barrier_text = f"正向能垒约 {fwd:.3f} eV，说明这条路径可以发生，但往往需要一定热激发。"
            else:
                barrier_text = f"正向能垒约 {fwd:.3f} eV，说明这条路径相对不容易发生。"
        else:
            barrier_text = "这次没有得到可靠的势垒高度，暂时无法判断反应或扩散是否容易发生。"

        if isinstance(rxn, float):
            if rxn < 0:
                reaction_text = f"反应能约 {rxn:.3f} eV，表示产物侧整体更低能、热力学上更有利。"
            elif rxn > 0:
                reaction_text = f"反应能约 +{rxn:.3f} eV，表示产物侧更高能、热力学上不如初态有利。"
            else:
                reaction_text = "反应能接近 0，说明初态和末态能量非常接近。"
        else:
            reaction_text = "这次没有得到可靠的初末态能量差。"

        html += f"""
        <div class="section">
          <h2>🧭 通俗解读</h2>
          <div class="summary-card">
            <p><strong>NEB 最重要的是看“翻山”有多难。</strong> 这里的势垒可以理解成体系从初态走到末态时需要跨过的能量山峰。</p>
            <ul>
              <li><strong>能垒越低</strong>，通常说明这条迁移或反应路径<strong>更容易发生</strong>；能垒越高，则往往更难发生。</li>
              <li>{barrier_text}</li>
              <li>{reaction_text}</li>
            </ul>
          </div>
        </div>
        """

        profile_img = self.data.get('visualizations', {}).get('neb_profile')
        if profile_img:
            png_dl = self._png_download_link('neb_profile.png', profile_img)
            html += f"""
        <div class="section">
          <h2>能量曲线图</h2>
          <img src="data:image/png;base64,{profile_img}" style="max-width:100%;">
          <div style="margin-top: 12px; padding: 12px 14px; background: #f7fafc; border-left: 4px solid #4299e1; border-radius: 6px; color: #2d3748;">
            <strong>怎么看这张图：</strong> 横轴是反应路径上的图像序号，纵轴是相对能量；曲线最高点通常就是最难跨过去的能量山峰，首尾的高度差则反映初态和末态谁更稳定。
          </div>
          {png_dl}
        </div>
            """

        rel_e = neb.get('relative_energies', [])
        if rel_e:
            html += self._csv_download_link(
                'neb_energies.csv',
                ['image_index', 'relative_energy_eV'],
                [[str(i), f'{e:.6f}'] for i, e in enumerate(rel_e)],
                label='下载能量数据 (CSV)',
            )

        return html

    def _generate_javascript(self) -> str:
        return ""


def generate_neb_report(work_dir: str, task_id: str = "neb",
                        output_dir: Optional[str] = None) -> str:
    """生成 NEB 分析 HTML 报告，返回 HTML 文件路径。"""
    import os

    analyzer = NEBAnalyzer(work_dir, task_id=task_id)
    data = analyzer.analyze()

    out_dir = output_dir or work_dir
    Path(out_dir).mkdir(parents=True, exist_ok=True)

    generator = NEBHTMLGenerator(data)
    html_path = os.path.join(out_dir, f"neb_report_{task_id}.html")
    generator.generate_html_report(html_path)
    return html_path
