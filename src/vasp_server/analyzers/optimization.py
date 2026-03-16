"""
VASP 结构优化分析器 — 从 optimization_analyzer.py 迁移。

包含:
- OUTCARAnalyzer: 解析 OUTCAR 并提取优化过程数据
- OptimizationHTMLGenerator: 生成 HTML 可视化报告
- generate_optimization_report: 入口便捷函数
"""

import re
import math
import os
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator


class OUTCARAnalyzer(BaseAnalyzer):
    """VASP 结构优化计算结果分析器。"""

    def _init_data(self) -> Dict[str, Any]:
        return {
            'file_info': {},
            'calculation_settings': {},
            'ionic_steps': [],
            'convergence_analysis': {},
            'final_results': {},
            'electronic_structure': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {},
        }

    def _do_analyze(self):
        content = self._get_outcar_content()
        if not content:
            return
        self._parse_ionic_steps(content)
        self._parse_electronic_structure(content)
        self._analyze_convergence()
        self._analyze_optimization_process()
        self._analyze_final_results()

    # ------------------------------------------------------------------ #
    #  离子步解析
    # ------------------------------------------------------------------ #
    def _parse_ionic_steps(self, content: str):
        ionic_steps: List[Dict[str, Any]] = []

        # 1) SCF 能量序列
        scf_energy_series: List[float] = []
        energy_candidates: List[Tuple[int, float, str]] = []
        energy_pattern = r'free\s+energy\s+TOTEN\s*=\s*([\-\d\.]+)\s*eV'
        energy_matches = list(re.finditer(energy_pattern, content))
        if len(energy_matches) == 0:
            alt_patterns = [
                r'energy\s*without\s*entropy[^=]*=\s*([\-\d\.]+)',
                r'TOTEN\s*=\s*([\-\d\.]+)',
                r'total\s*energy\s*=\s*([\-\d\.]+)',
            ]
            for pattern in alt_patterns:
                alt_matches = list(re.finditer(pattern, content, re.IGNORECASE))
                if alt_matches:
                    energy_matches = alt_matches
                    break

        for m in energy_matches:
            try:
                val = float(m.group(1))
                scf_energy_series.append(val)
                energy_candidates.append((m.start(), val, 'SCF'))
            except Exception:
                continue
        energy_candidates.sort(key=lambda x: x[0])

        # 2) TOTAL-FORCE 块计数离子步
        force_blocks = list(re.finditer(
            r'POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*\n\s*-{10,}\s*\n(.*?)\n\s*-{10,}',
            content, re.DOTALL,
        ))

        def find_prev_energy_value(before_pos: int) -> Optional[float]:
            left = [e for e in energy_candidates if e[0] < before_pos]
            return left[-1][1] if left else None

        step_num = 0
        for fm in force_blocks:
            block_start = fm.start()
            block = fm.group(1)
            lines = block.strip().split('\n') if block else []
            forces: List[List[float]] = []
            positions: List[List[float]] = []

            for raw in lines:
                line = raw.strip()
                if not line or line.startswith('-'):
                    continue
                if line.lower().startswith('total drift'):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    px, py, pz = float(parts[0]), float(parts[1]), float(parts[2])
                    fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                    positions.append([px, py, pz])
                    forces.append([fx, fy, fz])
                except Exception:
                    continue

            max_force_val: Optional[float] = None
            rms_force_val: Optional[float] = None
            if forces:
                mags = [math.sqrt(fx*fx + fy*fy + fz*fz) for fx, fy, fz in forces]
                if mags:
                    max_force_val = max(mags)
                    rms_force_val = math.sqrt(sum(m*m for m in mags) / len(mags))

            ref_energy = find_prev_energy_value(block_start) or 0.0

            step_num += 1
            ionic_steps.append({
                'step': step_num,
                'energy': float(ref_energy),
                'forces': forces,
                'max_force': max_force_val,
                'rms_force': rms_force_val,
                'positions': positions,
            })

        self.data['ionic_steps'] = ionic_steps
        self.data['scf_energies'] = scf_energy_series

    # ------------------------------------------------------------------ #
    #  电子结构
    # ------------------------------------------------------------------ #
    def _parse_electronic_structure(self, content: str):
        electronic_data: Dict[str, Any] = {}
        gap_match = re.search(r'E-fermi\s*:\s*([-\d\.]+).*?gap\s*=\s*([\d\.]+)', content)
        if gap_match:
            electronic_data['fermi_energy'] = float(gap_match.group(1))
            electronic_data['band_gap'] = float(gap_match.group(2))
        mag_match = re.search(r'number of electron\s+NELECT\s*=\s*([\d\.]+)', content)
        if mag_match:
            electronic_data['total_electrons'] = float(mag_match.group(1))
        self.data['electronic_structure'] = electronic_data

    # ------------------------------------------------------------------ #
    #  收敛性分析
    # ------------------------------------------------------------------ #
    def _analyze_convergence(self):
        convergence: Dict[str, Any] = {
            'force_convergence': {},
            'energy_convergence': {},
            'electronic_convergence': {},
            'tail_check': {},
            'overall_convergence': False,
        }

        # 尾部标志检查
        tail_info: Dict[str, Any] = {'matched': False, 'keywords': ['reached required accuracy', 'Voluntary']}
        try:
            if self.outcar_path:
                file_size = self.outcar_path.stat().st_size
                tail_bytes = min(1024, file_size)
                with open(self.outcar_path, 'rb') as f:
                    f.seek(-tail_bytes, os.SEEK_END)
                    last_lines = f.readlines()[-10:]
                    last_content = b''.join(last_lines).decode('utf-8', errors='ignore')
                    lc = last_content.lower()
                    tail_info['matched'] = ('reached required accuracy' in lc) or ('voluntary' in lc)
        except Exception as e:
            tail_info['exception'] = str(e)
        convergence['tail_check'] = tail_info

        overall_converged_by_tail = tail_info.get('matched', False)

        ionic_steps = self.data.get('ionic_steps', [])
        if not ionic_steps:
            convergence['overall_convergence'] = overall_converged_by_tail
            self.data['convergence_analysis'] = convergence
            return

        # 力收敛
        ediffg = self.data['calculation_settings'].get('EDIFFG', -0.01)
        force_threshold = abs(float(ediffg))
        max_forces = [s.get('max_force') for s in ionic_steps if s.get('max_force') is not None]
        if max_forces:
            final_max_force = max_forces[-1]
            force_converged = final_max_force < force_threshold
            convergence['force_convergence'] = {
                'converged': force_converged,
                'threshold': force_threshold,
                'final_max_force': final_max_force,
                'force_history': max_forces,
            }

        # 能量收敛
        scf_energies = self.data.get('scf_energies', [])
        if len(scf_energies) > 1:
            energy_changes = [abs(scf_energies[i] - scf_energies[i-1]) for i in range(1, len(scf_energies))]
            recent_changes = energy_changes[-5:] if len(energy_changes) >= 5 else energy_changes
            avg_recent_change = sum(recent_changes) / len(recent_changes) if recent_changes else 0.0
            energy_converged = avg_recent_change < 1e-4
            convergence['energy_convergence'] = {
                'converged': energy_converged,
                'final_energy': scf_energies[-1],
                'energy_history': scf_energies,
                'energy_changes': energy_changes,
                'avg_recent_change': avg_recent_change,
            }

        # 综合判断
        force_ok = convergence['force_convergence'].get('converged', False)
        energy_ok = convergence['energy_convergence'].get('converged', False)
        convergence['overall_convergence'] = overall_converged_by_tail or (force_ok and energy_ok)
        self.data['convergence_analysis'] = convergence

    # ------------------------------------------------------------------ #
    #  优化过程分析
    # ------------------------------------------------------------------ #
    def _analyze_optimization_process(self):
        if not self.data['ionic_steps']:
            return

        total_steps = len(self.data['ionic_steps'])
        max_steps = self.data['calculation_settings'].get('NSW', 500)
        is_converged = self.data.get('convergence_analysis', {}).get('overall_convergence', False)
        completion_ratio = 1.0 if is_converged else (total_steps / max_steps if max_steps > 0 else 0.0)

        process_data: Dict[str, Any] = {
            'total_steps': total_steps,
            'max_steps': max_steps,
            'completion_ratio': completion_ratio,
            'is_converged': is_converged,
            'energy_profile': [],
            'force_profile': [],
            'structure_changes': [],
        }

        for step in self.data['ionic_steps']:
            process_data['energy_profile'].append({'step': step['step'], 'energy': step['energy']})
            if step.get('max_force') is not None:
                process_data['force_profile'].append({
                    'step': step['step'],
                    'max_force': step.get('max_force'),
                    'rms_force': step.get('rms_force'),
                })

        # 结构变化
        for idx in range(1, len(self.data['ionic_steps'])):
            cur = self.data['ionic_steps'][idx].get('positions', [])
            prev = self.data['ionic_steps'][idx - 1].get('positions', [])
            if cur and prev and len(cur) == len(prev):
                displacements = [
                    math.sqrt(sum((c - p) ** 2 for c, p in zip(ca, pa)))
                    for ca, pa in zip(cur, prev)
                ]
                process_data['structure_changes'].append({
                    'step': self.data['ionic_steps'][idx]['step'],
                    'max_displacement': max(displacements) if displacements else 0.0,
                    'avg_displacement': sum(displacements) / len(displacements) if displacements else 0.0,
                })

        self.data['optimization_process'] = process_data

    # ------------------------------------------------------------------ #
    #  最终结果
    # ------------------------------------------------------------------ #
    def _analyze_final_results(self):
        if not self.data['ionic_steps']:
            return

        final_step = self.data['ionic_steps'][-1]
        scf_energies = self.data.get('scf_energies', [])
        final_energy = scf_energies[-1] if scf_energies else final_step['energy']

        final_results: Dict[str, Any] = {
            'final_energy': final_energy,
            'final_step': final_step['step'],
            'final_forces': final_step.get('forces', []),
            'final_positions': final_step.get('positions', []),
            'max_residual_force': final_step.get('max_force'),
            'rms_residual_force': final_step.get('rms_force'),
        }

        if final_step.get('forces'):
            mags = [math.sqrt(f[0]**2 + f[1]**2 + f[2]**2) for f in final_step['forces']]
            if mags:
                n = len(mags)
                mean_mag = sum(mags) / n
                variance = sum((m - mean_mag) ** 2 for m in mags) / n
                final_results['force_statistics'] = {
                    'mean_force_magnitude': mean_mag,
                    'std_force_magnitude': math.sqrt(variance),
                    'min_force_magnitude': min(mags),
                    'max_force_magnitude': max(mags),
                }

        self.data['final_results'] = final_results


# ====================================================================== #
#  HTML 报告生成器
# ====================================================================== #
class OptimizationHTMLGenerator(BaseHTMLGenerator):
    """结构优化 HTML 报告生成器。"""

    report_title = "VASP结构优化分析报告"
    report_emoji = "🔬"
    report_heading = "VASP结构优化分析报告"
    report_attribution = "由VASP API结构优化可视化分析模块生成"
    header_gradient = "linear-gradient(135deg, #667eea 0%, #764ba2 100%)"

    def _filter_contcar_velocities(self, lines: List[str]) -> List[str]:
        filtered: List[str] = []
        for i, line in enumerate(lines):
            if line.strip() == '' and i + 1 < len(lines):
                next_line = lines[i + 1].strip()
                if ('0.00000000E+00' in next_line or
                        next_line.replace(' ', '').replace('0', '').replace('.', '').replace('E', '').replace('+', '') == ''):
                    filtered.append(line)
                    break
                else:
                    filtered.append(line)
            else:
                filtered.append(line)
        return filtered

    def _generate_body_sections(self) -> str:
        return (
            self._generate_summary()
            + self._generate_convergence_section()
            + self._generate_optimization_process_section()
            + self._generate_final_results_section()
            + self._generate_structure_comparison_section()
        )

    def _generate_summary(self) -> str:
        convergence = self.data.get('convergence_analysis', {})
        process = self.data.get('optimization_process', {})
        final = self.data.get('final_results', {})

        overall_converged = convergence.get('overall_convergence', False)
        convergence_class = 'converged' if overall_converged else 'not-converged'
        convergence_text = '收敛' if overall_converged else '未收敛'

        return f"""
        <div class="section">
            <h2>📊 计算摘要</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>收敛状态</h3>
                    <div class="value">
                        <span class="convergence-status {convergence_class}">{convergence_text}</span>
                    </div>
                </div>
                <div class="summary-card">
                    <h3>离子步数</h3>
                    <div class="value">{process.get('total_steps', 0)} / {process.get('max_steps', 0)}</div>
                </div>
                <div class="summary-card">
                    <h3>最终能量</h3>
                    <div class="value">{final.get('final_energy', 0):.6f} eV</div>
                </div>
                <div class="summary-card">
                    <h3>最大剩余力</h3>
                    <div class="value">{(final.get('max_residual_force') or 0):.4f} eV/Å</div>
                </div>
            </div>
        </div>
        """

    def _generate_convergence_section(self) -> str:
        convergence = self.data.get('convergence_analysis', {})
        force_conv = convergence.get('force_convergence', {})
        energy_conv = convergence.get('energy_convergence', {})
        force_status = '✅ 收敛' if force_conv.get('converged', False) else '❌ 未收敛'
        energy_status = '✅ 收敛' if energy_conv.get('converged', False) else '❌ 未收敛'

        force_history = force_conv.get('force_history', [])
        energy_history = energy_conv.get('energy_history', self.data.get('scf_energies', []))

        force_dl = self._csv_download_link(
            'force_convergence.csv',
            ['离子步', '最大力(eV/Å)'],
            [[i + 1, v] for i, v in enumerate(force_history)],
        )
        energy_dl = self._csv_download_link(
            'scf_energy_convergence.csv',
            ['SCF步', '总能量(eV)'],
            [[i + 1, v] for i, v in enumerate(energy_history)],
        )

        return f"""
        <div class="section">
            <h2>🎯 收敛性检查</h2>
            <div class="grid-2">
                <div>
                    <h3>力收敛分析</h3>
                    <p><strong>状态:</strong> {force_status}</p>
                    <p><strong>收敛阈值:</strong> {(force_conv.get('threshold') or 0):.4f} eV/Å</p>
                    <p><strong>最终最大力:</strong> {(force_conv.get('final_max_force') or 0):.4f} eV/Å</p>
                    <div class="chart-container"><canvas id="forceChart"></canvas></div>
                    {force_dl}
                </div>
                <div>
                    <h3>能量收敛分析</h3>
                    <p><strong>状态:</strong> {energy_status}</p>
                    <p><strong>最终能量:</strong> {(energy_conv.get('final_energy') or 0):.6f} eV</p>
                    <p><strong>平均能量变化:</strong> {(energy_conv.get('avg_recent_change') or 0):.2e} eV</p>
                    <div class="chart-container"><canvas id="energyChart"></canvas></div>
                    {energy_dl}
                </div>
            </div>
        </div>
        """

    def _generate_optimization_process_section(self) -> str:
        process = self.data.get('optimization_process', {})
        is_converged = process.get('is_converged', False)
        completion_ratio = process.get('completion_ratio', 0)

        if is_converged:
            status_text, status_class = "已完成（收敛）", "converged"
        elif completion_ratio >= 0.95:
            status_text, status_class = f"{completion_ratio*100:.1f}% - 接近完成", "not-converged"
        else:
            status_text, status_class = f"{completion_ratio*100:.1f}% - 进行中", "not-converged"

        energy_profile = process.get('energy_profile', [])
        force_profile = process.get('force_profile', [])

        energy_evo_dl = self._csv_download_link(
            'energy_evolution.csv',
            ['离子步', '能量(eV)'],
            [[s['step'], s['energy']] for s in energy_profile],
        )
        force_evo_dl = self._csv_download_link(
            'force_evolution.csv',
            ['离子步', '最大力(eV/Å)', 'RMS力(eV/Å)'],
            [[s['step'], s.get('max_force'), s.get('rms_force')] for s in force_profile],
        )

        return f"""
        <div class="section">
            <h2>📈 优化过程监控</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>完成状态</h3>
                    <div class="value"><span class="convergence-status {status_class}">{status_text}</span></div>
                </div>
                <div class="summary-card">
                    <h3>步数进度</h3>
                    <div class="value">{process.get('total_steps', 0)} / {process.get('max_steps', 0)} 步</div>
                </div>
            </div>
            <div class="grid-2">
                <div>
                    <h3>能量演化曲线</h3>
                    <div class="chart-container"><canvas id="energyEvolutionChart"></canvas></div>
                    {energy_evo_dl}
                </div>
                <div>
                    <h3>力演化曲线</h3>
                    <div class="chart-container"><canvas id="forceEvolutionChart"></canvas></div>
                    {force_evo_dl}
                </div>
            </div>
        </div>
        """

    def _generate_final_results_section(self) -> str:
        final = self.data.get('final_results', {})
        force_stats = final.get('force_statistics', {})

        return f"""
        <div class="section">
            <h2>🏁 最终结果</h2>
            <div class="grid-2">
                <div>
                    <h3>能量信息</h3>
                    <table class="data-table">
                        <tr><td>最终总能量</td><td>{final.get('final_energy', 0):.6f} eV</td></tr>
                        <tr><td>完成步数</td><td>{final.get('final_step', 0)}</td></tr>
                    </table>
                </div>
                <div>
                    <h3>力统计信息</h3>
                    <table class="data-table">
                        <tr><td>最大剩余力</td><td>{(final.get('max_residual_force') or 0):.4f} eV/Å</td></tr>
                        <tr><td>RMS剩余力</td><td>{(final.get('rms_residual_force') or 0):.4f} eV/Å</td></tr>
                        <tr><td>平均力大小</td><td>{(force_stats.get('mean_force_magnitude') or 0):.4f} eV/Å</td></tr>
                        <tr><td>力标准差</td><td>{(force_stats.get('std_force_magnitude') or 0):.4f} eV/Å</td></tr>
                    </table>
                </div>
            </div>
        </div>
        """

    def _generate_structure_comparison_section(self) -> str:
        structure_files = self.data.get('structure_files', {})
        poscar_content = structure_files.get('poscar', '')
        contcar_content = structure_files.get('contcar', '')

        def format_content(content: str, title: str) -> str:
            if not content:
                return f"<p>❌ 未找到{title}文件</p>"
            lines = content.split('\n')
            if 'CONTCAR' in title:
                lines = self._filter_contcar_velocities(lines)
            if len(lines) > 50:
                display = '\n'.join(lines[:50]) + f'\n\n... (共{len(lines)}行，仅显示前50行)'
            else:
                display = '\n'.join(lines)
            return f"<pre><code>{display}</code></pre>"

        return f"""
        <div class="section">
            <h2>📋 结构文件对比</h2>
            <p>以下是初始结构(POSCAR)与优化后结构(CONTCAR)的对比：</p>
            <div class="grid-2">
                <div>
                    <h3>🏁 初始结构 (POSCAR)</h3>
                    <div style="background: #f8f9fa; padding: 15px; border-radius: 5px; overflow-x: auto; max-height: 600px; overflow-y: auto;">
                        {format_content(poscar_content, 'POSCAR')}
                    </div>
                </div>
                <div>
                    <h3>🎯 优化后结构 (CONTCAR)</h3>
                    <div style="background: #f8f9fa; padding: 15px; border-radius: 5px; overflow-x: auto; max-height: 600px; overflow-y: auto;">
                        {format_content(contcar_content, 'CONTCAR')}
                    </div>
                </div>
            </div>
            <div style="margin-top: 20px; padding: 15px; background: #e8f4f8; border-radius: 5px;">
                <p><strong>📝 说明：</strong></p>
                <ul>
                    <li>POSCAR包含计算开始时的初始原子坐标和晶格参数</li>
                    <li>CONTCAR包含结构优化完成后的最终原子坐标和晶格参数</li>
                    <li>通过对比可以观察到原子位置的变化和晶格的弛豫</li>
                    <li>如果显示行数过多，仅展示前50行内容</li>
                </ul>
            </div>
        </div>
        """

    def _generate_javascript(self) -> str:
        convergence = self.data.get('convergence_analysis', {})
        process = self.data.get('optimization_process', {})

        force_history = convergence.get('force_convergence', {}).get('force_history', [])
        energy_history = self.data.get('scf_energies', [])
        energy_profile = process.get('energy_profile', [])
        force_profile = process.get('force_profile', [])
        force_threshold = convergence.get('force_convergence', {}).get('threshold', 0.01)

        return f"""
        const forceCtx = document.getElementById('forceChart').getContext('2d');
        new Chart(forceCtx, {{
            type: 'line',
            data: {{
                labels: {list(range(1, len(force_history) + 1))},
                datasets: [{{
                    label: '最大力 (eV/Å)',
                    data: {force_history},
                    borderColor: 'rgb(255, 99, 132)',
                    backgroundColor: 'rgba(255, 99, 132, 0.2)',
                    tension: 0.1
                }}, {{
                    label: '收敛阈值',
                    data: Array({len(force_history)}).fill({force_threshold}),
                    borderColor: 'rgb(75, 192, 192)',
                    borderDash: [5, 5],
                    pointRadius: 0
                }}]
            }},
            options: {{
                responsive: true, maintainAspectRatio: false,
                scales: {{
                    y: {{ beginAtZero: false, title: {{ display: true, text: '力 (eV/Å)' }} }},
                    x: {{ title: {{ display: true, text: '离子步' }} }}
                }}
            }}
        }});

        const energyCtx = document.getElementById('energyChart').getContext('2d');
        new Chart(energyCtx, {{
            type: 'line',
            data: {{
                labels: {list(range(1, len(energy_history) + 1))},
                datasets: [{{
                    label: '总能量 (eV)',
                    data: {energy_history},
                    borderColor: 'rgb(54, 162, 235)',
                    backgroundColor: 'rgba(54, 162, 235, 0.2)',
                    tension: 0.1
                }}]
            }},
            options: {{
                responsive: true, maintainAspectRatio: false,
                scales: {{
                    y: {{ title: {{ display: true, text: '能量 (eV)' }} }},
                    x: {{ title: {{ display: true, text: 'SCF步' }} }}
                }}
            }}
        }});

        const energyEvolutionCtx = document.getElementById('energyEvolutionChart').getContext('2d');
        new Chart(energyEvolutionCtx, {{
            type: 'line',
            data: {{
                labels: {[step['step'] for step in energy_profile]},
                datasets: [{{
                    label: '总能量 (eV)',
                    data: {[step['energy'] for step in energy_profile]},
                    borderColor: 'rgb(153, 102, 255)',
                    backgroundColor: 'rgba(153, 102, 255, 0.2)',
                    tension: 0.1
                }}]
            }},
            options: {{
                responsive: true, maintainAspectRatio: false,
                scales: {{
                    y: {{ title: {{ display: true, text: '能量 (eV)' }} }},
                    x: {{ title: {{ display: true, text: '离子步（参考能量）' }} }}
                }}
            }}
        }});

        const forceEvolutionCtx = document.getElementById('forceEvolutionChart').getContext('2d');
        new Chart(forceEvolutionCtx, {{
            type: 'line',
            data: {{
                labels: {[step['step'] for step in force_profile]},
                datasets: [{{
                    label: '最大力 (eV/Å)',
                    data: {[step['max_force'] for step in force_profile]},
                    borderColor: 'rgb(255, 159, 64)',
                    backgroundColor: 'rgba(255, 159, 64, 0.2)',
                    tension: 0.1
                }}, {{
                    label: 'RMS力 (eV/Å)',
                    data: {[step['rms_force'] for step in force_profile]},
                    borderColor: 'rgb(201, 203, 207)',
                    backgroundColor: 'rgba(201, 203, 207, 0.2)',
                    tension: 0.1
                }}]
            }},
            options: {{
                responsive: true, maintainAspectRatio: false,
                scales: {{
                    y: {{ title: {{ display: true, text: '力 (eV/Å)' }} }},
                    x: {{ title: {{ display: true, text: '离子步' }} }}
                }}
            }}
        }});
        """


# ====================================================================== #
#  入口便捷函数
# ====================================================================== #
def generate_optimization_report(
    input_path: str,
    task_id: Optional[str] = None,
    output_dir: Optional[str] = None,
) -> str:
    """生成结构优化分析报告。"""
    try:
        analyzer = OUTCARAnalyzer(input_path, task_id)
        analysis_data = analyzer.analyze()

        output_dir_path = Path(output_dir) if output_dir else analyzer.work_dir
        output_file = output_dir_path / "optimization_analysis_report.html"

        generator = OptimizationHTMLGenerator(analysis_data)
        return generator.generate_html_report(str(output_file))
    except Exception as e:
        raise Exception(f"生成优化分析报告失败: {str(e)}")
