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
import csv
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator


class OUTCARAnalyzer(BaseAnalyzer):
    """VASP 结构优化计算结果分析器。"""

    @staticmethod
    def _coerce_int(value: Any, default: int) -> int:
        try:
            return int(float(value))
        except (TypeError, ValueError):
            return default

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
        tail_info: Dict[str, Any] = {'matched': False, 'keywords': ['reached required accuracy']}
        try:
            if self.outcar_path:
                file_size = self.outcar_path.stat().st_size
                tail_bytes = min(1024, file_size)
                with open(self.outcar_path, 'rb') as f:
                    f.seek(-tail_bytes, os.SEEK_END)
                    last_lines = f.readlines()[-10:]
                    last_content = b''.join(last_lines).decode('utf-8', errors='ignore')
                    lc = last_content.lower()
                    tail_info['matched'] = 'reached required accuracy' in lc
        except Exception as e:
            tail_info['exception'] = str(e)
        convergence['tail_check'] = tail_info

        overall_converged_by_tail = tail_info.get('matched', False)

        ionic_steps = self.data.get('ionic_steps', [])
        ediffg = self.data['calculation_settings'].get('EDIFFG', -0.01)
        force_threshold = abs(float(ediffg))

        if not ionic_steps:
            convergence['force_convergence'] = {
                'available': False,
                'converged': None,
                'threshold': force_threshold,
                'final_max_force': None,
                'force_history': [],
                'reason': 'missing_force_history',
            }
        # 力收敛
        max_forces = [s.get('max_force') for s in ionic_steps if s.get('max_force') is not None]
        if max_forces:
            final_max_force = max_forces[-1]
            force_converged = final_max_force < force_threshold
            convergence['force_convergence'] = {
                'available': True,
                'converged': force_converged,
                'threshold': force_threshold,
                'final_max_force': final_max_force,
                'force_history': max_forces,
            }
        elif ionic_steps:
            convergence['force_convergence'] = {
                'available': False,
                'converged': None,
                'threshold': force_threshold,
                'final_max_force': None,
                'force_history': [],
                'reason': 'missing_force_history',
            }

        # 能量收敛
        ediff = self.data['calculation_settings'].get('EDIFF', 1e-4)
        try:
            energy_threshold = abs(float(ediff))
        except (TypeError, ValueError):
            energy_threshold = 1e-4
        scf_energies = self.data.get('scf_energies', [])
        if len(scf_energies) > 1:
            energy_changes = [abs(scf_energies[i] - scf_energies[i-1]) for i in range(1, len(scf_energies))]
            recent_changes = energy_changes[-5:] if len(energy_changes) >= 5 else energy_changes
            avg_recent_change = sum(recent_changes) / len(recent_changes) if recent_changes else 0.0
            energy_converged = avg_recent_change < energy_threshold
            convergence['energy_convergence'] = {
                'available': True,
                'converged': energy_converged,
                'threshold': energy_threshold,
                'final_energy': scf_energies[-1],
                'energy_history': scf_energies,
                'energy_changes': energy_changes,
                'avg_recent_change': avg_recent_change,
            }
        elif len(scf_energies) == 1:
            convergence['energy_convergence'] = {
                'available': False,
                'converged': None,
                'reason': 'insufficient_energy_history',
                'threshold': energy_threshold,
                'final_energy': scf_energies[-1],
                'energy_history': scf_energies,
                'energy_changes': [],
                'avg_recent_change': None,
            }
        else:
            convergence['energy_convergence'] = {
                'available': False,
                'converged': None,
                'reason': 'missing_energy_history',
                'threshold': energy_threshold,
                'final_energy': None,
                'energy_history': [],
                'energy_changes': [],
                'avg_recent_change': None,
            }

        # 综合判断
        force_info = convergence['force_convergence']
        force_available = force_info.get('available', len(force_info.get('force_history', [])) > 0)
        force_ok = bool(force_info.get('converged', False)) if force_available else False
        energy_info = convergence['energy_convergence']
        energy_available = energy_info.get('available', len(energy_info.get('energy_history', [])) > 1)
        energy_ok = bool(energy_info.get('converged', False)) if energy_available else True
        convergence['overall_convergence'] = overall_converged_by_tail or (force_ok and energy_ok)
        self.data['convergence_analysis'] = convergence

    # ------------------------------------------------------------------ #
    #  优化过程分析
    # ------------------------------------------------------------------ #
    def _analyze_optimization_process(self):
        if not self.data['ionic_steps']:
            max_steps = self._coerce_int(self.data['calculation_settings'].get('NSW', 0), 0)
            self.data['optimization_process'] = {
                'total_steps': 0,
                'max_steps': max_steps,
                'completion_ratio': 0.0,
                'is_converged': self.data.get('convergence_analysis', {}).get('overall_convergence', False),
                'data_available': False,
                'energy_profile': [],
                'force_profile': [],
                'structure_changes': [],
            }
            return

        total_steps = len(self.data['ionic_steps'])
        max_steps = self._coerce_int(self.data['calculation_settings'].get('NSW', 500), 500)
        is_converged = self.data.get('convergence_analysis', {}).get('overall_convergence', False)
        completion_ratio = 1.0 if is_converged else (total_steps / max_steps if max_steps > 0 else 0.0)

        process_data: Dict[str, Any] = {
            'total_steps': total_steps,
            'max_steps': max_steps,
            'completion_ratio': completion_ratio,
            'is_converged': is_converged,
            'data_available': True,
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
            scf_energies = self.data.get('scf_energies', [])
            if not scf_energies:
                return
            self.data['final_results'] = {
                'final_energy': scf_energies[-1],
                'final_step': 0,
                'final_forces': [],
                'final_positions': [],
                'max_residual_force': None,
                'rms_residual_force': None,
            }
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
            + self._generate_plain_language_interpretation()
            + self._generate_convergence_section()
            + self._generate_optimization_process_section()
            + self._generate_final_results_section()
            + self._generate_structure_comparison_section()
        )

    def _chart_guide(self, text: str) -> str:
        return (
            '<div style="margin-top: 12px; padding: 12px 14px; '
            'background: #f7fafc; border-left: 4px solid #4299e1; '
            'border-radius: 6px; color: #2d3748;">'
            f'<strong>怎么看这张图：</strong> {text}'
            '</div>'
        )

    @staticmethod
    def _is_real_number(value: Any) -> bool:
        return (
            value is not None
            and not isinstance(value, bool)
            and isinstance(value, (int, float))
            and math.isfinite(float(value))
        )

    def _format_optional_float(self, value: Any, precision: int, *, scientific: bool = False) -> str:
        if not self._is_real_number(value):
            return "--"
        fmt = f".{precision}e" if scientific else f".{precision}f"
        return format(float(value), fmt)

    def _describe_energy_convergence(self, energy_conv: Dict[str, Any], energy_history: List[float]) -> Tuple[str, str]:
        available = energy_conv.get('available')
        if available is None:
            available = len(energy_history) > 1
        if not available:
            reason = energy_conv.get('reason')
            if reason == 'insufficient_energy_history':
                return '⚪ 数据不足', 'SCF 能量历史不足，无法判断尾部波动趋势。'
            return '⚪ 数据不足', '未解析到可用的 SCF 能量历史，无法判断能量尾部波动。'
        return (
            '✅ 已稳定' if energy_conv.get('converged', False) else '❌ 波动较大',
            '该指标描述 SCF 尾部能量波动，不直接代表 VASP 停止条件。'
        )

    def _describe_force_convergence(self, force_conv: Dict[str, Any], force_history: List[float]) -> Tuple[str, str]:
        available = force_conv.get('available')
        if available is None:
            available = len(force_history) > 0
        if not available:
            return '⚪ 数据不足', '未解析到有效的离子步力信息，无法判断力收敛。'
        return (
            '✅ 收敛' if force_conv.get('converged', False) else '❌ 未收敛',
            '',
        )

    def _generate_summary(self) -> str:
        process = self.data.get('optimization_process', {})
        final = self.data.get('final_results', {})

        return f"""
        <div class="section">
            <h2>📊 计算摘要</h2>
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>任务状态</h3>
                    <div class="value">
                        <span class="convergence-status converged">已完成</span>
                    </div>
                </div>
                <div class="summary-card">
                    <h3>实际离子步数</h3>
                    <div class="value">{process.get('total_steps', 0)} 步</div>
                </div>
                <div class="summary-card">
                    <h3>最终能量</h3>
                    <div class="value">{self._format_optional_float(final.get('final_energy'), 6)} eV</div>
                </div>
                <div class="summary-card">
                    <h3>最大剩余力</h3>
                    <div class="value">{self._format_optional_float(final.get('max_residual_force'), 4)} eV/Å</div>
                </div>
            </div>
        </div>
        """

    def _generate_plain_language_interpretation(self) -> str:
        """生成结构优化的通俗解读。"""
        convergence = self.data.get('convergence_analysis', {})
        final = self.data.get('final_results', {})
        force_conv = convergence.get('force_convergence', {})
        energy_conv = convergence.get('energy_convergence', {})

        force_ready = force_conv.get('available') and force_conv.get('converged')
        max_force = final.get('max_residual_force')
        gap_text = (
            "从结构角度看，当前原子位置已经基本定住。"
            if force_ready
            else "从结构角度看，原子位置可能还在继续调整。"
        )
        next_step = (
            "通常可以把这套结构拿去做后续 SCF、DOS、能带等电子结构计算。"
            if force_ready
            else "如果后续要做高精度电子结构分析，通常建议先把结构继续优化到更稳。"
        )
        energy_text = (
            "能量波动主要用来辅助理解 SCF 尾部是否平稳，它不是结构优化是否停止的唯一标准。"
            if energy_conv.get('available')
            else "这次没有足够的 SCF 能量历史，所以不建议只靠能量波动来判断优化质量。"
        )

        max_force_text = (
            f"当前最大剩余力约 {float(max_force):.4f} eV/Å。"
            if isinstance(max_force, (int, float))
            else "这次没有可靠的最终力统计。"
        )

        return f"""
        <div class="section">
            <h2>🧭 通俗解读</h2>
            <div class="summary-card">
                <p><strong>这份报告最重要的是看结构是否已经基本定住。</strong></p>
                <ul>
                    <li>{gap_text} {max_force_text}</li>
                    <li>{next_step}</li>
                    <li>{energy_text}</li>
                </ul>
            </div>
        </div>
        """

    def _generate_convergence_section(self) -> str:
        convergence = self.data.get('convergence_analysis', {})
        force_conv = convergence.get('force_convergence', {})
        energy_conv = convergence.get('energy_convergence', {})

        force_history = force_conv.get('force_history', [])
        energy_history = energy_conv.get('energy_history', self.data.get('scf_energies', []))
        force_status, force_note = self._describe_force_convergence(force_conv, force_history)
        energy_status, energy_note = self._describe_energy_convergence(energy_conv, energy_history)

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
                    <p><strong>收敛阈值:</strong> {self._format_optional_float(force_conv.get('threshold'), 4)} eV/Å</p>
                    <p><strong>最终最大力:</strong> {self._format_optional_float(force_conv.get('final_max_force'), 4)} eV/Å</p>
                    {f'<p><strong>说明:</strong> {force_note}</p>' if force_note else ''}
                    <div class="chart-container"><canvas id="forceChart"></canvas></div>
                    {self._chart_guide('先看红线是否逐步压到阈值线下方；如果最后几步稳定低于阈值，通常说明原子位置已经基本放松到位。')}
                    {force_dl}
                </div>
                <div>
                    <h3>能量波动分析</h3>
                    <p><strong>状态:</strong> {energy_status}</p>
                    <p><strong>收敛阈值:</strong> {self._format_optional_float(energy_conv.get('threshold'), 2, scientific=True)} eV</p>
                    <p><strong>最终能量:</strong> {self._format_optional_float(energy_conv.get('final_energy'), 6)} eV</p>
                    <p><strong>平均能量变化:</strong> {self._format_optional_float(energy_conv.get('avg_recent_change'), 2, scientific=True)} eV</p>
                    {f'<p><strong>说明:</strong> {energy_note}</p>' if energy_note else ''}
                    <div class="chart-container"><canvas id="energyChart"></canvas></div>
                    {self._chart_guide('看后几步是否只剩小幅波动；如果早期有明显跳变、后期逐渐变平，通常是结构弛豫里比较常见的稳定过程。')}
                    {energy_dl}
                </div>
            </div>
        </div>
        """

    def _generate_optimization_process_section(self) -> str:
        process = self.data.get('optimization_process', {})
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
            <h2>📈 迭代曲线参考</h2>
            <p><strong>说明:</strong> 以下曲线用于辅助理解本次计算的迭代行为，不表示当前任务状态。</p>
            {('<p><strong>补充:</strong> 未从 OUTCAR 中解析到完整离子步信息，以下曲线仅供参考。</p>' if not process.get('data_available', True) else '')}
            <div class="grid-2">
                <div>
                    <h3>能量演化曲线</h3>
                    <div class="chart-container"><canvas id="energyEvolutionChart"></canvas></div>
                    {self._chart_guide('看整体能量是否逐步下降并趋于平缓；中间偶尔的小回弹不一定异常，但如果长期大起大落，往往值得回头检查结构或参数设置。')}
                    {energy_evo_dl}
                </div>
                <div>
                    <h3>力演化曲线</h3>
                    <div class="chart-container"><canvas id="forceEvolutionChart"></canvas></div>
                    {self._chart_guide('看最大力和 RMS 力是否一起往下走；如果两条曲线后期仍然维持较高水平，通常说明结构还没有完全放松。')}
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
                        <tr><td>最终总能量</td><td>{self._format_optional_float(final.get('final_energy'), 6)} eV</td></tr>
                        <tr><td>完成步数</td><td>{final.get('final_step', 0)}</td></tr>
                    </table>
                </div>
                <div>
                    <h3>力统计信息</h3>
                    <table class="data-table">
                        <tr><td>最大剩余力</td><td>{self._format_optional_float(final.get('max_residual_force'), 4)} eV/Å</td></tr>
                        <tr><td>RMS剩余力</td><td>{self._format_optional_float(final.get('rms_residual_force'), 4)} eV/Å</td></tr>
                        <tr><td>平均力大小</td><td>{self._format_optional_float(force_stats.get('mean_force_magnitude'), 4)} eV/Å</td></tr>
                        <tr><td>力标准差</td><td>{self._format_optional_float(force_stats.get('std_force_magnitude'), 4)} eV/Å</td></tr>
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
        export_optimization_assets(analysis_data, output_dir_path)
        output_file = output_dir_path / "optimization_analysis_report.html"

        generator = OptimizationHTMLGenerator(analysis_data)
        return generator.generate_html_report(str(output_file))
    except Exception as e:
        raise Exception(f"生成优化分析报告失败: {str(e)}")


def export_optimization_assets(data: Dict[str, Any], output_dir: str | Path) -> None:
    output_dir_path = Path(output_dir)
    output_dir_path.mkdir(parents=True, exist_ok=True)

    convergence = data.get("convergence_analysis", {})
    process = data.get("optimization_process", {})

    force_history = convergence.get("force_convergence", {}).get("force_history", [])
    energy_history = convergence.get("energy_convergence", {}).get("energy_history", data.get("scf_energies", []))
    energy_profile = process.get("energy_profile", [])
    force_profile = process.get("force_profile", [])

    _write_csv(
        output_dir_path / "force_convergence.csv",
        ["ionic_step", "max_force_eV_per_A"],
        [[idx + 1, value] for idx, value in enumerate(force_history)],
    )
    _write_csv(
        output_dir_path / "scf_energy_convergence.csv",
        ["scf_step", "total_energy_eV"],
        [[idx + 1, value] for idx, value in enumerate(energy_history)],
    )
    _write_csv(
        output_dir_path / "energy_evolution.csv",
        ["ionic_step", "energy_eV"],
        [[item.get("step"), item.get("energy")] for item in energy_profile],
    )
    _write_csv(
        output_dir_path / "force_evolution.csv",
        ["ionic_step", "max_force_eV_per_A", "rms_force_eV_per_A"],
        [[item.get("step"), item.get("max_force"), item.get("rms_force")] for item in force_profile],
    )

    _save_line_plot(
        output_dir_path / "force_convergence.png",
        title="Force Convergence",
        x_label="Ionic Step",
        y_label="Max Force (eV/A)",
        x_values=list(range(1, len(force_history) + 1)),
        y_series=[("Max Force", force_history)],
    )
    _save_line_plot(
        output_dir_path / "scf_energy_convergence.png",
        title="SCF Energy Convergence",
        x_label="SCF Step",
        y_label="Total Energy (eV)",
        x_values=list(range(1, len(energy_history) + 1)),
        y_series=[("Total Energy", energy_history)],
    )
    _save_line_plot(
        output_dir_path / "energy_evolution.png",
        title="Energy Evolution",
        x_label="Ionic Step",
        y_label="Energy (eV)",
        x_values=[item.get("step") for item in energy_profile],
        y_series=[("Total Energy", [item.get("energy") for item in energy_profile])],
    )
    _save_line_plot(
        output_dir_path / "force_evolution.png",
        title="Force Evolution",
        x_label="Ionic Step",
        y_label="Force (eV/A)",
        x_values=[item.get("step") for item in force_profile],
        y_series=[
            ("Max Force", [item.get("max_force") for item in force_profile]),
            ("RMS Force", [item.get("rms_force") for item in force_profile]),
        ],
    )


def _write_csv(path: Path, headers: List[str], rows: List[List[Any]]) -> None:
    with path.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(headers)
        writer.writerows(rows)


def _save_line_plot(
    path: Path,
    *,
    title: str,
    x_label: str,
    y_label: str,
    x_values: List[Any],
    y_series: List[Tuple[str, List[Any]]],
) -> None:
    if not x_values or not any(series for _, series in y_series):
        return

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 5))
    for label, series in y_series:
        if not series:
            continue
        paired = [(x, y) for x, y in zip(x_values, series) if x is not None and y is not None]
        if not paired:
            continue
        xs, ys = zip(*paired)
        ax.plot(xs, ys, marker="o", linewidth=1.8, markersize=3, label=label)

    if not ax.lines:
        plt.close(fig)
        return

    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.grid(True, alpha=0.25)
    if len(ax.lines) > 1:
        ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)
