#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VASP结构优化可视化分析模块

功能：
1. 解析OUTCAR文件，提取优化过程中的关键数据
2. 分析收敛性（力收敛、能量收敛、电子步收敛）
3. 监控优化过程（离子步数、能量-步数曲线、结构变化、力的变化）
4. 分析最终结果（晶格参数、原子坐标、总能量、剩余力）
5. 生成HTML可视化报告

作者: VASP API Team
日期: 2025年
"""

import re
import json
import math
import numpy as np
import os
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime


class OUTCARAnalyzer:
    """VASP计算结果分析器"""
    
    def __init__(self, input_path: str, task_id: Optional[str] = None):
        """
        初始化分析器
        
        Args:
            input_path: OUTCAR文件路径或包含VASP文件的文件夹路径
            task_id: 任务ID
        """
        self.input_path = Path(input_path)
        self.task_id = task_id or "unknown"
        
        # 判断输入是文件还是文件夹
        if self.input_path.is_file():
            # 输入是OUTCAR文件
            self.work_dir = self.input_path.parent
            self.outcar_path = self.input_path
        elif self.input_path.is_dir():
            # 输入是文件夹
            self.work_dir = self.input_path
            # 兼容OUTCAR文件可能有.txt后缀的情况
            self.outcar_path = self._find_outcar_file(self.work_dir)
            if not self.outcar_path:
                tried_names = ["OUTCAR", "OUTCAR.txt", "outcar", "outcar.txt", "Outcar", "Outcar.txt"]
                raise FileNotFoundError(f"文件夹中未找到OUTCAR文件，尝试了: {', '.join(tried_names)}\n路径: {self.work_dir}")
        else:
            raise FileNotFoundError(f"输入路径不存在: {input_path}")
        
        # 检查其他VASP文件（兼容不同后缀）
        self.poscar_path = self._find_vasp_file(self.work_dir, "POSCAR")
        self.contcar_path = self._find_vasp_file(self.work_dir, "CONTCAR")
        
        self.data = {
            'file_info': {},
            'calculation_settings': {},
            'ionic_steps': [],
            'convergence_analysis': {},
            'final_results': {},
            'electronic_structure': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {}
        }
        
    def analyze(self) -> Dict[str, Any]:
        """
        执行完整分析
        
        Returns:
            包含所有分析结果的字典
        """
        print(f"🔍 开始分析OUTCAR文件: {self.outcar_path}")
        
        # 解析结构文件
        self._parse_structure_files()
        
        # 解析OUTCAR文件
        self._parse_outcar()
        
        # 分析收敛性
        self._analyze_convergence()
        
        # 分析优化过程
        self._analyze_optimization_process()
        
        # 分析最终结果
        self._analyze_final_results()
        
        print(f"✅ OUTCAR分析完成")
        return self.data
    
    def _find_outcar_file(self, work_dir: Path) -> Optional[Path]:
        """智能寻找OUTCAR文件，兼容不同后缀"""
        # 可能的OUTCAR文件名
        possible_names = [
            "OUTCAR",        # 标准名称
            "OUTCAR.txt",    # Mac下载可能带.txt后缀
            "outcar",        # 小写版本
            "outcar.txt",    # 小写+.txt后缀
            "Outcar",        # 首字母大写
            "Outcar.txt"     # 首字母大写+.txt后缀
        ]
        
        for name in possible_names:
            outcar_path = work_dir / name
            if outcar_path.exists() and outcar_path.is_file():
                print(f"   🔍 找到OUTCAR文件: {name}")
                return outcar_path
        
        return None
    
    def _find_vasp_file(self, work_dir: Path, base_name: str) -> Optional[Path]:
        """智能寻找VASP文件，兼容不同后缀和大小写"""
        # 生成可能的文件名变体
        possible_names = [
            base_name,                    # 原始名称，如 POSCAR
            f"{base_name}.txt",          # 带.txt后缀
            base_name.lower(),           # 小写版本
            f"{base_name.lower()}.txt",  # 小写+.txt后缀
            base_name.capitalize(),      # 首字母大写
            f"{base_name.capitalize()}.txt"  # 首字母大写+.txt后缀
        ]
        
        for name in possible_names:
            file_path = work_dir / name
            if file_path.exists() and file_path.is_file():
                print(f"   🔍 找到{base_name}文件: {name}")
                return file_path
        
        return None
    

    
    def _parse_structure_files(self):
        """解析结构文件POSCAR和CONTCAR"""
        print("🔍 解析结构文件...")
        
        # 解析POSCAR
        if self.poscar_path and self.poscar_path.exists():
            try:
                with open(self.poscar_path, 'r', encoding='utf-8', errors='ignore') as f:
                    poscar_content = f.read()
                self.data['structure_files']['poscar'] = poscar_content
                # 提取化学组成
                composition = self._extract_composition_from_poscar(poscar_content)
                self.data['task_info']['composition'] = composition
                print(f"   ✅ 化学组成: {composition}")
            except Exception as e:
                print(f"   ❌ 读取POSCAR失败: {e}")
        else:
            print("   ⚠️ 未找到POSCAR文件")
            self.data['task_info']['composition'] = "Unknown"
        
        # 解析CONTCAR
        if self.contcar_path and self.contcar_path.exists():
            try:
                with open(self.contcar_path, 'r', encoding='utf-8', errors='ignore') as f:
                    contcar_content = f.read()
                self.data['structure_files']['contcar'] = contcar_content
                print("   ✅ 找到CONTCAR文件")
            except Exception as e:
                print(f"   ❌ 读取CONTCAR失败: {e}")
        else:
            print("   ⚠️ 未找到CONTCAR文件")
    
    def _extract_composition_from_poscar(self, poscar_content: str) -> str:
        """从POSCAR文件提取化学组成"""
        try:
            lines = poscar_content.strip().split('\n')
            if len(lines) < 6:
                return "Unknown"
            
            # 第6行是元素符号（VASP 5.x格式）
            elements_line = lines[5].strip()
            # 第7行是原子数量
            counts_line = lines[6].strip()
            
            elements = elements_line.split()
            counts = [int(x) for x in counts_line.split()]
            
            if len(elements) != len(counts):
                return "Unknown"
            
            # 构建化学式
            composition = ""
            for elem, count in zip(elements, counts):
                if count > 1:
                    composition += f"{elem}{count}"
                else:
                    composition += elem
            
            return composition if composition else "Unknown"
        except Exception:
            return "Unknown"
    
    def _parse_outcar(self):
        """解析OUTCAR文件"""
        print("📖 正在解析OUTCAR文件...")
        
        if not self.outcar_path:
            raise FileNotFoundError("OUTCAR文件路径未设置")
        
        with open(self.outcar_path, 'r', encoding='utf-8', errors='ignore') as f:
            content = f.read()
        
        print(f"🔍 文件大小: {len(content)} 字符")
        print(f"🔍 文件前100字符: {content[:100]}")
        
        # 检查文件类型
        if 'vasp' not in content.lower() and 'POTCAR' not in content:
            print("⚠️ 警告: 这可能不是OUTCAR文件，检测不到VASP特征")
        
        # 解析基本信息
        self._parse_file_info(content)
        
        # 解析计算设置
        self._parse_calculation_settings(content)
        
        # 解析离子步
        self._parse_ionic_steps(content)
        
        # 解析电子结构信息
        self._parse_electronic_structure(content)
    
    def _parse_file_info(self, content: str):
        """解析文件基本信息"""
        print("🔍 解析文件基本信息...")
        
        # VASP版本
        version_match = re.search(r'vasp\.([\d\.]+)', content)
        if version_match:
            self.data['file_info']['vasp_version'] = version_match.group(1)
            print(f"   ✅ VASP版本: {version_match.group(1)}")
        else:
            print("   ❌ 未找到VASP版本信息")
        
        # 计算日期
        date_match = re.search(r'executed on.*date ([\d\.]+)', content)
        if date_match:
            self.data['file_info']['calculation_date'] = date_match.group(1)
            print(f"   ✅ 计算日期: {date_match.group(1)}")
        else:
            print("   ❌ 未找到计算日期信息")
        
        # 核数信息
        cores_match = re.search(r'running on\s+(\d+)\s+total cores', content)
        if cores_match:
            self.data['file_info']['total_cores'] = int(cores_match.group(1))
            print(f"   ✅ 核数: {cores_match.group(1)}")
        else:
            print("   ❌ 未找到核数信息")
    
    def _parse_calculation_settings(self, content: str):
        """解析计算设置"""
        print("🔍 解析计算设置...")
        settings = {}
        
        # 提取INCAR参数
        incar_section = re.search(r'INCAR:(.*?)POTCAR:', content, re.DOTALL)
        if incar_section:
            incar_text = incar_section.group(1)
            print(f"   ✅ 找到INCAR部分，长度: {len(incar_text)} 字符")
            
            # 关键参数
            params = {
                'EDIFFG': r'EDIFFG\s*=\s*([-\d\.E]+)',
                'NSW': r'NSW\s*=\s*(\d+)',
                'POTIM': r'POTIM\s*=\s*([\d\.]+)',
                'IBRION': r'IBRION\s*=\s*(\d+)',
                'ISIF': r'ISIF\s*=\s*(\d+)',
                'ENCUT': r'ENCUT\s*=\s*([\d\.]+)',
                'EDIFF': r'EDIFF\s*=\s*([\d\.E-]+)',
                'NELM': r'NELM\s*=\s*(\d+)',
                'PREC': r'PREC\s*=\s*(\w+)'
            }
            
            for param, pattern in params.items():
                match = re.search(pattern, incar_text)
                if match:
                    value = match.group(1)
                    # 尝试转换为数字
                    try:
                        if '.' in value or 'E' in value or 'e' in value:
                            settings[param] = float(value)
                        else:
                            settings[param] = int(value)
                    except ValueError:
                        settings[param] = value
                    print(f"   ✅ {param}: {settings[param]}")
                else:
                    print(f"   ❌ 未找到 {param}")
        else:
            print("   ❌ 未找到INCAR部分")
        
        self.data['calculation_settings'] = settings
        print(f"   📊 总共解析了 {len(settings)} 个参数")
    
    def _parse_ionic_steps(self, content: str):
        """解析离子步信息（以力块为主；SCF能量序列用于画图，非一一对应）"""
        print("🔍 解析离子步信息...")
        ionic_steps: List[Dict[str, Any]] = []
        
        # 1) 收集SCF能量序列（来自TOTEN等行）
        scf_energy_series: List[float] = []
        energy_candidates: List[Tuple[int, float, str]] = []  # (pos, value, tag)
        # 先尝试标准格式
        energy_pattern = r'free\s+energy\s+TOTEN\s*=\s*([\-\d\.]+)\s*eV'
        energy_matches = list(re.finditer(energy_pattern, content))
        print("----")
        print(len(energy_matches))
        # 若未命中，再按旧的备选顺序回退
        if len(energy_matches) == 0:
            print("   🔍 尝试其他能量格式用于SCF序列...")
            alt_patterns = [
                r'energy\s*without\s*entropy[^=]*=\s*([\-\d\.]+)',
                r'TOTEN\s*=\s*([\-\d\.]+)',
                r'total\s*energy\s*=\s*([\-\d\.]+)'
            ]
            for i, pattern in enumerate(alt_patterns):
                alt_matches = list(re.finditer(pattern, content, re.IGNORECASE))
                print(f"   🔍 SCF替代格式 {i+1}: 找到 {len(alt_matches)} 个匹配")
                if len(alt_matches) > 0:
                    energy_matches = alt_matches
                    break
        # 用选中的匹配构建SCF能量序列与候选位点
        for m in energy_matches:
            try:
                val = float(m.group(1))
                scf_energy_series.append(val)
                energy_candidates.append((m.start(), val, 'SCF'))
            except Exception:
                continue
        energy_candidates.sort(key=lambda x: x[0])
        print(f"   🔍 SCF能量点数: {len(scf_energy_series)}")
        if len(scf_energy_series) < 2:
            print("   ⚠️ SCF能量序列过短，能量收敛与能量曲线信息有限")
        
        # 2) 按TOTAL-FORCE块计数离子步
        force_blocks = list(re.finditer(r'POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*\n\s*-{10,}\s*\n(.*?)\n\s*-{10,}', content, re.DOTALL))
        print(f"   🔍 找到 {len(force_blocks)} 个力信息部分（按此计数离子步）")
    
        def find_prev_energy_value(before_pos: int) -> Optional[float]:
            left = [e for e in energy_candidates if e[0] < before_pos]
            return left[-1][1] if left else None
        
        step_num = 0
        for j, fm in enumerate(force_blocks):
            block_start = fm.start()
            block = fm.group(1)
            lines = block.strip().split('\n') if block else []
            forces: List[List[float]] = []
            positions: List[List[float]] = []
            parsed_lines = 0
            
            for raw in lines:
                line = raw.strip()
                if not line or line.startswith('-'):
                    continue
                ls = line.lower()
                if ls.startswith('total drift'):
                    continue
                parts = line.split()
                if len(parts) < 6:
                    print(f"   ⚠️ 力行列数不足(跳过): '{line}'")
                    continue
                try:
                    px, py, pz = float(parts[0]), float(parts[1]), float(parts[2])
                    fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                    positions.append([px, py, pz])
                    forces.append([fx, fy, fz])
                    parsed_lines += 1
                except Exception:
                    print(f"   ⚠️ 力行解析失败(跳过): '{line}'")
                    continue
            
            if j < 2:
                print(f"   ✅ 力块 {j+1}: 行数={len(lines)}, 解析成功={parsed_lines}")
                if lines:
                    print(f"      示例首行: {lines[0]}")
                    print(f"      示例尾行: {lines[-1]}")
            
            max_force_val: Optional[float] = None
            rms_force_val: Optional[float] = None
            if forces:
                mags = [math.sqrt(fx*fx + fy*fy + fz*fz) for fx, fy, fz in forces]
                if mags:
                    max_force_val = max(mags)
                    rms_force_val = math.sqrt(sum(m*m for m in mags) / len(mags))
            else:
                print(f"   ❌ 力块 {j+1}: 未解析到有效的原子力行")
            
            # 为该离子步附带一个“就近回溯”的参考能量（仅用于展示，不代表离子步能量）
            ref_energy = find_prev_energy_value(block_start)
            if ref_energy is None:
                print(f"   ⚠️ 力块 {j+1}: 前无能量记录，参考能量置为 0.0")
                ref_energy = 0.0
            
            step_num += 1
            ionic_steps.append({
                'step': step_num,
                'energy': float(ref_energy),
                'forces': forces,
                'max_force': max_force_val,
                'rms_force': rms_force_val,
                'positions': positions
            })
            if j < 3:
                print(f"   ▶️ 步骤 {step_num}: 参考能量={ref_energy:.6f} eV, maxF={(max_force_val or 0):.4f} eV/Å")
        
        self.data['ionic_steps'] = ionic_steps
        # 存储SCF能量序列
        if 'electronic_structure' not in self.data:
            self.data['electronic_structure'] = {}
        self.data['scf_energies'] = scf_energy_series
        print(f"📊 解析离子步完成: {len(ionic_steps)} 步；SCF能量点: {len(scf_energy_series)}")
    
    def _parse_electronic_structure(self, content: str):
        """解析电子结构信息"""
        electronic_data = {}
        
        # 查找带隙信息
        gap_pattern = r'E-fermi\s*:\s*([-\d\.]+).*?gap\s*=\s*([\d\.]+)'
        gap_match = re.search(gap_pattern, content)
        if gap_match:
            electronic_data['fermi_energy'] = float(gap_match.group(1))
            electronic_data['band_gap'] = float(gap_match.group(2))
        
        # 查找磁矩信息
        mag_pattern = r'number of electron\s+NELECT\s*=\s*([\d\.]+)'
        mag_match = re.search(mag_pattern, content)
        if mag_match:
            electronic_data['total_electrons'] = float(mag_match.group(1))
        
        self.data['electronic_structure'] = electronic_data
    
    def _analyze_convergence(self):
        """分析收敛性（参考OUTCAR尾部标志，并提供详细调试信息）"""
        print("🧪 正在进行收敛性分析...")
        convergence = {
            'force_convergence': {},
            'energy_convergence': {},
            'electronic_convergence': {},
            'tail_check': {},
            'overall_convergence': False
        }
        
        # 1) 优先进行尾部收敛性检查（与vasp_worker一致）
        tail_info = {
            'matched': False,
            'keywords': ['reached required accuracy', 'Voluntary'],
            'tail_bytes': 1024,
            'exception': None
        }
        try:
            if not self.outcar_path:
                raise Exception("OUTCAR文件路径未设置")
            file_size = self.outcar_path.stat().st_size
            tail_bytes = min(1024, file_size)
            with open(self.outcar_path, 'rb') as f:
                f.seek(-tail_bytes, os.SEEK_END)
                last_lines = f.readlines()[-10:]
                last_content = b''.join(last_lines).decode('utf-8', errors='ignore')
                lc = last_content.lower()
                matched = ('reached required accuracy' in lc) or ('voluntary' in lc)
                tail_info['matched'] = matched
                tail_info['tail_bytes'] = tail_bytes
                if not matched:
                    # 提供一个简要的尾部内容片段，辅助调试
                    snippet = last_content[-300:].replace('\n', ' | ')
                    print(f"   🔎 尾部未匹配到收敛关键词，最后片段: {snippet}")
                else:
                    print("   ✅ 尾部标志显示已收敛（匹配关键词）")
        except Exception as e:
            tail_info['exception'] = str(e)
            print(f"   ⚠️ 读取OUTCAR尾部失败: {e}")
        finally:
            convergence['tail_check'] = tail_info
        
        # 如果尾部标志明确收敛，则直接认为整体收敛
        overall_converged_by_tail = tail_info.get('matched', False)
        
        # 2) 解析型收敛（力/能量）分析与调试信息
        ionic_steps = self.data.get('ionic_steps', [])
        print(f"   📊 离子步统计: steps = {len(ionic_steps)}")
        if not ionic_steps:
            print("   ❌ 未找到离子步数据（ionic_steps为空）。将仅依据尾部标志判断。")
            convergence['overall_convergence'] = overall_converged_by_tail
            self.data['convergence_analysis'] = convergence
            return
        
        # 力收敛分析
        ediffg = self.data['calculation_settings'].get('EDIFFG', -0.01)
        if 'EDIFFG' not in self.data['calculation_settings']:
            print("   ⚠️ 未在INCAR中找到EDIFFG，使用默认值 -0.01 eV/Å")
        force_threshold = abs(ediffg)
        
        max_forces = [step.get('max_force', None) for step in ionic_steps if step.get('max_force') is not None]
        if not max_forces:
            print("   ❌ 未能提取到任何max_force（可能未解析到TOTAL-FORCE块）。")
        else:
            final_max_force = max_forces[-1]
            force_converged = final_max_force < force_threshold
            print(f"   🔧 力分析: 阈值={force_threshold:.4f} eV/Å, 最终max={final_max_force:.4f} eV/Å, 收敛={force_converged}")
            convergence['force_convergence'] = {
                'converged': force_converged,
                'threshold': force_threshold,
                'final_max_force': final_max_force,
                'force_history': max_forces
            }
        
        # 能量收敛分析
        scf_energies = self.data.get('scf_energies', [])
        if len(scf_energies) <= 1:
            print(f"   ❌ SCF能量点不足以评估能量收敛（共 {len(scf_energies)} 个）。")
        else:
            energy_changes = [abs(scf_energies[i] - scf_energies[i-1]) for i in range(1, len(scf_energies))]
            recent_changes = energy_changes[-5:] if len(energy_changes) >= 5 else energy_changes
            avg_recent_change = sum(recent_changes) / len(recent_changes) if recent_changes else 0.0
            energy_converged = avg_recent_change < 1e-4  # 0.1 meV 阈值
            print(f"   🔧 能量分析(SCF): 最终能量={scf_energies[-1]:.6f} eV, 最近ΔE均值={avg_recent_change:.2e} eV, 收敛={energy_converged}")
            convergence['energy_convergence'] = {
                'converged': energy_converged,
                'final_energy': scf_energies[-1],
                'energy_history': scf_energies,
                'energy_changes': energy_changes,
                'avg_recent_change': avg_recent_change
            }
        
        # 3) 组合总体收敛判断
        force_ok = convergence['force_convergence'].get('converged', False)
        energy_ok = convergence['energy_convergence'].get('converged', False)
        if overall_converged_by_tail:
            print("   ✅ 使用尾部标志认定整体收敛。")
            overall = True
        else:
            overall = force_ok and energy_ok
            reason = []
            if not force_ok:
                reason.append("力未收敛")
            if not energy_ok:
                reason.append("能量未收敛")
            if overall:
                print("   ✅ 力与能量均满足阈值，判定收敛。")
            else:
                if reason:
                    print("   ❌ 综合判定未收敛: " + ", ".join(reason))
                else:
                    print("   ❌ 综合判定未收敛（缺少充分信息）。")
        
        convergence['overall_convergence'] = overall
        self.data['convergence_analysis'] = convergence
    
    def _analyze_optimization_process(self):
        """分析优化过程"""
        if not self.data['ionic_steps']:
            return
        
        total_steps = len(self.data['ionic_steps'])
        max_steps = self.data['calculation_settings'].get('NSW', 500)
        
        # 修复完成度计算：如果已收敛，则完成度为100%
        convergence_data = self.data.get('convergence_analysis', {})
        is_converged = convergence_data.get('overall_convergence', False)
        if is_converged:
            completion_ratio = 1.0  # 已收敛则100%完成
        else:
            completion_ratio = total_steps / max_steps if max_steps > 0 else 0.0
        
        process_data = {
            'total_steps': total_steps,
            'max_steps': max_steps,
            'completion_ratio': completion_ratio,
            'is_converged': is_converged,
            'energy_profile': [],
            'force_profile': [],
            'structure_changes': []
        }
        
        # 能量和力的演化
        for step in self.data['ionic_steps']:
            process_data['energy_profile'].append({
                'step': step['step'],
                'energy': step['energy']
            })
            
            if step.get('max_force') is not None:
                process_data['force_profile'].append({
                    'step': step['step'],
                    'max_force': step.get('max_force'),
                    'rms_force': step.get('rms_force')
                })
        
        # 计算结构变化（如果有位置信息）
        if len(self.data['ionic_steps']) > 1:
            for step_idx in range(1, len(self.data['ionic_steps'])):
                current_positions = self.data['ionic_steps'][step_idx].get('positions', [])
                previous_positions = self.data['ionic_steps'][step_idx-1].get('positions', [])
                
                if current_positions and previous_positions and len(current_positions) == len(previous_positions):
                    # 计算原子位移
                    displacements = []
                    for atom_idx in range(len(current_positions)):
                        curr = current_positions[atom_idx]
                        prev = previous_positions[atom_idx]
                        displacement = math.sqrt((curr[0] - prev[0])**2 + (curr[1] - prev[1])**2 + (curr[2] - prev[2])**2)
                        displacements.append(displacement)
                    
                    max_displacement = max(displacements) if displacements else 0.0
                    avg_displacement = sum(displacements) / len(displacements) if displacements else 0.0
                    
                    process_data['structure_changes'].append({
                        'step': self.data['ionic_steps'][step_idx]['step'],
                        'max_displacement': max_displacement,
                        'avg_displacement': avg_displacement
                    })
        
        self.data['optimization_process'] = process_data
    
    def _analyze_final_results(self):
        """分析最终结果"""
        if not self.data['ionic_steps']:
            return
        
        final_step = self.data['ionic_steps'][-1]
        
        # 优先取SCF序列末值作为最终能量
        scf_energies = self.data.get('scf_energies', [])
        final_energy_value = final_step['energy']
        if scf_energies:
            final_energy_value = scf_energies[-1]

        final_results = {
            'final_energy': final_energy_value,
            'final_step': final_step['step'],
            'final_forces': final_step.get('forces', []),
            'final_positions': final_step.get('positions', []),
            'max_residual_force': final_step.get('max_force'),
            'rms_residual_force': final_step.get('rms_force')
        }
        
        # 计算一些统计信息
        if final_step.get('forces'):
            force_magnitudes = []
            for force in final_step['forces']:
                magnitude = math.sqrt(force[0]**2 + force[1]**2 + force[2]**2)
                force_magnitudes.append(magnitude)
            
            if force_magnitudes:
                n = len(force_magnitudes)
                mean_mag = sum(force_magnitudes) / n
                variance = sum((mag - mean_mag)**2 for mag in force_magnitudes) / n
                std_mag = math.sqrt(variance)
                
                final_results['force_statistics'] = {
                    'mean_force_magnitude': mean_mag,
                    'std_force_magnitude': std_mag,
                    'min_force_magnitude': min(force_magnitudes),
                    'max_force_magnitude': max(force_magnitudes)
                }
        
        self.data['final_results'] = final_results


class OptimizationHTMLGenerator:
    """HTML报告生成器"""
    
    def __init__(self, analysis_data: Dict[str, Any]):
        """
        初始化HTML生成器
        
        Args:
            analysis_data: 分析结果数据
        """
        self.data = analysis_data
    
    def _filter_contcar_velocities(self, lines: List[str]) -> List[str]:
        """过滤CONTCAR文件中的原子速度部分"""
        filtered_lines = []
        in_velocity_section = False
        
        for i, line in enumerate(lines):
            # 检测是否到达原子速度部分
            # 通常在原子坐标后有一个空行，然后是速度数据
            if not in_velocity_section:
                # 如果当前行是空行，检查下一行是否是速度数据
                if line.strip() == '' and i + 1 < len(lines):
                    next_line = lines[i + 1].strip()
                    # 检查下一行是否包含科学计数法的零值（速度数据特征）
                    if ('0.00000000E+00' in next_line or 
                        next_line.replace(' ', '').replace('0', '').replace('.', '').replace('E', '').replace('+', '') == ''):
                        in_velocity_section = True
                        filtered_lines.append(line)  # 保留空行
                        break  # 不再添加后续的速度数据
                    else:
                        filtered_lines.append(line)
                else:
                    filtered_lines.append(line)
            # 如果已经进入速度部分，就不再添加任何行
        
        return filtered_lines
    
    def generate_html_report(self, output_path: str) -> str:
        """
        生成HTML报告
        
        Args:
            output_path: 输出HTML文件路径
            
        Returns:
            生成的HTML文件路径
        """
        print(f"📄 正在生成HTML报告: {output_path}")
        
        html_content = self._generate_html_content()
        
        output_file = Path(output_path)
        output_file.parent.mkdir(parents=True, exist_ok=True)
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        print(f"✅ HTML报告已生成: {output_file}")
        return str(output_file)
    
    def _generate_html_content(self) -> str:
        """生成HTML内容"""
        html = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>VASP结构优化分析报告</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
    <style>
        {self._get_css_styles()}
    </style>
</head>
<body>
    <div class="container">
        {self._generate_header()}
        {self._generate_summary()}
        {self._generate_convergence_section()}
        {self._generate_optimization_process_section()}
        {self._generate_final_results_section()}
        {self._generate_structure_comparison_section()}
        {self._generate_footer()}
    </div>
    
    <script>
        {self._generate_javascript()}
    </script>
</body>
</html>
"""
        return html
    
    def _get_css_styles(self) -> str:
        """获取CSS样式"""
        return """
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background-color: #f5f5f5;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
        }
        
        .header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
            border-radius: 10px;
            margin-bottom: 30px;
            text-align: center;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .section {
            background: white;
            padding: 25px;
            margin-bottom: 20px;
            border-radius: 10px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }
        
        .section h2 {
            color: #4a5568;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 2px solid #e2e8f0;
            font-size: 1.5em;
        }
        
        .summary-grid {
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 20px;
        }
        
        .summary-card {
            background: #f7fafc;
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #4299e1;
        }
        
        .summary-card h3 {
            color: #2d3748;
            margin-bottom: 10px;
        }
        
        .summary-card .value {
            font-size: 1.5em;
            font-weight: bold;
            color: #4299e1;
        }
        
        .convergence-status {
            display: inline-block;
            padding: 5px 15px;
            border-radius: 20px;
            font-weight: bold;
            color: white;
        }
        
        .converged {
            background-color: #48bb78;
        }
        
        .not-converged {
            background-color: #f56565;
        }
        
        .chart-container {
            position: relative;
            height: 400px;
            margin: 20px 0;
        }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
        }
        
        .data-table th,
        .data-table td {
            padding: 12px;
            text-align: left;
            border-bottom: 1px solid #e2e8f0;
        }
        
        .data-table th {
            background-color: #edf2f7;
            font-weight: bold;
            color: #4a5568;
        }
        
        .data-table tr:hover {
            background-color: #f7fafc;
        }
        
        .footer {
            text-align: center;
            margin-top: 40px;
            padding: 20px;
            color: #718096;
            font-size: 0.9em;
        }
        
        .grid-2 {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
        }
        
        @media (max-width: 768px) {
            .grid-2 {
                grid-template-columns: 1fr;
            }
            
            .summary-grid {
                grid-template-columns: 1fr;
            }
        }
        """
    
    def _generate_header(self) -> str:
        """生成页面头部"""
        file_info = self.data.get('file_info', {})
        task_info = self.data.get('task_info', {})
        return f"""
        <div class="header">
            <h1>🔬 VASP结构优化分析报告</h1>
            <p>任务ID: <strong>{task_info.get('task_id', 'Unknown')}</strong> | 
               材料组成: <strong>{task_info.get('composition', 'Unknown')}</strong></p>
            <p>生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>VASP版本: {file_info.get('vasp_version', 'Unknown')} | 
               计算日期: {file_info.get('calculation_date', 'Unknown')} | 
               核数: {file_info.get('total_cores', 'Unknown')}</p>
        </div>
        """
    
    def _generate_summary(self) -> str:
        """生成摘要部分"""
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
                        <span class="convergence-status {convergence_class}">
                            {convergence_text}
                        </span>
                    </div>
                </div>
                
                <div class="summary-card">
                    <h3>离子步数</h3>
                    <div class="value">
                        {process.get('total_steps', 0)} / {process.get('max_steps', 0)}
                    </div>
                </div>
                
                <div class="summary-card">
                    <h3>最终能量</h3>
                    <div class="value">
                        {final.get('final_energy', 0):.6f} eV
                    </div>
                </div>
                
                <div class="summary-card">
                    <h3>最大剩余力</h3>
                    <div class="value">
                        {(final.get('max_residual_force') or 0):.4f} eV/Å
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_convergence_section(self) -> str:
        """生成收敛性分析部分"""
        convergence = self.data.get('convergence_analysis', {})
        force_conv = convergence.get('force_convergence', {})
        energy_conv = convergence.get('energy_convergence', {})
        
        force_status = '✅ 收敛' if force_conv.get('converged', False) else '❌ 未收敛'
        energy_status = '✅ 收敛' if energy_conv.get('converged', False) else '❌ 未收敛'
        
        return f"""
        <div class="section">
            <h2>🎯 收敛性检查</h2>
            
            <div class="grid-2">
                <div>
                    <h3>力收敛分析</h3>
                    <p><strong>状态:</strong> {force_status}</p>
                    <p><strong>收敛阈值:</strong> {(force_conv.get('threshold') or 0):.4f} eV/Å</p>
                    <p><strong>最终最大力:</strong> {(force_conv.get('final_max_force') or 0):.4f} eV/Å</p>
                    
                    <div class="chart-container">
                        <canvas id="forceChart"></canvas>
                    </div>
                </div>
                
                <div>
                    <h3>能量收敛分析</h3>
                    <p><strong>状态:</strong> {energy_status}</p>
                    <p><strong>最终能量:</strong> {(energy_conv.get('final_energy') or 0):.6f} eV</p>
                    <p><strong>平均能量变化:</strong> {(energy_conv.get('avg_recent_change') or 0):.2e} eV</p>
                    
                    <div class="chart-container">
                        <canvas id="energyChart"></canvas>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_optimization_process_section(self) -> str:
        """生成优化过程监控部分"""
        process = self.data.get('optimization_process', {})
        is_converged = process.get('is_converged', False)
        completion_ratio = process.get('completion_ratio', 0)
        
        # 确定状态描述
        if is_converged:
            status_text = "已完成（收敛）"
            status_class = "converged"
        elif completion_ratio >= 0.95:
            status_text = f"{completion_ratio*100:.1f}% - 接近完成"
            status_class = "not-converged"
        else:
            status_text = f"{completion_ratio*100:.1f}% - 进行中"
            status_class = "not-converged"
        
        return f"""
        <div class="section">
            <h2>📈 优化过程监控</h2>
            
            <div class="summary-grid">
                <div class="summary-card">
                    <h3>完成状态</h3>
                    <div class="value">
                        <span class="convergence-status {status_class}">
                            {status_text}
                        </span>
                    </div>
                </div>
                <div class="summary-card">
                    <h3>步数进度</h3>
                    <div class="value">
                        {process.get('total_steps', 0)} / {process.get('max_steps', 0)} 步
                    </div>
                </div>
            </div>
            
            <div class="grid-2">
                <div>
                    <h3>能量演化曲线</h3>
                    <div class="chart-container">
                        <canvas id="energyEvolutionChart"></canvas>
                    </div>
                </div>
                
                <div>
                    <h3>力演化曲线</h3>
                    <div class="chart-container">
                        <canvas id="forceEvolutionChart"></canvas>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_final_results_section(self) -> str:
        """生成最终结果部分"""
        final = self.data.get('final_results', {})
        force_stats = final.get('force_statistics', {})
        
        return f"""
        <div class="section">
            <h2>🏁 最终结果</h2>
            
            <div class="grid-2">
                <div>
                    <h3>能量信息</h3>
                    <table class="data-table">
                        <tr>
                            <td>最终总能量</td>
                            <td>{final.get('final_energy', 0):.6f} eV</td>
                        </tr>
                        <tr>
                            <td>完成步数</td>
                            <td>{final.get('final_step', 0)}</td>
                        </tr>
                    </table>
                </div>
                
                <div>
                    <h3>力统计信息</h3>
                    <table class="data-table">
                        <tr>
                            <td>最大剩余力</td>
                            <td>{(final.get('max_residual_force') or 0):.4f} eV/Å</td>
                        </tr>
                        <tr>
                            <td>RMS剩余力</td>
                            <td>{(final.get('rms_residual_force') or 0):.4f} eV/Å</td>
                        </tr>
                        <tr>
                            <td>平均力大小</td>
                            <td>{(force_stats.get('mean_force_magnitude') or 0):.4f} eV/Å</td>
                        </tr>
                        <tr>
                            <td>力标准差</td>
                            <td>{(force_stats.get('std_force_magnitude') or 0):.4f} eV/Å</td>
                        </tr>
                    </table>
                </div>
            </div>
        </div>
        """
    
    def _generate_structure_comparison_section(self) -> str:
        """生成结构文件对比部分"""
        structure_files = self.data.get('structure_files', {})
        poscar_content = structure_files.get('poscar', '')
        contcar_content = structure_files.get('contcar', '')
        
        # 处理文件内容以便在HTML中显示
        def format_content(content: str, title: str) -> str:
            if not content:
                return f"<p>❌ 未找到{title}文件</p>"
            
            lines = content.split('\n')
            
            # 如果是CONTCAR文件，需要过滤掉原子速度部分
            if 'CONTCAR' in title:
                lines = self._filter_contcar_velocities(lines)
            
            # 限制显示行数避免页面过长
            if len(lines) > 50:
                display_content = '\n'.join(lines[:50]) + f'\n\n... (共{len(lines)}行，仅显示前50行)'
            else:
                display_content = '\n'.join(lines)
            return f"<pre><code>{display_content}</code></pre>"
        
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
    
    def _generate_footer(self) -> str:
        """生成页面底部"""
        return f"""
        <div class="footer">
            <p>📊 VASP结构优化分析报告 | 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
            <p>由VASP API结构优化可视化分析模块生成</p>
        </div>
        """
    
    def _generate_javascript(self) -> str:
        """生成JavaScript代码"""
        # 准备图表数据
        convergence = self.data.get('convergence_analysis', {})
        process = self.data.get('optimization_process', {})
        electronic = self.data.get('electronic_structure', {})
        
        force_history = convergence.get('force_convergence', {}).get('force_history', [])
        # 使用SCF能量序列作为能量图表的数据源
        energy_history = self.data.get('scf_energies', [])
        energy_profile = process.get('energy_profile', [])
        force_profile = process.get('force_profile', [])
        
        force_threshold = convergence.get('force_convergence', {}).get('threshold', 0.01)
        
        return f"""
        // 力收敛图表
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
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        beginAtZero: false,
                        title: {{
                            display: true,
                            text: '力 (eV/Å)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '离子步'
                        }}
                    }}
                }}
            }}
        }});
        
        // 能量收敛图表（横坐标为SCF步）
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
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        title: {{
                            display: true,
                            text: '能量 (eV)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: 'SCF步'
                        }}
                    }}
                }}
            }}
        }});
        
        // 能量演化图表（仍展示离子步对应的参考能量）
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
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        title: {{
                            display: true,
                            text: '能量 (eV)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '离子步（参考能量）'
                        }}
                    }}
                }}
            }}
        }});
        
        // 力演化图表
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
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        title: {{
                            display: true,
                            text: '力 (eV/Å)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '离子步'
                        }}
                    }}
                }}
            }}
        }});
        """


def generate_optimization_report(input_path: str, task_id: Optional[str] = None, output_dir: Optional[str] = None) -> str:
    """
    生成结构优化分析报告
    
    Args:
        input_path: OUTCAR文件路径或包含VASP文件的文件夹路径
        task_id: 任务ID（可选）
        output_dir: 输出目录，默认与输入路径同目录
        
    Returns:
        生成的HTML报告文件路径
    """
    try:
        # 分析VASP计算结果
        analyzer = OUTCARAnalyzer(input_path, task_id)
        analysis_data = analyzer.analyze()
        
        # 确定输出路径
        if output_dir is None:
            output_dir_path = analyzer.work_dir
        else:
            output_dir_path = Path(output_dir)
        
        # 生成HTML报告
        output_file = output_dir_path / "optimization_analysis_report.html"
        
        generator = OptimizationHTMLGenerator(analysis_data)
        html_path = generator.generate_html_report(str(output_file))
        
        return html_path
        
    except Exception as e:
        raise Exception(f"生成优化分析报告失败: {str(e)}")


if __name__ == "__main__":
    # 测试代码
    import sys
    
    # 测试文件夹路径（假设包含OUTCAR、POSCAR、CONTCAR）
    test_path = "/Users/ysl/Desktop/Code/MCP_IC_Tool"
    test_task_id = "test_optimization_001"
    print(f"🔍 测试路径: {test_path}")
    try:
        html_report = generate_optimization_report(test_path, test_task_id)
        print(f"✅ HTML报告已生成: {html_report}")
    except Exception as e:
        print(f"❌ 错误: {e}")
        import traceback
        traceback.print_exc()
