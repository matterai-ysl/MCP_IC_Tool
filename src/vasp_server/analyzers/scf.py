#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
VASP 自洽场 (SCF) 计算分析器 — 从 scf_analyzer.py 迁移。

包含:
- SCFAnalyzer: 解析 OUTCAR 并提取 SCF 计算数据
- SCFHTMLGenerator: 生成 HTML 可视化报告
- generate_scf_report: 入口便捷函数
"""

import re
import json
import math
import numpy as np
import os
import base64
import io
from pathlib import Path
from typing import Dict, List, Tuple, Any, Optional
from datetime import datetime

from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator
from .file_utils import find_vasp_file

# matplotlib导入（用于ELF可视化）
try:
    import matplotlib
    matplotlib.use('Agg')  # 使用非GUI后端
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# pymatgen导入（可选，用于ELFCAR高级分析）
try:
    from pymatgen.io.vasp import VolumetricData
    from pymatgen.core import Structure, Lattice
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    from pymatgen.analysis.local_env import CrystalNN
    from pymatgen.analysis.bond_valence import BVAnalyzer
    HAS_PYMATGEN = True
except ImportError:
    HAS_PYMATGEN = False


class SCFAnalyzer(BaseAnalyzer):
    """VASP自洽场计算结果分析器"""

    def __init__(self, input_path: str, task_id: Optional[str] = None):
        super().__init__(input_path, task_id)
        # 额外的VASP文件路径（SCF特有）
        self.chgcar_path = find_vasp_file(self.work_dir, "CHGCAR")
        self.doscar_path = find_vasp_file(self.work_dir, "DOSCAR")
        self.potcar_path = find_vasp_file(self.work_dir, "POTCAR")
        self.acf_path = find_vasp_file(self.work_dir, "ACF")
        self.elfcar_path = find_vasp_file(self.work_dir, "ELFCAR")

    def _init_data(self) -> Dict[str, Any]:
        return {
            'file_info': {},
            'calculation_settings': {},
            'electronic_convergence': {},
            'electronic_structure': {},
            'forces_and_stress': {},
            'magnetic_properties': {},
            'bader_analysis': {},
            'elfcar_analysis': {},
            'final_results': {},
            'task_info': {'task_id': self.task_id},
            'structure_files': {}
        }

    def _do_analyze(self):
        content = self._get_outcar_content()
        if not content:
            return

        # Parse SCF-specific fields from OUTCAR that the base class doesn't cover
        self._parse_scf_specific_file_info(content)

        # 解析电子步信息
        self._parse_electronic_steps(content)

        # 解析力和应力信息
        self._parse_forces_stress(content)

        # 解析磁性信息
        self._parse_magnetic_info(content)

        # 分析电子收敛性
        self._analyze_electronic_convergence()

        # 分析电子结构
        self._analyze_electronic_structure(content)

        # 分析Bader电荷
        self._analyze_bader_charges()

        # 分析ELFCAR电子局域函数
        self._analyze_elfcar()

        # 分析最终结果
        self._analyze_final_results()

    # ------------------------------------------------------------------ #
    #  SCF-specific file info (calculation_date, total_cores)
    # ------------------------------------------------------------------ #
    def _parse_scf_specific_file_info(self, content: str):
        """Parse SCF-specific file info fields that the base class doesn't extract."""
        info = self.data['file_info']

        # calculation_date (base class stores as run_date; SCF HTML expects calculation_date)
        date_match = re.search(r'executed on.*date ([\d\.]+)', content)
        if date_match:
            info['calculation_date'] = date_match.group(1)

        # total_cores
        cores_match = re.search(r'running on\s+(\d+)\s+total cores', content)
        if cores_match:
            info['total_cores'] = int(cores_match.group(1))

    # ------------------------------------------------------------------ #
    #  电子步解析
    # ------------------------------------------------------------------ #
    def _parse_electronic_steps(self, content: str):
        """解析电子步信息"""
        print("🔍 解析电子步信息...")

        # 查找所有电子步的能量
        electronic_steps = []

        lines = content.split('\n')
        current_iteration = 0
        current_step = 1

        # 记录当前迭代的信息
        current_toten = None
        current_energy_change = None
        current_magnetization = None

        for i, line in enumerate(lines):
            # 匹配迭代标志行
            iteration_match = re.match(r'---+\s+Iteration\s+\d+\(\s*(\d+)\s*\)\s+---+', line)
            if iteration_match:
                current_step = int(iteration_match.group(1))
                current_iteration += 1
                continue

            # 匹配能量变化行
            if 'total energy-change (2. order)' in line:
                # 直接从冒号后提取，处理各种格式： : 0.123E+04 或 :-0.123E+04
                colon_idx = line.find(':')
                if colon_idx != -1:
                    try:
                        # 提取冒号后的第一个数值部分
                        energy_part = line[colon_idx+1:].strip().split()[0]
                        # 处理科学计数法格式转换
                        energy_part = energy_part.replace('E+', 'e+').replace('E-', 'e-')
                        # 处理只有E的情况（如1.23E变成1.23e+0）
                        if energy_part.endswith('E'):
                            energy_part = energy_part.replace('E', 'e+0')
                        current_energy_change = float(energy_part)
                    except (ValueError, IndexError):
                        print(f"   ⚠️ 无法解析能量变化值: {line.strip()}")
                        continue
                continue

            # 匹配磁矩信息
            magnetization_match = re.match(r'\s*number of electron\s+[\d\.]+\s+magnetization\s+([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)', line)
            if magnetization_match:
                try:
                    current_magnetization = float(magnetization_match.group(1))
                except ValueError:
                    print(f"   ⚠️ 无法解析磁矩值: {line.strip()}")
                continue

            # 匹配总能量行
            toten_match = re.match(r'\s*free energy\s+TOTEN\s*=\s*([-+]?\d+\.?\d*(?:[eE][-+]?\d+)?)\s+eV', line)
            if toten_match:
                try:
                    current_toten = float(toten_match.group(1))

                    # 如果有完整的信息，添加到电子步列表
                    if current_toten is not None and current_energy_change is not None:
                        electronic_steps.append({
                            'step': current_step,
                            'iteration': current_iteration,
                            'free_energy': current_toten,
                            'energy_change': current_energy_change,
                            'magnetization': current_magnetization,
                            'energy_sigma0': current_toten  # TOTEN 就是 energy(sigma->0)
                        })

                    # 重置变量
                    current_energy_change = None
                    current_magnetization = None
                except ValueError:
                    print(f"   ⚠️ 无法解析总能量值: {line.strip()}")
                continue

        self.data['electronic_convergence']['electronic_steps'] = electronic_steps
        print(f"   📊 解析了 {len(electronic_steps)} 个电子步")

        # 提取最终的SCF能量序列
        self._extract_scf_energies(content)

    def _extract_scf_energies(self, content: str):
        """提取SCF能量序列"""
        scf_energies = []

        # 查找 "free energy    TOTEN" 模式
        energy_pattern = r'free\s+energy\s+TOTEN\s*=\s*([-\d\.E]+)\s*eV'
        energy_matches = re.findall(energy_pattern, content)

        for match in energy_matches:
            try:
                energy = float(match)
                scf_energies.append(energy)
            except ValueError:
                continue

        self.data['electronic_convergence']['scf_energies'] = scf_energies
        print(f"   📊 提取了 {len(scf_energies)} 个SCF能量点")

    # ------------------------------------------------------------------ #
    #  力和应力解析
    # ------------------------------------------------------------------ #
    def _parse_forces_stress(self, content: str):
        """解析力和应力信息"""
        print("🔍 解析力和应力信息...")

        # 解析最终的原子力
        forces = []
        force_section = re.search(r'POSITION\s+TOTAL-FORCE \(eV/Angst\)\s*\n\s*-+\s*\n(.*?)\n\s*-+', content, re.DOTALL)
        if force_section:
            force_lines = force_section.group(1).strip().split('\n')
            for line in force_lines:
                parts = line.strip().split()
                if len(parts) >= 6:
                    try:
                        px, py, pz = float(parts[0]), float(parts[1]), float(parts[2])
                        fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                        force_magnitude = math.sqrt(fx*fx + fy*fy + fz*fz)
                        forces.append({
                            'position': [px, py, pz],
                            'force': [fx, fy, fz],
                            'magnitude': force_magnitude
                        })
                    except ValueError:
                        continue

        # 解析应力张量
        stress_tensor = None
        stress_match = re.search(r'in kB\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)', content)
        if stress_match:
            stress_values = [float(x) for x in stress_match.groups()]
            stress_tensor = {
                'xx': stress_values[0],
                'yy': stress_values[1],
                'zz': stress_values[2],
                'xy': stress_values[3],
                'yz': stress_values[4],
                'zx': stress_values[5]
            }

        self.data['forces_and_stress'] = {
            'forces': forces,
            'stress_tensor': stress_tensor,
            'max_force': max([f['magnitude'] for f in forces]) if forces else None,
            'rms_force': math.sqrt(sum(f['magnitude']**2 for f in forces) / len(forces)) if forces else None
        }

        print(f"   📊 解析了 {len(forces)} 个原子的力信息")
        if stress_tensor:
            print("   ✅ 解析了应力张量")

    # ------------------------------------------------------------------ #
    #  磁性信息解析
    # ------------------------------------------------------------------ #
    def _parse_magnetic_info(self, content: str):
        """解析磁性信息"""
        print("🔍 解析磁性信息...")

        magnetic_data = {}

        # 检查是否是自旋极化计算
        ispin = self.data['calculation_settings'].get('ISPIN', 1)
        # handle string vs int
        try:
            ispin = int(ispin)
        except (ValueError, TypeError):
            ispin = 1
        is_spin_polarized = ispin == 2

        if is_spin_polarized:
            # 提取总磁矩
            total_mag_match = re.search(r'number of electron\s+(\d+\.\d+)\s+magnetization\s+([-\d\.]+)', content)
            if total_mag_match:
                magnetic_data['total_magnetization'] = float(total_mag_match.group(2))

            # 提取原子磁矩
            atom_mag_pattern = r'magnetization \(x\)\s*\n(.*?)total magnetization'
            atom_mag_match = re.search(atom_mag_pattern, content, re.DOTALL)
            if atom_mag_match:
                atom_mag_lines = atom_mag_match.group(1).strip().split('\n')
                atom_magnetizations = []
                for line in atom_mag_lines:
                    if line.strip() and not line.startswith('#'):
                        parts = line.strip().split()
                        if len(parts) >= 4:
                            try:
                                atom_idx = int(parts[0]) - 1  # 转为0索引
                                mag_x = float(parts[1])
                                mag_y = float(parts[2])
                                mag_z = float(parts[3])
                                total_mag = math.sqrt(mag_x**2 + mag_y**2 + mag_z**2)
                                atom_magnetizations.append({
                                    'atom_index': atom_idx,
                                    'magnetization': [mag_x, mag_y, mag_z],
                                    'total_magnitude': total_mag
                                })
                            except (ValueError, IndexError):
                                continue

                magnetic_data['atom_magnetizations'] = atom_magnetizations
                print(f"   📊 解析了 {len(atom_magnetizations)} 个原子的磁矩")

        magnetic_data['is_spin_polarized'] = is_spin_polarized
        self.data['magnetic_properties'] = magnetic_data

    # ------------------------------------------------------------------ #
    #  电子收敛性分析
    # ------------------------------------------------------------------ #
    def _analyze_electronic_convergence(self):
        """分析电子收敛性"""
        print("🧪 分析电子收敛性...")

        electronic_steps = self.data['electronic_convergence'].get('electronic_steps', [])
        scf_energies = self.data['electronic_convergence'].get('scf_energies', [])

        convergence_analysis = {
            'total_electronic_steps': len(electronic_steps),
            'energy_convergence': {},
            'converged': False
        }

        if scf_energies:
            # 分析能量收敛
            if len(scf_energies) > 1:
                energy_changes = [abs(scf_energies[i] - scf_energies[i-1]) for i in range(1, len(scf_energies))]
                final_energy_change = energy_changes[-1] if energy_changes else 0

                # 检查收敛标准
                ediff = self.data['calculation_settings'].get('EDIFF', 1e-4)
                try:
                    ediff = float(ediff)
                except (ValueError, TypeError):
                    ediff = 1e-4
                converged = final_energy_change < ediff

                convergence_analysis['energy_convergence'] = {
                    'final_energy': scf_energies[-1],
                    'final_energy_change': final_energy_change,
                    'ediff_threshold': ediff,
                    'converged': converged,
                    'energy_changes': energy_changes
                }

                convergence_analysis['converged'] = converged
                print(f"   📊 能量收敛: {converged}, 最终变化: {final_energy_change:.2e} eV")

        self.data['electronic_convergence']['analysis'] = convergence_analysis

    # ------------------------------------------------------------------ #
    #  电子结构分析
    # ------------------------------------------------------------------ #
    def _analyze_electronic_structure(self, content: str):
        """分析电子结构"""
        print("🔬 分析电子结构...")

        electronic_structure = {}

        # 提取费米能级
        fermi_pattern = r'E-fermi\s*:\s*([-\d\.]+)'
        fermi_match = re.search(fermi_pattern, content)
        if fermi_match:
            electronic_structure['fermi_energy'] = float(fermi_match.group(1))
            print(f"   ✅ 费米能级: {fermi_match.group(1)} eV")

        # 提取带隙信息
        gap_pattern = r'E-fermi\s*:\s*([-\d\.]+).*?gap\s*=\s*([\d\.]+)'
        gap_match = re.search(gap_pattern, content, re.DOTALL)
        if gap_match:
            electronic_structure['band_gap'] = float(gap_match.group(2))
            print(f"   ✅ 带隙: {gap_match.group(2)} eV")

        # 提取总电子数
        nelect_pattern = r'NELECT\s*=\s*([\d\.]+)'
        nelect_match = re.search(nelect_pattern, content)
        if nelect_match:
            electronic_structure['total_electrons'] = float(nelect_match.group(1))
            print(f"   ✅ 总电子数: {nelect_match.group(1)}")

        # 提取K点信息
        kpoint_pattern = r'Found\s+(\d+)\s+irreducible\s+k-points'
        kpoint_match = re.search(kpoint_pattern, content)
        if kpoint_match:
            electronic_structure['irreducible_kpoints'] = int(kpoint_match.group(1))
            print(f"   ✅ 不可约K点数: {kpoint_match.group(1)}")

        self.data['electronic_structure'] = electronic_structure

    # ------------------------------------------------------------------ #
    #  Bader 电荷分析
    # ------------------------------------------------------------------ #
    def _analyze_bader_charges(self):
        """分析Bader电荷"""
        print("🔋 分析Bader电荷...")

        bader_data = {
            'available': False,
            'atom_charges': [],
            'summary': {},
            'error_message': None
        }

        # 检查所需文件是否存在
        required_files = {
            'POSCAR': self.poscar_path,
            'POTCAR': self.potcar_path,
            'ACF.dat': self.acf_path
        }

        missing_files = []
        for file_name, file_path in required_files.items():
            if not file_path or not file_path.exists():
                missing_files.append(file_name)

        if missing_files:
            error_msg = f"缺少Bader分析所需文件: {', '.join(missing_files)}"
            print(f"   ⚠️ {error_msg}")
            bader_data['error_message'] = error_msg
            self.data['bader_analysis'] = bader_data
            return

        try:
            # 解析POSCAR获取原子信息
            elements, atom_list = self._parse_poscar_elements(self.poscar_path)  # type: ignore
            if not elements:
                raise Exception("无法解析POSCAR中的元素信息")

            # 解析POTCAR获取价电子数
            zval_map = self._parse_potcar_zval(elements, self.potcar_path)  # type: ignore
            if not zval_map:
                raise Exception("无法解析POTCAR中的价电子数信息")

            # 解析ACF.dat文件
            acf_data = self._parse_acf_file(self.acf_path)  # type: ignore
            if not acf_data:
                raise Exception("无法解析ACF.dat文件")

            # 计算电荷转移
            atom_charges = []
            total_charge_transfer = 0.0

            for acf_entry in acf_data:
                atom_index = acf_entry['atom_index']
                if atom_index <= len(atom_list):
                    element = atom_list[atom_index - 1]
                    zval = zval_map[element]
                    bader_charge = acf_entry['charge']
                    charge_transfer = zval - bader_charge

                    atom_charges.append({
                        'atom_index': atom_index,
                        'element': element,
                        'position': acf_entry['position'],
                        'zval': zval,
                        'bader_charge': bader_charge,
                        'charge_transfer': charge_transfer
                    })

                    total_charge_transfer += abs(charge_transfer)

            # 计算统计信息
            if atom_charges:
                charge_transfers = [atom['charge_transfer'] for atom in atom_charges]
                bader_charges = [atom['bader_charge'] for atom in atom_charges]

                summary = {
                    'total_atoms': len(atom_charges),
                    'total_charge_transfer': total_charge_transfer,
                    'avg_charge_transfer': sum(charge_transfers) / len(charge_transfers),
                    'max_charge_transfer': max(charge_transfers),
                    'min_charge_transfer': min(charge_transfers),
                    'avg_bader_charge': sum(bader_charges) / len(bader_charges),
                    'element_summary': self._calculate_element_summary(atom_charges)
                }

                bader_data.update({
                    'available': True,
                    'atom_charges': atom_charges,
                    'summary': summary
                })

                print(f"   ✅ 成功分析 {len(atom_charges)} 个原子的Bader电荷")
                print(f"   📊 总电荷转移: {total_charge_transfer:.4f}")

            else:
                bader_data['error_message'] = "未找到有效的原子电荷数据"

        except Exception as e:
            error_msg = f"Bader分析失败: {str(e)}"
            print(f"   ❌ {error_msg}")
            bader_data['error_message'] = error_msg

        self.data['bader_analysis'] = bader_data

    def _parse_poscar_elements(self, poscar_path: Path) -> Tuple[List[str], List[str]]:
        """解析POSCAR文件获取元素信息"""
        try:
            with open(poscar_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            if len(lines) < 7:
                return [], []

            # 第6行是元素符号（VASP 5.x格式）
            elements_line = lines[5].strip()
            # 第7行是原子数量
            counts_line = lines[6].strip()

            elements = elements_line.split()
            counts = [int(x) for x in counts_line.split()]

            if len(elements) != len(counts):
                # 尝试VASP 4.x格式（没有元素行）
                counts = [int(x) for x in lines[5].strip().split()]
                elements = [f"Elem{i+1}" for i in range(len(counts))]

            # 构建原子列表
            atom_list = []
            for elem, count in zip(elements, counts):
                atom_list.extend([elem] * count)

            return elements, atom_list

        except Exception as e:
            print(f"   ❌ 解析POSCAR失败: {e}")
            return [], []

    def _parse_potcar_zval(self, elements: List[str], potcar_path: Path) -> Dict[str, float]:
        """解析POTCAR文件获取价电子数"""
        zval_map = {}
        try:
            with open(potcar_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()

            # 按空行分割POTCAR的不同元素块
            blocks = re.split(r'\n\s*\n', content.strip())

            # 创建元素到ZVAL的映射
            potcar_zvals = {}
            for i, block in enumerate(blocks):
                if not block.strip():
                    continue

                # 查找VRHFIN或TITEL行获取元素名
                element_name = None

                # 先尝试VRHFIN行
                vrhfin_match = re.search(r'VRHFIN\s*=\s*([A-Za-z_]+)', block, re.IGNORECASE)
                if vrhfin_match:
                    element_name = vrhfin_match.group(1).strip()
                else:
                    # 如果没有VRHFIN，尝试TITEL行
                    titel_match = re.search(r'TITEL\s*=\s*\S*\s+([A-Za-z_]+)', block, re.IGNORECASE)
                    if titel_match:
                        element_name = titel_match.group(1).strip()

                if element_name:
                    # 提取元素符号（第一个大写字母 + 可能的小写字母）
                    element_match = re.match(r'([A-Z][a-z]?)', element_name)
                    if element_match:
                        element_symbol = element_match.group(1)

                        # 在同一块中查找ZVAL（可能在POMASS行中或独立行）
                        # 匹配格式：ZVAL = 3.000 或 ; ZVAL = 3.000
                        zval_match = re.search(r'[;\s]*ZVAL\s*=\s*([\d\.]+)', block, re.IGNORECASE)
                        if zval_match:
                            potcar_zvals[element_symbol] = float(zval_match.group(1))
                        else:
                            # 如果找不到ZVAL，记录到日志但不打印详细信息
                            pass
                    else:
                        print(f"   ⚠️ 无法从 '{element_name}' 中提取元素符号")

            # 为每个需要的元素获取ZVAL
            for element in elements:
                if element in potcar_zvals:
                    zval_map[element] = potcar_zvals[element]
                else:
                    # 如果找不到，使用常见元素的默认ZVAL
                    default_zvals = {
                        'H': 1.0, 'He': 2.0, 'Li': 1.0, 'Be': 2.0, 'B': 3.0, 'C': 4.0,
                        'N': 5.0, 'O': 6.0, 'F': 7.0, 'Ne': 8.0, 'Na': 1.0, 'Mg': 2.0,
                        'Al': 3.0, 'Si': 4.0, 'P': 5.0, 'S': 6.0, 'Cl': 7.0, 'Ar': 8.0,
                        'K': 1.0, 'Ca': 2.0, 'Sc': 3.0, 'Ti': 4.0, 'V': 5.0, 'Cr': 6.0,
                        'Mn': 7.0, 'Fe': 8.0, 'Co': 9.0, 'Ni': 10.0, 'Cu': 11.0, 'Zn': 12.0,
                        'Ga': 3.0, 'Ge': 4.0, 'As': 5.0, 'Se': 6.0, 'Br': 7.0, 'Kr': 8.0,
                        'Rb': 1.0, 'Sr': 2.0, 'Y': 3.0, 'Zr': 4.0, 'Nb': 5.0, 'Mo': 6.0,
                        'Tc': 7.0, 'Ru': 8.0, 'Rh': 9.0, 'Pd': 10.0, 'Ag': 11.0, 'Cd': 12.0,
                        'In': 3.0, 'Sn': 4.0, 'Sb': 5.0, 'Te': 6.0, 'I': 7.0, 'Xe': 8.0,
                        'Cs': 1.0, 'Ba': 2.0, 'La': 3.0, 'Ce': 4.0, 'Pr': 3.0, 'Nd': 3.0,
                        'Pm': 3.0, 'Sm': 3.0, 'Eu': 2.0, 'Gd': 3.0, 'Tb': 3.0, 'Dy': 3.0,
                        'Ho': 3.0, 'Er': 3.0, 'Tm': 3.0, 'Yb': 2.0, 'Lu': 3.0, 'Hf': 4.0,
                        'Ta': 5.0, 'W': 6.0, 'Re': 7.0, 'Os': 8.0, 'Ir': 9.0, 'Pt': 10.0,
                        'Au': 11.0, 'Hg': 2.0, 'Tl': 3.0, 'Pb': 4.0, 'Bi': 5.0
                    }
                    zval_map[element] = default_zvals.get(element, 4.0)
                    print(f"   ⚠️ 未在POTCAR中找到{element}的ZVAL，使用默认值: {zval_map[element]}")

            return zval_map

        except Exception as e:
            print(f"   ❌ 解析POTCAR失败: {e}")
            # 发生错误时返回默认值
            default_zvals = {
                'H': 1.0, 'He': 2.0, 'Li': 1.0, 'Be': 2.0, 'B': 3.0, 'C': 4.0,
                'N': 5.0, 'O': 6.0, 'F': 7.0, 'Ne': 8.0, 'Na': 1.0, 'Mg': 2.0,
                'Al': 3.0, 'Si': 4.0, 'P': 5.0, 'S': 6.0, 'Cl': 7.0, 'Ar': 8.0,
                'K': 1.0, 'Ca': 2.0, 'Ti': 4.0, 'V': 5.0, 'Cr': 6.0, 'Mn': 7.0,
                'Fe': 8.0, 'Co': 9.0, 'Ni': 10.0, 'Cu': 11.0, 'Zn': 12.0
            }
            return {element: default_zvals.get(element, 4.0) for element in elements}

    def _parse_acf_file(self, acf_path: Path) -> List[Dict[str, Any]]:
        """解析ACF.dat文件"""
        acf_data = []
        try:
            with open(acf_path, 'r', encoding='utf-8', errors='ignore') as f:
                lines = f.readlines()

            # 跳过头部，寻找数据行
            for line in lines:
                line = line.strip()
                if not line or line.startswith('#') or 'ATOM' in line.upper():
                    continue

                parts = line.split()
                if len(parts) >= 5:
                    try:
                        atom_index = int(parts[0])
                        x = float(parts[1])
                        y = float(parts[2])
                        z = float(parts[3])
                        charge = float(parts[4])

                        acf_data.append({
                            'atom_index': atom_index,
                            'position': [x, y, z],
                            'charge': charge
                        })
                    except (ValueError, IndexError):
                        continue

            return acf_data

        except Exception as e:
            print(f"   ❌ 解析ACF.dat失败: {e}")
            return []

    def _calculate_element_summary(self, atom_charges: List[Dict[str, Any]]) -> Dict[str, Dict[str, float]]:
        """计算各元素的电荷统计"""
        element_data = {}

        for atom in atom_charges:
            element = atom['element']
            if element not in element_data:
                element_data[element] = {
                    'count': 0,
                    'total_charge_transfer': 0.0,
                    'avg_charge_transfer': 0.0,
                    'avg_bader_charge': 0.0,
                    'charges': [],
                    'transfers': []
                }

            element_data[element]['count'] += 1
            element_data[element]['total_charge_transfer'] += abs(atom['charge_transfer'])
            element_data[element]['charges'].append(atom['bader_charge'])
            element_data[element]['transfers'].append(atom['charge_transfer'])

        # 计算平均值
        for element, data in element_data.items():
            data['avg_charge_transfer'] = sum(data['transfers']) / data['count']
            data['avg_bader_charge'] = sum(data['charges']) / data['count']
            # 清理临时列表
            del data['charges']
            del data['transfers']

        return element_data

    # ------------------------------------------------------------------ #
    #  最终结果分析
    # ------------------------------------------------------------------ #
    def _analyze_final_results(self):
        """分析最终结果"""
        print("📊 分析最终结果...")

        scf_energies = self.data['electronic_convergence'].get('scf_energies', [])
        convergence = self.data['electronic_convergence'].get('analysis', {})
        forces_data = self.data['forces_and_stress']
        electronic_structure = self.data['electronic_structure']

        final_results = {
            'final_energy': scf_energies[-1] if scf_energies else None,
            'converged': convergence.get('converged', False),
            'total_electronic_steps': convergence.get('total_electronic_steps', 0),
            'fermi_energy': electronic_structure.get('fermi_energy'),
            'band_gap': electronic_structure.get('band_gap'),
            'max_force': forces_data.get('max_force'),
            'rms_force': forces_data.get('rms_force'),
            'total_magnetization': self.data['magnetic_properties'].get('total_magnetization')
        }

        self.data['final_results'] = final_results

    # ------------------------------------------------------------------ #
    #  ELFCAR 分析
    # ------------------------------------------------------------------ #
    def _analyze_elfcar(self):
        """分析ELFCAR电子局域函数"""
        print("🧮 分析ELFCAR电子局域函数...")

        elfcar_data = {
            'available': False,
            'pymatgen_available': HAS_PYMATGEN,
            'matplotlib_available': HAS_MATPLOTLIB,
            'lattice_vectors': [],
            'grid_dimensions': [],
            'statistics': {},
            'crystal_analysis': {},
            'material_properties': {},
            'chemical_bonding': {},
            'visualization': {},
            'error_message': None
        }

        if not self.elfcar_path or not self.elfcar_path.exists():
            error_msg = "未找到ELFCAR文件"
            print(f"   ⚠️ {error_msg}")
            elfcar_data['error_message'] = error_msg
            self.data['elfcar_analysis'] = elfcar_data
            return

        try:
            if HAS_PYMATGEN:
                # 使用pymatgen进行高级分析
                elfcar_data = self._analyze_elfcar_with_pymatgen()
            else:
                # 基础ELFCAR分析
                elfcar_data = self._analyze_elfcar_basic()

            print(f"   ✅ 成功分析ELFCAR文件")

        except Exception as e:
            error_msg = f"ELFCAR分析失败: {str(e)}"
            print(f"   ❌ {error_msg}")
            elfcar_data['error_message'] = error_msg

        self.data['elfcar_analysis'] = elfcar_data

    def _analyze_elfcar_with_pymatgen(self) -> Dict[str, Any]:
        """使用pymatgen进行ELFCAR高级分析"""
        print("   🔧 使用pymatgen进行高级ELFCAR分析...")

        # 解析ELFCAR文件
        poscar, data_dict, data_aug = VolumetricData.parse_file(str(self.elfcar_path))
        structure = poscar.structure

        # 获取ELF数据
        if 'total' in data_dict:
            elf_3d = data_dict['total']
        else:
            first_key = list(data_dict.keys())[0]
            elf_3d = data_dict[first_key]

        grid_dims = elf_3d.shape

        # 基础统计分析
        statistics = {
            'min': float(np.min(elf_3d)),
            'max': float(np.max(elf_3d)),
            'mean': float(np.mean(elf_3d)),
            'std': float(np.std(elf_3d)),
            'median': float(np.median(elf_3d))
        }

        # 晶体学分析
        crystal_analysis = self._analyze_crystal_structure_pymatgen(structure)

        # 材料属性分析
        material_properties = self._analyze_material_properties_pymatgen(structure, elf_3d)

        # 化学键分析
        chemical_bonding = self._analyze_chemical_bonding(elf_3d, statistics)

        # 生成2D可视化
        visualization = self._generate_elf_visualization(elf_3d, statistics)

        return {
            'available': True,
            'pymatgen_available': True,
            'matplotlib_available': HAS_MATPLOTLIB,
            'lattice_vectors': structure.lattice.matrix.tolist(),
            'grid_dimensions': list(grid_dims),
            'statistics': statistics,
            'crystal_analysis': crystal_analysis,
            'material_properties': material_properties,
            'chemical_bonding': chemical_bonding,
            'visualization': visualization,
            'error_message': None
        }

    def _analyze_elfcar_basic(self) -> Dict[str, Any]:
        """基础ELFCAR分析（不使用pymatgen）"""
        print("   🛠️ 进行基础ELFCAR分析...")

        if not self.elfcar_path:
            raise Exception("ELFCAR文件路径未设置")

        with open(self.elfcar_path, 'r', encoding='utf-8', errors='ignore') as f:
            lines = f.readlines()

        # 解析基础信息
        scale = float(lines[1].strip())

        # 晶格向量
        lattice_vectors = []
        for i in range(2, 5):
            lattice_vector = [float(x) * scale for x in lines[i].strip().split()]
            lattice_vectors.append(lattice_vector)

        # 找到网格维度
        grid_line_idx = None
        for i, line in enumerate(lines):
            line = line.strip()
            if line and len(line.split()) == 3:
                try:
                    nx, ny, nz = int(line.split()[0]), int(line.split()[1]), int(line.split()[2])
                    if 10 <= nx <= 1000 and 10 <= ny <= 1000 and 10 <= nz <= 1000:
                        grid_line_idx = i
                        break
                except ValueError:
                    continue

        if grid_line_idx is None:
            raise Exception("无法找到网格维度信息")

        nx, ny, nz = [int(x) for x in lines[grid_line_idx].strip().split()]

        # 读取ELF数据
        elf_data = []
        for i in range(grid_line_idx + 1, len(lines)):
            line = lines[i].strip()
            if line:
                values = line.split()
                for value_str in values:
                    try:
                        value = float(value_str)
                        elf_data.append(value)
                    except ValueError:
                        continue

        # 统计分析
        if len(elf_data) >= nx * ny * nz:
            elf_data = elf_data[:nx * ny * nz]
            elf_array = np.array(elf_data)

            statistics = {
                'min': float(np.min(elf_array)),
                'max': float(np.max(elf_array)),
                'mean': float(np.mean(elf_array)),
                'std': float(np.std(elf_array)),
                'median': float(np.median(elf_array))
            }

            # 基础化学键分析
            elf_3d = elf_array.reshape(nx, ny, nz)
            chemical_bonding = self._analyze_chemical_bonding(elf_3d, statistics)

            # 生成2D可视化
            visualization = self._generate_elf_visualization(elf_3d, statistics)
        else:
            statistics = {'error': 'ELF数据不完整'}
            chemical_bonding = {}
            visualization = {}

        return {
            'available': True,
            'pymatgen_available': False,
            'matplotlib_available': HAS_MATPLOTLIB,
            'lattice_vectors': lattice_vectors,
            'grid_dimensions': [nx, ny, nz],
            'statistics': statistics,
            'crystal_analysis': {'basic_info': '需要pymatgen进行详细分析'},
            'material_properties': {'basic_info': '需要pymatgen进行详细分析'},
            'chemical_bonding': chemical_bonding,
            'visualization': visualization,
            'error_message': None
        }

    def _analyze_crystal_structure_pymatgen(self, structure) -> Dict[str, Any]:
        """使用pymatgen进行晶体结构分析"""
        analysis = {}

        try:
            # 基础信息
            analysis['formula'] = str(structure.composition)
            analysis['reduced_formula'] = structure.composition.reduced_formula
            analysis['density'] = structure.density
            analysis['volume'] = structure.volume
            analysis['num_sites'] = len(structure)

            # 晶格参数
            lattice = structure.lattice
            analysis['lattice_params'] = {
                'a': lattice.a,
                'b': lattice.b,
                'c': lattice.c,
                'alpha': lattice.alpha,
                'beta': lattice.beta,
                'gamma': lattice.gamma,
                'volume': lattice.volume
            }

            # 对称性分析
            try:
                sga = SpacegroupAnalyzer(structure)
                analysis['symmetry'] = {
                    'spacegroup_symbol': sga.get_space_group_symbol(),
                    'spacegroup_number': sga.get_space_group_number(),
                    'crystal_system': sga.get_crystal_system(),
                    'lattice_type': sga.get_lattice_type(),
                    'point_group': sga.get_point_group_symbol()
                }
                print(f"   🔍 空间群: {analysis['symmetry']['spacegroup_symbol']} ({analysis['symmetry']['spacegroup_number']})")
            except Exception as e:
                print(f"   ⚠️ 对称性分析失败: {e}")
                analysis['symmetry'] = {'error': str(e)}

            # 配位环境分析
            try:
                cn_analyzer = CrystalNN()
                coord_numbers = []
                for i, site in enumerate(structure):
                    try:
                        cn = cn_analyzer.get_cn(structure, i)
                        coord_numbers.append(cn)
                    except:
                        coord_numbers.append(None)

                valid_cn = [cn for cn in coord_numbers if cn is not None]
                if valid_cn:
                    analysis['coordination'] = {
                        'coordination_numbers': coord_numbers,
                        'average_cn': np.mean(valid_cn),
                        'cn_range': [min(valid_cn), max(valid_cn)]
                    }
                    print(f"   🔗 平均配位数: {analysis['coordination']['average_cn']:.2f}")
            except Exception as e:
                print(f"   ⚠️ 配位数分析失败: {e}")

        except Exception as e:
            print(f"   ❌ 晶体结构分析失败: {e}")
            analysis['error'] = str(e)

        return analysis

    def _analyze_material_properties_pymatgen(self, structure, elf_3d) -> Dict[str, Any]:
        """使用pymatgen进行材料属性分析"""
        properties = {}

        try:
            composition = structure.composition

            # 基础物理性质
            properties['molecular_weight'] = composition.weight
            properties['electrons_per_formula'] = composition.total_electrons

            # 磁性分析
            magnetic_elements = ['Fe', 'Co', 'Ni', 'Mn', 'Cr', 'V', 'Ti']
            has_magnetic = any(str(el) in magnetic_elements for el in composition.elements)

            if has_magnetic:
                properties['magnetic_elements'] = [str(el) for el in composition.elements
                                                   if str(el) in magnetic_elements]
                print(f"   🧲 发现磁性元素: {properties['magnetic_elements']}")
            else:
                properties['magnetic_analysis'] = "非磁性材料"

            # 电负性分析
            total_electronegativity = sum(el.X * amt for el, amt in composition.items())
            avg_electronegativity = total_electronegativity / composition.num_atoms

            electronegativities = [el.X for el in composition.elements]
            electronegativity_diff = max(electronegativities) - min(electronegativities)

            properties['electronegativity_analysis'] = {
                'average': avg_electronegativity,
                'range': electronegativity_diff,
                'ionic_character': min(electronegativity_diff / 3.0, 1.0)
            }

            print(f"   ⚡ 平均电负性: {avg_electronegativity:.2f}")
            print(f"   ⚡ 离子性: {properties['electronegativity_analysis']['ionic_character']:.2f}")

        except Exception as e:
            print(f"   ❌ 材料属性分析失败: {e}")
            properties['error'] = str(e)

        return properties

    def _analyze_chemical_bonding(self, elf_data, statistics) -> Dict[str, Any]:
        """分析化学键特征"""
        bonding = {}

        try:
            # ELF值分布分析
            if isinstance(elf_data, np.ndarray):
                high_localization = np.sum(elf_data > 0.8) / elf_data.size
                metallic_character = np.sum(elf_data < 0.2) / elf_data.size
                covalent_character = np.sum((elf_data >= 0.4) & (elf_data <= 0.8)) / elf_data.size
            else:
                # 如果是1D数组
                elf_array = np.array(elf_data)
                high_localization = np.sum(elf_array > 0.8) / elf_array.size
                metallic_character = np.sum(elf_array < 0.2) / elf_array.size
                covalent_character = np.sum((elf_array >= 0.4) & (elf_array <= 0.8)) / elf_array.size

            bonding['elf_distribution'] = {
                'high_localization_fraction': high_localization,
                'metallic_character_fraction': metallic_character,
                'covalent_character_fraction': covalent_character,
                'average_elf': statistics.get('mean', 0)
            }

            # 键合类型预测
            if high_localization > 0.15:
                bond_type = "强共价键特征"
                bond_description = "高度局域化的电子，典型的共价键体系"
            elif metallic_character > 0.6:
                bond_type = "金属键特征"
                bond_description = "低局域化电子，典型的金属键体系"
            elif covalent_character > 0.4:
                bond_type = "共价键特征"
                bond_description = "中等局域化电子，共价键为主"
            else:
                bond_type = "离子键特征"
                bond_description = "电子局域化分布不均，可能为离子键"

            bonding['predicted_bonding'] = {
                'type': bond_type,
                'description': bond_description,
                'confidence': max(float(high_localization), float(metallic_character), float(covalent_character))
            }

            print(f"   🔗 高局域化区域: {high_localization*100:.1f}%")
            print(f"   🔗 金属性特征: {metallic_character*100:.1f}%")
            print(f"   🔗 共价键特征: {covalent_character*100:.1f}%")
            print(f"   🎯 预测键合类型: {bond_type}")

        except Exception as e:
            print(f"   ❌ 化学键分析失败: {e}")
            bonding['error'] = str(e)

        return bonding

    def _generate_elf_visualization(self, elf_3d, statistics) -> Dict[str, Any]:
        """生成ELF 2D可视化图像"""
        visualization = {
            'available': False,
            'slices': {},
            'error_message': None
        }

        if not HAS_MATPLOTLIB:
            visualization['error_message'] = "matplotlib不可用，无法生成可视化图像"
            return visualization

        try:
            print("   🎨 生成ELF 2D可视化图像...")

            # 设置中文字体（如果需要的话）
            plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
            plt.rcParams['axes.unicode_minus'] = False

            # 获取维度
            nx, ny, nz = elf_3d.shape

            # 生成三个方向的切面图像
            slices = {}

            # XY切面 (Z方向中间)
            z_mid = nz // 2
            xy_slice = elf_3d[:, :, z_mid]
            xy_image = self._create_slice_image(xy_slice, f'XY Plane (Z = {z_mid}/{nz-1})', 'viridis', statistics)
            if xy_image:
                slices['xy'] = {
                    'title': f'XY Plane (Z = {z_mid})',
                    'image_base64': xy_image,
                    'slice_index': z_mid,
                    'total_slices': nz
                }

            # XZ切面 (Y方向中间)
            y_mid = ny // 2
            xz_slice = elf_3d[:, y_mid, :]
            xz_image = self._create_slice_image(xz_slice, f'XZ Plane (Y = {y_mid}/{ny-1})', 'plasma', statistics)
            if xz_image:
                slices['xz'] = {
                    'title': f'XZ Plane (Y = {y_mid})',
                    'image_base64': xz_image,
                    'slice_index': y_mid,
                    'total_slices': ny
                }

            # YZ切面 (X方向中间)
            x_mid = nx // 2
            yz_slice = elf_3d[x_mid, :, :]
            yz_image = self._create_slice_image(yz_slice, f'YZ Plane (X = {x_mid}/{nx-1})', 'hot', statistics)
            if yz_image:
                slices['yz'] = {
                    'title': f'YZ Plane (X = {x_mid})',
                    'image_base64': yz_image,
                    'slice_index': x_mid,
                    'total_slices': nx
                }

            visualization['available'] = True
            visualization['slices'] = slices

            print(f"   ✅ 成功生成 {len(slices)} 个ELF切面图像")

        except Exception as e:
            error_msg = f"ELF可视化生成失败: {str(e)}"
            print(f"   ❌ {error_msg}")
            visualization['error_message'] = error_msg

        return visualization

    def _create_slice_image(self, slice_data, title, colormap, statistics) -> str:
        """创建单个切面图像并返回base64编码"""
        try:
            fig, ax = plt.subplots(figsize=(8, 6))

            # 创建图像
            im = ax.imshow(slice_data, origin='lower', cmap=colormap, aspect='auto')
            ax.set_title(title, fontsize=14, fontweight='bold')
            ax.set_xlabel('Grid Index', fontsize=12)
            ax.set_ylabel('Grid Index', fontsize=12)

            # 添加颜色条
            cbar = plt.colorbar(im, ax=ax, shrink=0.8)
            cbar.set_label('ELF Value', fontsize=12, rotation=270, labelpad=20)

            # 添加统计信息
            stats_text = f"Min: {statistics.get('min', 0):.3f}\nMax: {statistics.get('max', 1):.3f}\nMean: {statistics.get('mean', 0.5):.3f}"
            ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=10,
                    verticalalignment='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

            plt.tight_layout()

            # 保存为base64
            buffer = io.BytesIO()
            plt.savefig(buffer, format='png', dpi=150, bbox_inches='tight')
            buffer.seek(0)
            image_base64 = base64.b64encode(buffer.getvalue()).decode()

            plt.close(fig)
            return image_base64

        except Exception as e:
            print(f"   ⚠️ 创建切面图像失败: {e}")
            return ""


class SCFHTMLGenerator(BaseHTMLGenerator):
    """SCF计算HTML报告生成器"""

    report_title = "VASP自洽场(SCF)计算分析报告"
    report_emoji = "⚛️"
    report_heading = "VASP自洽场(SCF)计算分析报告"
    report_attribution = "由VASP API自洽场可视化分析模块生成"
    header_gradient = "linear-gradient(135deg, #4a90e2 0%, #7b68ee 100%)"

    def _generate_body_sections(self) -> str:
        return (
            self._generate_summary()
            + self._generate_electronic_convergence_section()
            + self._generate_electronic_structure_section()
            + self._generate_forces_stress_section()
            + self._generate_magnetic_section()
            + self._generate_bader_section()
            + self._generate_elfcar_section()
            + self._generate_final_results_section()
        )

    # ------------------------------------------------------------------ #
    #  摘要
    # ------------------------------------------------------------------ #
    def _generate_summary(self) -> str:
        """生成摘要部分"""
        final = self.data.get('final_results', {})
        convergence = self.data['electronic_convergence'].get('analysis', {})

        overall_converged = final.get('converged', False)
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
                    <h3>电子步数</h3>
                    <div class="value">
                        {final.get('total_electronic_steps', 0)}
                    </div>
                </div>

                <div class="summary-card">
                    <h3>总能量</h3>
                    <div class="value">
                        {(final.get('final_energy') or 0):.6f} eV
                    </div>
                </div>

                <div class="summary-card">
                    <h3>费米能级</h3>
                    <div class="value">
                        {(final.get('fermi_energy') or 0):.4f} eV
                    </div>
                </div>
            </div>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  电子收敛性
    # ------------------------------------------------------------------ #
    def _generate_electronic_convergence_section(self) -> str:
        """生成电子收敛性分析部分"""
        convergence = self.data['electronic_convergence'].get('analysis', {})
        energy_conv = convergence.get('energy_convergence', {})

        energy_status = '✅ 收敛' if energy_conv.get('converged', False) else '❌ 未收敛'

        electronic_steps = self.data['electronic_convergence'].get('electronic_steps', [])
        scf_energies = self.data['electronic_convergence'].get('scf_energies', [])

        elec_dl = self._csv_download_link(
            'electronic_steps.csv',
            ['电子步', '自由能(eV)', '能量变化(eV)'],
            [[s['step'], s.get('free_energy'), s.get('energy_change')] for s in electronic_steps],
        )
        scf_dl = self._csv_download_link(
            'scf_energy.csv',
            ['SCF步', '总能量(eV)'],
            [[i + 1, v] for i, v in enumerate(scf_energies)],
        )

        return f"""
        <div class="section">
            <h2>⚡ 电子收敛性分析</h2>

            <div class="grid-2">
                <div>
                    <h3>能量收敛详情</h3>
                    <p><strong>状态:</strong> {energy_status}</p>
                    <p><strong>收敛阈值(EDIFF):</strong> {(energy_conv.get('ediff_threshold') or 0):.2e} eV</p>
                    <p><strong>最终能量变化:</strong> {(energy_conv.get('final_energy_change') or 0):.2e} eV</p>
                    <p><strong>最终能量:</strong> {(energy_conv.get('final_energy') or 0):.6f} eV</p>
                </div>

                <div>
                    <h3>收敛统计</h3>
                    <p><strong>总电子步数:</strong> {convergence.get('total_electronic_steps', 0)}</p>
                    <p><strong>SCF能量点数:</strong> {len(scf_energies)}</p>
                </div>
            </div>

            <div class="grid-2">
                <div>
                    <h3>电子步能量变化</h3>
                    <p style="color: #718096; font-size: 0.9em;">
                        🔴 红色：能量增加 &nbsp;&nbsp; 🟢 绿色：能量降低
                    </p>
                    <div class="chart-container">
                        <canvas id="electronicStepsChart"></canvas>
                    </div>
                    {elec_dl}
                </div>

                <div>
                    <h3>SCF能量收敛</h3>
                    <div class="chart-container">
                        <canvas id="scfEnergyChart"></canvas>
                    </div>
                    {scf_dl}
                </div>
            </div>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  电子结构
    # ------------------------------------------------------------------ #
    def _generate_electronic_structure_section(self) -> str:
        """生成电子结构分析部分"""
        electronic = self.data.get('electronic_structure', {})

        return f"""
        <div class="section">
            <h2>🔬 电子结构分析</h2>

            <div class="grid-3">
                <div class="summary-card">
                    <h3>费米能级位置</h3>
                    <div class="value">
                        {(electronic.get('fermi_energy') or 0):.4f} eV
                    </div>
                </div>

                <div class="summary-card">
                    <h3>带隙</h3>
                    <div class="value">
                        {(electronic.get('band_gap') or 0):.4f} eV
                    </div>
                </div>

                <div class="summary-card">
                    <h3>总电子数</h3>
                    <div class="value">
                        {(electronic.get('total_electrons') or 0):.2f}
                    </div>
                </div>
            </div>

            <table class="data-table">
                <tr>
                    <th>属性</th>
                    <th>值</th>
                    <th>说明</th>
                </tr>
                <tr>
                    <td>费米能级</td>
                    <td>{(electronic.get('fermi_energy') or 0):.4f} eV</td>
                    <td>电子填充的化学势</td>
                </tr>
                <tr>
                    <td>带隙</td>
                    <td>{(electronic.get('band_gap') or 0):.4f} eV</td>
                    <td>价带顶到导带底的能量差</td>
                </tr>
                <tr>
                    <td>总电子数</td>
                    <td>{(electronic.get('total_electrons') or 0):.2f}</td>
                    <td>系统中的总电子数</td>
                </tr>
                <tr>
                    <td>不可约K点数</td>
                    <td>{electronic.get('irreducible_kpoints', 'N/A')}</td>
                    <td>布里渊区采样点数</td>
                </tr>
            </table>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  力和应力
    # ------------------------------------------------------------------ #
    def _generate_forces_stress_section(self) -> str:
        """生成力和应力分析部分"""
        forces_data = self.data.get('forces_and_stress', {})
        forces = forces_data.get('forces', [])
        stress_tensor = forces_data.get('stress_tensor')

        return f"""
        <div class="section">
            <h2>⚖️ 原子受力与应力分析</h2>

            <div class="grid-2">
                <div>
                    <h3>力统计信息</h3>
                    <table class="data-table">
                        <tr>
                            <td>原子数量</td>
                            <td>{len(forces)}</td>
                        </tr>
                        <tr>
                            <td>最大力</td>
                            <td>{(forces_data.get('max_force') or 0):.4f} eV/Å</td>
                        </tr>
                        <tr>
                            <td>RMS力</td>
                            <td>{(forces_data.get('rms_force') or 0):.4f} eV/Å</td>
                        </tr>
                    </table>

                    <div class="chart-container">
                        <canvas id="forcesChart"></canvas>
                    </div>
                    {self._csv_download_link(
                        'atomic_forces.csv',
                        ['原子索引', '力大小(eV/Å)', 'Fx(eV/Å)', 'Fy(eV/Å)', 'Fz(eV/Å)'],
                        [[i + 1, f.get('magnitude'), f.get('fx'), f.get('fy'), f.get('fz')] for i, f in enumerate(forces)],
                    )}
                </div>

                <div>
                    <h3>应力张量 (kB)</h3>
                    {self._generate_stress_table(stress_tensor)}
                </div>
            </div>
        </div>
        """

    def _generate_stress_table(self, stress_tensor) -> str:
        """生成应力张量表格"""
        if not stress_tensor:
            return "<p>❌ 未找到应力张量信息</p>"

        return f"""
        <table class="data-table">
            <tr>
                <th>分量</th>
                <th>值 (kB)</th>
            </tr>
            <tr>
                <td>σ<sub>xx</sub></td>
                <td>{stress_tensor.get('xx', 0):.4f}</td>
            </tr>
            <tr>
                <td>σ<sub>yy</sub></td>
                <td>{stress_tensor.get('yy', 0):.4f}</td>
            </tr>
            <tr>
                <td>σ<sub>zz</sub></td>
                <td>{stress_tensor.get('zz', 0):.4f}</td>
            </tr>
            <tr>
                <td>σ<sub>xy</sub></td>
                <td>{stress_tensor.get('xy', 0):.4f}</td>
            </tr>
            <tr>
                <td>σ<sub>yz</sub></td>
                <td>{stress_tensor.get('yz', 0):.4f}</td>
            </tr>
            <tr>
                <td>σ<sub>zx</sub></td>
                <td>{stress_tensor.get('zx', 0):.4f}</td>
            </tr>
        </table>
        """

    # ------------------------------------------------------------------ #
    #  磁性分析
    # ------------------------------------------------------------------ #
    def _generate_magnetic_section(self) -> str:
        """生成磁性分析部分"""
        magnetic = self.data.get('magnetic_properties', {})
        is_spin_polarized = magnetic.get('is_spin_polarized', False)

        if not is_spin_polarized:
            return f"""
            <div class="section">
                <h2>🧲 磁性分析</h2>
                <p>⚠️ 该计算未考虑自旋极化 (ISPIN=1)</p>
            </div>
            """

        total_mag = magnetic.get('total_magnetization', 0)
        atom_mags = magnetic.get('atom_magnetizations', [])

        return f"""
        <div class="section">
            <h2>🧲 磁性分析</h2>

            <div class="grid-2">
                <div>
                    <h3>总磁矩</h3>
                    <div class="summary-card">
                        <div class="value">{total_mag:.4f} μ<sub>B</sub></div>
                    </div>

                    <h3>磁矩统计</h3>
                    <table class="data-table">
                        <tr>
                            <td>自旋极化</td>
                            <td>是 (ISPIN=2)</td>
                        </tr>
                        <tr>
                            <td>总磁矩</td>
                            <td>{total_mag:.4f} μ<sub>B</sub></td>
                        </tr>
                        <tr>
                            <td>原子磁矩数</td>
                            <td>{len(atom_mags)}</td>
                        </tr>
                    </table>
                </div>

                <div>
                    <h3>原子磁矩分布</h3>
                    <div class="chart-container">
                        <canvas id="magnetizationChart"></canvas>
                    </div>
                    {self._csv_download_link(
                        'atom_magnetization.csv',
                        ['原子索引', '磁矩大小(μB)'],
                        [[mag['atom_index'] + 1, mag['total_magnitude']] for mag in atom_mags],
                    )}
                </div>
            </div>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  Bader 电荷分析
    # ------------------------------------------------------------------ #
    def _generate_bader_section(self) -> str:
        """生成Bader电荷分析部分"""
        bader_data = self.data.get('bader_analysis', {})

        if not bader_data.get('available', False):
            error_msg = bader_data.get('error_message', '未进行Bader分析')
            return f"""
            <div class="section">
                <h2>🔋 Bader电荷分析</h2>
                <p>⚠️ {error_msg}</p>
                <p>📝 提示：Bader电荷分析需要POSCAR、POTCAR和ACF.dat文件</p>
            </div>
            """

        summary = bader_data.get('summary', {})
        atom_charges = bader_data.get('atom_charges', [])
        element_summary = summary.get('element_summary', {})

        # 生成元素统计表格
        element_table_rows = ""
        for element, data in element_summary.items():
            element_table_rows += f"""
            <tr>
                <td>{element}</td>
                <td>{data['count']}</td>
                <td>{data['avg_bader_charge']:.4f}</td>
                <td>{data['avg_charge_transfer']:.4f}</td>
                <td>{data['total_charge_transfer']:.4f}</td>
            </tr>
            """

        return f"""
        <div class="section">
            <h2>🔋 Bader电荷分析</h2>

            <div class="grid-2">
                <div>
                    <h3>Bader电荷统计</h3>
                    <table class="data-table">
                        <tr>
                            <td>分析原子总数</td>
                            <td>{summary.get('total_atoms', 0)}</td>
                        </tr>
                        <tr>
                            <td>总电荷转移</td>
                            <td>{summary.get('total_charge_transfer', 0):.4f} e</td>
                        </tr>
                        <tr>
                            <td>平均电荷转移</td>
                            <td>{summary.get('avg_charge_transfer', 0):.4f} e</td>
                        </tr>
                        <tr>
                            <td>最大电荷转移</td>
                            <td>{summary.get('max_charge_transfer', 0):.4f} e</td>
                        </tr>
                        <tr>
                            <td>最小电荷转移</td>
                            <td>{summary.get('min_charge_transfer', 0):.4f} e</td>
                        </tr>
                        <tr>
                            <td>平均Bader电荷</td>
                            <td>{summary.get('avg_bader_charge', 0):.4f} e</td>
                        </tr>
                    </table>
                </div>

                <div>
                    <h3>Bader电荷分布</h3>
                    <div class="chart-container">
                        <canvas id="baderChargeChart"></canvas>
                    </div>
                    {self._csv_download_link(
                        'bader_charges.csv',
                        ['原子索引', '元素', 'Bader电荷(e)', '电荷转移(e)'],
                        [[a['atom_index'], a.get('element'), a.get('bader_charge'), a.get('charge_transfer')] for a in atom_charges],
                    )}
                </div>
            </div>

            <div class="grid-2">
                <div>
                    <h3>各元素电荷统计</h3>
                    <table class="data-table">
                        <tr>
                            <th>元素</th>
                            <th>原子数</th>
                            <th>平均Bader电荷</th>
                            <th>平均电荷转移</th>
                            <th>总电荷转移</th>
                        </tr>
                        {element_table_rows}
                    </table>
                </div>

                <div>
                    <h3>电荷转移分布</h3>
                    <div class="chart-container">
                        <canvas id="chargeTransferChart"></canvas>
                    </div>
                    {self._csv_download_link(
                        'charge_transfer.csv',
                        ['原子索引', '元素', '电荷转移(e)'],
                        [[a['atom_index'], a.get('element'), a.get('charge_transfer')] for a in atom_charges],
                    )}
                </div>
            </div>

            <h3>详细原子电荷信息</h3>
            <div style="max-height: 400px; overflow-y: auto;">
                <table class="data-table">
                    <tr>
                        <th>原子索引</th>
                        <th>元素</th>
                        <th>位置 (x, y, z)</th>
                        <th>价电子数</th>
                        <th>Bader电荷</th>
                        <th>电荷转移</th>
                    </tr>
                    {self._generate_atom_charge_rows(atom_charges)}
                </table>
            </div>
        </div>
        """

    def _generate_atom_charge_rows(self, atom_charges: List[Dict[str, Any]]) -> str:
        """生成原子电荷表格行"""
        rows = ""
        for atom in atom_charges:
            position_str = f"({atom['position'][0]:.3f}, {atom['position'][1]:.3f}, {atom['position'][2]:.3f})"

            # 根据电荷转移的正负设置颜色
            charge_transfer = atom['charge_transfer']
            if charge_transfer > 0.1:
                color_class = "style='color: #e53e3e;'"  # 红色：失去电子
            elif charge_transfer < -0.1:
                color_class = "style='color: #3182ce;'"  # 蓝色：获得电子
            else:
                color_class = ""

            rows += f"""
            <tr>
                <td>{atom['atom_index']}</td>
                <td><strong>{atom['element']}</strong></td>
                <td>{position_str}</td>
                <td>{atom['zval']:.1f}</td>
                <td>{atom['bader_charge']:.4f}</td>
                <td {color_class}>{charge_transfer:+.4f}</td>
            </tr>
            """

        return rows

    # ------------------------------------------------------------------ #
    #  ELFCAR 分析
    # ------------------------------------------------------------------ #
    def _generate_elfcar_section(self) -> str:
        """生成ELFCAR电子局域函数分析部分"""
        elfcar_data = self.data.get('elfcar_analysis', {})

        if not elfcar_data.get('available', False):
            error_msg = elfcar_data.get('error_message', '未进行ELFCAR分析')
            pymatgen_status = "✅ Available" if elfcar_data.get('pymatgen_available', False) else "❌ Not Available"
            matplotlib_status = "✅ Available" if elfcar_data.get('matplotlib_available', False) else "❌ Not Available"
            return f"""
            <div class="section">
                <h2>🧮 Electron Localization Function (ELF) Analysis</h2>
                <p>⚠️ {error_msg}</p>
                <p>📝 Note: ELF analysis requires ELFCAR file (generated by setting LELF=.TRUE. in VASP)</p>
                <p>🔧 pymatgen status: {pymatgen_status}</p>
                <p>🎨 matplotlib status: {matplotlib_status}</p>
                <div style="background: #fff3cd; padding: 15px; border-radius: 5px; margin-top: 15px; border-left: 4px solid #ffc107;">
                    <p style="margin: 0; font-size: 0.9em;">
                        <strong>🔬 For detailed ELF analysis, please use professional tools:</strong><br>
                        • <strong>VESTA</strong> - 3D visualization and analysis<br>
                        • <strong>VMD</strong> - Advanced molecular visualization<br>
                        • <strong>ChemCraft</strong> - Chemical data analysis<br>
                        • <strong>pymatgen</strong> - Python materials analysis library<br>
                        • <strong>ASE</strong> - Atomic Simulation Environment
                    </p>
                </div>
            </div>
            """

        statistics = elfcar_data.get('statistics', {})
        crystal_analysis = elfcar_data.get('crystal_analysis', {})
        material_properties = elfcar_data.get('material_properties', {})
        chemical_bonding = elfcar_data.get('chemical_bonding', {})
        visualization = elfcar_data.get('visualization', {})
        grid_dims = elfcar_data.get('grid_dimensions', [])

        # 基础统计信息
        stats_html = self._generate_elf_statistics_table(statistics, grid_dims)

        # 2D可视化
        visualization_html = self._generate_elf_visualization_html(visualization)

        # 晶体学分析
        crystal_html = self._generate_crystal_analysis_table(crystal_analysis)

        # 材料属性分析
        material_html = self._generate_material_properties_table(material_properties)

        # 化学键分析
        bonding_html = self._generate_chemical_bonding_analysis(chemical_bonding)

        return f"""
        <div class="section">
            <h2>🧮 Electron Localization Function (ELF) Analysis</h2>

            <div style="background: #fff3cd; padding: 15px; border-radius: 5px; margin-bottom: 20px; border-left: 4px solid #ffc107;">
                <p style="margin: 0; font-size: 0.9em;">
                    <strong>🔬 For detailed ELF analysis and advanced visualization, please use professional tools:</strong><br>
                    • <strong>VESTA</strong> - 3D visualization and isosurface analysis<br>
                    • <strong>VMD</strong> - Advanced molecular visualization and scripting<br>
                    • <strong>ChemCraft</strong> - Chemical data analysis and orbital visualization<br>
                    • <strong>pymatgen</strong> - Python materials analysis library<br>
                    • <strong>ASE</strong> - Atomic Simulation Environment
                </p>
            </div>

            <div class="grid-2">
                <div>
                    <h3>ELF Statistics</h3>
                    {stats_html}
                </div>
                <div>
                    <h3>ELF Theory Overview</h3>
                    <div style="background: #f7fafc; padding: 15px; border-radius: 5px; font-size: 0.9em;">
                        <p><strong>Electron Localization Function (ELF) Features:</strong></p>
                        <ul style="margin: 10px 0; padding-left: 20px;">
                            <li><strong>ELF ≈ 1:</strong> High localization (covalent bonds, lone pairs)</li>
                            <li><strong>ELF ≈ 0.5:</strong> Similar to uniform electron gas</li>
                            <li><strong>ELF ≈ 0:</strong> Low electron density regions</li>
                        </ul>
                        <p style="color: #666; font-size: 0.85em;">
                            💡 Tip: ELF helps identify chemical bond types and electronic structure features
                        </p>
                    </div>
                </div>
            </div>

            {visualization_html}
            {crystal_html}
            {material_html}
            {bonding_html}
        </div>
        """

    def _generate_elf_statistics_table(self, statistics, grid_dims) -> str:
        """生成ELF统计信息表格"""
        if 'error' in statistics:
            return f"<p>❌ {statistics['error']}</p>"

        total_points = np.prod(grid_dims) if grid_dims else 0

        return f"""
        <table class="data-table">
            <tr>
                <td>Grid Dimension</td>
                <td>{' × '.join(map(str, grid_dims)) if grid_dims else 'N/A'}</td>
            </tr>
            <tr>
                <td>Total Grid Points</td>
                <td>{total_points:,}</td>
            </tr>
            <tr>
                <td>ELF Minimum</td>
                <td>{statistics.get('min', 0):.6f}</td>
            </tr>
            <tr>
                <td>ELF Maximum</td>
                <td>{statistics.get('max', 0):.6f}</td>
            </tr>
            <tr>
                <td>ELF Mean</td>
                <td>{statistics.get('mean', 0):.6f}</td>
            </tr>
            <tr>
                <td>ELF Std Dev</td>
                <td>{statistics.get('std', 0):.6f}</td>
            </tr>
            <tr>
                <td>ELF Median</td>
                <td>{statistics.get('median', 0):.6f}</td>
            </tr>
        </table>
        """

    def _generate_elf_visualization_html(self, visualization) -> str:
        """生成ELF 2D可视化HTML内容"""
        if not visualization or not visualization.get('available', False):
            error_msg = visualization.get('error_message', 'Visualization not available')
            return f"""
            <div style="margin-top: 20px;">
                <h3>📊 ELF 2D Slice Visualization</h3>
                <p>⚠️ {error_msg}</p>
                <p>Note: 2D visualization requires matplotlib library</p>
            </div>
            """

        slices = visualization.get('slices', {})
        if not slices:
            return f"""
            <div style="margin-top: 20px;">
                <h3>📊 ELF 2D Slice Visualization</h3>
                <p>⚠️ No slice images generated</p>
            </div>
            """

        # 生成图像HTML
        images_html = ""
        for plane, slice_info in slices.items():
            title = slice_info.get('title', f'{plane.upper()} Plane')
            image_base64 = slice_info.get('image_base64', '')
            slice_index = slice_info.get('slice_index', 0)
            total_slices = slice_info.get('total_slices', 1)

            if image_base64:
                images_html += f"""
                <div style="text-align: center; margin-bottom: 20px;">
                    <h4>{title}</h4>
                    <img src="data:image/png;base64,{image_base64}"
                         style="max-width: 100%; height: auto; border: 1px solid #ddd; border-radius: 5px;"
                         alt="{title}">
                    <p style="font-size: 0.9em; color: #666; margin-top: 5px;">
                        Slice {slice_index + 1} of {total_slices} | Middle slice of the {plane.upper()} plane
                    </p>
                </div>
                """

        return f"""
        <div style="margin-top: 30px;">
            <h3>📊 ELF 2D Slice Visualization</h3>
            <p style="color: #666; margin-bottom: 20px;">
                Representative 2D slices showing ELF distribution in different crystallographic planes.
                Colors represent ELF values from low (dark) to high (bright).
            </p>
            <div class="grid-1" style="max-width: 800px; margin: 0 auto;">
                {images_html}
            </div>
            <div style="background: #e3f2fd; padding: 15px; border-radius: 5px; margin-top: 20px; border-left: 4px solid #2196f3;">
                <p style="margin: 0; font-size: 0.9em;">
                    <strong>📌 Note:</strong> These are representative middle slices. For complete 3D analysis,
                    interactive visualization, and custom slice selection, please use specialized software
                    such as <strong>VESTA</strong>, <strong>VMD</strong>, or <strong>pymatgen</strong>.
                </p>
            </div>
        </div>
        """

    def _generate_crystal_analysis_table(self, crystal_analysis) -> str:
        """生成晶体学分析表格"""
        if not crystal_analysis or 'error' in crystal_analysis:
            return ""

        if 'basic_info' in crystal_analysis:
            return f"""
            <div style="margin-top: 20px;">
                <h3>Crystallographic Analysis</h3>
                <p>⚠️ {crystal_analysis['basic_info']}</p>
            </div>
            """

        # 详细晶体学分析
        symmetry = crystal_analysis.get('symmetry', {})
        lattice_params = crystal_analysis.get('lattice_params', {})
        coordination = crystal_analysis.get('coordination', {})

        symmetry_html = ""
        if 'spacegroup_symbol' in symmetry:
            symmetry_html = f"""
            <h4>Symmetry Information</h4>
            <table class="data-table">
                <tr><td>Space Group</td><td>{symmetry['spacegroup_symbol']} (No. {symmetry.get('spacegroup_number', 'N/A')})</td></tr>
                <tr><td>Crystal System</td><td>{symmetry.get('crystal_system', 'N/A')}</td></tr>
                <tr><td>Lattice Type</td><td>{symmetry.get('lattice_type', 'N/A')}</td></tr>
                <tr><td>Point Group</td><td>{symmetry.get('point_group', 'N/A')}</td></tr>
            </table>
            """

        lattice_html = ""
        if lattice_params:
            lattice_html = f"""
            <h4>Lattice Parameters</h4>
            <table class="data-table">
                <tr><td>a</td><td>{lattice_params.get('a', 0):.4f} Å</td></tr>
                <tr><td>b</td><td>{lattice_params.get('b', 0):.4f} Å</td></tr>
                <tr><td>c</td><td>{lattice_params.get('c', 0):.4f} Å</td></tr>
                <tr><td>α</td><td>{lattice_params.get('alpha', 0):.2f}°</td></tr>
                <tr><td>β</td><td>{lattice_params.get('beta', 0):.2f}°</td></tr>
                <tr><td>γ</td><td>{lattice_params.get('gamma', 0):.2f}°</td></tr>
                <tr><td>Volume</td><td>{lattice_params.get('volume', 0):.4f} Å³</td></tr>
            </table>
            """

        coordination_html = ""
        if coordination:
            coordination_html = f"""
            <h4>Coordination Environment</h4>
            <table class="data-table">
                <tr><td>Average Coordination Number</td><td>{coordination.get('average_cn', 0):.2f}</td></tr>
                <tr><td>CN Range</td><td>{coordination.get('cn_range', [0, 0])[0]:.1f} - {coordination.get('cn_range', [0, 0])[1]:.1f}</td></tr>
            </table>
            """

        basic_info_html = ""
        if 'formula' in crystal_analysis:
            basic_info_html = f"""
            <h4>Basic Information</h4>
            <table class="data-table">
                <tr><td>Chemical Formula</td><td>{crystal_analysis.get('formula', 'N/A')}</td></tr>
                <tr><td>Reduced Formula</td><td>{crystal_analysis.get('reduced_formula', 'N/A')}</td></tr>
                <tr><td>Density</td><td>{crystal_analysis.get('density', 0):.4f} g/cm³</td></tr>
                <tr><td>Volume</td><td>{crystal_analysis.get('volume', 0):.4f} Å³</td></tr>
                <tr><td>Number of Sites</td><td>{crystal_analysis.get('num_sites', 0)}</td></tr>
            </table>
            """

        if basic_info_html or symmetry_html or lattice_html or coordination_html:
            return f"""
            <div style="margin-top: 20px;">
                <h3>Pymatgen Crystallographic Analysis</h3>
                <div class="grid-2">
                    <div>
                        {basic_info_html}
                        {symmetry_html}
                    </div>
                    <div>
                        {lattice_html}
                        {coordination_html}
                    </div>
                </div>
            </div>
            """

        return ""

    def _generate_material_properties_table(self, material_properties) -> str:
        """生成材料属性分析表格"""
        if not material_properties or 'error' in material_properties:
            return ""

        if 'basic_info' in material_properties:
            return f"""
            <div style="margin-top: 20px;">
                <h3>Material Properties Analysis</h3>
                <p>⚠️ {material_properties['basic_info']}</p>
            </div>
            """

        # 基础物理性质
        basic_html = ""
        if 'molecular_weight' in material_properties:
            basic_html = f"""
            <h4>Basic Physical Properties</h4>
            <table class="data-table">
                <tr><td>Molecular Weight</td><td>{material_properties.get('molecular_weight', 0):.4f} g/mol</td></tr>
                <tr><td>Total Electrons</td><td>{material_properties.get('electrons_per_formula', 0):.0f}</td></tr>
            </table>
            """

        # 磁性分析
        magnetic_html = ""
        if 'magnetic_elements' in material_properties:
            magnetic_html = f"""
            <h4>Magnetic Analysis</h4>
            <table class="data-table">
                <tr><td>Magnetic Elements</td><td>{', '.join(material_properties['magnetic_elements'])}</td></tr>
            </table>
            """
        elif 'magnetic_analysis' in material_properties:
            magnetic_html = f"""
            <h4>Magnetic Analysis</h4>
            <table class="data-table">
                <tr><td>Magnetic Character</td><td>{material_properties['magnetic_analysis']}</td></tr>
            </table>
            """

        # 电负性分析
        electronegativity_html = ""
        if 'electronegativity_analysis' in material_properties:
            en_analysis = material_properties['electronegativity_analysis']
            electronegativity_html = f"""
            <h4>Electronegativity Analysis</h4>
            <table class="data-table">
                <tr><td>Average Electronegativity</td><td>{en_analysis.get('average', 0):.3f}</td></tr>
                <tr><td>Electronegativity Range</td><td>{en_analysis.get('range', 0):.3f}</td></tr>
                <tr><td>Ionic Character</td><td>{en_analysis.get('ionic_character', 0):.3f}</td></tr>
            </table>
            """

        if basic_html or magnetic_html or electronegativity_html:
            return f"""
            <div style="margin-top: 20px;">
                <h3>Material Properties Analysis</h3>
                <div class="grid-2">
                    <div>
                        {basic_html}
                        {magnetic_html}
                    </div>
                    <div>
                        {electronegativity_html}
                    </div>
                </div>
            </div>
            """

        return ""

    def _generate_chemical_bonding_analysis(self, chemical_bonding) -> str:
        """生成化学键分析"""
        if not chemical_bonding or 'error' in chemical_bonding:
            return ""

        elf_distribution = chemical_bonding.get('elf_distribution', {})
        predicted_bonding = chemical_bonding.get('predicted_bonding', {})

        if not elf_distribution and not predicted_bonding:
            return ""

        # ELF分布分析
        distribution_html = ""
        if elf_distribution:
            high_loc = elf_distribution.get('high_localization_fraction', 0) * 100
            metallic = elf_distribution.get('metallic_character_fraction', 0) * 100
            covalent = elf_distribution.get('covalent_character_fraction', 0) * 100
            avg_elf = elf_distribution.get('average_elf', 0)

            distribution_html = f"""
            <div>
                <h4>ELF Distribution Features</h4>
                <table class="data-table">
                    <tr><td>High Localization Regions (ELF > 0.8)</td><td>{high_loc:.1f}%</td></tr>
                    <tr><td>Metallic Character (ELF < 0.2)</td><td>{metallic:.1f}%</td></tr>
                    <tr><td>Covalent Character (0.4 ≤ ELF ≤ 0.8)</td><td>{covalent:.1f}%</td></tr>
                    <tr><td>Average ELF Value</td><td>{avg_elf:.4f}</td></tr>
                </table>
            </div>
            """

        # 键合类型预测
        bonding_html = ""
        if predicted_bonding:
            bond_type = predicted_bonding.get('type', 'Unknown')
            description = predicted_bonding.get('description', '')
            confidence = predicted_bonding.get('confidence', 0) * 100

            # 根据键合类型设置颜色
            if '共价键' in bond_type:
                color_class = "style='color: #e53e3e; font-weight: bold;'"
            elif '金属键' in bond_type:
                color_class = "style='color: #3182ce; font-weight: bold;'"
            elif '离子键' in bond_type:
                color_class = "style='color: #38a169; font-weight: bold;'"
            else:
                color_class = "style='font-weight: bold;'"

            bonding_html = f"""
            <div>
                <h4>Bonding Type Prediction</h4>
                <table class="data-table">
                    <tr><td>Predicted Bonding Type</td><td {color_class}>{bond_type}</td></tr>
                    <tr><td>Confidence</td><td>{confidence:.1f}%</td></tr>
                </table>
                <div style="background: #f7fafc; padding: 10px; border-radius: 5px; margin-top: 10px;">
                    <p style="font-size: 0.9em; color: #4a5568;"><strong>Analysis Description:</strong> {description}</p>
                </div>
            </div>
            """

        return f"""
        <div style="margin-top: 20px;">
            <h3>🔗 Chemical Bonding Analysis</h3>
            <div class="grid-2">
                {distribution_html}
                {bonding_html}
            </div>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  最终结果
    # ------------------------------------------------------------------ #
    def _generate_final_results_section(self) -> str:
        """生成最终结果部分"""
        final = self.data.get('final_results', {})

        return f"""
        <div class="section">
            <h2>🏁 最终结果汇总</h2>

            <div class="grid-2">
                <div>
                    <h3>能量与收敛</h3>
                    <table class="data-table">
                        <tr>
                            <td>最终总能量</td>
                            <td>{(final.get('final_energy') or 0):.6f} eV</td>
                        </tr>
                        <tr>
                            <td>收敛状态</td>
                            <td>{'✅ 收敛' if final.get('converged', False) else '❌ 未收敛'}</td>
                        </tr>
                        <tr>
                            <td>总电子步数</td>
                            <td>{final.get('total_electronic_steps', 0)}</td>
                        </tr>
                    </table>
                </div>

                <div>
                    <h3>电子结构</h3>
                    <table class="data-table">
                        <tr>
                            <td>费米能级</td>
                            <td>{(final.get('fermi_energy') or 0):.4f} eV</td>
                        </tr>
                        <tr>
                            <td>带隙</td>
                            <td>{(final.get('band_gap') or 0):.4f} eV</td>
                        </tr>
                        <tr>
                            <td>总磁矩</td>
                            <td>{(final.get('total_magnetization') or 0):.4f} μ<sub>B</sub></td>
                        </tr>
                    </table>
                </div>
            </div>
        </div>
        """

    # ------------------------------------------------------------------ #
    #  JavaScript
    # ------------------------------------------------------------------ #
    def _generate_javascript(self) -> str:
        """生成JavaScript代码"""
        # 准备图表数据
        electronic_steps = self.data['electronic_convergence'].get('electronic_steps', [])
        scf_energies = self.data['electronic_convergence'].get('scf_energies', [])
        forces = self.data['forces_and_stress'].get('forces', [])
        atom_mags = self.data['magnetic_properties'].get('atom_magnetizations', [])

        # 电子步数据
        electronic_step_nums = [step['step'] for step in electronic_steps]
        electronic_energies = [step['free_energy'] for step in electronic_steps]
        energy_changes = [step['energy_change'] for step in electronic_steps]

        # 力数据
        force_magnitudes = [f['magnitude'] for f in forces]
        atom_indices = list(range(1, len(forces) + 1))

        # 磁矩数据
        mag_atom_indices = [mag['atom_index'] + 1 for mag in atom_mags]
        mag_magnitudes = [mag['total_magnitude'] for mag in atom_mags]

        # Bader电荷数据
        bader_data = self.data.get('bader_analysis', {})
        atom_charges = bader_data.get('atom_charges', [])
        bader_atom_indices = [atom['atom_index'] for atom in atom_charges]
        bader_charges = [atom['bader_charge'] for atom in atom_charges]
        charge_transfers = [atom['charge_transfer'] for atom in atom_charges]

        return f"""
        // 电子步能量变化图表
        const electronicCtx = document.getElementById('electronicStepsChart').getContext('2d');
        new Chart(electronicCtx, {{
            type: 'bar',
            data: {{
                labels: {electronic_step_nums},
                datasets: [{{
                    label: '能量变化 (eV)',
                    data: {energy_changes},
                    backgroundColor: function(context) {{
                        const value = context.parsed.y;
                        return value >= 0 ? 'rgba(220, 53, 69, 0.6)' : 'rgba(40, 167, 69, 0.6)';
                    }},
                    borderColor: function(context) {{
                        const value = context.parsed.y;
                        return value >= 0 ? 'rgb(220, 53, 69)' : 'rgb(40, 167, 69)';
                    }},
                    borderWidth: 1
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        type: 'linear',
                        display: true,
                        title: {{
                            display: true,
                            text: '能量变化 (eV)'
                        }},
                        beginAtZero: true,
                        grid: {{
                            color: function(context) {{
                                if (context.tick.value === 0) {{
                                    return 'rgba(0, 0, 0, 0.8)';
                                }}
                                return 'rgba(0, 0, 0, 0.1)';
                            }},
                            lineWidth: function(context) {{
                                if (context.tick.value === 0) {{
                                    return 2;
                                }}
                                return 1;
                            }}
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '电子步'
                        }}
                    }}
                }},
                plugins: {{
                    legend: {{
                        display: true
                    }},
                    tooltip: {{
                        callbacks: {{
                            afterLabel: function(context) {{
                                const value = context.parsed.y;
                                return value >= 0 ? '能量增加' : '能量降低';
                            }}
                        }}
                    }}
                }}
            }}
        }});

        // SCF能量收敛图表
        const scfCtx = document.getElementById('scfEnergyChart').getContext('2d');
        new Chart(scfCtx, {{
            type: 'line',
            data: {{
                labels: {list(range(1, len(scf_energies) + 1))},
                datasets: [{{
                    label: '总能量 (eV)',
                    data: {scf_energies},
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
                            text: 'SCF步'
                        }}
                    }}
                }}
            }}
        }});

        // 原子受力图表
        const forcesCtx = document.getElementById('forcesChart').getContext('2d');
        new Chart(forcesCtx, {{
            type: 'bar',
            data: {{
                labels: {atom_indices},
                datasets: [{{
                    label: '力大小 (eV/Å)',
                    data: {force_magnitudes},
                    backgroundColor: 'rgba(255, 159, 64, 0.6)',
                    borderColor: 'rgb(255, 159, 64)',
                    borderWidth: 1
                }}]
            }},
            options: {{
                responsive: true,
                maintainAspectRatio: false,
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{
                            display: true,
                            text: '力大小 (eV/Å)'
                        }}
                    }},
                    x: {{
                        title: {{
                            display: true,
                            text: '原子索引'
                        }}
                    }}
                }}
            }}
        }});

        // 磁矩分布图表
        if ({len(atom_mags)} > 0) {{
            const magCtx = document.getElementById('magnetizationChart').getContext('2d');
            new Chart(magCtx, {{
                type: 'bar',
                data: {{
                    labels: {mag_atom_indices},
                    datasets: [{{
                        label: '磁矩大小 (μB)',
                        data: {mag_magnitudes},
                        backgroundColor: 'rgba(75, 192, 192, 0.6)',
                        borderColor: 'rgb(75, 192, 192)',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        y: {{
                            title: {{
                                display: true,
                                text: '磁矩大小 (μB)'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: '原子索引'
                            }}
                        }}
                    }}
                }}
            }});
        }}

        // Bader电荷分布图表
        if ({len(atom_charges)} > 0) {{
            const baderCtx = document.getElementById('baderChargeChart').getContext('2d');
            new Chart(baderCtx, {{
                type: 'scatter',
                data: {{
                    labels: {bader_atom_indices},
                    datasets: [{{
                        label: 'Bader电荷 (e)',
                        data: {[{'x': idx, 'y': charge} for idx, charge in zip(bader_atom_indices, bader_charges)]},
                        backgroundColor: 'rgba(255, 193, 7, 0.6)',
                        borderColor: 'rgb(255, 193, 7)',
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        y: {{
                            title: {{
                                display: true,
                                text: 'Bader电荷 (e)'
                            }}
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: '原子索引'
                            }}
                        }}
                    }}
                }}
            }});

            // 电荷转移分布图表
            const transferCtx = document.getElementById('chargeTransferChart').getContext('2d');
            new Chart(transferCtx, {{
                type: 'bar',
                data: {{
                    labels: {bader_atom_indices},
                    datasets: [{{
                        label: '电荷转移 (e)',
                        data: {charge_transfers},
                        backgroundColor: function(context) {{
                            const value = context.parsed.y;
                            return value > 0 ? 'rgba(220, 53, 69, 0.6)' : 'rgba(40, 167, 69, 0.6)';
                        }},
                        borderColor: function(context) {{
                            const value = context.parsed.y;
                            return value > 0 ? 'rgb(220, 53, 69)' : 'rgb(40, 167, 69)';
                        }},
                        borderWidth: 1
                    }}]
                }},
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    scales: {{
                        y: {{
                            title: {{
                                display: true,
                                text: '电荷转移 (e)'
                            }},
                            beginAtZero: true
                        }},
                        x: {{
                            title: {{
                                display: true,
                                text: '原子索引'
                            }}
                        }}
                    }},
                    plugins: {{
                        legend: {{
                            display: true
                        }},
                        tooltip: {{
                            callbacks: {{
                                afterLabel: function(context) {{
                                    const value = context.parsed.y;
                                    return value > 0 ? '失去电子' : '获得电子';
                                }}
                            }}
                        }}
                    }}
                }}
            }});
        }}
        """


def generate_scf_report(input_path: str, task_id: Optional[str] = None, output_dir: Optional[str] = None) -> str:
    """
    生成自洽场(SCF)计算分析报告

    Args:
        input_path: OUTCAR文件路径或包含VASP文件的文件夹路径
        task_id: 任务ID（可选）
        output_dir: 输出目录，默认与输入路径同目录

    Returns:
        生成的HTML报告文件路径
    """
    try:
        # 分析VASP计算结果
        analyzer = SCFAnalyzer(input_path, task_id)
        analysis_data = analyzer.analyze()

        # 确定输出路径
        if output_dir is None:
            output_dir_path = analyzer.work_dir
        else:
            output_dir_path = Path(output_dir)

        # 生成HTML报告
        output_file = output_dir_path / "scf_analysis_report.html"

        generator = SCFHTMLGenerator(analysis_data)
        html_path = generator.generate_html_report(str(output_file))

        return html_path

    except Exception as e:
        raise Exception(f"生成SCF分析报告失败: {str(e)}")
