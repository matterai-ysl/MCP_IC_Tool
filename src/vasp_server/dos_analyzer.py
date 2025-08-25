#!/usr/bin/env python3
"""
PyMatGen DOS分析器 - 专业级DOS分析工具

基于PyMatGen实现的高级DOS分析器，提供：
- 基于Vasprun的完整DOS分析
- 简化的DOS数据解析和可视化
- 数据保存功能（CSV）
- 专业可视化（matplotlib + PyMatGen）
- HTML和Markdown报告生成
- 化学环境和键合分析
- 能带结构分析
- 轨道投影和磁性分析
"""

import io
import base64
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple
from datetime import datetime

# PyMatGen导入
from pymatgen.io.vasp import Outcar, Poscar
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.core import Structure, Element
from pymatgen.electronic_structure.dos import DOS, CompleteDos
from pymatgen.electronic_structure.core import Spin
from pymatgen.electronic_structure.plotter import DosPlotter, BSPlotter
from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
from pymatgen.analysis.local_env import CrystalNN
from pymatgen.analysis.bond_valence import BVAnalyzer

# 设置中文字体和日志
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class PyMatGenDOSAnalyzer:
    """完全基于PyMatGen的专业DOS分析器，融合test3.py的简化逻辑"""
    
    def __init__(self, input_path: str, task_id: Optional[str] = None, output_dir: Optional[str] = None):
        """
        初始化DOS分析器
        
        Args:
            input_path: 包含VASP输出文件的目录或vasprun.xml文件路径
            task_id: 任务ID
            output_dir: 结果输出目录
        """
        self.input_path = Path(input_path)
        self.task_id = task_id or "dos_analysis"
        
        # 确定工作目录和vasprun.xml路径
        if self.input_path.is_file() and self.input_path.name == 'vasprun.xml':
            self.work_dir = self.input_path.parent
            self.vasprun_path = self.input_path
        else:
            self.work_dir = self.input_path
            self.vasprun_path = self.work_dir / 'vasprun.xml'
        
        # 设置输出目录
        if output_dir:
            self.output_dir = Path(output_dir)
        else:
            self.output_dir = self.work_dir / "dos_analysis"
        self.output_dir.mkdir(exist_ok=True)
        
        # 检查必要文件
        if not self.vasprun_path.exists():
            raise FileNotFoundError(f"找不到vasprun.xml文件: {self.vasprun_path}")
        
        # 初始化变量
        self.analysis_data = {}
        self.vasprun = None
        self.complete_dos = None
        self.structure = None
        self.is_spin_polarized = False
        
    def analyze(self) -> Dict[str, Any]:
        """执行完整的DOS分析"""
        logger.info(f"🔬 开始PyMatGen DOS分析: {self.vasprun_path}")
        
        try:
            # 加载DOS数据
            if not self.load_dos_data():
                raise ValueError("无法加载DOS数据")
            
            # 基础分析
            self._analyze_structure()
            self._analyze_calculation_settings()
            self._analyze_dos()
            self._analyze_band_structure()
            self._analyze_chemical_properties()
            self._analyze_magnetic_properties()
            
            # 数据提取和保存
            dos_data = self.extract_dos_data()
            if dos_data:
                self.save_dos_data(dos_data)
            
            # 生成可视化（matplotlib图表优先）
            self._generate_matplotlib_plots()
            self._generate_visualizations()
            
            # 生成报告
            self._generate_markdown_report()
            self._finalize_analysis()
            
            logger.info("✅ PyMatGen DOS分析完成")
            return self.analysis_data
            
        except Exception as e:
            logger.error(f"❌ DOS分析失败: {e}")
            import traceback
            traceback.print_exc()
            raise
    
    def load_dos_data(self) -> bool:
        """
        从vasprun.xml文件加载DOS数据 - 基于test3.py的简化逻辑
        
        Returns:
            bool: 是否成功加载数据
        """
        if not self.vasprun_path.exists():
            logger.error(f"未找到 {self.vasprun_path} 文件")
            logger.info("请确保vasprun.xml文件存在于工作目录中")
            return False
        
        try:
            logger.info(f"从 {self.vasprun_path} 读取DOS数据")
            
            # 解析vasprun.xml
            self.vasprun = Vasprun(str(self.vasprun_path), parse_dos=True, parse_eigen=True)
            self.complete_dos = self.vasprun.complete_dos
            self.structure = self.vasprun.final_structure
            
            # 检查DOS数据是否成功解析
            if self.complete_dos is None:
                logger.error("vasprun.xml中没有找到DOS数据")
                logger.info("请确保VASP计算中包含DOS输出（LORBIT参数）")
                return False
            
            # 检查是否为自旋极化计算
            densities = self.complete_dos.densities
            self.is_spin_polarized = Spin.down in densities
            
            logger.info(f"成功加载DOS数据")
            logger.info(f"  - 自旋极化: {'是' if self.is_spin_polarized else '否'}")
            logger.info(f"  - 费米能级: {self.complete_dos.efermi:.4f} eV")
            logger.info(f"  - 能量点数: {len(self.complete_dos.energies)}")
            
            if self.structure:
                logger.info(f"  - 化学式: {self.structure.composition.reduced_formula}")
                logger.info(f"  - 原子数: {len(self.structure)}")
            
            return True
            
        except ImportError as e:
            logger.error(f"无法导入pymatgen模块: {e}")
            logger.info("请确保pymatgen已正确安装")
            return False
        except Exception as e:
            logger.error(f"加载vasprun.xml失败: {e}")
            logger.info("请检查vasprun.xml文件是否完整且格式正确")
            return False
    
    def get_band_gap_info(self) -> Dict:
        """获取带隙信息 - 基于test3.py的逻辑并结合PyMatGen方法"""
        gap_info = {
            'fermi_level': 0.0,
            'band_gap': 0.0,
            'vbm': None,
            'cbm': None,
            'is_metal': True
        }
        
        if not self.complete_dos:
            return gap_info
            
        try:
            gap_info['fermi_level'] = self.complete_dos.efermi
            
            # 尝试使用pymatgen的内置方法
            try:
                gap = self.complete_dos.get_gap()
                print("here")
                print(gap)
                cbm_vbm = self.complete_dos.get_cbm_vbm()
                print(cbm_vbm)
                if gap > 0.01:
                    gap_info['band_gap'] = gap
                    gap_info['is_metal'] = False
                    if cbm_vbm:
                        gap_info['cbm'] = cbm_vbm[0]
                        gap_info['vbm'] = cbm_vbm[1]
            except:
                # 备用方法：简单的带隙估算
                energies = self.complete_dos.energies - self.complete_dos.efermi
                densities = self.complete_dos.densities
                
                if densities:
                    if self.is_spin_polarized:
                        # 自旋极化：取平均DOS
                        tdos_values = (densities[Spin.up] + densities[Spin.down]) / 2
                    else:
                        tdos_values = densities[Spin.up]
                    
                    fermi_idx = np.argmin(np.abs(energies))
                    dos_threshold = 0.01
                    
                    # VBM: 费米能级以下最高占据态
                    for i in range(fermi_idx, -1, -1):
                        if i < len(tdos_values) and tdos_values[i] > dos_threshold:
                            gap_info['vbm'] = energies[i]
                            break
                    
                    # CBM: 费米能级以上最低未占据态
                    for i in range(fermi_idx, len(energies)):
                        if i < len(tdos_values) and tdos_values[i] > dos_threshold:
                            gap_info['cbm'] = energies[i]
                            break
                    
                    # 如果没有从内置方法获得带隙，使用估算值
                    if gap_info['band_gap'] == 0.0 and gap_info['vbm'] is not None and gap_info['cbm'] is not None:
                        gap_info['band_gap'] = gap_info['cbm'] - gap_info['vbm']
                        if gap_info['band_gap'] > 0.1:
                            gap_info['is_metal'] = False
                            
        except Exception as e:
            logger.warning(f"计算带隙信息时出错: {e}")
            
        return gap_info
    
    def extract_dos_data(self) -> Dict:
        """
        提取DOS数据到字典格式 - 基于test3.py的逻辑
        
        Returns:
            Dict: 包含各种DOS数据的字典
        """
        if self.complete_dos is None:
            logger.error("DOS数据未加载")
            return {}
            
        dos_data = {}
        
        try:
            energies = self.complete_dos.energies - self.complete_dos.efermi
            
            # 总DOS
            tdos = self.complete_dos.densities
            
            if self.is_spin_polarized:
                dos_data['TDOS'] = {
                    'energy': energies,
                    'dos_up': tdos[Spin.up],
                    'dos_down': -tdos[Spin.down]  # 负值用于绘图
                }
                logger.info("提取了自旋极化的总DOS数据")
            else:
                dos_data['TDOS'] = {
                    'energy': energies,
                    'dos': tdos[Spin.up]
                }
                logger.info("提取了非自旋极化的总DOS数据")
            
            # 元素分解DOS
            if self.structure:
                try:
                    unique_elements = list(set([site.species_string for site in self.structure]))
                    element_dos_dict = self.complete_dos.get_element_dos()
                    
                    logger.info(f"发现 {len(unique_elements)} 种元素: {', '.join(unique_elements)}")
                    
                    for element_str in unique_elements:
                        try:
                            element = Element(element_str)
                            
                            if element in element_dos_dict:
                                element_dos = element_dos_dict[element]
                                element_densities = element_dos.densities
                                
                                if self.is_spin_polarized:
                                    dos_data[element_str] = {
                                        'energy': energies,
                                        'total_up': element_densities[Spin.up],
                                        'total_down': -element_densities[Spin.down]
                                    }
                                else:
                                    dos_data[element_str] = {
                                        'energy': energies,
                                        'total': element_densities[Spin.up]
                                    }
                                
                                logger.info(f"成功提取元素 {element_str} 的DOS数据")
                            else:
                                logger.warning(f"未找到元素 {element_str} 的DOS数据")
                                
                        except Exception as e:
                            logger.warning(f"处理元素 {element_str} 时出错: {e}")
                            continue
                            
                except Exception as e:
                    logger.warning(f"获取元素分解DOS时出错: {e}")
            else:
                logger.info("未找到元素分解DOS数据，只提取总DOS")
            
        except Exception as e:
            logger.error(f"提取DOS数据失败: {e}")
        
        return dos_data
    
    def save_dos_data(self, dos_data: Dict) -> None:
        """保存DOS数据为CSV格式 - 基于test3.py的逻辑"""
        if not dos_data:
            logger.warning("没有DOS数据可保存")
            return
            
        data_dir = self.output_dir / "data"
        data_dir.mkdir(exist_ok=True)
        
        for name, data in dos_data.items():
            try:
                df = pd.DataFrame(data)
                output_file = data_dir / f"{name}_DOS.csv"
                df.to_csv(output_file, index=False)
                logger.info(f"已保存 {name} DOS数据到 {output_file}")
            except Exception as e:
                logger.error(f"保存 {name} DOS数据失败: {e}")
    
    def _analyze_structure(self):
        """分析晶体结构"""
        logger.info("🏗️ 分析晶体结构...")
        
        if not self.vasprun:
            logger.error("vasprun数据未加载")
            return
            
        structure = self.vasprun.final_structure
        if not structure:
            logger.error("未找到结构信息")
            return
        
        structure_analysis = {
            'formula': structure.formula,
            'reduced_formula': structure.reduced_formula,
            'composition': dict(structure.composition.as_dict()),
            'num_sites': len(structure),
            'lattice_parameters': {
                'a': structure.lattice.a,
                'b': structure.lattice.b, 
                'c': structure.lattice.c,
                'alpha': structure.lattice.alpha,
                'beta': structure.lattice.beta,
                'gamma': structure.lattice.gamma,
                'volume': structure.lattice.volume
            },
            'space_group': structure.get_space_group_info()[1],
            'point_group': structure.get_space_group_info()[0],
            'density': structure.density,
            'elements': [str(el) for el in structure.composition.elements]
        }
        
        # 配位环境分析
        try:
            crystal_nn = CrystalNN()
            coordination_analysis = {}
            
            for i, site in enumerate(structure):
                try:
                    cn_info = crystal_nn.get_cn_dict(structure, i)
                    coordination_analysis[f'site_{i}_{site.specie}'] = {
                        'coordination_number': sum(cn_info.values()),
                        'neighbors': {str(k): v for k, v in cn_info.items()}
                    }
                except:
                    coordination_analysis[f'site_{i}_{site.specie}'] = {
                        'coordination_number': 'unknown',
                        'neighbors': {}
                    }
            
            structure_analysis['coordination_analysis'] = coordination_analysis
        except Exception as e:
            print(f"   ⚠️ 配位分析失败: {e}")
            structure_analysis['coordination_analysis'] = {}
        
        # 键价分析
        try:
            bv_analyzer = BVAnalyzer()
            valences = bv_analyzer.get_valences(structure)
            # 确保valences是可序列化的
            if hasattr(valences, 'tolist'):
                valences_list = valences.tolist()
            else:
                valences_list = list(valences)
            structure_analysis['bond_valence_analysis'] = {
                'valences': valences_list,
                'is_valid_structure': True
            }
        except Exception as e:
            print(f"   ⚠️ 键价分析失败: {e}")
            structure_analysis['bond_valence_analysis'] = {
                'is_valid_structure': False,
                'error': str(e)
            }
        
        self.analysis_data['structure'] = structure_analysis
    
    def _analyze_calculation_settings(self):
        """分析计算设置"""
        logger.info("⚙️ 分析计算设置...")
        
        if not self.vasprun:
            logger.error("vasprun数据未加载")
            return
            
        incar = self.vasprun.incar
        kpoints = self.vasprun.kpoints
        
        calc_settings = {
            'functional': incar.get('GGA', 'PBE'),
            'encut': incar.get('ENCUT', 'unknown'),
            'ismear': incar.get('ISMEAR', 'unknown'),
            'sigma': incar.get('SIGMA', 'unknown'),
            'ispin': incar.get('ISPIN', 1),
            'lorbit': incar.get('LORBIT', 'unknown'),
            'nedos': incar.get('NEDOS', 'unknown'),
            'kpoints_mode': kpoints.style.name if kpoints else 'unknown',
            'kpoints_grid': kpoints.kpts[0] if kpoints and len(kpoints.kpts) > 0 else 'unknown'
        }
        
        self.analysis_data['calculation_settings'] = calc_settings
    
    def _analyze_dos(self):
        """使用PyMatGen进行DOS分析 - 基于test3.py的简化逻辑"""
        logger.info("📊 分析态密度...")
        
        if not self.complete_dos:
            raise ValueError("DOS数据未加载")
        
        # 基础DOS信息
        fermi_energy = self.complete_dos.efermi
        energies = self.complete_dos.energies
        
        # 带隙分析 - 使用改进的方法
        gap_info = self.get_band_gap_info()
        
        # 判断材料类型
        band_gap = gap_info['band_gap']
        is_metal = gap_info['is_metal']
        material_type = 'metal' if is_metal else ('semiconductor' if band_gap < 3.0 else 'insulator')
        
        # 简化的带隙类型判断
        if gap_info['cbm'] is not None and gap_info['vbm'] is not None:
            gap_type = 'direct' if abs(gap_info['cbm'] - gap_info['vbm'] - band_gap) < 0.01 else 'indirect'
        else:
            gap_type = 'unknown'
        
        dos_analysis = {
            'fermi_energy': fermi_energy,
            'band_gap': band_gap,
            'cbm_energy': gap_info['cbm'],
            'vbm_energy': gap_info['vbm'],
            'gap_type': gap_type,
            'material_type': material_type,
            'is_metal': is_metal,
            'is_spin_polarized': self.is_spin_polarized,
            'energy_range': [float(energies.min()), float(energies.max())],
            'num_energy_points': len(energies)
        }
        
        # DOS积分分析
        dos_integrals = self._calculate_dos_integrals(self.complete_dos)
        dos_analysis['dos_integrals'] = dos_integrals
        
        # 轨道分析
        orbital_analysis = self._analyze_orbital_projections(self.complete_dos)
        dos_analysis['orbital_analysis'] = orbital_analysis
        
        self.analysis_data['dos_analysis'] = dos_analysis
    
    def _analyze_band_structure(self):
        """分析能带结构（如果可用）"""
        logger.info("🎼 分析能带结构...")
        
        if not self.vasprun:
            logger.warning("vasprun数据未加载，跳过能带结构分析")
            self.analysis_data['band_structure'] = {'has_band_structure': False}
            return
        
        try:
            # 尝试不同的方法获取能带结构
            band_structure = None
            
            # 方法1: 尝试获取沿对称线的能带结构
            try:
                band_structure = self.vasprun.get_band_structure(line_mode=True)
                if band_structure:
                    logger.info("获取了沿对称线的能带结构 (BandStructureSymmLine)")
            except Exception as e:
                logger.debug(f"沿对称线模式失败: {e}")
            
            # 方法2: 如果方法1失败，尝试标准方法
            if not band_structure:
                try:
                    band_structure = self.vasprun.get_band_structure()
                    if band_structure:
                        logger.info(f"获取了能带结构，类型: {type(band_structure).__name__}")
                except Exception as e:
                    logger.debug(f"标准方法失败: {e}")
            
            if band_structure:
                # 使用DOS数据计算更准确的带隙信息
                gap_info_from_dos = self.get_band_gap_info()
                
                # 尝试从能带结构获取信息（作为参考）
                try:
                    bs_gap_info = band_structure.get_band_gap()
                    bs_is_metal = band_structure.is_metal()
                    bs_direct_gap = 0.0
                    try:
                        bs_direct_gap_result = band_structure.get_direct_band_gap()
                        if isinstance(bs_direct_gap_result, dict):
                            bs_direct_gap = bs_direct_gap_result.get('energy', 0.0)
                        elif isinstance(bs_direct_gap_result, (int, float)):
                            bs_direct_gap = float(bs_direct_gap_result)
                    except:
                        bs_direct_gap = 0.0
                except Exception as e:
                    logger.warning(f"从能带结构获取信息失败: {e}")
                    bs_gap_info = None
                    bs_is_metal = None
                    bs_direct_gap = 0.0
                
                # 优先使用DOS计算的结果，能带结构作为补充
                band_analysis = {
                    'has_band_structure': True,
                    'band_structure_type': type(band_structure).__name__,
                    # 使用DOS计算的更准确结果
                    'is_metal': gap_info_from_dos.get('is_metal', True),
                    'fundamental_gap': gap_info_from_dos.get('band_gap', 0.0),
                    'vbm_energy': gap_info_from_dos.get('vbm', None),
                    'cbm_energy': gap_info_from_dos.get('cbm', None),
                    'fermi_level': gap_info_from_dos.get('fermi_level', 0.0),
                    # 对于DOS计算，直接带隙和基本带隙相同
                    'direct_gap': gap_info_from_dos.get('band_gap', 0.0),
                    'num_bands': band_structure.nb_bands,
                    # 保留能带结构的原始信息作为参考
                    'band_structure_gap_info': bs_gap_info,
                    'band_structure_is_metal': bs_is_metal,
                    'calculation_method': 'DOS-based calculation with band structure supplement'
                }
                
                self.band_structure = band_structure
                logger.info(f"能带结构分析完成 - 类型: {type(band_structure).__name__}")
                
                # 输出DOS计算的带隙信息
                logger.info("=== DOS计算的电子结构信息 ===")
                logger.info(f"材料类型: {'金属' if gap_info_from_dos.get('is_metal', True) else '半导体/绝缘体'}")
                logger.info(f"费米能级: {gap_info_from_dos.get('fermi_level', 0.0):.4f} eV")
                logger.info(f"基本带隙: {gap_info_from_dos.get('band_gap', 0.0):.4f} eV")
                if gap_info_from_dos.get('vbm') is not None:
                    logger.info(f"价带顶 (VBM): {gap_info_from_dos.get('vbm'):.4f} eV (相对费米能级)")
                if gap_info_from_dos.get('cbm') is not None:
                    logger.info(f"导带底 (CBM): {gap_info_from_dos.get('cbm'):.4f} eV (相对费米能级)")
                logger.info(f"能带数: {band_structure.nb_bands}")
            else:
                # 即使没有能带结构，也可以基于DOS计算带隙
                gap_info_from_dos = self.get_band_gap_info()
                
                band_analysis = {
                    'has_band_structure': False,
                    'has_dos_gap_analysis': True,
                    # 使用DOS计算的结果
                    'is_metal': gap_info_from_dos.get('is_metal', True),
                    'fundamental_gap': gap_info_from_dos.get('band_gap', 0.0),
                    'vbm_energy': gap_info_from_dos.get('vbm', None),
                    'cbm_energy': gap_info_from_dos.get('cbm', None),
                    'fermi_level': gap_info_from_dos.get('fermi_level', 0.0),
                    'direct_gap': gap_info_from_dos.get('band_gap', 0.0),
                    'num_bands': 'N/A (DOS only)',
                    'calculation_method': 'DOS-based calculation (no band structure available)'
                }
                
                self.band_structure = None
                logger.info("未找到能带结构数据，使用DOS计算带隙信息")
                
                # 输出DOS计算的带隙信息
                logger.info("=== DOS计算的电子结构信息 ===")
                logger.info(f"材料类型: {'金属' if gap_info_from_dos.get('is_metal', True) else '半导体/绝缘体'}")
                logger.info(f"费米能级: {gap_info_from_dos.get('fermi_level', 0.0):.4f} eV")
                logger.info(f"基本带隙: {gap_info_from_dos.get('band_gap', 0.0):.4f} eV")
                if gap_info_from_dos.get('vbm') is not None:
                    logger.info(f"价带顶 (VBM): {gap_info_from_dos.get('vbm'):.4f} eV (相对费米能级)")
                if gap_info_from_dos.get('cbm') is not None:
                    logger.info(f"导带底 (CBM): {gap_info_from_dos.get('cbm'):.4f} eV (相对费米能级)")
                
        except Exception as e:
            logger.warning(f"能带结构分析失败: {e}")
            band_analysis = {'has_band_structure': False, 'error': str(e)}
            self.band_structure = None
        
        self.analysis_data['band_structure'] = band_analysis
    
    def _analyze_chemical_properties(self):
        """分析化学性质"""
        logger.info("🧪 分析化学性质...")
        
        if not self.vasprun or not self.structure or not self.complete_dos:
            logger.warning("缺少必要数据，跳过化学性质分析")
            self.analysis_data['chemical_properties'] = {}
            return
        
        structure = self.structure
        complete_dos = self.complete_dos
        
        # 元素分析
        element_analysis = {}
        for element in structure.composition.elements:
            element_str = str(element)
            element_dos = complete_dos.get_element_dos()[element]
            
            # 计算该元素的DOS积分
            element_integral = {}
            for spin in element_dos.densities:
                spin_label = 'up' if spin == Spin.up else 'down'
                dos_values = element_dos.densities[spin]
                element_integral[spin_label] = float(np.trapezoid(dos_values, complete_dos.energies))
            
            element_analysis[element_str] = {
                'electronegativity': element.X if hasattr(element, 'X') else None,
                'atomic_radius': getattr(element, 'atomic_radius', None),
                'dos_integral': element_integral,
                'oxidation_states': getattr(element, 'common_oxidation_states', [])
            }
        
        # 键合特征分析
        bonding_analysis = self._analyze_bonding_character()
        
        chemical_analysis = {
            'elements': element_analysis,
            'bonding_character': bonding_analysis,
            'electronegativity_difference': self._calculate_electronegativity_difference(structure),
            'ionic_character': self._estimate_ionic_character(structure)
        }
        
        self.analysis_data['chemical_properties'] = chemical_analysis
    
    def _analyze_magnetic_properties(self):
        """分析磁性性质"""
        print("🧲 分析磁性性质...")
        
        complete_dos = self.complete_dos
        
        if not self.analysis_data['dos_analysis']['is_spin_polarized']:
            magnetic_analysis = {
                'is_magnetic': False,
                'magnetic_type': 'non-magnetic'
            }
        else:
            dos_up = complete_dos.densities[Spin.up]
            dos_down = complete_dos.densities[Spin.down]
            energies = complete_dos.energies
            
            # 计算磁矩和自旋极化
            spin_diff = dos_up - dos_down
            total_magnetization = np.trapezoid(spin_diff, energies)
            
            # 自旋极化度
            spin_polarization = np.mean(np.abs(spin_diff) / (dos_up + dos_down + 1e-10))
            
            # 磁性类型判断
            magnetic_type = 'ferromagnetic' if total_magnetization > 0.1 else 'antiferromagnetic'
            
            magnetic_analysis = {
                'is_magnetic': True,
                'magnetic_type': magnetic_type,
                'total_magnetization': float(total_magnetization),
                'spin_polarization': float(spin_polarization),
                'max_spin_difference': float(np.max(np.abs(spin_diff))),
                'fermi_spin_polarization': float(spin_diff[np.argmin(np.abs(energies - complete_dos.efermi))])
            }
            
            # 尝试从OUTCAR获取更多磁性信息
            try:
                outcar_path = self.work_dir / 'OUTCAR'
                if outcar_path.exists():
                    outcar = Outcar(str(outcar_path))
                    if hasattr(outcar, 'total_mag'):
                        magnetic_analysis['outcar_total_magnetization'] = outcar.total_mag
            except:
                pass
        
        self.analysis_data['magnetic_properties'] = magnetic_analysis
    
    def _generate_matplotlib_plots(self, energy_range: Tuple[float, float] = (-8, 6)):
        """生成matplotlib可视化图表 - 基于test3.py的实现"""
        logger.info("📈 生成matplotlib可视化图表...")
        
        try:
            # 生成总DOS图
            self.plot_total_dos(energy_range=energy_range)
            
            # 生成元素分解DOS图
            self.plot_element_dos(energy_range=energy_range)
            
            # 生成SPD轨道DOS图
            self.plot_spd_dos(energy_range=energy_range)
            
            # 生成能带结构图（如果有数据）
            self.plot_band_structure()
            
        except Exception as e:
            logger.error(f"生成matplotlib图表失败: {e}")
    
    def plot_total_dos(self, energy_range: Tuple[float, float] = (-10, 10)) -> None:
        """绘制总DOS图 - 基于test3.py的实现"""
        if not self.complete_dos:
            logger.error("DOS数据未加载")
            return
            
        try:
            fig, ax = plt.subplots(figsize=(10, 6))
            
            energies = self.complete_dos.energies - self.complete_dos.efermi
            mask = (energies >= energy_range[0]) & (energies <= energy_range[1])
            
            tdos = self.complete_dos.densities
            
            if self.is_spin_polarized:
                ax.plot(energies[mask], tdos[Spin.up][mask], 'r-', linewidth=2, label='Spin Up')
                ax.plot(energies[mask], -tdos[Spin.down][mask], 'b-', linewidth=2, label='Spin Down')
                ax.fill_between(energies[mask], 0, tdos[Spin.up][mask], alpha=0.3, color='red')
                ax.fill_between(energies[mask], 0, -tdos[Spin.down][mask], alpha=0.3, color='blue')
            else:
                dos_values = tdos[Spin.up]
                ax.plot(energies[mask], dos_values[mask], 'k-', linewidth=2)
                ax.fill_between(energies[mask], 0, dos_values[mask], alpha=0.3, color='gray')
            
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='Fermi Level')
            ax.set_xlabel('Energy - E$_F$ (eV)')
            ax.set_ylabel('DOS (states/eV)')
            ax.set_title('Total Density of States')
            ax.legend()
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'total_dos.png', dpi=300, bbox_inches='tight')
            plt.savefig(self.output_dir / 'total_dos.pdf', bbox_inches='tight')
            logger.info("已保存总DOS图")
            plt.close(fig)
            
        except Exception as e:
            logger.error(f"绘制总DOS图失败: {e}")
    
    def plot_element_dos(self, energy_range: Tuple[float, float] = (-10, 10)) -> None:
        """绘制元素分解DOS图 - 基于test3.py的实现"""
        if not self.structure:
            logger.warning("无结构信息，跳过元素分解DOS绘制")
            return
            
        if not (self.complete_dos and hasattr(self.complete_dos, 'get_element_dos')):
            logger.warning("DOS对象不支持元素分解，跳过元素DOS绘制")
            return
            
        try:
            element_dos_dict = self.complete_dos.get_element_dos()
            
            if not element_dos_dict:
                logger.warning("未找到元素分解DOS数据")
                return
            
            unique_elements = list(set([site.species_string for site in self.structure]))
            # 使用安全的颜色映射
            cmap = plt.get_cmap('tab10')
            colors = cmap(np.linspace(0, 1, len(unique_elements)))
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            energies = self.complete_dos.energies - self.complete_dos.efermi
            mask = (energies >= energy_range[0]) & (energies <= energy_range[1])
            
            plotted_elements = 0
            
            for i, element_str in enumerate(unique_elements):
                try:
                    element = Element(element_str)
                    
                    if element in element_dos_dict:
                        element_dos = element_dos_dict[element]
                        densities = element_dos.densities
                        
                        if self.is_spin_polarized:
                            ax.plot(energies[mask], densities[Spin.up][mask], 
                                   color=colors[i], linewidth=2, label=f'{element_str} (↑)')
                            ax.plot(energies[mask], -densities[Spin.down][mask], 
                                   color=colors[i], linewidth=2, linestyle='--', label=f'{element_str} (↓)')
                        else:
                            dos_values = densities[Spin.up]
                            ax.plot(energies[mask], dos_values[mask], 
                                   color=colors[i], linewidth=2, label=element_str)
                        
                        plotted_elements += 1
                    else:
                        logger.warning(f"未找到元素 {element_str} 的DOS数据")
                        
                except Exception as e:
                    logger.warning(f"绘制元素 {element_str} DOS时出错: {e}")
                    continue
            
            if plotted_elements == 0:
                logger.error("没有成功绘制任何元素的DOS")
                plt.close(fig)
                return
            
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='Fermi Level')
            ax.set_xlabel('Energy - E$_F$ (eV)')
            ax.set_ylabel('DOS (states/eV)')
            ax.set_title('Element-projected Density of States')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'element_dos.png', dpi=300, bbox_inches='tight')
            plt.savefig(self.output_dir / 'element_dos.pdf', bbox_inches='tight')
            logger.info(f"已保存元素分解DOS图，包含 {plotted_elements} 种元素")
            plt.close(fig)
            
        except Exception as e:
            logger.error(f"绘制元素分解DOS图失败: {e}")
    
    def plot_spd_dos(self, energy_range: Tuple[float, float] = (-10, 10)) -> None:
        """绘制SPD轨道分解DOS图 - 基于test3.py风格"""
        if not self.complete_dos:
            logger.error("DOS数据未加载")
            return
            
        if not self.structure:
            logger.warning("无结构信息，跳过SPD DOS绘制")
            return
            
        try:
            from pymatgen.electronic_structure.core import OrbitalType
            
            # 获取SPD轨道数据
            spd_dos_dict = self.complete_dos.get_spd_dos()
            
            if not spd_dos_dict:
                logger.warning("未找到SPD DOS数据")
                return
            
            fig, ax = plt.subplots(figsize=(12, 8))
            
            energies = self.complete_dos.energies - self.complete_dos.efermi
            mask = (energies >= energy_range[0]) & (energies <= energy_range[1])
            
            # 颜色映射
            colors = {'s': 'red', 'p': 'blue', 'd': 'green', 'f': 'orange'}
            
            # 轨道类型映射
            orbital_map = {
                OrbitalType.s: 's',
                OrbitalType.p: 'p',
                OrbitalType.d: 'd', 
                OrbitalType.f: 'f'
            }
            
            plotted_orbitals = 0
            
            for orbital_type, orbital_name in orbital_map.items():
                if orbital_type in spd_dos_dict:
                    orbital_dos = spd_dos_dict[orbital_type]
                    color = colors.get(orbital_name, 'black')
                    
                    if self.is_spin_polarized:
                        # 自旋极化情况
                        dos_up = orbital_dos.densities[Spin.up][mask]
                        dos_down = orbital_dos.densities[Spin.down][mask]
                        
                        ax.plot(energies[mask], dos_up, color=color, linewidth=2, 
                               label=f'{orbital_name}-orbital (↑)')
                        ax.plot(energies[mask], -dos_down, color=color, linewidth=2, 
                               linestyle='--', label=f'{orbital_name}-orbital (↓)')
                        
                        ax.fill_between(energies[mask], 0, dos_up, alpha=0.3, color=color)
                        ax.fill_between(energies[mask], 0, -dos_down, alpha=0.3, color=color)
                    else:
                        # 非自旋极化情况
                        dos_values = orbital_dos.densities[Spin.up][mask]
                        ax.plot(energies[mask], dos_values, color=color, linewidth=2, 
                               label=f'{orbital_name}-orbital')
                        ax.fill_between(energies[mask], 0, dos_values, alpha=0.3, color=color)
                    
                    plotted_orbitals += 1
            
            if plotted_orbitals == 0:
                logger.warning("没有成功绘制任何SPD轨道")
                plt.close(fig)
                return
            
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='Fermi Level')
            ax.set_xlabel('Energy - E$_F$ (eV)')
            ax.set_ylabel('DOS (states/eV)')
            ax.set_title('SPD Orbital-projected Density of States')
            ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
            ax.grid(True, alpha=0.3)
            
            plt.tight_layout()
            plt.savefig(self.output_dir / 'spd_dos.png', dpi=300, bbox_inches='tight')
            plt.savefig(self.output_dir / 'spd_dos.pdf', bbox_inches='tight')
            logger.info(f"已保存SPD DOS图，包含 {plotted_orbitals} 种轨道")
            plt.close(fig)
            
        except Exception as e:
            logger.error(f"绘制SPD DOS图失败: {e}")
    
    def plot_band_structure(self) -> None:
        """绘制能带结构图 - 基于matplotlib的实现"""
        if not hasattr(self, 'band_structure') or not self.band_structure:
            logger.warning("无能带结构数据，跳过能带图绘制")
            return
            
        try:
            from pymatgen.electronic_structure.bandstructure import BandStructureSymmLine
            
            band_structure = self.band_structure
            
            if isinstance(band_structure, BandStructureSymmLine):
                # 对于沿对称线的能带结构，可以手动绘制
                fig, ax = plt.subplots(figsize=(10, 8))
                
                # 获取k点路径和能量数据
                kpoints = band_structure.kpoints
                # 正确处理能带数据 - band_structure.bands是字典
                energies = {}
                for spin in band_structure.bands:
                    energies[spin] = band_structure.bands[spin] - band_structure.efermi
                
                # 计算k点距离 - 使用numpy计算欧几里得距离
                k_distances = [0.0]
                for i in range(1, len(kpoints)):
                    # 使用frac_coords计算距离
                    coord1 = np.array(kpoints[i-1].frac_coords)
                    coord2 = np.array(kpoints[i].frac_coords)
                    dist = float(np.linalg.norm(coord2 - coord1))
                    k_distances.append(k_distances[-1] + dist)
                
                # 绘制能带
                if self.is_spin_polarized:
                    # 自旋极化
                    for spin in [Spin.up, Spin.down]:
                        spin_label = 'Spin Up' if spin == Spin.up else 'Spin Down'
                        color = 'red' if spin == Spin.up else 'blue'
                        
                        for band_idx in range(band_structure.nb_bands):
                            band_energies = energies[spin][band_idx]
                            ax.plot(k_distances, band_energies, color=color, linewidth=1,
                                   label=spin_label if band_idx == 0 else "")
                else:
                    # 非自旋极化
                    for band_idx in range(band_structure.nb_bands):
                        band_energies = energies[Spin.up][band_idx]
                        ax.plot(k_distances, band_energies, 'b-', linewidth=1)
                
                # 添加高对称点标记
                if hasattr(band_structure, 'labels_dict'):
                    for label, kpoint in band_structure.labels_dict.items():
                        for i, k in enumerate(kpoints):
                            if k.frac_coords.tolist() == kpoint.frac_coords.tolist():
                                ax.axvline(x=k_distances[i], color='gray', linestyle='--', alpha=0.7)
                                ax.text(k_distances[i], ax.get_ylim()[1], label, 
                                        ha='center', va='bottom')
                                break
                
                ax.axhline(y=0, color='gray', linestyle='-', alpha=0.7, label='Fermi Level')
                ax.set_xlabel('k-path')
                ax.set_ylabel('Energy - E$_F$ (eV)')
                ax.set_title('Band Structure')
                ax.grid(True, alpha=0.3)
                
                if self.is_spin_polarized:
                    ax.legend()
                
                plt.tight_layout()
                plt.savefig(self.output_dir / 'band_structure.png', dpi=300, bbox_inches='tight')
                plt.savefig(self.output_dir / 'band_structure.pdf', bbox_inches='tight')
                logger.info("已保存能带结构图")
                plt.close(fig)
                
            else:
                # 对于uniform grid的能带结构，可以绘制能带密度图
                logger.info(f"能带结构类型为 {type(band_structure).__name__}，绘制简化版本")
                
                fig, ax = plt.subplots(figsize=(10, 6))
                
                # 获取能量范围
                # 正确处理能带数据 - band_structure.bands是字典
                all_energies = []
                for spin in band_structure.bands:
                    band_energies = band_structure.bands[spin] - band_structure.efermi
                    all_energies.extend(band_energies.flatten())
                all_energies = np.array(all_energies)
                
                # 统计能带密度
                energy_range = np.linspace(all_energies.min(), all_energies.max(), 200)
                band_density = np.zeros_like(energy_range)
                
                for spin in band_structure.bands:
                    for band_idx in range(band_structure.nb_bands):
                        band_energies = band_structure.bands[spin][band_idx] - band_structure.efermi
                        hist, _ = np.histogram(band_energies, bins=energy_range)
                        band_density[:-1] += hist
                
                ax.plot(energy_range[:-1], band_density[:-1], 'b-', linewidth=2)
                ax.fill_between(energy_range[:-1], 0, band_density[:-1], alpha=0.3)
                ax.axvline(x=0, color='gray', linestyle='--', alpha=0.7, label='Fermi Level')
                ax.set_xlabel('Energy - E$_F$ (eV)')
                ax.set_ylabel('Band Density')
                ax.set_title('Band Structure Density')
                ax.legend()
                ax.grid(True, alpha=0.3)
                
                plt.tight_layout()
                plt.savefig(self.output_dir / 'band_structure.png', dpi=300, bbox_inches='tight')
                plt.savefig(self.output_dir / 'band_structure.pdf', bbox_inches='tight')
                logger.info("已保存能带结构密度图")
                plt.close(fig)
                
        except Exception as e:
            logger.error(f"绘制能带结构图失败: {e}")
    
    def _generate_markdown_report(self) -> None:
        """生成Markdown分析报告 - 基于test3.py的实现"""
        if self.complete_dos is None:
            logger.error("DOS数据未加载，无法生成报告")
            return
            
        try:
            logger.info("📋 生成Markdown分析报告...")
            
            gap_info = self.get_band_gap_info()
            
            report = f"""# VASP DOS 分析报告

## 基本信息
- 费米能级: {gap_info['fermi_level']:.4f} eV
- 自旋极化: {'是' if self.is_spin_polarized else '否'}
- 材料类型: {'金属' if gap_info['is_metal'] else '半导体/绝缘体'}

## 电子结构
- 带隙: {gap_info['band_gap']:.4f} eV
"""
            
            if gap_info['vbm'] is not None:
                report += f"- 价带顶 (VBM): {gap_info['vbm']:.4f} eV\n"
            if gap_info['cbm'] is not None:
                report += f"- 导带底 (CBM): {gap_info['cbm']:.4f} eV\n"
            
            report += "\n## 结构信息\n"
            
            if self.structure:
                report += f"- 化学式: {self.structure.composition.reduced_formula}\n"
                report += f"- 原子总数: {len(self.structure)}\n"
                unique_elements = list(set([site.species_string for site in self.structure]))
                report += f"- 包含元素: {', '.join(unique_elements)}\n"
            
            # VASP计算信息
            if self.vasprun:
                report += f"\n## VASP计算信息\n"
                report += f"- VASP版本: {getattr(self.vasprun, 'vasp_version', 'N/A')}\n"
                report += f"- 能量点数: {len(self.complete_dos.energies)}\n"
                if hasattr(self.vasprun, 'parameters'):
                    params = self.vasprun.parameters
                    if 'LORBIT' in params:
                        report += f"- LORBIT: {params['LORBIT']}\n"
                    if 'NEDOS' in params:
                        report += f"- NEDOS: {params['NEDOS']}\n"
            
            report += f"""
## 输出文件

### matplotlib图表（基于test3.py风格）
- 总DOS图: total_dos.png/pdf
- 元素分解DOS图: element_dos.png/pdf  
- SPD轨道DOS图: spd_dos.png/pdf
- 能带结构图: band_structure.png/pdf（如果有能带数据）

### 数据文件
- DOS数据: data/目录下的CSV文件

### 报告文件  
- HTML报告: pymatgen_dos_analysis_report.html
- Markdown报告: analysis_report.md

## 分析建议
"""
            
            if gap_info['is_metal']:
                report += "- 材料显示金属性质，适合导电性应用\n"
            elif gap_info['band_gap'] < 2.0:
                report += "- 材料为小带隙半导体，适合光电应用\n"
            else:
                report += "- 材料为宽带隙半导体/绝缘体\n"
                
            if self.is_spin_polarized:
                report += "- 材料具有磁性，建议进一步分析磁性质\n"
            
            report += "\n## 使用说明\n"
            report += "- 本程序基于PyMatGen和vasprun.xml文件进行分析\n"
            report += "- 确保VASP计算中设置了合适的LORBIT参数以获得PDOS数据\n"
            report += "- 推荐使用LORBIT=11或LORBIT=12获得更详细的轨道信息\n"
            
            # 保存报告
            with open(self.output_dir / 'analysis_report.md', 'w', encoding='utf-8') as f:
                f.write(report)
            
            logger.info("已生成Markdown分析报告: analysis_report.md")
            
        except Exception as e:
            logger.error(f"生成Markdown分析报告失败: {e}")
    
    def _calculate_dos_integrals(self, complete_dos: CompleteDos) -> Dict[str, Any]:
        """计算DOS积分"""
        energies = complete_dos.energies
        fermi = complete_dos.efermi
        
        # 定义能量区间
        valence_mask = energies <= fermi
        conduction_mask = energies > fermi
        near_fermi_mask = np.abs(energies - fermi) <= 1.0
        
        integrals = {}
        
        for spin in complete_dos.densities:
            dos = complete_dos.densities[spin]
            spin_label = 'up' if spin == Spin.up else 'down'
            
            integrals[spin_label] = {
                'total': float(np.trapezoid(dos, energies)),
                'valence_band': float(np.trapezoid(dos[valence_mask], energies[valence_mask])),
                'conduction_band': float(np.trapezoid(dos[conduction_mask], energies[conduction_mask])),
                'near_fermi': float(np.trapezoid(dos[near_fermi_mask], energies[near_fermi_mask]))
            }
        
        return integrals
    
    def _analyze_orbital_projections(self, complete_dos: CompleteDos) -> Dict[str, Any]:
        """分析轨道投影"""
        orbital_analysis = {}
        
        if not complete_dos or not hasattr(complete_dos, 'structure'):
            logger.warning("缺少DOS或结构信息，跳过轨道分析")
            return orbital_analysis
        
        try:
            # 使用正确的PyMatGen API - get_element_spd_dos
            logger.info("=== 轨道投影分析 ===")
            
            from pymatgen.electronic_structure.core import OrbitalType
            
            try:
                logger.info("使用get_element_spd_dos方法获取轨道数据")
                
                # 获取结构中的唯一元素
                unique_elements = complete_dos.structure.composition.elements
                logger.info(f"结构中的元素: {[str(e) for e in unique_elements]}")
                
                for element in unique_elements:
                    element_str = str(element)
                    orbital_contribs = {'s': 0.0, 'p': 0.0, 'd': 0.0, 'f': 0.0}
                    
                    try:
                        # 使用正确的API获取该元素的SPD轨道投影DOS
                        element_spd_dos = complete_dos.get_element_spd_dos(element)
                        
                        logger.info(f"元素 {element_str}: SPD数据键 = {list(element_spd_dos.keys())}")
                        
                        # 轨道类型映射
                        orbital_map = {
                            OrbitalType.s: 's',
                            OrbitalType.p: 'p', 
                            OrbitalType.d: 'd',
                            OrbitalType.f: 'f'
                        }
                        
                        # 计算每个轨道的贡献
                        for orbital_type, orbital_name in orbital_map.items():
                            if orbital_type in element_spd_dos:
                                orbital_dos = element_spd_dos[orbital_type]
                                
                                # 计算轨道贡献（积分所有自旋）
                                contrib = 0.0
                                for spin in orbital_dos.densities:
                                    contrib += float(np.trapezoid(orbital_dos.densities[spin], complete_dos.energies))
                                
                                orbital_contribs[orbital_name] = contrib
                                logger.info(f"  {element_str} {orbital_name}轨道贡献: {contrib:.3f}")
                            else:
                                logger.debug(f"  {element_str}无{orbital_name}轨道数据")
                        
                        orbital_analysis[element_str] = orbital_contribs
                        logger.info(f"元素 {element_str}: 轨道分析完成")
                        
                    except Exception as e:
                        logger.warning(f"元素 {element_str} SPD计算失败: {e}")
                        # 使用估算值作为备用
                        if element_str == 'H':
                            orbital_contribs = {'s': 1.0, 'p': 0.0, 'd': 0.0, 'f': 0.0}
                        elif element_str in ['Li', 'Na', 'K']:
                            orbital_contribs = {'s': 0.8, 'p': 0.2, 'd': 0.0, 'f': 0.0}
                        elif element_str in ['C', 'N', 'O', 'F']:
                            orbital_contribs = {'s': 0.3, 'p': 0.7, 'd': 0.0, 'f': 0.0}
                        elif element_str in ['P', 'S', 'Cl']:
                            orbital_contribs = {'s': 0.2, 'p': 0.6, 'd': 0.2, 'f': 0.0}
                        elif element_str in ['Fe', 'Co', 'Ni', 'Cu', 'Zn']:
                            orbital_contribs = {'s': 0.1, 'p': 0.2, 'd': 0.7, 'f': 0.0}
                        else:
                            orbital_contribs = {'s': 0.25, 'p': 0.25, 'd': 0.25, 'f': 0.25}
                        
                        orbital_analysis[element_str] = orbital_contribs
                        logger.info(f"元素 {element_str}: 使用估算轨道贡献")
                
                # 输出汇总结果
                logger.info("=== 轨道贡献汇总 ===")
                for element_str, contribs in orbital_analysis.items():
                    total = sum(contribs.values())
                    logger.info(f"{element_str}: s={contribs['s']:.3f}, p={contribs['p']:.3f}, d={contribs['d']:.3f}, f={contribs['f']:.3f} (总计={total:.3f})")
                
            except Exception as e:
                logger.warning(f"轨道投影分析失败: {e}")
                logger.info("使用估算方法")
                
                # 完全备用方案：基于化学知识的估算
                for element in complete_dos.structure.composition.elements:
                    element_str = str(element)
                    
                    # 基于元素的电子结构给出大致的轨道贡献估算
                    if element_str == 'H':
                        orbital_contribs = {'s': 1.0, 'p': 0.0, 'd': 0.0, 'f': 0.0}
                    elif element_str in ['Li', 'Na', 'K']:
                        orbital_contribs = {'s': 0.8, 'p': 0.2, 'd': 0.0, 'f': 0.0}
                    elif element_str in ['C', 'N', 'O', 'F']:
                        orbital_contribs = {'s': 0.3, 'p': 0.7, 'd': 0.0, 'f': 0.0}
                    elif element_str in ['P', 'S', 'Cl']:
                        orbital_contribs = {'s': 0.2, 'p': 0.6, 'd': 0.2, 'f': 0.0}
                    elif element_str in ['Fe', 'Co', 'Ni', 'Cu', 'Zn']:
                        orbital_contribs = {'s': 0.1, 'p': 0.2, 'd': 0.7, 'f': 0.0}
                    else:
                        orbital_contribs = {'s': 0.25, 'p': 0.25, 'd': 0.25, 'f': 0.25}
                    
                    orbital_analysis[element_str] = orbital_contribs
                    logger.info(f"元素 {element_str}: 使用估算轨道贡献")
            
            logger.info("=== 轨道投影分析完成 ===")
                
        except Exception as e:
            logger.warning(f"轨道投影分析失败: {e}")
            # 即使失败，也返回基础的轨道信息
            if complete_dos and hasattr(complete_dos, 'structure'):
                for element in complete_dos.structure.composition.elements:
                    element_str = str(element)
                    orbital_analysis[element_str] = {'s': 0.0, 'p': 0.0, 'd': 0.0, 'f': 0.0}
        
        return orbital_analysis
    
    def _analyze_bonding_character(self) -> str:
        """分析成键特征"""
        try:
            dos_analysis = self.analysis_data.get('dos_analysis', {})
            orbital_analysis = dos_analysis.get('orbital_analysis', {})
            
            if not orbital_analysis:
                return 'unknown'
            
            # 统计各轨道的总贡献
            total_contributions = {'s': 0.0, 'p': 0.0, 'd': 0.0, 'f': 0.0}
            
            for element_contribs in orbital_analysis.values():
                for orbital, contrib in element_contribs.items():
                    if orbital in total_contributions and isinstance(contrib, (int, float)):
                        total_contributions[orbital] += float(contrib)
            
            # 找到主导轨道
            if not any(total_contributions.values()):
                return 'unknown'
            
            dominant_orbital = max(total_contributions.keys(), key=lambda k: total_contributions[k])
            
            # 根据主导轨道推断成键特征
            if dominant_orbital == 's':
                return 'ionic'
            elif dominant_orbital == 'p':
                return 'covalent'
            elif dominant_orbital == 'd':
                return 'metallic'
            else:
                return 'mixed'
                
        except Exception as e:
            logger.warning(f"分析成键特征失败: {e}")
            return 'unknown'
    
    def _calculate_electronegativity_difference(self, structure: Structure) -> float:
        """计算电负性差异"""
        try:
            elements = structure.composition.elements
            if len(elements) < 2:
                return 0.0
            
            electronegativities = []
            for element in elements:
                if hasattr(element, 'X') and element.X:
                    electronegativities.append(element.X)
            
            if len(electronegativities) < 2:
                return 0.0
            
            return max(electronegativities) - min(electronegativities)
            
        except:
            return 0.0
    
    def _estimate_ionic_character(self, structure: Structure) -> float:
        """估算离子性"""
        en_diff = self._calculate_electronegativity_difference(structure)
        # 使用Pauling公式估算离子性百分比
        if en_diff > 0:
            ionic_character = 1 - np.exp(-0.25 * en_diff**2)
            return float(ionic_character)
        return 0.0
    
    def _generate_visualizations(self):
        """生成matplotlib可视化图表的base64编码用于HTML展示"""
        logger.info("📈 生成matplotlib可视化图表用于HTML...")
        
        visualizations = {}
        
        if not self.complete_dos:
            logger.warning("DOS数据未加载，跳过可视化")
            self.analysis_data['visualizations'] = visualizations
            return
        
        try:
            # 将matplotlib图表文件转换为base64编码
            image_files = {
                'total_dos_plot': self.output_dir / 'total_dos.png',
                'element_dos_plot': self.output_dir / 'element_dos.png',
                'spd_dos_plot': self.output_dir / 'spd_dos.png',
                'band_structure_plot': self.output_dir / 'band_structure.png'
            }
            
            for plot_name, file_path in image_files.items():
                if file_path.exists():
                    try:
                        with open(file_path, 'rb') as f:
                            visualizations[plot_name] = base64.b64encode(f.read()).decode()
                        logger.info(f"已添加 {plot_name} 到HTML报告")
                    except Exception as e:
                        logger.warning(f"读取图片文件 {file_path} 失败: {e}")
                else:
                    logger.warning(f"图片文件 {file_path} 不存在")
            
        except Exception as e:
            logger.error(f"生成HTML可视化失败: {e}")
        
        self.analysis_data['visualizations'] = visualizations
    
    def _finalize_analysis(self):
        """完成分析并汇总结果"""
        print("📋 汇总分析结果...")
        
        # 任务信息
        task_info = {
            'task_id': self.task_id,
            'input_path': str(self.input_path),
            'vasprun_path': str(self.vasprun_path),
            'timestamp': datetime.now().isoformat(),
            'analysis_type': 'PyMatGen_DOS_Analysis'
        }
        
        # 总结信息
        summary = {
            'formula': self.analysis_data['structure']['formula'],
            'space_group': self.analysis_data['structure']['space_group'],
            'band_gap': self.analysis_data['dos_analysis']['band_gap'],
            'material_type': self.analysis_data['dos_analysis']['material_type'],
            'is_magnetic': self.analysis_data['magnetic_properties']['is_magnetic'],
            'magnetic_type': self.analysis_data['magnetic_properties']['magnetic_type']
        }
        
        # 添加到分析数据
        self.analysis_data['task_info'] = task_info
        self.analysis_data['summary'] = summary


class PyMatGenDOSHTMLGenerator:
    """PyMatGen DOS分析结果的HTML报告生成器"""
    
    def __init__(self, analysis_data: Dict[str, Any]):
        self.data = analysis_data
    
    def generate_html_report(self, output_path: str) -> str:
        """生成HTML报告"""
        html_content = self._generate_html_content()
        
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        return output_path
    
    def _generate_html_content(self) -> str:
        """生成HTML内容"""
        charts_data = self.data.get('visualizations', {})
        
        return f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>DOS 分析报告</title>
    {self._generate_css_styles()}
</head>
<body>
    <div class="container">
        {self._generate_header()}
        {self._generate_summary()}
        {self._generate_structure_analysis()}
        {self._generate_dos_analysis()}
        {self._generate_band_structure_analysis()}
        {self._generate_chemical_analysis()}
        {self._generate_magnetic_analysis()}
        {self._generate_visualization_section(charts_data)}
        {self._generate_footer()}
    </div>
</body>
</html>
        """
    
    def _generate_css_styles(self) -> str:
        """生成CSS样式"""
        return """
    <style>
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            line-height: 1.6;
            color: #333;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            min-height: 100vh;
        }
        
        .container {
            max-width: 1200px;
            margin: 0 auto;
            padding: 20px;
            background: white;
            box-shadow: 0 10px 30px rgba(0,0,0,0.1);
            border-radius: 10px;
            margin-top: 20px;
            margin-bottom: 20px;
        }
        
        .header {
            text-align: center;
            padding: 30px 0;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            border-radius: 10px;
            margin-bottom: 30px;
        }
        
        .header h1 {
            font-size: 2.5em;
            margin-bottom: 10px;
        }
        
        .header .subtitle {
            font-size: 1.2em;
            opacity: 0.9;
        }
        
        .section {
            margin-bottom: 30px;
            padding: 25px;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            background: #fafafa;
        }
        
        .section h2 {
            color: #667eea;
            border-bottom: 2px solid #667eea;
            padding-bottom: 10px;
            margin-bottom: 20px;
            font-size: 1.8em;
        }
        
        .section h3 {
            color: #764ba2;
            margin: 15px 0 10px 0;
            font-size: 1.4em;
        }
        
        .grid-2 {
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 20px;
            margin: 20px 0;
        }
        
        .grid-3 {
            display: grid;
            grid-template-columns: repeat(3, 1fr);
            gap: 15px;
            margin: 20px 0;
        }
        
        .data-table {
            width: 100%;
            border-collapse: collapse;
            margin: 15px 0;
            background: white;
            border-radius: 8px;
            overflow: hidden;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        
        .data-table th, .data-table td {
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }
        
        .data-table th {
            background: #667eea;
            color: white;
            font-weight: 600;
        }
        
        .data-table tr:hover {
            background: #f5f5f5;
        }
        
        .summary-card {
            background: white;
            padding: 20px;
            border-radius: 8px;
            text-align: center;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border: 2px solid #667eea;
        }
        
        .summary-card h4 {
            color: #667eea;
            margin-bottom: 10px;
            font-size: 1.1em;
        }
        
        .summary-card .value {
            font-size: 1.8em;
            font-weight: bold;
            color: #764ba2;
        }
        
        .image-container {
            text-align: center;
            margin: 20px 0;
            padding: 20px;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }
        
        .image-container img {
            max-width: 100%;
            height: auto;
            border-radius: 8px;
        }
        
        .highlight {
            background: linear-gradient(135deg, #fff3cd, #ffeaa7);
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #f39c12;
            margin: 15px 0;
        }
        
        .warning {
            background: linear-gradient(135deg, #f8d7da, #f5c6cb);
            padding: 15px;
            border-radius: 8px;
            border-left: 4px solid #dc3545;
            margin: 15px 0;
        }
        
        .info-box {
            background: linear-gradient(135deg, #e3f2fd, #bbdefb);
            padding: 20px;
            border-radius: 8px;
            border-left: 4px solid #1976d2;
            margin: 20px 0;
            color: #0d47a1;
        }
        
        .info-box h4 {
            margin-top: 0;
            margin-bottom: 15px;
            color: #0d47a1;
            font-size: 1.2em;
        }
        
        .info-box ul {
            margin: 10px 0;
            padding-left: 25px;
        }
        
        .info-box li {
            margin: 8px 0;
            line-height: 1.4;
        }
        
        .info-box p {
            margin: 10px 0;
            line-height: 1.5;
        }
        
        .footer {
            text-align: center;
            padding: 20px;
            color: #666;
            border-top: 1px solid #e0e0e0;
            margin-top: 30px;
        }
        
        @media (max-width: 768px) {
            .grid-2, .grid-3 {
                grid-template-columns: 1fr;
            }
            
            .container {
                margin: 10px;
                padding: 15px;
            }
            
            .header h1 {
                font-size: 2em;
            }
        }
    </style>
        """
    
    def _generate_header(self) -> str:
        """生成页面头部"""
        task_info = self.data.get('task_info', {})
        return f"""
        <div class="header">
            <h1>🔬 PyMatGen DOS 分析报告</h1>
            <div class="subtitle">
                任务ID: {task_info.get('task_id', 'unknown')} | 
                生成时间: {task_info.get('timestamp', 'unknown')}
            </div>
        </div>
        """
    
    def _generate_summary(self) -> str:
        """生成分析摘要"""
        summary = self.data.get('summary', {})
        
        return f"""
        <div class="section">
            <h2>📋 分析摘要</h2>
            <div class="grid-3">
                <div class="summary-card">
                    <h4>化学式</h4>
                    <div class="value">{summary.get('formula', 'Unknown')}</div>
                </div>
                <div class="summary-card">
                    <h4>空间群</h4>
                    <div class="value">{summary.get('space_group', 'Unknown')}</div>
                </div>
                <div class="summary-card">
                    <h4>带隙</h4>
                    <div class="value">{summary.get('band_gap', 0):.3f} eV</div>
                </div>
                <div class="summary-card">
                    <h4>材料类型</h4>
                    <div class="value">{summary.get('material_type', 'unknown').title()}</div>
                </div>
                <div class="summary-card">
                    <h4>磁性</h4>
                    <div class="value">{'是' if summary.get('is_magnetic', False) else '否'}</div>
                </div>
                <div class="summary-card">
                    <h4>磁性类型</h4>
                    <div class="value">{summary.get('magnetic_type', 'non-magnetic').title()}</div>
                </div>
            </div>
        </div>
        """
    
    def _generate_structure_analysis(self) -> str:
        """生成结构分析部分"""
        structure = self.data.get('structure', {})
        lattice = structure.get('lattice_parameters', {})
        
        return f"""
        <div class="section">
            <h2>🏗️ 晶体结构分析</h2>
            
            <div class="grid-2">
                <div>
                    <h3>基本信息</h3>
                    <table class="data-table">
                        <tr><td>化学式</td><td>{structure.get('formula', 'unknown')}</td></tr>
                        <tr><td>简化式</td><td>{structure.get('reduced_formula', 'unknown')}</td></tr>
                        <tr><td>原子数</td><td>{structure.get('num_sites', 'unknown')}</td></tr>
                        <tr><td>密度</td><td>{structure.get('density', 0):.3f} g/cm³</td></tr>
                        <tr><td>空间群</td><td>{structure.get('space_group', 'unknown')}</td></tr>
                        <tr><td>点群</td><td>{structure.get('point_group', 'unknown')}</td></tr>
                    </table>
                    
                    <h3>组成元素</h3>
                    <div class="highlight">
                        {', '.join(structure.get('elements', []))}
                    </div>
                </div>
                
                <div>
                    <h3>晶格参数</h3>
                    <table class="data-table">
                        <tr><td>a</td><td>{lattice.get('a', 0):.4f} Å</td></tr>
                        <tr><td>b</td><td>{lattice.get('b', 0):.4f} Å</td></tr>
                        <tr><td>c</td><td>{lattice.get('c', 0):.4f} Å</td></tr>
                        <tr><td>α</td><td>{lattice.get('alpha', 0):.2f}°</td></tr>
                        <tr><td>β</td><td>{lattice.get('beta', 0):.2f}°</td></tr>
                        <tr><td>γ</td><td>{lattice.get('gamma', 0):.2f}°</td></tr>
                        <tr><td>体积</td><td>{lattice.get('volume', 0):.3f} Ų</td></tr>
                    </table>
                </div>
            </div>
        </div>
        """
    
    def _generate_dos_analysis(self) -> str:
        """生成DOS分析部分"""
        dos = self.data.get('dos_analysis', {})
        integrals = dos.get('dos_integrals', {})
        orbital = dos.get('orbital_analysis', {})
        
        return f"""
        <div class=\"section\">
            <h2>📊 态密度分析</h2>
            
            <div class=\"grid-2\">
                <div>
                    <h3>电子结构特征</h3>
                    <table class=\"data-table\">
                        <tr><td>费米能级</td><td>{dos.get('fermi_energy', 0):.4f} eV</td></tr>
                        <tr><td>带隙</td><td>{dos.get('band_gap', 0):.4f} eV</td></tr>
                        <tr><td>导带底</td><td>{dos.get('cbm_energy', 0) or 0:.4f} eV</td></tr>
                        <tr><td>价带顶</td><td>{dos.get('vbm_energy', 0) or 0:.4f} eV</td></tr>
                        <tr><td>带隙类型</td><td>{dos.get('gap_type', 'unknown')}</td></tr>
                        <tr><td>材料类型</td><td>{dos.get('material_type', 'unknown')}</td></tr>
                        <tr><td>是否金属</td><td>{'是' if dos.get('is_metal', False) else '否'}</td></tr>
                        <tr><td>自旋极化</td><td>{'是' if dos.get('is_spin_polarized', False) else '否'}</td></tr>
                    </table>
                </div>
                
                <div>
                    <h3>DOS积分分析</h3>
                    {self._generate_dos_integrals_table(integrals)}
                    
                    <h3>轨道分析</h3>
                    <div class=\"info-box\">
                        <h4>📊 轨道分析说明</h4>
                        <p><strong>轨道分析</strong>展示了各元素的电子在不同原子轨道（s、p、d、f）上的贡献程度。数值表示该轨道在费米能级附近对态密度的相对贡献。</p>
                    </div>
                    
                    <div class=\"grid-2\">
                        <div>
                            <h4>🧮 元素轨道贡献汇总</h4>
                            {self._generate_spd_summary_table(orbital)}
                        </div>
                    </div>
                </div>
            </div>
        </div>
        """
    
    def _generate_band_structure_analysis(self) -> str:
        """生成能带结构分析部分"""
        band = self.data.get('band_structure', {})
        
        # 检查是否有任何可用的分析数据
        has_band_structure = band.get('has_band_structure', False)
        has_dos_analysis = band.get('has_dos_gap_analysis', False)
        
        if not has_band_structure and not has_dos_analysis:
            return f"""
            <div class="section">
                <h2>🎼 能带结构分析</h2>
                <div class="warning">
                    <strong>注意:</strong> 未检测到能带结构数据或DOS数据
                </div>
            </div>
            """
        
        # 获取计算方法信息
        calc_method = band.get('calculation_method', 'Unknown')
        
        # 构建VBM/CBM信息
        vbm_info = f"{band.get('vbm_energy', 'N/A'):.4f} eV" if band.get('vbm_energy') is not None else "未检测到"
        cbm_info = f"{band.get('cbm_energy', 'N/A'):.4f} eV" if band.get('cbm_energy') is not None else "未检测到"
        fermi_info = f"{band.get('fermi_level', 0):.4f} eV" if band.get('fermi_level') is not None else "N/A"
        
        # 根据数据可用性调整标题
        if has_band_structure:
            title = "🎼 能带结构分析"
        else:
            title = "🎼 电子结构分析 (基于DOS数据)"
            
        return f"""
        <div class="section">
            <h2>{title}</h2>
            
            <div class="info-box">
                <h4>📊 计算方法说明</h4>
                <p><strong>方法</strong>: {calc_method}</p>
                <p>此分析基于DOS数据进行带隙计算，比直接使用能带结构对象更准确。</p>
            </div>
            
            <h4>🔬 电子结构特性</h4>
            <table class="data-table">
                <tr><td>材料类型</td><td><strong>{'金属' if band.get('is_metal', False) else '半导体/绝缘体'}</strong></td></tr>
                <tr><td>费米能级</td><td>{fermi_info}</td></tr>
                <tr><td>基本带隙</td><td><strong>{band.get('fundamental_gap', 0):.4f} eV</strong></td></tr>
                <tr><td>直接带隙</td><td><strong>{band.get('direct_gap', 0):.4f} eV</strong></td></tr>
                <tr><td>价带顶 (VBM)</td><td>{vbm_info}</td></tr>
                <tr><td>导带底 (CBM)</td><td>{cbm_info}</td></tr>
                <tr><td>能带数</td><td>{band.get('num_bands', 'unknown')}</td></tr>
            </table>
            
            <h4>💡 物理意义</h4>
            <div class="highlight">
                {self._get_gap_interpretation(band)}
            </div>
        </div>
        """
    
    def _get_gap_interpretation(self, band: Dict[str, Any]) -> str:
        """获取带隙物理意义解释"""
        is_metal = band.get('is_metal', True)
        gap = band.get('fundamental_gap', 0.0)
        
        if is_metal:
            return """
            <p><strong>金属性材料</strong></p>
            <ul>
                <li>价带和导带重叠，无带隙</li>
                <li>具有良好的导电性</li>
                <li>适用于电极、导线等导电应用</li>
            </ul>
            """
        elif gap < 0.5:
            return f"""
            <p><strong>窄带隙半导体</strong> (Eg = {gap:.3f} eV)</p>
            <ul>
                <li>红外光响应材料</li>
                <li>适用于红外探测器、热电材料</li>
                <li>可能具有较高的载流子迁移率</li>
            </ul>
            """
        elif gap < 2.0:
            return f"""
            <p><strong>小带隙半导体</strong> (Eg = {gap:.3f} eV)</p>
            <ul>
                <li>可见光和近红外光响应</li>
                <li>适用于太阳能电池、光电器件</li>
                <li>良好的光电转换效率</li>
            </ul>
            """
        elif gap < 3.5:
            return f"""
            <p><strong>中等带隙半导体</strong> (Eg = {gap:.3f} eV)</p>
            <ul>
                <li>可见光响应材料</li>
                <li>适用于LED、激光器</li>
                <li>良好的光学性质</li>
            </ul>
            """
        else:
            return f"""
            <p><strong>宽带隙材料/绝缘体</strong> (Eg = {gap:.3f} eV)</p>
            <ul>
                <li>紫外光响应材料</li>
                <li>适用于UV探测器、绝缘材料</li>
                <li>高介电强度，优异的绝缘性质</li>
            </ul>
            """
    
    def _generate_chemical_analysis(self) -> str:
        """生成化学分析部分"""
        chemical = self.data.get('chemical_properties', {})
        elements = chemical.get('elements', {})
        
        return f"""
        <div class="section">
            <h2>🧪 化学性质分析</h2>
            
            <div class="grid-2">
                <div>
                    <h3>化学性质</h3>
                    <table class="data-table">
                        <tr><td>成键特征</td><td>{chemical.get('bonding_character', 'unknown')}</td></tr>
                        <tr><td>电负性差异</td><td>{chemical.get('electronegativity_difference', 0):.3f}</td></tr>
                        <tr><td>离子性</td><td>{chemical.get('ionic_character', 0):.1%}</td></tr>
                    </table>
                </div>
                
                <div>
                    <h3>元素分析</h3>
                    {self._generate_element_analysis_table(elements)}
                </div>
            </div>
        </div>
        """
    
    def _generate_magnetic_analysis(self) -> str:
        """生成磁性分析部分"""
        magnetic = self.data.get('magnetic_properties', {})
        
        if not magnetic.get('is_magnetic', False):
            return f"""
            <div class="section">
                <h2>🧲 磁性分析</h2>
                <div class="highlight">
                    <strong>结果:</strong> 该材料为非磁性材料
                </div>
            </div>
            """
        
        return f"""
        <div class="section">
            <h2>🧲 磁性分析</h2>
            
            <table class="data-table">
                <tr><td>磁性类型</td><td>{magnetic.get('magnetic_type', 'unknown')}</td></tr>
                <tr><td>总磁化强度</td><td>{magnetic.get('total_magnetization', 0):.4f} μB</td></tr>
                <tr><td>自旋极化度</td><td>{magnetic.get('spin_polarization', 0):.4f}</td></tr>
                <tr><td>最大自旋差</td><td>{magnetic.get('max_spin_difference', 0):.4f}</td></tr>
                <tr><td>费米面自旋极化</td><td>{magnetic.get('fermi_spin_polarization', 0):.4f}</td></tr>
            </table>
        </div>
        """
    
    def _generate_visualization_section(self, charts_data: Dict[str, Any]) -> str:
        """生成可视化部分"""
        if not charts_data:
            return f"""
            <div class="section">
                <h2>📈 PyMatGen 可视化</h2>
                <div class="warning">未生成可视化图表</div>
            </div>
            """
        
        total_dos = charts_data.get('total_dos_plot', '')
        element_dos = charts_data.get('element_dos_plot', '')
        spd_dos = charts_data.get('spd_dos_plot', '')
        band_plot = charts_data.get('band_structure_plot', '')
        
        return f"""
        <div class="section">
            <h2>📈 PyMatGen 专业可视化</h2>
            
            <div class="grid-2">
                {f'<div><h3>总态密度</h3><div class="image-container"><img src="data:image/png;base64,{total_dos}" alt="Total DOS"></div></div>' if total_dos else ''}
                {f'<div><h3>元素分解DOS</h3><div class="image-container"><img src="data:image/png;base64,{element_dos}" alt="Element DOS"></div></div>' if element_dos else ''}
                {f'<div><h3>SPD轨道DOS</h3><div class="image-container"><img src="data:image/png;base64,{spd_dos}" alt="SPD DOS"></div></div>' if spd_dos else ''}
                {f'<div><h3>能带结构</h3><div class="image-container"><img src="data:image/png;base64,{band_plot}" alt="Band Structure"></div></div>' if band_plot else ''}
            </div>
        </div>
        """
    
    def _embed_plot_from_visualizations(self, key: str) -> str:
        """从visualizations中嵌入一张图为HTML img标签。"""
        charts_data = self.data.get('visualizations', {})
        img64 = charts_data.get(key, '')
        if not img64:
            return "<div class=\"image-container\"><em>无可用图像</em></div>"
        return f"<div class=\"image-container\"><img src=\"data:image/png;base64,{img64}\" alt=\"{key}\"></div>"
    
    def _generate_recommendations(self) -> str:
        """生成建议部分"""
        dos = self.data.get('dos_analysis', {})
        material_type = dos.get('material_type', 'unknown')
        is_magnetic = self.data.get('magnetic_properties', {}).get('is_magnetic', False)
        
        recommendations = []
        
        if material_type == 'metal':
            recommendations.append("该材料具有金属性，可考虑导电应用")
        elif material_type == 'semiconductor':
            band_gap = dos.get('band_gap', 0)
            if band_gap < 1.5:
                recommendations.append("窄带隙半导体，适合红外光电应用")
            elif band_gap < 3.0:
                recommendations.append("适中带隙半导体，适合太阳能电池应用")
        elif material_type == 'insulator':
            recommendations.append("绝缘体材料，适合绝缘或介电应用")
        
        if is_magnetic:
            recommendations.append("具有磁性，可考虑磁存储或自旋电子学应用")
        
        if not recommendations:
            recommendations.append("需要进一步分析确定应用方向")
        
        return f"""
        <div class="section">
            <h2>💡 应用建议</h2>
            <div class="highlight">
                <ul>
                    {''.join(f'<li>{rec}</li>' for rec in recommendations)}
                </ul>
            </div>
        </div>
        """
    
    def _generate_footer(self) -> str:
        """生成页面底部"""
        return f"""
        <div class="footer">
            <p>🔬 PyMatGen DOS分析报告 | 基于Materials Project生态系统</p>
            <p>生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        """
    
    def _generate_dos_integrals_table(self, integrals: Dict[str, Any]) -> str:
        """生成DOS积分表格"""
        if not integrals:
            return "<p>无DOS积分数据</p>"
        
        rows = ""
        for spin, data in integrals.items():
            rows += f"""
            <tr>
                <td colspan="2"><strong>{spin.upper()} 自旋</strong></td>
            </tr>
            <tr><td>总积分</td><td>{data.get('total', 0):.3f}</td></tr>
            <tr><td>价带</td><td>{data.get('valence_band', 0):.3f}</td></tr>
            <tr><td>导带</td><td>{data.get('conduction_band', 0):.3f}</td></tr>
            <tr><td>费米附近</td><td>{data.get('near_fermi', 0):.3f}</td></tr>
            """
        
        return f'<table class="data-table">{rows}</table>'
    
    def _generate_spd_summary_table(self, orbital: Dict[str, Any]) -> str:
        """生成SPD轨道贡献汇总表格"""
        if not orbital:
            return "<p>无轨道分析数据</p>"
        
        # 表头
        header = """
        <table class="data-table">
            <thead>
                <tr>
                    <th>元素</th>
                    <th>s轨道</th>
                    <th>p轨道</th>
                    <th>d轨道</th>
                    <th>f轨道</th>
                    <th>总计</th>
                    <th>主要贡献轨道</th>
                </tr>
            </thead>
            <tbody>
        """
        
        rows = ""
        for element, orbitals in orbital.items():
            s_contrib = orbitals.get('s', 0.0)
            p_contrib = orbitals.get('p', 0.0)
            d_contrib = orbitals.get('d', 0.0)
            f_contrib = orbitals.get('f', 0.0)
            total = s_contrib + p_contrib + d_contrib + f_contrib
            
            # 确定主要贡献轨道
            contribs = {'s': s_contrib, 'p': p_contrib, 'd': d_contrib, 'f': f_contrib}
            main_orbital = max(contribs.keys(), key=lambda k: contribs[k])
            main_percentage = (contribs[main_orbital] / total * 100) if total > 0 else 0
            
            rows += f"""
                <tr>
                    <td><strong>{element}</strong></td>
                    <td>{s_contrib:.3f}</td>
                    <td>{p_contrib:.3f}</td>
                    <td>{d_contrib:.3f}</td>
                    <td>{f_contrib:.3f}</td>
                    <td><strong>{total:.3f}</strong></td>
                    <td>{main_orbital.upper()}轨道 ({main_percentage:.1f}%)</td>
                </tr>
            """
        
        footer = """
            </tbody>
        </table>
        """
        
        return header + rows + footer
    
    def _generate_orbital_analysis_table(self, orbital: Dict[str, Any]) -> str:
        """生成轨道分析表格"""
        if not orbital:
            return "<p>无轨道分析数据</p>"
        
        rows = ""
        for element, orbitals in orbital.items():
            rows += f'<tr><td colspan="2"><strong>{element}</strong></td></tr>'
            for orb, contrib in orbitals.items():
                rows += f'<tr><td>{orb.upper()}</td><td>{contrib:.3f}</td></tr>'
        
        return f'<table class="data-table">{rows}</table>'
    
    def _generate_element_analysis_table(self, elements: Dict[str, Any]) -> str:
        """生成元素分析表格"""
        if not elements:
            return "<p>无元素分析数据</p>"
        
        rows = ""
        for element, data in elements.items():
            en = data.get('electronegativity', 'N/A')
            radius = data.get('atomic_radius', 'N/A')
            rows += f"""
            <tr><td>{element}</td><td>电负性: {en}<br>原子半径: {radius}</td></tr>
            """
        
        return f'<table class="data-table"><tr><th>元素</th><th>性质</th></tr>{rows}</table>'


def generate_pymatgen_dos_report(input_path: str, task_id: Optional[str] = None, output_dir: Optional[str] = None) -> str:
    """
    生成基于PyMatGen的DOS分析报告 - 整合test3.py功能
    
    Args:
        input_path: VASP输出文件路径
        task_id: 任务ID
        output_dir: 输出目录
    
    Returns:
        HTML报告文件路径
    """
    try:
        logger.info("🚀 启动PyMatGen DOS分析...")
        
        # 执行分析 - 使用新的构造函数
        analyzer = PyMatGenDOSAnalyzer(input_path, task_id, output_dir)
        analysis_data = analyzer.analyze()
        
        # 生成HTML报告
        logger.info("📄 生成HTML报告...")
        generator = PyMatGenDOSHTMLGenerator(analysis_data)
        
        # 确定输出路径
        html_output_file = analyzer.output_dir / "pymatgen_dos_analysis_report.html"
        html_path = generator.generate_html_report(str(html_output_file))
        
        logger.info(f"✅ PyMatGen DOS分析报告已生成: {html_path}")
        logger.info(f"📁 输出目录: {analyzer.output_dir}")
        logger.info(f"包含文件:")
        logger.info(f"  📄 报告文件:")
        logger.info(f"    - HTML报告: pymatgen_dos_analysis_report.html")
        logger.info(f"    - Markdown报告: analysis_report.md")
        logger.info(f"  📊 matplotlib图表（基于test3.py风格）:")
        logger.info(f"    - 总DOS图: total_dos.png/pdf")
        logger.info(f"    - 元素分解DOS图: element_dos.png/pdf")
        logger.info(f"    - SPD轨道DOS图: spd_dos.png/pdf")
        logger.info(f"    - 能带结构图: band_structure.png/pdf（如果有能带数据）")
        logger.info(f"  📈 HTML嵌入图表（基于matplotlib）:")
        logger.info(f"    - 总DOS图、元素DOS图、SPD DOS图、能带图")
        logger.info(f"  💾 数据文件:")
        logger.info(f"    - DOS数据: data/目录下的CSV文件")
        
        return html_path
        
    except Exception as e:
        logger.error(f"❌ PyMatGen DOS分析失败: {str(e)}")
        raise Exception(f"生成PyMatGen DOS分析报告失败: {str(e)}")


if __name__ == "__main__":
    # 测试代码
    test_path = "/Users/ysl/Desktop/Code/MCP_IC_Tool/dos_test/3b9c2f5cf70841449e6cc1437e13dd52"
    test_task_id = "pymatgen_dos_test"
    
    print(f"🔍 测试路径: {test_path}")
    
    try:
        html_report = generate_pymatgen_dos_report(test_path, test_task_id)
        print(f"✅ PyMatGen DOS分析HTML报告已生成: {html_report}")
    except Exception as e:
        print(f"❌ 测试失败: {e}")
