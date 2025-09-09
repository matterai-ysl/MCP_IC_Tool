
"""
VASP 分子动力学智能分析器 - 基于 PyMatGen 的MD分析与可视化
功能概览（与 `dos_analyzer.py` 风格一致）：
- 扩散性质分析：
  - MSD 计算（按元素分类）
  - 扩散系数自动计算（线性拟合尾段）
  - 离子电导率计算（Nernst–Einstein）
  - Arrhenius 图生成（多温数据自动聚合）
  - 激活能自动拟合（ln D vs 1/T）
- 径向分布函数（RDF）分析：
  - 全体/分元素 RDF 计算
  - 配位数统计（积分至首极小值）
  - 峰位识别与分析（简易峰值检测）
  - 结构演化追踪（首峰随时间演化）
  - 与实验 XRD/中子散射对比（提供覆盖接口）
- 系统稳定性监控：
  - 能量/温度/压力演化
  - 晶格参数变化追踪（a, b, c, 体积）
  - 密度波动分析
  - 平衡态识别（滚动斜率阈值）
  - 异常状态预警（简单Z分数阈值）

输入建议：包含 XDATCAR（轨迹）、OUTCAR（热力学信息）、vasprun.xml（补充）、INCAR（步长POTIM等）
输出：在工作目录下生成 MD_output/ 图表与数据，并生成 HTML 报告。
"""

import os
import io
import json
import math
import base64
import logging
import importlib
from pathlib import Path
from typing import Dict, Any, Optional, List, Tuple, Iterable

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from datetime import datetime

# PyMatGen
from pymatgen.io.vasp import Outcar
from pymatgen.io.vasp.outputs import Xdatcar, Vasprun
from pymatgen.io.vasp.inputs import Incar
from pymatgen.core import Structure, Element
from pymatgen.analysis.diffraction.xrd import XRDCalculator

from pymatgen.analysis.diffusion.analyzer import DiffusionAnalyzer

try:
    # 某些版本提供中子衍射计算器
    from pymatgen.analysis.diffraction.neutron import NDCalculator as NeutronDiffractionCalculator  # type: ignore
except Exception:  # pragma: no cover - 可选依赖
    NeutronDiffractionCalculator = None

try:
    # pymatgen的RDF分析工具
    from pymatgen.analysis.diffusion.aimd.rdf import RadialDistributionFunction
except ImportError:
    try:
        # 备用导入路径
        from pymatgen.analysis.diffusion.rdf import RadialDistributionFunction  # type: ignore
    except ImportError:
        RadialDistributionFunction = None


# 自定义的基于pymatgen核心功能的RDF分析器
class PyMatGenRDFAnalyzer:
    """基于PyMatGen核心功能的RDF分析器"""
    
    def __init__(self, structures: List[Structure], rmax: float = 10.0, nbins: int = 100):
        self.structures = structures
        self.rmax = rmax
        self.nbins = nbins
        self.dr = rmax / nbins
        self.r = np.linspace(0, rmax, nbins)
    
    def compute_rdf(self, indices_a: List[int], indices_b: List[int]) -> Tuple[np.ndarray, np.ndarray]:
        """计算指定原子索引间的RDF"""
        rdf_sum = np.zeros(self.nbins)
        
        for structure in self.structures:
            # 获取原子位置
            coords_a = np.array([structure[i].coords for i in indices_a])  # type: ignore
            coords_b = np.array([structure[i].coords for i in indices_b])  # type: ignore
            
            # 计算距离矩阵
            distances = []
            for coord_a in coords_a:
                for j, coord_b in enumerate(coords_b):
                    if indices_a != indices_b or coord_a is not coord_b:  # 避免自相关
                        # 使用pymatgen的最短距离计算（考虑周期性边界条件）
                        dist = structure.lattice.get_distance_and_image(coord_a, coord_b)[0]
                        if dist <= self.rmax:
                            distances.append(dist)
            
            # 统计直方图
            hist, _ = np.histogram(distances, bins=self.nbins, range=(0, self.rmax))
            rdf_sum += hist
        
        # 归一化
        n_a = len(indices_a)
        n_b = len(indices_b)
        volume = np.mean([s.lattice.volume for s in self.structures])
        density = n_b / volume
        
        # RDF归一化：g(r) = (V/(N_a*N_b)) * (1/(4πr²dr)) * hist
        normalization = np.zeros_like(self.r)
        for i in range(len(self.r)):
            r_val = self.r[i]
            if r_val > 0:
                shell_volume = 4 * np.pi * r_val**2 * self.dr
                normalization[i] = volume / (n_a * shell_volume * len(self.structures))
        
        # 避免除零错误
        normalization[0] = 0
        rdf_sum[0] = 0
        
        rdf = rdf_sum * normalization
        return self.r, rdf



# 图形/日志
plt.rcParams['font.sans-serif'] = ['SimHei', 'Arial Unicode MS', 'DejaVu Sans']
plt.rcParams['axes.unicode_minus'] = False
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


# 物理常数
KB_J_PER_K = 1.380649e-23  # J/K
E_CHARGE_C = 1.602176634e-19  # C
ANG3_TO_M3 = 1e-30
PS_TO_S = 1e-12


class VASP_MDAnalyzer:
    """VASP 分子动力学智能分析器

    参考 `dos_analyzer.py` 的面向对象结构与输出风格。
    """

    def __init__(
        self,
        input_path: str,
        task_id: Optional[str] = None,
        output_dir: Optional[str] = None,
        mobile_species: Optional[List[str]] = None,
        rdf_rmax: float = 10.0,
        rdf_bin_width: float = 0.1,
        rdf_stride: int = 5,
    ) -> None:
        """初始化

        Args:
            input_path: 单次MD目录，或包含多温度子目录的上级目录
            task_id: 任务ID
            output_dir: 输出目录；默认 <work_dir>/MD_output
            mobile_species: 待评估电导率的迁移离子（默认自动识别 Li/Na/K/Mg...）
            rdf_rmax: RDF 截止半径 (Å)
            rdf_bin_width: RDF 直方图 bin 宽度 (Å)
            rdf_stride: RDF 采样步长（降低计算量）
        """
        self.input_path = Path(input_path)
        self.task_id = task_id or 'md_analysis'
        self.work_dir = self.input_path
        self.output_dir = Path(output_dir) if output_dir else self.work_dir / 'MD_output'
        self.output_dir.mkdir(exist_ok=True)

        # 运行内状态
        self.mobile_species = mobile_species
        self.rdf_rmax = float(rdf_rmax)
        self.rdf_bin_width = float(rdf_bin_width)
        self.rdf_stride = int(max(1, rdf_stride))

        # 文件路径
        self.xdatcar_path = self.work_dir / 'XDATCAR'
        self.outcar_path = self.work_dir / 'OUTCAR'
        self.vasprun_path = self.work_dir / 'vasprun.xml'
        self.incar_path = self.work_dir / 'INCAR'

        # 数据容器
        self.structures: List[Structure] = []
        self.lattice_volumes: List[float] = []
        self.time_step_ps: Optional[float] = None
        self.temperature_K: Optional[float] = None
        self.temperature_series: List[float] = []
        self.pressure_series: List[float] = []
        self.energy_series: List[float] = []
        self.lattice_series: Dict[str, List[float]] = {'a': [], 'b': [], 'c': [], 'vol': []}

        self.analysis_data: Dict[str, Any] = {}

    # =============================
    # 主入口
    # =============================
    def analyze(self) -> Dict[str, Any]:
        """执行完整MD分析（单温度目录）；若目录内含子目录，则自动尝试做 Arrhenius 聚合。

        Returns:
            Dict[str, Any]: 分析数据汇总
        """
        logger.info(f"🚀 开始 VASP MD 分析: {self.work_dir}")

        try:
            # 检测是否多运行（用于Arrhenius）
            md_subdirs = self._discover_md_subdirs(self.work_dir)
            if md_subdirs:
                logger.info(f"检测到 {len(md_subdirs)} 个子运行，执行 Arrhenius 聚合分析…")
                arrh = self._analyze_arrhenius_across_subdirs(md_subdirs)
                self.analysis_data['arrhenius'] = arrh

            # 单运行常规分析
            self._load_md_data()
            if not self.structures:
                raise FileNotFoundError('未能加载MD结构数据（XDATCAR）')

            self._analyze_stability()
            diffusion = self._analyze_diffusion()
            self.analysis_data['diffusion'] = diffusion

            rdf = self._analyze_rdf()
            self.analysis_data['rdf'] = rdf

            self._generate_visualizations()
            self._generate_markdown_report()

            # HTML 报告
            html_path = self._generate_html_report()
            self.analysis_data['report_html'] = html_path

            logger.info('✅ VASP MD 分析完成')
            return self.analysis_data

        except Exception as e:
            logger.error(f"❌ MD 分析失败: {e}")
            import traceback
            traceback.print_exc()
            raise

    # =============================
    # 数据加载与通用工具
    # =============================
    def _discover_md_subdirs(self, root: Path) -> List[Path]:
        """扫描 root 下的子目录，若存在 XDATCAR 则认为是一个独立MD运行。

        仅用于 Arrhenius 聚合，不影响当前目录分析。
        """
        subdirs: List[Path] = []
        try:
            for p in root.iterdir():
                if p.is_dir() and (p / 'XDATCAR').exists():
                    subdirs.append(p)
        except Exception:
            pass
        return subdirs

    def _load_md_data(self) -> None:
        """加载 XDATCAR 轨迹与 OUTCAR/INCAR/vasprun 的辅助信息"""
        logger.info('📦 加载MD数据…')

        # 1) XDATCAR
        if not self.xdatcar_path.exists():
            raise FileNotFoundError(f'未找到 XDATCAR: {self.xdatcar_path}')
        xdatcar = Xdatcar(str(self.xdatcar_path))
        self.structures = xdatcar.structures
        logger.info(f'加载结构帧数: {len(self.structures)}')

        # 2) 采样晶格参数
        self.lattice_volumes = [float(s.lattice.volume) for s in self.structures]
        self.lattice_series['a'] = [float(s.lattice.a) for s in self.structures]
        self.lattice_series['b'] = [float(s.lattice.b) for s in self.structures]
        self.lattice_series['c'] = [float(s.lattice.c) for s in self.structures]
        self.lattice_series['vol'] = self.lattice_volumes

        # 3) INCAR: POTIM 估算时间步长（fs→ps）
        if self.incar_path.exists():
            try:
                incar = Incar.from_file(str(self.incar_path))
                potim = float(incar.get('POTIM', 1.0))  # fs
                self.time_step_fs = potim
                self.time_step_ps = potim * 1e-3
                logger.info(f'时间步长 POTIM = {potim} fs ({self.time_step_ps:.6f} ps)')
            except Exception as e:
                logger.warning(f'读取 INCAR 失败: {e}')
        # 4) vasprun.xml 补充（如果有）
        if self.vasprun_path.exists():
            try:
                v = Vasprun(str(self.vasprun_path), exception_on_bad_xml=False)
                if self.time_step_ps is None:
                    # 有些版本提供 ionic_step_time（fs）或时序，可回退估算
                    try:
                        ts_fs = float(getattr(v, 'ionic_step_time', None) or 1.0)
                        self.time_step_fs = ts_fs
                        self.time_step_ps = ts_fs * 1e-3
                    except Exception:
                        pass
                # 估计温度（若是恒温器MD，可读取 TEBEG/TEEND；此处优先 OUTCAR）
            except Exception as e:
                logger.warning(f'读取 vasprun.xml 失败: {e}')

        # 5) OUTCAR: 温度、压力、能量随时间
        if self.outcar_path.exists():
            try:
                outcar = Outcar(str(self.outcar_path))
                # 尝试模式匹配温度/压力/能量
                # 温度模式 - 尝试多种可能的格式
                temp_patterns = [
                    r"temperature\s+([0-9\.]+)\s+K",      # temperature 300.0 K
                    # r"TEBEG\s*=\s*([0-9\.]+)",            # TEBEG = 300.0 (从INCAR部分)
                    # r"T=\s*([0-9\.]+)",                   # T= 300.0
                    # r"TEMP\s*=\s*([0-9\.]+)",             # TEMP = 300.0
                ]
                
                temps = []
                for pattern in temp_patterns:
                    try_temps = self._read_outcar_pattern(outcar, pattern)
                    if try_temps:
                        temps = try_temps
                        logger.info(f"使用温度模式: {pattern}")
                        break
                
                if temps:
                    self.temperature_series = [float(t[0]) for t in temps]
                    # 平均温度
                    self.temperature_K = float(np.mean(self.temperature_series)) if self.temperature_series else None
                else:
                    self.temperature_K = None

                # 压力（kB）或 bar，VASP OUTCAR 常见行："external pressure =   xxx kB"
                pressure_patterns = [
                    r"external\s+pressure\s*=\s*(-?[\d\.\-]+)\s*kB",    # external pressure = -0.12 kB
                    r"pressure\s*=\s*(-?[\d\.\-]+)",                    # pressure = -0.12
                    r"PRESS\s*=\s*(-?[\d\.\-]+)",                       # PRESS = -0.12
                ]
                
                pressures = []
                for pattern in pressure_patterns:
                    try_pressures = self._read_outcar_pattern(outcar, pattern)
                    if try_pressures:
                        pressures = try_pressures
                        logger.info(f"使用压力模式: {pattern}")
                        break
                
                if pressures:
                    self.pressure_series = [float(p[0]) for p in pressures]

                # 能量 TOTEN（eV）- 修正正则表达式匹配实际OUTCAR格式
                print("outcar",outcar)
                # 尝试多种能量模式
                energy_patterns = [
                    # r"energy\s+without\s+entropy\s*=\s*(-?[0-9\.\-]+)",  # energy without entropy = -189.32447661
                    # r"energy\(sigma->0\)\s*=\s*(-?[0-9\.\-]+)",          # energy(sigma->0) = -189.32688170  
                    r"free\s+energy\s+TOTEN\s*=\s*(-?[0-9\.\-]+)",       # free energy TOTEN = xxx
                    # r"TOTEN\s*=\s*(-?[0-9\.\-]+)",                       # TOTEN = xxx
                ]
                
                energies = []
                for pattern in energy_patterns:
                    try_energies = self._read_outcar_pattern(outcar, pattern)

                    if try_energies:
                        energies = try_energies
                        logger.info(f"使用能量模式: {pattern}")
                        break
                
                if energies:
                    energy_values = [float(e[0]) for e in energies]
                    # 只取与结构帧数相同的能量点（每个MD步的最后一个能量）
                    n_structures = len(self.structures)
                    if len(energy_values) > n_structures:
                        # 每n个点取最后一个，或者简单地按步长采样
                        stride = len(energy_values) // n_structures
                        self.energy_series = energy_values[stride-1::stride][:n_structures]
                        logger.info(f"能量数据降采样: {len(energy_values)} -> {len(self.energy_series)} 点")
                    else:
                        self.energy_series = energy_values

                logger.info(
                    f"OUTCAR 统计: 温度点={len(self.temperature_series)}, 压力点={len(self.pressure_series)}, 能量点={len(self.energy_series)}"
                )
            except Exception as e:
                logger.warning(f'解析 OUTCAR 失败: {e}')

        # 默认时间步长兜底
        if self.time_step_ps is None:
            self.time_step_ps = 1.0  # 1 ps 兜底，避免为 None
            logger.warning('未能从 INCAR/vasprun 获取时间步长，使用默认 1.0 ps')

        # 默认温度兜底
        if self.temperature_K is None:
            # 若无系列，尝试 INCAR 的 TEBEG
            try:
                if self.incar_path.exists():
                    incar = Incar.from_file(str(self.incar_path))
                    self.temperature_K = float(incar.get('TEBEG', 300))
            except Exception:
                self.temperature_K = 300.0
            logger.warning(f'未能从 OUTCAR 提取温度，使用 {self.temperature_K} K')

        # 自动识别迁移离子
        if not self.mobile_species:
            self.mobile_species = self._auto_mobile_species(self.structures[0])
            logger.info(f"自动识别迁移离子: {', '.join(self.mobile_species) if self.mobile_species else '无'}")

    def _read_outcar_pattern(self, outcar: Outcar, regex: str) -> List[List[str]]:
        """直接从OUTCAR文件读取并使用正则表达式匹配"""
        import re
        try:
            with open(self.outcar_path, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read()
            
            pattern = re.compile(regex)
            matches = pattern.findall(content)
            # 将匹配结果转换为List[List[str]]格式
            if matches:
                return [[match] if isinstance(match, str) else list(match) for match in matches]
            return []
        except Exception as e:
            logger.warning(f"读取OUTCAR模式匹配失败: {e}")
            return []

    def _auto_mobile_species(self, structure: Structure) -> List[str]:
        """根据组成自动推断可能的迁移离子（优先 Li/Na/K/Mg/Ag/H）"""
        candidates = ['Li', 'Na', 'K', 'Mg', 'Ag', 'H']
        present = set([str(el) for el in structure.composition.elements])
        return [s for s in candidates if s in present]

    # =============================
    # 扩散性质分析
    # =============================
    def _analyze_diffusion(self) -> Dict[str, Any]:
        logger.info('📈 扩散性质分析（MSD/扩散系数/电导率/Arrhenius）…')
        # 兜底确保非 None
        time_step_ps_val = float(self.time_step_ps if self.time_step_ps is not None else 1.0)

        time_step_s = time_step_ps_val * PS_TO_S

        # MSD 按元素
        msd_by_elem = self._compute_msd_by_element(time_step_s=time_step_s)

        # 使用 DiffusionAnalyzer 直接计算扩散系数和电导率
        diffusion_by_elem: Dict[str, Dict[str, float]] = {}
        conductivity: Dict[str, float] = {}
        
        structures = self.structures
        species_labels = [site.species_string for site in structures[0]]
        unique_elements = sorted(set(species_labels))
        print("time_step_fs",self.time_step_fs)

        for elem in unique_elements:
            try:
                tempK = float(self.temperature_K or 300.0)
                da = DiffusionAnalyzer.from_structures( # type: ignore
                    structures,
                    specie=elem,
                    temperature=tempK,
                    time_step=self.time_step_fs,
                    step_skip=1,
                    smoothed=False,
                )
                
                # 直接从 DiffusionAnalyzer 获取扩散系数
                D_cm2_per_s = da.diffusivity  # cm²/s
                D_m2_per_s = float(D_cm2_per_s * 1e-4)  # 转换为 m²/s
                # 获取拟合信息（如果可用）
                fit_slope = getattr(da, 'slope', 0.0)
                fit_intercept = getattr(da, 'intercept', 0.0) 
                fit_r2 = getattr(da, 'r2', 0.0)
                
                diffusion_by_elem[elem] = {
                    'D_m2_per_s': D_m2_per_s,
                    'D_cm2_per_s': float(D_cm2_per_s),
                    'fit_slope': float(fit_slope),
                    'fit_intercept': float(fit_intercept),
                    'fit_r2': float(fit_r2)
                }
                
                # 计算电导率（仅对迁移离子）
                if self.mobile_species and elem in self.mobile_species:
                    try:
                        # 使用 DiffusionAnalyzer 的内置电导率计算
                        if hasattr(da, 'conductivity'):
                            sigma_S_per_cm = da.conductivity  # S/cm
                            conductivity[elem] = float(sigma_S_per_cm * 100)  # 转换为 S/m
                            print("sigma_S_per_cm",sigma_S_per_cm)
                        else:
                            # 手工计算 Nernst-Einstein 电导率
                            avg_vol_m3 = float(np.mean(self.lattice_volumes)) * ANG3_TO_M3
                            comp = self.structures[0].composition
                            n_species = comp[Element(elem)]  # 原胞中原子数
                            number_density = (n_species / avg_vol_m3)  # 1/m^3
                            z = self._guess_charge_number(elem)
                            sigma = number_density * (z * E_CHARGE_C) ** 2 * D_m2_per_s / (KB_J_PER_K * tempK)
                            conductivity[elem] = float(sigma)
                    except Exception as e:
                        logger.warning(f"计算元素 {elem} 电导率失败: {e}")
                        
            except Exception as e:
                logger.warning(f"DiffusionAnalyzer 分析元素 {elem} 失败: {e}")
                # 回退到手工MSD方法
                if elem in msd_by_elem:
                    times = msd_by_elem[elem]['time_s']
                    values = msd_by_elem[elem]['msd']
                    slope, intercept, r2 = self._fit_linear_tail(times, values)
                    # Einstein relation: <r^2> = 2*d*D*t（3D）
                    D = max(slope / (2 * 3), 0.0)
                    diffusion_by_elem[elem] = {
                        'D_m2_per_s': float(D),
                        'D_cm2_per_s': float(D * 1e4),
                        'fit_slope': float(slope),
                        'fit_intercept': float(intercept),
                        'fit_r2': float(r2)
                    }
                    print(elem,diffusion_by_elem[elem])
        # Arrhenius（若单目录，仅返回当前点；多目录分析由 _analyze_arrhenius_across_subdirs 提供）
        arrhenius = {
            'T_list_K': [float(self.temperature_K or 300.0)],
            'D_list_m2_per_s': [float(np.mean([v['D_m2_per_s'] for v in diffusion_by_elem.values()]) if diffusion_by_elem else 0.0)]
        }

        # 导出 CSV/JSON
        self._save_msd_results(msd_by_elem, diffusion_by_elem, conductivity)

        return {
            'msd_by_element': msd_by_elem,
            'diffusion_by_element': diffusion_by_elem,
            'ionic_conductivity_S_per_m': conductivity,
            'arrhenius_single': arrhenius
        }

    def _compute_msd_by_element(self, time_step_s: float) -> Dict[str, Dict[str, np.ndarray]]:
        """按元素计算 MSD（优先使用 PyMatGen DiffusionAnalyzer；否则手工PBC解缠+累积）。"""
        structures = self.structures
        if len(structures) < 2:
            return {}

        species_labels = [site.species_string for site in structures[0]]
        unique_elements = sorted(set(species_labels))

        # 优先使用 DiffusionAnalyzer（如可用）
        if DiffusionAnalyzer is not None:
            try:
                tempK = float(self.temperature_K or 300.0)
                msd_by_elem_da: Dict[str, Dict[str, np.ndarray]] = {}
                for elem in unique_elements:
                    try:
                        print("DiffusionAnalyzer is available")
                        print("elem", elem)
                        da = DiffusionAnalyzer.from_structures(
                            structures,
                            specie=elem,
                            temperature=tempK,
                            time_step=time_step_s,
                            step_skip=1,
                            smoothed=False,
                        )
                        times_arr = None
                        msd_arr = None
                        if hasattr(da, 'msd') and getattr(da, 'msd') is not None:
                            msd_raw = np.array(getattr(da, 'msd'))
                            if msd_raw.ndim == 1:
                                msd_arr = msd_raw
                            else:
                                # 若为 (n, 3/4) 等多列，取各向平均
                                msd_arr = np.mean(msd_raw, axis=1)
                            dt_val = float(getattr(da, 'dt', time_step_s))
                            times_arr = np.arange(len(msd_arr), dtype=float) * dt_val
                        if times_arr is not None and msd_arr is not None:
                            msd_by_elem_da[elem] = {'time_s': np.array(times_arr), 'msd': np.array(msd_arr)}
                    except Exception:
                        continue
                if msd_by_elem_da:
                    return msd_by_elem_da
            except Exception as e:
                logger.warning(f'DiffusionAnalyzer 不可用或发生错误，使用手工MSD: {e}')

        # 手工：基于PBC最短像解缠，累计位移
        frac_0 = np.array([site.frac_coords for site in structures[0]])
        num_atoms = frac_0.shape[0]
        cum_disp_cart = np.zeros((num_atoms, 3), dtype=float)

        # 时间序列
        times = np.arange(len(structures), dtype=float) * time_step_s
        msd_all_frames = []  # (frame, per-atom |r|^2)

        prev_frac = frac_0
        for k in range(1, len(structures)):
            s_prev = structures[k - 1]
            s_curr = structures[k]
            frac_curr = np.array([site.frac_coords for site in s_curr])
            df = frac_curr - prev_frac
            df -= np.round(df)  # unwrap to [-0.5, 0.5]
            # 使用上一帧晶格将分数位移变换到笛卡尔
            dcart = s_prev.lattice.get_cartesian_coords(df)
            cum_disp_cart += dcart
            # 相对初态位移的平方
            msd_all_frames.append(np.sum(cum_disp_cart ** 2, axis=1))
            prev_frac = frac_curr

        # 汇总到元素
        msd_by_elem: Dict[str, Dict[str, np.ndarray]] = {}
        msd_all_frames_arr = np.array(msd_all_frames)  # shape (n-1, N)
        times_eff = times[1:]

        for elem in unique_elements:
            mask = np.array([sp == elem for sp in species_labels])
            if not np.any(mask):
                continue
            msd_elem = np.mean(msd_all_frames_arr[:, mask], axis=1)
            msd_by_elem[elem] = {
                'time_s': times_eff,
                'msd': msd_elem
            }

        return msd_by_elem

    def _fit_linear_tail(self, x: np.ndarray, y: np.ndarray, start_frac: float = 0.2, end_frac: float = 0.8) -> Tuple[float, float, float]:
        """对尾段进行线性拟合，返回 slope/intercept/R^2。"""
        if len(x) < 5 or len(y) < 5:
            return 0.0, 0.0, 0.0
        n = len(x)
        i0 = int(max(0, math.floor(n * start_frac)))
        i1 = int(min(n, math.ceil(n * end_frac)))
        if i1 - i0 < 3:
            i0 = 0
            i1 = n
        xx = x[i0:i1]
        yy = y[i0:i1]
        A = np.vstack([xx, np.ones_like(xx)]).T
        slope, intercept = np.linalg.lstsq(A, yy, rcond=None)[0]
        # R^2
        yhat = slope * xx + intercept
        ss_res = float(np.sum((yy - yhat) ** 2))
        ss_tot = float(np.sum((yy - np.mean(yy)) ** 2)) + 1e-20
        r2 = 1.0 - ss_res / ss_tot
        return float(slope), float(intercept), float(r2)

    def _guess_charge_number(self, elem: str) -> int:
        """基于常见价态猜测电荷数 |z|（仅用于估算电导率）"""
        mapping = {
            'Li': 1, 'Na': 1, 'K': 1, 'Ag': 1,
            'Mg': 2, 'Ca': 2, 'Al': 3,
            'H': 1
        }
        return mapping.get(elem, 1)

    def _save_msd_results(
        self,
        msd_by_elem: Dict[str, Dict[str, np.ndarray]],
        diffusion_by_elem: Dict[str, Dict[str, float]],
        conductivity: Dict[str, float],
    ) -> None:
        data_dir = self.output_dir / 'data'
        data_dir.mkdir(exist_ok=True)
        # MSD CSV
        for elem, data in msd_by_elem.items():
            df = pd.DataFrame({'time_s': data['time_s'], 'msd_m2': data['msd']})
            df.to_csv(data_dir / f'MSD_{elem}.csv', index=False)
        # 扩散 & 电导率 JSON
        results = {
            'diffusion_by_element': diffusion_by_elem,
            'ionic_conductivity_S_per_m': conductivity,
            'temperature_K': float(self.temperature_K or 300.0),
            'mobile_species': self.mobile_species or [],
            'analysis_method': 'DiffusionAnalyzer + Manual_fallback',
            'units': {
                'diffusivity': 'm²/s and cm²/s',
                'conductivity': 'S/m',
                'temperature': 'K'
            }
        }
        with open(data_dir / 'diffusion_results.json', 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, ensure_ascii=False)

    # =============================
    # RDF 分析
    # =============================
    def _analyze_rdf(self) -> Dict[str, Any]:
        logger.info('🔍 RDF/配位/峰位/结构演化分析（基于PyMatGen）…')
        if not self.structures:
            return {}

        # 优先使用pymatgen的RadialDistributionFunction
        if RadialDistributionFunction is not None:
            logger.info('使用 PyMatGen RadialDistributionFunction 进行RDF分析')
            return self._analyze_rdf_pymatgen()
        else:
            logger.info('使用自定义 PyMatGen RDF 分析器（基于pymatgen核心功能）')
            return self._analyze_rdf_pymatgen_custom()

    def _analyze_rdf_pymatgen(self) -> Dict[str, Any]:
        """使用pymatgen的RadialDistributionFunction进行RDF分析"""
        structures = self.structures[:: self.rdf_stride]
        ngrid = int(self.rdf_rmax / self.rdf_bin_width)
        
        # 获取所有元素
        species = [str(el) for el in structures[0].composition.elements]
        logger.info(f'检测到的元素: {species}')
        
        # 生成所有元素对
        pairs: List[Tuple[str, str]] = []
        for i, a in enumerate(species):
            for j, b in enumerate(species):
                if j < i:
                    continue
                pairs.append((a, b))
        
        g_pairs: Dict[str, Dict[str, Any]] = {}
        r_array = None
        
        for a, b in pairs:
            try:
                # 获取原子索引
                indices_a = self._get_species_indices(structures[0], a)
                indices_b = self._get_species_indices(structures[0], b)
                
                if not indices_a or not indices_b:
                    logger.warning(f'元素对 {a}-{b} 索引为空，跳过')
                    continue
                
                # 使用pymatgen的RadialDistributionFunction
                rdf_analyzer = RadialDistributionFunction(  # type: ignore
                    structures=structures,
                    indices=indices_a,
                    reference_indices=indices_b,
                    rmax=self.rdf_rmax,
                    ngrid=ngrid,
                    sigma=0.1  # 高斯平滑参数
                )
                
                # 提取结果 - 使用正确的属性名
                r_array = rdf_analyzer.interval  # r轴数组
                g_rdf = rdf_analyzer.rdf  # g(r)数组
                
                # 分析峰位和配位数
                peaks, mins, cn = self._analyze_rdf_features(r_array, g_rdf, structures[0], a, b)
                
                g_pairs[f'{a}-{b}'] = {
                    'g_r': g_rdf,
                    'peaks': peaks,
                    'mins': mins,
                    'coordination_number': cn,
                    'method': 'pymatgen_RadialDistributionFunction'
                }
                
                logger.info(f'✅ {a}-{b} RDF 分析完成 (CN={cn:.2f})')
                
            except Exception as e:
                logger.warning(f'元素对 {a}-{b} RDF 分析失败: {e}')
                continue
        
        # 总RDF（所有原子对的平均）
        if r_array is not None and g_pairs:
            g_all = np.mean([val['g_r'] for val in g_pairs.values()], axis=0)
        else:
            # 回退到手工计算
            logger.warning('PyMatGen RDF分析失败，回退到手工计算')
            return self._analyze_rdf_manual()
        
        # 结构演化分析
        evolution = self._track_first_peak_evolution_pymatgen(structures, species[0] if species else 'all')
        
        # 保存数据
        self._save_rdf_results(r_array, g_all, g_pairs)
        
        return {
            'r_A': r_array,
            'g_all': g_all,
            'g_pairs': g_pairs,
            'evolution': evolution,
            'method': 'pymatgen'
        }

    def _analyze_rdf_pymatgen_custom(self) -> Dict[str, Any]:
        """使用自定义PyMatGen RDF分析器"""
        structures = self.structures[:: self.rdf_stride]
        nbins = int(self.rdf_rmax / self.rdf_bin_width)
        
        # 创建自定义RDF分析器
        rdf_analyzer = PyMatGenRDFAnalyzer(structures, rmax=self.rdf_rmax, nbins=nbins)
        
        # 获取所有元素
        species = [str(el) for el in structures[0].composition.elements]
        logger.info(f'检测到的元素: {species}')
        
        # 生成所有元素对
        pairs: List[Tuple[str, str]] = []
        for i, a in enumerate(species):
            for j, b in enumerate(species):
                if j < i:
                    continue
                pairs.append((a, b))
        
        g_pairs: Dict[str, Dict[str, Any]] = {}
        r_array = None
        
        for a, b in pairs:
            try:
                # 获取原子索引
                indices_a = self._get_species_indices(structures[0], a)
                indices_b = self._get_species_indices(structures[0], b)
                
                if not indices_a or not indices_b:
                    logger.warning(f'元素对 {a}-{b} 索引为空，跳过')
                    continue
                
                # 使用自定义RDF分析器
                r_array, g_rdf = rdf_analyzer.compute_rdf(indices_a, indices_b)
                
                # 分析峰位和配位数
                peaks, mins, cn = self._analyze_rdf_features(r_array, g_rdf, structures[0], a, b)
                
                g_pairs[f'{a}-{b}'] = {
                    'g_r': g_rdf,
                    'peaks': peaks,
                    'mins': mins,
                    'coordination_number': cn,
                    'method': 'pymatgen_custom_rdf'
                }
                
                logger.info(f'✅ {a}-{b} RDF 分析完成 (CN={cn:.2f})')
                
            except Exception as e:
                logger.warning(f'元素对 {a}-{b} RDF 分析失败: {e}')
                continue
        
        # 总RDF（所有原子对的平均）
        if r_array is not None and g_pairs:
            g_all = np.mean([val['g_r'] for val in g_pairs.values()], axis=0)
        else:
            # 回退到手工计算
            logger.warning('自定义 PyMatGen RDF分析失败，回退到手工计算')
            return self._analyze_rdf_manual()
        
        # 结构演化分析
        evolution = self._track_first_peak_evolution_custom(structures, species[0] if species else 'all', rdf_analyzer)
        
        # 保存数据
        self._save_rdf_results(r_array, g_all, g_pairs)
        
        return {
            'r_A': r_array,
            'g_all': g_all,
            'g_pairs': g_pairs,
            'evolution': evolution,
            'method': 'pymatgen_custom'
        }

    def _track_first_peak_evolution_custom(self, structures: List[Structure], reference_species: str, 
                                         rdf_analyzer: PyMatGenRDFAnalyzer) -> Dict[str, Any]:
        """使用自定义RDF分析器追踪首峰演化"""
        if not structures:
            return {}
        
        first_peak_rs: List[float] = []
        
        # 获取参考元素的索引
        ref_indices = self._get_species_indices(structures[0], reference_species)
        all_indices = list(range(len(structures[0])))
        
        if not ref_indices:
            return {}
        
        # 创建单帧RDF分析器
        for i in range(0, len(structures), max(1, len(structures) // 50)):  # 采样50个点
            try:
                single_rdf = PyMatGenRDFAnalyzer([structures[i]], rmax=self.rdf_rmax, nbins=rdf_analyzer.nbins)
                r, g = single_rdf.compute_rdf(ref_indices, all_indices)
                
                # 找到首峰
                if len(g) > 1:
                    max_idx = np.argmax(g[1:]) + 1  # 跳过r=0点
                    first_peak_rs.append(float(r[max_idx]))
                else:
                    first_peak_rs.append(0.0)
                    
            except Exception:
                # 如果单帧分析失败
                if first_peak_rs:
                    first_peak_rs.append(first_peak_rs[-1])  # 复制上一个值
                else:
                    first_peak_rs.append(0.0)
        
        return {
            'first_peak_r_series_A': first_peak_rs,
            'stride': max(1, len(structures) // 50),
            'reference_species': reference_species
        }

    def _analyze_rdf_manual(self) -> Dict[str, Any]:
        """手工RDF计算（原有实现）"""
        structures = self.structures[:: self.rdf_stride]
        r_edges = np.arange(0.0, self.rdf_rmax + self.rdf_bin_width, self.rdf_bin_width)
        r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])

        # 所有原子 RDF
        g_all = self._compute_rdf_all(structures, r_edges)

        # 按元素对（αβ）
        species = [str(el) for el in structures[0].composition.elements]
        pairs: List[Tuple[str, str]] = []
        for i, a in enumerate(species):
            for j, b in enumerate(species):
                if j < i:
                    continue
                pairs.append((a, b))

        g_pairs: Dict[str, Dict[str, Any]] = {}
        for a, b in pairs:
            g = self._compute_rdf_pair(structures, r_edges, a, b)
            peaks, mins, cn = self._analyze_rdf_features(r_centers, g, structures[0], a, b)
            g_pairs[f'{a}-{b}'] = {
                'g_r': g,
                'peaks': peaks,
                'mins': mins,
                'coordination_number': cn,
                'method': 'manual_calculation'
            }

        # 结构演化：首峰位置随时间
        evolution = self._track_first_peak_evolution(structures)

        # 保存数据
        self._save_rdf_results(r_centers, g_all, g_pairs)

        return {
            'r_A': r_centers,
            'g_all': g_all,
            'g_pairs': g_pairs,
            'evolution': evolution,
            'method': 'manual'
        }

    def _get_species_indices(self, structure: Structure, species: str) -> List[int]:
        """获取指定元素的原子索引"""
        indices = []
        for i, site in enumerate(structure):
            if str(site.specie) == species:
                indices.append(i)
        return indices

    def _track_first_peak_evolution_pymatgen(self, structures: List[Structure], reference_species: str) -> Dict[str, Any]:
        """使用pymatgen追踪首峰演化"""
        if not structures or RadialDistributionFunction is None:
            return {}
        
        first_peak_rs: List[float] = []
        ngrid = int(self.rdf_rmax / self.rdf_bin_width)
        
        # 获取参考元素的索引
        ref_indices = self._get_species_indices(structures[0], reference_species)
        if not ref_indices:
            return {}
        
        for i in range(0, len(structures), max(1, len(structures) // 50)):  # 采样50个点
            try:
                # 单帧RDF分析
                rdf_analyzer = RadialDistributionFunction(  # type: ignore
                    structures=[structures[i]],
                    indices=ref_indices,
                    reference_indices=list(range(len(structures[i]))),  # 与所有原子的RDF
                    rmax=self.rdf_rmax,
                    ngrid=ngrid,
                    sigma=0.1
                )
                
                # 找到首峰 - 使用正确的属性名
                r = rdf_analyzer.interval  # r轴数组
                g = rdf_analyzer.rdf  # g(r)数组
                
                # 简单峰检测
                max_idx = np.argmax(g[1:]) + 1  # 跳过r=0点
                first_peak_rs.append(float(r[max_idx]))
                
            except Exception:
                # 如果单帧分析失败，使用手工方法
                if first_peak_rs:
                    first_peak_rs.append(first_peak_rs[-1])  # 复制上一个值
                else:
                    first_peak_rs.append(0.0)
        
        return {
            'first_peak_r_series_A': first_peak_rs,
            'stride': max(1, len(structures) // 50),
            'reference_species': reference_species
        }

    def _save_rdf_results(self, r_array: np.ndarray, g_all: np.ndarray, g_pairs: Dict[str, Dict[str, Any]]) -> None:
        """保存RDF结果"""
        data_dir = self.output_dir / 'data'
        data_dir.mkdir(exist_ok=True)
        
        # 保存总RDF
        pd.DataFrame({'r_A': r_array, 'g_all': g_all}).to_csv(data_dir / 'RDF_all.csv', index=False)
        
        # 保存各元素对RDF
        for key, val in g_pairs.items():
            df = pd.DataFrame({'r_A': r_array, 'g': val['g_r']})
            df.to_csv(data_dir / f'RDF_{key}.csv', index=False)
        
        # 保存RDF分析汇总
        rdf_summary = {
            'analysis_method': g_pairs[list(g_pairs.keys())[0]]['method'] if g_pairs else 'unknown',
            'rmax_A': float(self.rdf_rmax),
            'bin_width_A': float(self.rdf_bin_width),
            'stride': self.rdf_stride,
            'pairs_analyzed': list(g_pairs.keys()),
            'coordination_numbers': {k: v['coordination_number'] for k, v in g_pairs.items()},
            'peak_positions': {k: [p['r'] for p in v['peaks'][:3]] for k, v in g_pairs.items()}  # 前3个峰
        }
        
        with open(data_dir / 'rdf_summary.json', 'w', encoding='utf-8') as f:
            json.dump(rdf_summary, f, indent=2, ensure_ascii=False)

    def _compute_rdf_all(self, structures: List[Structure], r_edges: np.ndarray) -> np.ndarray:
        counts = np.zeros(len(r_edges) - 1, dtype=float)
        n_frames = len(structures)
        s0 = structures[0]
        N = len(s0)
        V = float(np.mean([s.lattice.volume for s in structures]))
        rho = N / V

        for s in structures:
            # 利用距离矩阵（最短像）
            dm = s.distance_matrix
            # 只取上三角 i<j
            iu = np.triu_indices(N, k=1)
            dists = dm[iu]
            hist, _ = np.histogram(dists, bins=r_edges)
            counts += hist

        # 归一化: g(r) = V / (N(N-1)) / (4π r^2 dr) * counts
        r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
        shell_vol = 4.0 * math.pi * (r_centers ** 2) * np.diff(r_edges)
        norm = (V / (N * (N - 1))) / shell_vol
        g = counts * norm / n_frames
        return g

    def _compute_rdf_pair(self, structures: List[Structure], r_edges: np.ndarray, a: str, b: str) -> np.ndarray:
        counts = np.zeros(len(r_edges) - 1, dtype=float)
        n_frames = len(structures)
        s0 = structures[0]
        species_labels = [site.species_string for site in s0]
        idx_a = [i for i, sp in enumerate(species_labels) if sp == a]
        idx_b = [i for i, sp in enumerate(species_labels) if sp == b]
        N_a = len(idx_a)
        N_b = len(idx_b)
        if N_a == 0 or N_b == 0:
            return np.zeros(len(r_edges) - 1)

        V = float(np.mean([s.lattice.volume for s in structures]))

        for s in structures:
            dm = s.distance_matrix
            # 取 a-b 对（若 a==b，避免重复/自对）
            if a == b:
                iu = np.triu_indices(N_a, k=1)
                dists = dm[np.ix_(idx_a, idx_b)][iu]
            else:
                dists = dm[np.ix_(idx_a, idx_b)].ravel()
            hist, _ = np.histogram(dists, bins=r_edges)
            counts += hist

        r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
        shell_vol = 4.0 * math.pi * (r_centers ** 2) * np.diff(r_edges)
        # g_ab(r) = V/(N_a N_b) / shell * counts；若 a==b 用 N_a (N_a-1)
        denom = (N_a * (N_b - 1)) if a == b else (N_a * N_b)
        denom = max(denom, 1)
        norm = (V / denom) / shell_vol
        g = counts * norm / n_frames
        return g

    def _analyze_rdf_features(
        self,
        r_centers: np.ndarray,
        g: np.ndarray,
        structure: Structure,
        a: Optional[str] = None,
        b: Optional[str] = None,
    ) -> Tuple[List[Dict[str, float]], List[Dict[str, float]], float]:
        """峰位/极小值/配位数（至首极小值）"""
        if len(g) < 5:
            return [], [], 0.0

        # 平滑（简单移动平均）
        k = 3
        g_pad = np.pad(g, (k, k), mode='edge')
        g_smooth = np.convolve(g_pad, np.ones(2 * k + 1) / (2 * k + 1), mode='valid')

        peaks = []
        mins = []
        for i in range(1, len(g_smooth) - 1):
            if g_smooth[i] > g_smooth[i - 1] and g_smooth[i] > g_smooth[i + 1] and g_smooth[i] > 1.05:
                peaks.append({'r': float(r_centers[i]), 'g': float(g_smooth[i])})
            if g_smooth[i] < g_smooth[i - 1] and g_smooth[i] < g_smooth[i + 1]:
                mins.append({'r': float(r_centers[i]), 'g': float(g_smooth[i])})

        # 找首峰与其后的首极小
        cn = 0.0
        if peaks:
            first_peak_r = peaks[0]['r']
            r_min = None
            for m in mins:
                if m['r'] > first_peak_r:
                    r_min = m['r']
                    break
            if r_min is None:
                r_min = float(r_centers[-1])
            # 配位数 CN = 4π ∫ g_ab(r) r^2 ρ_b dr
            comp = structure.composition
            V = float(structure.lattice.volume)
            rho_b = 0.0
            if b is not None:
                try:
                    N_b = comp[Element(b)]
                    rho_b = N_b / V
                except Exception:
                    rho_b = len(structure) / V
            else:
                rho_b = len(structure) / V

            mask = r_centers <= r_min
            dr = np.gradient(r_centers)
            cn = float(np.trapz(4.0 * math.pi * rho_b * g[mask] * (r_centers[mask] ** 2), x=r_centers[mask]))

        return peaks[:5], mins[:5], float(cn)

    def _track_first_peak_evolution(self, structures: List[Structure]) -> Dict[str, Any]:
        r_edges = np.arange(0.0, self.rdf_rmax + self.rdf_bin_width, self.rdf_bin_width)
        r_centers = 0.5 * (r_edges[:-1] + r_edges[1:])
        first_peak_rs: List[float] = []
        for s in structures:
            g = self._compute_rdf_all([s], r_edges)  # 单帧近似
            # 简易峰检测
            idx = np.argmax(g)
            first_peak_rs.append(float(r_centers[idx]))
        return {
            'first_peak_r_series_A': first_peak_rs,
            'stride': self.rdf_stride
        }

    # =============================
    # 稳定性监控
    # =============================
    def _analyze_stability(self) -> None:
        logger.info('🛡️ 系统稳定性监控…')
        time_step_ps_val = float(self.time_step_ps if self.time_step_ps is not None else 1.0)
        time_step_s = time_step_ps_val * PS_TO_S
        n_frames = len(self.structures)
        times = np.arange(n_frames, dtype=float) * time_step_s

        # 密度（kg/m^3）
        mass_list: List[float] = []
        for site in self.structures[0]:
            am = getattr(site.specie, 'atomic_mass', None)
            try:
                mass_list.append(float(am) if am is not None else 0.0)
            except Exception:
                mass_list.append(0.0)
        mass_amu = float(sum(mass_list))
        # 1 amu = 1.66053906660e-27 kg
        mass_kg = mass_amu * 1.66053906660e-27
        vols_m3 = np.array(self.lattice_volumes) * ANG3_TO_M3
        density_series = mass_kg / vols_m3

        # 平衡识别：能量/温度滚动斜率阈值
        eq_idx_energy = self._detect_equilibrium_index(times[: len(self.energy_series)], np.array(self.energy_series)) if self.energy_series else None
        eq_idx_temp = self._detect_equilibrium_index(times[: len(self.temperature_series)], np.array(self.temperature_series)) if self.temperature_series else None

        # 异常预警：Z分数
        anomalies = {
            'energy': self._detect_anomaly(np.array(self.energy_series)),
            'temperature': self._detect_anomaly(np.array(self.temperature_series)),
            'pressure': self._detect_anomaly(np.array(self.pressure_series)),
        }

        self.analysis_data['stability'] = {
            'times_s': times.tolist(),
            'energy_eV': self.energy_series,
            'temperature_K': self.temperature_series,
            'pressure_kB': self.pressure_series,
            'lattice': self.lattice_series,
            'density_kg_per_m3': density_series.tolist(),
            'equilibrium_index_energy': int(eq_idx_energy) if eq_idx_energy is not None else None,
            'equilibrium_index_temperature': int(eq_idx_temp) if eq_idx_temp is not None else None,
            'anomalies': anomalies
        }

        # 保存时间序列 CSV
        data_dir = self.output_dir / 'data'
        data_dir.mkdir(exist_ok=True)
        df = pd.DataFrame({
            'time_s': times,
            'energy_eV': self.energy_series[: len(times)] if self.energy_series else [np.nan] * len(times),
            'temperature_K': self.temperature_series[: len(times)] if self.temperature_series else [np.nan] * len(times),
            'pressure_kB': self.pressure_series[: len(times)] if self.pressure_series else [np.nan] * len(times),
            'a_A': self.lattice_series['a'][: len(times)],
            'b_A': self.lattice_series['b'][: len(times)],
            'c_A': self.lattice_series['c'][: len(times)],
            'volume_A3': self.lattice_series['vol'][: len(times)],
            'density_kg_m3': density_series[: len(times)]
        })
        df.to_csv(data_dir / 'time_series.csv', index=False)

    def _detect_equilibrium_index(self, t: np.ndarray, y: np.ndarray, window: int = 50, slope_thr: float = 1e-4) -> Optional[int]:
        if len(t) < window + 2:
            return None
        # 简易：移动窗口线性拟合斜率绝对值小于阈值即视为平衡
        for i in range(len(t) - window):
            xx = t[i:i + window]
            yy = y[i:i + window]
            A = np.vstack([xx, np.ones_like(xx)]).T
            slope, _ = np.linalg.lstsq(A, yy, rcond=None)[0]
            if abs(slope) < slope_thr:
                return i + window
        return None

    def _detect_anomaly(self, y: np.ndarray, z_thr: float = 3.0) -> List[int]:
        if y is None or len(y) == 0:
            return []
        m = float(np.nanmean(y))
        s = float(np.nanstd(y) + 1e-12)
        z = np.abs((y - m) / s)
        return [int(i) for i in np.where(z > z_thr)[0]]

    # =============================
    # Arrhenius 聚合（多温度子目录）
    # =============================
    def _analyze_arrhenius_across_subdirs(self, subdirs: List[Path]) -> Dict[str, Any]:
        T_list: List[float] = []
        D_list: List[float] = []

        for p in subdirs:
            try:
                sub = VASP_MDAnalyzer(str(p), mobile_species=self.mobile_species, rdf_rmax=self.rdf_rmax, rdf_bin_width=self.rdf_bin_width, rdf_stride=self.rdf_stride)
                sub._load_md_data()
                diff = sub._analyze_diffusion()
                T = float(sub.temperature_K or 300.0)
                D_mean = float(np.mean([v['D_m2_per_s'] for v in diff['diffusion_by_element'].values()])) if diff['diffusion_by_element'] else 0.0
                if D_mean > 0:
                    T_list.append(T)
                    D_list.append(D_mean)
            except Exception as e:
                logger.warning(f'子目录 {p} Arrhenius 数据收集失败: {e}')

        Ea_eV = None
        slope = None
        intercept = None
        r2 = None
        if len(T_list) >= 2:
            x = 1.0 / np.array(T_list)  # 1/K
            y = np.log(np.array(D_list))
            A = np.vstack([x, np.ones_like(x)]).T
            slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
            yhat = slope * x + intercept
            ss_res = float(np.sum((y - yhat) ** 2))
            ss_tot = float(np.sum((y - np.mean(y)) ** 2)) + 1e-20
            r2 = 1.0 - ss_res / ss_tot
            # D = D0 * exp(-Ea/(kB T)) => ln D = ln D0 - Ea/kB * 1/T
            Ea_J = -slope * KB_J_PER_K
            Ea_eV = float(Ea_J / E_CHARGE_C)

        # 保存图（若有）
        if len(T_list) >= 2:
            fig, ax = plt.subplots(figsize=(8, 6))
            x = 1.0 / np.array(T_list)
            y = np.log(np.array(D_list))
            ax.scatter(x, y, c='k', label='数据点')
            if slope is not None and intercept is not None:
                xx = np.linspace(x.min(), x.max(), 100)
                ax.plot(xx, slope * xx + intercept, 'r-', label=f'拟合: ln D = a (1/T) + b\nEa={Ea_eV:.3f} eV, R^2={r2:.3f}')
            ax.set_xlabel('1/T (1/K)')
            ax.set_ylabel('ln D (m^2/s)')
            ax.set_title('Arrhenius 图')
            ax.grid(True, alpha=0.3)
            ax.legend()
            plt.tight_layout()
            plt.savefig(self.output_dir / 'arrhenius.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

        return {
            'T_list_K': T_list,
            'D_list_m2_per_s': D_list,
            'fit': {
                'slope': float(slope) if slope is not None else None,
                'intercept': float(intercept) if intercept is not None else None,
                'R2': float(r2) if r2 is not None else None,
                'Ea_eV': Ea_eV
            }
        }

    # =============================
    # 可视化
    # =============================
    def _generate_visualizations(self) -> None:
        logger.info('🖼️ 生成可视化图…')
        self._plot_msd()
        self._plot_stability()
        self._plot_rdf()

    def _plot_msd(self) -> None:
        diff = self.analysis_data.get('diffusion', {})
        msd_by_elem = diff.get('msd_by_element', {})
        diffusion_by_elem = diff.get('diffusion_by_element', {})
        if not msd_by_elem:
            return
        fig, ax = plt.subplots(figsize=(10, 6))
        for elem, data in msd_by_elem.items():
            t = np.array(data['time_s'])
            y = np.array(data['msd'])
            ax.plot(t, y, linewidth=2, label=f'{elem}')
            if elem in diffusion_by_elem:
                slope = diffusion_by_elem[elem]['fit_slope']
                intercept = diffusion_by_elem[elem]['fit_intercept']
                ax.plot(t, slope * t + intercept, linestyle='--', alpha=0.6)
        ax.set_xlabel('时间 (s)')
        ax.set_ylabel('MSD (m^2)')
        ax.set_title('按元素 MSD 及线性拟合')
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.savefig(self.output_dir / 'msd.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

    def _plot_stability(self) -> None:
        stab = self.analysis_data.get('stability', {})
        if not stab:
            return
        t = np.array(stab.get('times_s', []))
        if len(t) == 0:
            return
        fig, axs = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
        # 能量
        e = np.array(stab.get('energy_eV', []))
        axs[0].plot(t[: len(e)], e, 'k-', linewidth=1)
        axs[0].set_ylabel('能量 (eV)')
        axs[0].grid(True, alpha=0.3)
        # 温度
        temp = np.array(stab.get('temperature_K', []))
        axs[1].plot(t[: len(temp)], temp, 'r-', linewidth=1)
        axs[1].set_ylabel('温度 (K)')
        axs[1].grid(True, alpha=0.3)
        # 压力
        p = np.array(stab.get('pressure_kB', []))
        axs[2].plot(t[: len(p)], p, 'b-', linewidth=1)
        axs[2].set_ylabel('压力 (kB)')
        axs[2].set_xlabel('时间 (s)')
        axs[2].grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'stability_energy_temp_pressure.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 晶格与密度
        fig, ax = plt.subplots(figsize=(12, 6))
        a = np.array(stab['lattice'].get('a', []))
        b = np.array(stab['lattice'].get('b', []))
        c = np.array(stab['lattice'].get('c', []))
        vol = np.array(stab['lattice'].get('vol', []))
        ax.plot(t[: len(a)], a, label='a (Å)')
        ax.plot(t[: len(b)], b, label='b (Å)')
        ax.plot(t[: len(c)], c, label='c (Å)')
        ax2 = ax.twinx()
        dens = np.array(stab.get('density_kg_per_m3', []))
        ax2.plot(t[: len(dens)], dens, 'k--', alpha=0.6, label='密度 (kg/m^3)')
        ax.set_xlabel('时间 (s)')
        ax.set_ylabel('晶格参数 (Å)')
        ax2.set_ylabel('密度 (kg/m^3)')
        ax.legend(loc='upper left')
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(self.output_dir / 'lattice_density.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

    def _plot_rdf(self) -> None:
        rdf = self.analysis_data.get('rdf', {})
        if not rdf:
            return
        r = np.array(rdf.get('r_A', []))
        if len(r) == 0:
            return
        # 总 RDF
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(r, rdf.get('g_all', []), 'k-', linewidth=2, label='所有原子')
        ax.set_xlabel('r (Å)')
        ax.set_ylabel('g(r)')
        ax.set_title('总 RDF')
        ax.grid(True, alpha=0.3)
        ax.legend()
        plt.tight_layout()
        plt.savefig(self.output_dir / 'rdf_all.png', dpi=300, bbox_inches='tight')
        plt.close(fig)

        # 分元素对 RDF（保存所有元素对）
        g_pairs = rdf.get('g_pairs', {})
        for key, val in g_pairs.items():
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.plot(r, val.get('g_r', []), linewidth=2, label=key)
            # 标注峰
            for pk in val.get('peaks', [])[:3]:
                ax.axvline(x=pk['r'], color='r', linestyle='--', alpha=0.5)
            ax.set_xlabel('r (Å)')
            ax.set_ylabel('g(r)')
            ax.set_title(f'RDF: {key} (CN≈{val.get("coordination_number", 0):.2f})')
            ax.grid(True, alpha=0.3)
            ax.legend()
            plt.tight_layout()
            plt.savefig(self.output_dir / f'rdf_{key}.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

        # 首峰演化
        evo = rdf.get('evolution', {})
        r1 = np.array(evo.get('first_peak_r_series_A', []))
        if len(r1) > 0:
            time_step_ps_val = float(self.time_step_ps if self.time_step_ps is not None else 1.0)
            t = np.arange(len(r1), dtype=float) * time_step_ps_val * PS_TO_S * int(evo.get('stride', 1))
            fig, ax = plt.subplots(figsize=(10, 5))
            ax.plot(t, r1, 'b-', linewidth=1.5)
            ax.set_xlabel('时间 (s)')
            ax.set_ylabel('首峰位置 r1 (Å)')
            ax.set_title('RDF 首峰位置演化')
            ax.grid(True, alpha=0.3)
            plt.tight_layout()
            plt.savefig(self.output_dir / 'rdf_first_peak_evolution.png', dpi=300, bbox_inches='tight')
            plt.close(fig)

    # =============================
    # 报告
    # =============================
    def _generate_markdown_report(self) -> None:
        logger.info('📝 生成 Markdown 报告…')
        lines: List[str] = []
        lines.append('# VASP MD 分析报告')
        lines.append('')
        lines.append('## 基本信息')
        lines.append(f'- 任务ID: {self.task_id}')
        lines.append(f'- 帧数: {len(self.structures)}')
        lines.append(f'- 时间步长: {self.time_step_ps:.6f} ps')
        lines.append(f'- 平均温度: {self.temperature_K:.1f} K')
        lines.append('')
        lines.append('## 扩散性质 (DiffusionAnalyzer)')
        diff = self.analysis_data.get('diffusion', {})
        for elem, v in (diff.get('diffusion_by_element', {}) or {}).items():
            D_m2 = v.get("D_m2_per_s", 0)
            D_cm2 = v.get("D_cm2_per_s", 0)
            r2 = v.get("fit_r2", 0)
            lines.append(f'- {elem}: D = {D_m2:.3e} m²/s = {D_cm2:.3e} cm²/s (R²={r2:.3f})')
        if diff.get('ionic_conductivity_S_per_m'):
            lines.append('')
            lines.append('### 离子电导率 (基于 DiffusionAnalyzer)')
            for elem, sig in diff['ionic_conductivity_S_per_m'].items():
                lines.append(f'- {elem}: σ = {sig:.3e} S/m')
            mobile_species = diff.get('mobile_species', [])
            if mobile_species:
                lines.append(f'- 识别的迁移离子: {", ".join(mobile_species)}')
        lines.append('')
        if 'arrhenius' in self.analysis_data:
            arrh = self.analysis_data['arrhenius']
            fit = arrh.get('fit', {})
            if fit and fit.get('Ea_eV') is not None:
                lines.append(f'- Arrhenius 拟合激活能: {fit["Ea_eV"]:.3f} eV (R^2={fit.get("R2", 0):.3f})')
        lines.append('')
        lines.append('## RDF/配位 (PyMatGen优化)')
        rdf = self.analysis_data.get('rdf', {})
        if rdf:
            method = rdf.get('method', 'unknown')
            lines.append(f'- 分析方法: {method}')
            lines.append(f'- 计算范围: 0 - {self.rdf_rmax} Å')
            
            # 元素对RDF结果
            g_pairs = rdf.get('g_pairs', {})
            if g_pairs:
                lines.append('- 元素对 RDF 分析结果:')
                for pair, data in g_pairs.items():
                    cn = data.get('coordination_number', 0)
                    peaks = data.get('peaks', [])
                    first_peak = peaks[0]['r'] if peaks else 'N/A'
                    lines.append(f'  - {pair}: 配位数={cn:.2f}, 首峰位置={first_peak} Å')
            
            # 结构演化
            evolution = rdf.get('evolution', {})
            if evolution and evolution.get('first_peak_r_series_A'):
                ref_species = evolution.get('reference_species', 'unknown')
                lines.append(f'- 首峰演化追踪: {ref_species} 原子（{len(evolution["first_peak_r_series_A"])} 个时间点）')
        lines.append('')
        lines.append('## 稳定性监控')
        stab = self.analysis_data.get('stability', {})
        if stab:
            lines.append(f'- 平衡能量索引: {stab.get("equilibrium_index_energy")}')
            lines.append(f'- 平衡温度索引: {stab.get("equilibrium_index_temperature")}')
            anom = stab.get('anomalies', {})
            lines.append(f'- 异常预警: 能量{len(anom.get("energy", []))} 温度{len(anom.get("temperature", []))} 压力{len(anom.get("pressure", []))}')
        with open(self.output_dir / 'analysis_report.md', 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines))
        logger.info('已生成 Markdown 报告: analysis_report.md')

    def _generate_html_report(self) -> str:
        logger.info('🌐 生成 HTML 报告…')
        generator = MDHTMLReportGenerator(self.analysis_data)
        html_path = str(self.output_dir / 'md_analysis_report.html')
        generator.generate_html_report(html_path)
        logger.info(f'HTML 报告: {html_path}')
        return html_path


class MDHTMLReportGenerator:
    """MD 分析 HTML 报告生成器（风格参考 PyMatGenDOSHTMLGenerator）"""
    def __init__(self, analysis_data: Dict[str, Any]) -> None:
        self.data = analysis_data

    def generate_html_report(self, output_path: str) -> str:
        html = self._render()
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html)
        return output_path

    def _render(self) -> str:
        charts = self._collect_charts()
        return f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="UTF-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1.0" />
  <title>VASP MD 分析报告</title>
  {self._css()}
  <style>img{{max-width:100%;height:auto;}}</style>
</head>
<body>
  <div class="container">
    {self._header()}
    {self._summary()}
    {self._section_diffusion(charts)}
    {self._section_rdf(charts)}
    {self._section_stability(charts)}
    {self._footer()}
  </div>
  
</body>
</html>
        """

    def _collect_charts(self) -> Dict[str, str]:
        out = {}
        img_files = {
            'msd': 'msd.png',
            'arrhenius': 'arrhenius.png',
            'rdf_all': 'rdf_all.png',
            'rdf_first_peak_evolution': 'rdf_first_peak_evolution.png',
            'stability': 'stability_energy_temp_pressure.png',
            'lattice_density': 'lattice_density.png',
        }
        for k, fn in img_files.items():
            fp = Path(self.data.get('task_info', {}).get('output_dir', ''))
            if not fp:
                fp = Path('.')
            # 优先使用分析输出目录
            # 由于此生成器可能并不知晓 output_dir，尝试当前工作目录下的 MD_output
            try_candidates = []
            if 'task_info' in self.data and 'output_dir' in self.data['task_info']:
                try_candidates.append(Path(self.data['task_info']['output_dir']) / fn)
            if 'report_html' in self.data:
                try_candidates.append(Path(self.data['report_html']).parent / fn)
            try_candidates.append(Path('MD_output') / fn)

            for p in try_candidates:
                if p.exists():
                    with open(p, 'rb') as f:
                        out[k] = base64.b64encode(f.read()).decode()
                    break
        
        # 额外：收集所有元素对 RDF 图像（文件名形如 rdf_*.png）
        try:
            rdf_imgs = {}
            search_dirs = []
            if 'task_info' in self.data and 'output_dir' in self.data['task_info']:
                search_dirs.append(Path(self.data['task_info']['output_dir']))
            if 'report_html' in self.data:
                search_dirs.append(Path(self.data['report_html']).parent)
            search_dirs.append(Path('MD_output'))

            for d in search_dirs:
                if d and d.exists():
                    for p in d.glob('rdf_*.png'):
                        key = p.stem  # 如 rdf_Li-O
                        if key not in out:
                            with open(p, 'rb') as f:
                                rdf_imgs[key] = base64.b64encode(f.read()).decode()
            # 合并到 out（不覆盖已存在键）
            out.update({k: v for k, v in rdf_imgs.items() if k not in out})
        except Exception:
            pass
        return out

    def _css(self) -> str:
        return """
<style>
  body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, 'Noto Sans', 'PingFang SC', 'Microsoft YaHei', sans-serif; background:#f5f7fb; }
  .container { max-width: 1200px; margin: 20px auto; background: #fff; border-radius: 12px; box-shadow: 0 10px 30px rgba(0,0,0,.06); padding: 24px; }
  h1 { margin: 0 0 10px; }
  .subtitle { color: #666; margin-bottom: 20px; }
  .section { border:1px solid #eee; border-radius: 10px; padding: 16px; margin: 16px 0; background: #fafbff; }
  .section h2 { border-left: 4px solid #667eea; padding-left: 10px; color:#333 }
  .grid-2 { display:grid; grid-template-columns: 1fr 1fr; gap: 16px; }
  .grid-3 { display:grid; grid-template-columns: 1fr 1fr 1fr; gap: 16px; }
  table { width:100%; border-collapse: collapse; }
  th, td { padding: 8px 10px; border-bottom: 1px solid #eee; text-align:left; }
  th { background: #667eea; color: #fff; }
  .note { background: #fff7e6; padding: 10px; border-radius: 8px; border-left: 4px solid #faad14; }
  
  /* RDF专用样式 */
  .rdf-table { margin: 10px 0; }
  .rdf-table th { background: #4ECDC4; color: #fff; }
  .rdf-table td { vertical-align: middle; }
  
  /* 元素对颜色标识 */
  .color-li { background: linear-gradient(90deg, #ff9a9e 0%, #fecfef 100%); }
  .color-o { background: linear-gradient(90deg, #a8edea 0%, #fed6e3 100%); }
  .color-p { background: linear-gradient(90deg, #ffecd2 0%, #fcb69f 100%); }
  .color-fe { background: linear-gradient(90deg, #667eea 0%, #764ba2 100%); color: white; }
  .color-default { background: #f8f9fa; }
  
  .method-tag { 
    background: #e3f2fd; 
    color: #1565c0; 
    padding: 2px 6px; 
    border-radius: 4px; 
    font-size: 0.8em; 
    font-weight: bold;
  }
  
  /* RDF卡片网格 */
  .rdf-grid { 
    display: grid; 
    grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); 
    gap: 16px; 
    margin: 16px 0; 
  }
  
  .rdf-card { 
    border: 1px solid #e0e0e0; 
    border-radius: 8px; 
    overflow: hidden; 
    background: #fff; 
    box-shadow: 0 2px 8px rgba(0,0,0,0.1);
    transition: transform 0.2s ease, box-shadow 0.2s ease;
  }
  
  .rdf-card:hover { 
    transform: translateY(-2px); 
    box-shadow: 0 4px 12px rgba(0,0,0,0.15); 
  }
  
  .rdf-card-header { 
    background: linear-gradient(135deg, #667eea 0%, #764ba2 100%); 
    color: white; 
    padding: 12px 16px; 
    display: flex; 
    justify-content: space-between; 
    align-items: center; 
  }
  
  .rdf-card-header h5 { 
    margin: 0; 
    font-size: 1.1em; 
    font-weight: bold; 
  }
  
  .cn-badge { 
    background: rgba(255,255,255,0.2); 
    padding: 4px 8px; 
    border-radius: 12px; 
    font-size: 0.9em; 
    font-weight: bold; 
  }
  
  .rdf-card img { 
    width: 100%; 
    height: auto; 
    display: block; 
  }
  
  /* 交互式图表控制按钮 */
  .chart-controls { 
    text-align: center; 
  }
  
  .btn { 
    background: #667eea; 
    color: white; 
    border: none; 
    padding: 8px 16px; 
    border-radius: 4px; 
    cursor: pointer; 
    margin: 0 4px; 
    transition: background 0.2s ease;
  }
  
  .btn:hover { 
    background: #5a67d8; 
  }
  
  .btn-sm { 
    padding: 6px 12px; 
    font-size: 0.9em; 
  }
  
  /* 响应式设计 */
  @media (max-width: 768px) {
    .grid-2, .grid-3 { grid-template-columns: 1fr; }
    .rdf-grid { grid-template-columns: 1fr; }
    .container { margin: 10px; padding: 16px; }
  }
</style>
        """

    def _header(self) -> str:
        task = self.data.get('task_info', {})
        return f"""
  <div>
    <h1>🧪 VASP MD 分析报告</h1>
    <div class="subtitle">任务ID: {task.get('task_id', 'unknown')} | 生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
  </div>
        """

    def _summary(self) -> str:
        stab = self.data.get('stability', {})
        rdf = self.data.get('rdf', {})
        diff = self.data.get('diffusion', {})
        # 选取一个代表扩散系数
        D_mean = 0.0
        if diff.get('diffusion_by_element'):
            D_mean = float(np.mean([v['D_m2_per_s'] for v in diff['diffusion_by_element'].values()]))
        return f"""
  <div class="section">
    <h2>📋 摘要</h2>
    <div class="grid-3">
      <div>
        <table>
          <tr><th>属性</th><th>值</th></tr>
          <tr><td>帧数</td><td>{len(self.data.get('stability', {}).get('lattice', {}).get('a', []))}</td></tr>
          <tr><td>平均温度 (K)</td><td>{self.data.get('diffusion', {}).get('arrhenius_single', {}).get('T_list_K', [0])[0]:.1f}</td></tr>
          <tr><td>平均扩散系数 (m^2/s)</td><td>{D_mean:.3e}</td></tr>
        </table>
      </div>
      <div>
        <table>
          <tr><th>扩散 (示例)</th><th>值</th></tr>
          {self._diffusion_rows()}
        </table>
      </div>
      <div>
        <table>
          <tr><th>RDF/配位</th><th>信息</th></tr>
          <tr><td>RDF 半径</td><td>{self.data.get('rdf', {}).get('r_A', [0])[-1] if rdf else 'N/A'} Å</td></tr>
          <tr><td>首峰追踪点数</td><td>{len(self.data.get('rdf', {}).get('evolution', {}).get('first_peak_r_series_A', [])) if rdf else 0}</td></tr>
        </table>
      </div>
    </div>
  </div>
        """

    def _diffusion_rows(self) -> str:
        diff = self.data.get('diffusion', {})
        rows = []
        if not diff or not diff.get('diffusion_by_element'):
            return ''
        for elem, v in diff['diffusion_by_element'].items():
            rows.append(f"<tr><td>{elem}</td><td>{v['D_m2_per_s']:.3e} m^2/s</td></tr>")
        return '\n'.join(rows[:6])

    def _img(self, img64: str, title: str) -> str:
        if not img64:
            return '<div class="note">无图像</div>'
        return f'<img alt="{title}" src="data:image/png;base64,{img64}" />'
    
    def _generate_rdf_table(self, g_pairs: Dict[str, Dict[str, Any]]) -> str:
        """生成RDF配位数汇总表格"""
        if not g_pairs:
            return '<div class="note">无RDF数据</div>'
        
        rows = []
        for pair, data in g_pairs.items():
            cn = data.get('coordination_number', 0)
            peaks = data.get('peaks', [])
            first_peak = f"{peaks[0]['r']:.2f}" if peaks else 'N/A'
            method = data.get('method', 'unknown')
            
            # 为不同元素对添加颜色标识
            color_class = self._get_pair_color_class(pair)
            
            rows.append(f"""
                <tr class="{color_class}">
                    <td><strong>{pair}</strong></td>
                    <td>{cn:.2f}</td>
                    <td>{first_peak} Å</td>
                    <td>{len(peaks)}</td>
                    <td><span class="method-tag">{method}</span></td>
                </tr>
            """)
        
        return f"""
        <table class="rdf-table">
            <thead>
                <tr>
                    <th>元素对</th>
                    <th>配位数</th>
                    <th>首峰位置</th>
                    <th>峰数量</th>
                    <th>分析方法</th>
                </tr>
            </thead>
            <tbody>
                {''.join(rows)}
            </tbody>
        </table>
        """
    
    def _get_pair_color_class(self, pair: str) -> str:
        """为不同元素对分配颜色类"""
        color_map = {
            'Li': 'color-li',
            'O': 'color-o', 
            'P': 'color-p',
            'Fe': 'color-fe'
        }
        
        elements = pair.split('-')
        if len(elements) == 2:
            elem1, elem2 = elements
            if elem1 == elem2:
                return color_map.get(elem1, 'color-default')
            else:
                # 混合元素对使用渐变色
                return f"color-pair-{elem1.lower()}-{elem2.lower()}"
        return 'color-default'
    
    def _generate_rdf_pair_images(self, charts: Dict[str, str], g_pairs: Dict[str, Dict[str, Any]]) -> str:
        """生成元素对RDF图像网格"""
        if not g_pairs:
            return '<div class="note">无元素对RDF图像</div>'
        
        # 按元素类型分组显示
        pair_groups = self._group_pairs_by_type(g_pairs.keys())
        
        html_parts = []
        for group_name, pairs in pair_groups.items():
            if not pairs:
                continue
                
            html_parts.append(f'<h4>{group_name}</h4>')
            html_parts.append('<div class="rdf-grid">')
            
            for pair in pairs:
                chart_key = f'rdf_{pair}'
                img64 = charts.get(chart_key, '')
                cn = g_pairs[pair].get('coordination_number', 0)
                
                html_parts.append(f'''
                    <div class="rdf-card">
                        <div class="rdf-card-header">
                            <h5>{pair}</h5>
                            <span class="cn-badge">CN: {cn:.2f}</span>
                        </div>
                        {self._img(img64, f'rdf_{pair}')}
                    </div>
                ''')
            
            html_parts.append('</div>')
        
        return ''.join(html_parts)
    
    def _group_pairs_by_type(self, pairs: Iterable[str]) -> Dict[str, List[str]]:
        """按元素类型对元素对分组"""
        groups = {
            '同元素对': [],
            '锂相关': [],
            '氧相关': [], 
            '磷相关': [],
            '铁相关': [],
            '其他': []
        }
        
        for pair in pairs:
            elements = pair.split('-')
            if len(elements) == 2:
                elem1, elem2 = elements
                if elem1 == elem2:
                    groups['同元素对'].append(pair)
                elif 'Li' in elements:
                    groups['锂相关'].append(pair)
                elif 'O' in elements:
                    groups['氧相关'].append(pair)
                elif 'P' in elements:
                    groups['磷相关'].append(pair)
                elif 'Fe' in elements:
                    groups['铁相关'].append(pair)
                else:
                    groups['其他'].append(pair)
        
        # 移除空组
        return {k: v for k, v in groups.items() if v}
    
    def _generate_interactive_rdf_chart(self, r_array: List[float], g_pairs: Dict[str, Dict[str, Any]]) -> str:
        """生成交互式RDF对比图表（使用Chart.js）"""
        if len(r_array) == 0 or not g_pairs:
            return ""
        
        # 准备数据
        datasets = []
        colors = ['#FF6384', '#36A2EB', '#FFCE56', '#4BC0C0', '#9966FF', '#FF9F40', '#FF6384', '#C9CBCF']
        
        for i, (pair, data) in enumerate(g_pairs.items()):
            g_r = data.get('g_r', [])
            if len(g_r) != len(r_array):
                continue
                
            color = colors[i % len(colors)]
            datasets.append({
                'label': pair,
                'data': [{'x': float(r), 'y': float(g)} for r, g in zip(r_array, g_r)],
                'borderColor': color,
                'backgroundColor': color + '20',  # 20% 透明度
                'borderWidth': 2,
                'fill': False,
                'pointRadius': 0,
                'pointHoverRadius': 4
            })
        
        chart_data = {
            'datasets': datasets
        }
        
        return f"""
        <div style="margin: 20px 0;">
            <h3>交互式RDF对比图</h3>
            <div style="position: relative; height: 400px; width: 100%;">
                <canvas id="rdfChart"></canvas>
            </div>
            <div class="chart-controls" style="margin-top: 10px;">
                <button onclick="toggleAllSeries()" class="btn btn-sm">显示/隐藏全部</button>
                <button onclick="resetZoom()" class="btn btn-sm">重置缩放</button>
            </div>
        </div>
        
        <script src="https://cdn.jsdelivr.net/npm/chart.js"></script>
        <script>
            const chartData = {json.dumps(chart_data)};
            const ctx = document.getElementById('rdfChart').getContext('2d');
            
            const rdfChart = new Chart(ctx, {{
                type: 'line',
                data: chartData,
                options: {{
                    responsive: true,
                    maintainAspectRatio: false,
                    interaction: {{
                        intersect: false,
                        mode: 'index'
                    }},
                    plugins: {{
                        legend: {{
                            position: 'top',
                            onClick: function(e, legendItem) {{
                                const index = legendItem.datasetIndex;
                                const chart = this.chart;
                                const meta = chart.getDatasetMeta(index);
                                meta.hidden = meta.hidden === null ? !chart.data.datasets[index].hidden : null;
                                chart.update();
                            }}
                        }},
                        title: {{
                            display: true,
                            text: '径向分布函数 g(r) 对比'
                        }},
                        tooltip: {{
                            callbacks: {{
                                title: function(context) {{
                                    return 'r = ' + context[0].parsed.x.toFixed(2) + ' Å';
                                }},
                                label: function(context) {{
                                    return context.dataset.label + ': g(r) = ' + context.parsed.y.toFixed(3);
                                }}
                            }}
                        }}
                    }},
                    scales: {{
                        x: {{
                            type: 'linear',
                            position: 'bottom',
                            title: {{
                                display: true,
                                text: 'r (Å)'
                            }},
                            min: 0,
                            max: 10
                        }},
                        y: {{
                            title: {{
                                display: true,
                                text: 'g(r)'
                            }},
                            min: 0
                        }}
                    }},
                    elements: {{
                        line: {{
                            tension: 0.1
                        }}
                    }}
                }}
            }});
            
            function toggleAllSeries() {{
                const chart = rdfChart;
                const allHidden = chart.data.datasets.every(dataset => {{
                    const meta = chart.getDatasetMeta(chart.data.datasets.indexOf(dataset));
                    return meta.hidden === true;
                }});
                
                chart.data.datasets.forEach((dataset, index) => {{
                    const meta = chart.getDatasetMeta(index);
                    meta.hidden = !allHidden;
                }});
                chart.update();
            }}
            
            function resetZoom() {{
                rdfChart.resetZoom();
            }}
        </script>
        """

    def _section_diffusion(self, charts: Dict[str, str]) -> str:
        # 电导率表格
        diff = self.data.get('diffusion', {})
        conductivity = diff.get('ionic_conductivity_S_per_m', {}) or {}
        cond_rows = ''
        if conductivity:
            cond_rows = '\n'.join([
                f"<tr><td>{elem}</td><td>{val:.3e} S/m</td></tr>" for elem, val in conductivity.items()
            ])
        conductivity_table = f"""
            <div>
              <h3>离子电导率（Nernst–Einstein 或 DiffusionAnalyzer）</h3>
              <table>
                <tr><th>元素</th><th>σ (S/m)</th></tr>
                {cond_rows if cond_rows else '<tr><td colspan=2>无数据</td></tr>'}
              </table>
            </div>
        """

        return f"""
  <div class="section">
    <h2>📈 扩散性质</h2>
    <div class="grid-2">
      <div>
        <h3>MSD（按元素）</h3>
        {self._img(charts.get('msd', ''), 'msd')}
      </div>
      <div>
        <h3>Arrhenius（多温聚合时显示）</h3>
        {self._img(charts.get('arrhenius', ''), 'arrhenius')}
      </div>
    </div>
    {conductivity_table}
  </div>
        """

    def _section_rdf(self, charts: Dict[str, str]) -> str:
        # 获取RDF数据用于交互式可视化
        rdf_data = self.data.get('rdf', {})
        g_pairs = rdf_data.get('g_pairs', {})
        r_array = rdf_data.get('r_A', [])
        
        # 生成元素对RDF表格
        rdf_table = self._generate_rdf_table(g_pairs)
        
        # 生成元素对RDF图表网格
        rdf_pair_images = self._generate_rdf_pair_images(charts, g_pairs)
        
        # 交互式RDF图表（如果有数据）
        interactive_chart = self._generate_interactive_rdf_chart(r_array, g_pairs) if len(r_array) > 0 and g_pairs else ""
        
        return f"""
  <div class="section">
    <h2>🔬 RDF/配位分析</h2>
    
    <!-- 总览图 -->
    <div class="grid-2">
      <div>
        <h3>总 RDF</h3>
        {self._img(charts.get('rdf_all', ''), 'rdf_all')}
      </div>
      <div>
        <h3>首峰演化</h3>
        {self._img(charts.get('rdf_first_peak_evolution', ''), 'rdf_first_peak_evolution')}
      </div>
    </div>
    
    <!-- 配位数汇总表 -->
    <div style="margin: 20px 0;">
      <h3>配位数与峰位汇总</h3>
      {rdf_table}
    </div>
    
    <!-- 交互式RDF对比图 -->
    {interactive_chart}
    
    <!-- 元素对RDF详细图表 -->
    <div style="margin: 20px 0;">
      <h3>元素对RDF详细分析</h3>
      {rdf_pair_images}
    </div>
  </div>
        """

    def _section_stability(self, charts: Dict[str, str]) -> str:
        return f"""
  <div class="section">
    <h2>🛡️ 稳定性监控</h2>
    <div class="grid-2">
      <div>
        <h3>能量/温度/压力</h3>
        {self._img(charts.get('stability', ''), 'stability')}
      </div>
      <div>
        <h3>晶格/密度</h3>
        {self._img(charts.get('lattice_density', ''), 'lattice_density')}
      </div>
    </div>
  </div>
        """

    def _footer(self) -> str:
        return f"""
  <div style="text-align:center;color:#888;margin-top:20px;">
    <div>VASP MD 分析报告 | 生成时间 {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</div>
  </div>
        """


def generate_md_analysis_report(input_path: str, task_id: Optional[str] = None, output_dir: Optional[str] = None) -> str:
    """便捷函数：执行完整 MD 分析并生成 HTML 报告。

    Args:
        input_path: MD 输出目录（含 XDATCAR/OUTCAR/vasprun.xml）或其上级（含多个子运行）
        task_id: 任务ID
        output_dir: 输出目录

    Returns:
        HTML 报告路径
    """
    analyzer = VASP_MDAnalyzer(input_path, task_id=task_id, output_dir=output_dir)

    # 任务信息（供HTML收集图像路径）
    analyzer.analysis_data['task_info'] = {
        'task_id': analyzer.task_id,
        'input_path': str(analyzer.input_path),
        'output_dir': str(analyzer.output_dir),
        'timestamp': datetime.now().isoformat(),
        'analysis_type': 'VASP_MD_Analysis'
    }
    analyzer.analyze()
    return str(analyzer.output_dir / 'md_analysis_report.html')


if __name__ == '__main__':
    # 简易测试入口（请根据实际MD数据路径调整）
    test_dir = '/Users/ysl/Desktop/Code/MCP_IC_Tool/md_test/25e88a44d8274fad84caf4e6c5561868'
    output_dir = '/Users/ysl/Desktop/Code/MCP_IC_Tool/md_test/25e88a44d8274fad84caf4e6c5561868/output'
    try:
        html = generate_md_analysis_report(test_dir, task_id='md_test')
        print(f'✅ MD 分析报告已生成: {html}')
    except Exception as e:
        print(f'❌ 测试失败: {e}')
