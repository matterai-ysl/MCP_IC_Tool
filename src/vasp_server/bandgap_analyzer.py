import numpy as np
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin

class BandGapAnalyzer:
    def __init__(self, vasprun_path="vasprun.xml"):
        """
        初始化带隙分析器
        
        Args:
            vasprun_path: vasprun.xml文件路径
        """
        self.vasprun = Vasprun(vasprun_path)
        self.eigenvalues = self.vasprun.eigenvalues
        self.kpoints = self.vasprun.actual_kpoints
        self.fermi_level = self.vasprun.efermi
        self.is_spin_polarized = self.vasprun.is_spin
    
    def find_band_edges(self, spin=Spin.up):
        """
        找到价带顶(VBM)和导带底(CBM)
        
        Args:
            spin: 自旋通道
            
        Returns:
            dict: 包含VBM和CBM信息的字典
        """
        eigenvals = self.eigenvalues[spin] #type: ignore
        
        vbm_energy = -float('inf')
        cbm_energy = float('inf')
        vbm_kpoint_idx = -1
        cbm_kpoint_idx = -1
        vbm_band_idx = -1
        cbm_band_idx = -1
        
        # 遍历所有k点和能带
        for k_idx, kpoint_data in enumerate(eigenvals):
            for band_idx, (energy, occupation) in enumerate(kpoint_data):
                if occupation > 0.5:  # 占据态
                    if energy > vbm_energy:
                        vbm_energy = energy
                        vbm_kpoint_idx = k_idx
                        vbm_band_idx = band_idx
                else:  # 空态
                    if energy < cbm_energy:
                        cbm_energy = energy
                        cbm_kpoint_idx = k_idx
                        cbm_band_idx = band_idx
        
        return {
            'vbm_energy': vbm_energy,
            'cbm_energy': cbm_energy,
            'vbm_kpoint_idx': vbm_kpoint_idx,
            'cbm_kpoint_idx': cbm_kpoint_idx,
            'vbm_kpoint': self.kpoints[vbm_kpoint_idx],
            'cbm_kpoint': self.kpoints[cbm_kpoint_idx],
            'vbm_band_idx': vbm_band_idx,
            'cbm_band_idx': cbm_band_idx
        }
    
    def calculate_global_bandgap(self):
        """
        计算考虑所有自旋通道的全局带隙，包括直接、间接和基本带隙
        
        Returns:
            dict: 全局带隙信息
        """
        global_vbm = -float('inf')
        global_cbm = float('inf')
        global_vbm_info = {}
        global_cbm_info = {}
        all_direct_gaps = []
        
        spins_to_check = [Spin.up, Spin.down] if self.is_spin_polarized else [Spin.up]
        
        # 找全局VBM和CBM
        for spin in spins_to_check:
            edges = self.find_band_edges(spin)
            spin_label = "up" if spin == Spin.up else "down"
            
            if edges['vbm_energy'] > global_vbm:
                global_vbm = edges['vbm_energy']
                global_vbm_info = {**edges, 'spin': spin_label}
                
            if edges['cbm_energy'] < global_cbm:
                global_cbm = edges['cbm_energy']
                global_cbm_info = {**edges, 'spin': spin_label}
            
            # 收集所有自旋通道的直接带隙
            direct_gaps = self.calculate_direct_gaps(spin)
            for gap_info in direct_gaps:
                gap_info['spin'] = spin_label
            all_direct_gaps.extend(direct_gaps)
        
        # 基本带隙（全局最小带隙）
        fundamental_gap = global_cbm - global_vbm
        
        # 判断基本带隙是直接还是间接
        same_kpoint = np.allclose(global_vbm_info['vbm_kpoint'], 
                                 global_cbm_info['cbm_kpoint'], atol=1e-6)
        fundamental_type = "direct" if same_kpoint else "indirect"
        
        # 最小直接带隙（所有k点中最小的直接带隙）
        min_direct_gap_info = min(all_direct_gaps, key=lambda x: x['gap']) if all_direct_gaps else None
        min_direct_gap = min_direct_gap_info['gap'] if min_direct_gap_info else None
        
        # 间接带隙（如果基本带隙是直接的，则间接带隙需要从其他k点组合计算）
        if fundamental_type == "direct":
            # 基本带隙就是直接带隙，寻找最小的间接带隙
            indirect_gap = None
            indirect_gap_info = {}
            
            # 尝试不同k点组合寻找最小间接带隙
            for spin in spins_to_check:
                edges = self.find_band_edges(spin)
                spin_label = "up" if spin == Spin.up else "down"
                
                # 如果这个自旋通道的VBM或CBM与全局的不在同一k点
                vbm_diff_k = not np.allclose(edges['vbm_kpoint'], global_cbm_info['cbm_kpoint'], atol=1e-6)
                cbm_diff_k = not np.allclose(edges['cbm_kpoint'], global_vbm_info['vbm_kpoint'], atol=1e-6)
                
                if vbm_diff_k:
                    gap = global_cbm - edges['vbm_energy']
                    if indirect_gap is None or gap < indirect_gap:
                        indirect_gap = gap
                        indirect_gap_info = {
                            'gap': gap,
                            'vbm_energy': edges['vbm_energy'],
                            'cbm_energy': global_cbm,
                            'vbm_kpoint': edges['vbm_kpoint'].tolist(),
                            'cbm_kpoint': global_cbm_info['cbm_kpoint'],
                            'vbm_spin': spin_label,
                            'cbm_spin': global_cbm_info['spin']
                        }
                
                if cbm_diff_k:
                    gap = edges['cbm_energy'] - global_vbm
                    if indirect_gap is None or gap < indirect_gap:
                        indirect_gap = gap
                        indirect_gap_info = {
                            'gap': gap,
                            'vbm_energy': global_vbm,
                            'cbm_energy': edges['cbm_energy'],
                            'vbm_kpoint': global_vbm_info['vbm_kpoint'],
                            'cbm_kpoint': edges['cbm_kpoint'].tolist(),
                            'vbm_spin': global_vbm_info['spin'],
                            'cbm_spin': spin_label
                        }
        else:
            # 基本带隙就是间接带隙
            indirect_gap = fundamental_gap
            indirect_gap_info = {
                'gap': indirect_gap,
                'vbm_energy': global_vbm,
                'cbm_energy': global_cbm,
                'vbm_kpoint': global_vbm_info['vbm_kpoint'],
                'cbm_kpoint': global_cbm_info['cbm_kpoint'],
                'vbm_spin': global_vbm_info['spin'],
                'cbm_spin': global_cbm_info['spin']
            }
        
        return {
            # 基本带隙（全局最小）
            'fundamental_gap': fundamental_gap,
            'fundamental_type': fundamental_type,
            
            # 直接带隙（最小的直接带隙）
            'direct_gap': min_direct_gap,
            'direct_gap_kpoint': min_direct_gap_info['kpoint'] if min_direct_gap_info else None,
            'direct_gap_kpoint_idx': min_direct_gap_info['kpoint_idx'] if min_direct_gap_info else None,
            'direct_gap_spin': min_direct_gap_info['spin'] if min_direct_gap_info else None,
            
            # 间接带隙
            'indirect_gap': indirect_gap,
            'indirect_gap_info': indirect_gap_info if indirect_gap else None,
            
            # 系统性质
            'is_metal': fundamental_gap <= 0.0,
            
            # 全局带边信息
            'global_vbm_energy': global_vbm,
            'global_cbm_energy': global_cbm,
            'vbm_spin': global_vbm_info['spin'],
            'cbm_spin': global_cbm_info['spin'],
            'vbm_kpoint_idx': global_vbm_info['vbm_kpoint_idx'],
            'cbm_kpoint_idx': global_cbm_info['cbm_kpoint_idx'],
            'vbm_kpoint': global_vbm_info['vbm_kpoint'],
            'cbm_kpoint': global_cbm_info['cbm_kpoint'],
            'vbm_band_idx': global_vbm_info['vbm_band_idx'],
            'cbm_band_idx': global_cbm_info['cbm_band_idx']
        }
    
    def calculate_direct_gaps(self, spin=Spin.up):
        """
        计算所有k点的直接带隙
        
        Args:
            spin: 自旋通道
            
        Returns:
            list: 每个k点的直接带隙信息
        """
        eigenvals = self.eigenvalues[spin] #type: ignore
        direct_gaps = []
        
        for k_idx, kpoint_data in enumerate(eigenvals):
            # 找到该k点的最高占据态和最低空态
            occupied_energies = []
            unoccupied_energies = []
            
            for energy, occupation in kpoint_data:
                if occupation > 0.5:
                    occupied_energies.append(energy)
                else:
                    unoccupied_energies.append(energy)
            
            if occupied_energies and unoccupied_energies:
                local_vbm = max(occupied_energies)
                local_cbm = min(unoccupied_energies)
                direct_gap = local_cbm - local_vbm
                direct_gaps.append({
                    'gap': direct_gap,
                    'kpoint_idx': k_idx,
                    'kpoint': self.kpoints[k_idx],
                    'vbm': local_vbm,
                    'cbm': local_cbm
                })
        
        return direct_gaps
    
    def calculate_spin_resolved_gaps(self):
        """
        计算每个自旋通道的带隙信息
        
        Returns:
            dict: 每个自旋通道的详细信息
        """
        results = {}
        
        spins_to_check = [Spin.up, Spin.down] if self.is_spin_polarized else [Spin.up]
        
        for spin in spins_to_check:
            spin_label = "up" if spin == Spin.up else "down"
            edges = self.find_band_edges(spin)
            direct_gaps = self.calculate_direct_gaps(spin)
            
            # 找最小的直接带隙
            min_direct_gap = min(direct_gaps, key=lambda x: x['gap']) if direct_gaps else None
            
            # 基本带隙
            fundamental_gap = edges['cbm_energy'] - edges['vbm_energy']
            
            # 判断直接还是间接
            same_kpoint = np.allclose(edges['vbm_kpoint'], edges['cbm_kpoint'], atol=1e-6)
            gap_type = "direct" if same_kpoint else "indirect"
            
            results[spin_label] = {
                'fundamental_gap': fundamental_gap,
                'gap_type': gap_type,
                'is_metal': fundamental_gap <= 0.0,
                'vbm_energy': edges['vbm_energy'],
                'cbm_energy': edges['cbm_energy'],
                'vbm_kpoint': edges['vbm_kpoint'],
                'cbm_kpoint': edges['cbm_kpoint'],
                'vbm_kpoint_idx': edges['vbm_kpoint_idx'],
                'cbm_kpoint_idx': edges['cbm_kpoint_idx'],
                'vbm_band_idx': edges['vbm_band_idx'],
                'cbm_band_idx': edges['cbm_band_idx'],
                'min_direct_gap': min_direct_gap['gap'] if min_direct_gap else None,
                'min_direct_gap_kpoint': min_direct_gap['kpoint'] if min_direct_gap else None,
                'min_direct_gap_kpoint_idx': min_direct_gap['kpoint_idx'] if min_direct_gap else None,
                'all_direct_gaps': direct_gaps
            }
        
        return results
    
    def analyze(self):
        """
        执行完整的带隙分析
        
        Returns:
            dict: 完整的带隙分析结果
        """
        # 全局带隙分析
        global_info = self.calculate_global_bandgap()
        
        # 自旋分辨带隙分析
        spin_resolved = self.calculate_spin_resolved_gaps()
        
        # 系统基本信息
        system_info = {
            'fermi_level': self.fermi_level,
            'is_spin_polarized': self.is_spin_polarized,
            'n_kpoints': len(self.kpoints),
            'n_bands': len(self.eigenvalues[Spin.up][0]) #type: ignore
        }
        
        return {
            'system_info': system_info,
            'global_bandgap': global_info,
            'spin_resolved': spin_resolved
        }

# 使用示例
def analyze_bandgap(vasprun_path="vasprun.xml"):
    """
    便捷函数：分析带隙并返回结果字典
    
    Args:
        vasprun_path: vasprun.xml文件路径
        
    Returns:
        dict: 完整的带隙分析结果
    """
    analyzer = BandGapAnalyzer(vasprun_path)
    return analyzer.analyze()
