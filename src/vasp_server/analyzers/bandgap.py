"""
VASP 带隙分析器 — 从 bandgap_analyzer.py 迁移。

包含:
- BandGapAnalyzer: 基于 vasprun.xml 的带隙分析
- analyze_bandgap: 入口便捷函数
"""

import numpy as np
from pymatgen.io.vasp import Vasprun
from pymatgen.electronic_structure.core import Spin


class BandGapAnalyzer:
    def __init__(self, vasprun_path="vasprun.xml"):
        self.vasprun = Vasprun(vasprun_path)
        self.eigenvalues = self.vasprun.eigenvalues
        self.kpoints = self.vasprun.actual_kpoints
        self.fermi_level = self.vasprun.efermi
        self.is_spin_polarized = self.vasprun.is_spin

    def find_band_edges(self, spin=Spin.up):
        eigenvals = self.eigenvalues[spin]  # type: ignore
        vbm_energy = -float('inf')
        cbm_energy = float('inf')
        vbm_kpoint_idx = cbm_kpoint_idx = -1
        vbm_band_idx = cbm_band_idx = -1

        for k_idx, kpoint_data in enumerate(eigenvals):
            for band_idx, (energy, occupation) in enumerate(kpoint_data):
                if occupation > 0.5:
                    if energy > vbm_energy:
                        vbm_energy = energy
                        vbm_kpoint_idx = k_idx
                        vbm_band_idx = band_idx
                else:
                    if energy < cbm_energy:
                        cbm_energy = energy
                        cbm_kpoint_idx = k_idx
                        cbm_band_idx = band_idx

        return {
            'vbm_energy': vbm_energy, 'cbm_energy': cbm_energy,
            'vbm_kpoint_idx': vbm_kpoint_idx, 'cbm_kpoint_idx': cbm_kpoint_idx,
            'vbm_kpoint': self.kpoints[vbm_kpoint_idx],
            'cbm_kpoint': self.kpoints[cbm_kpoint_idx],
            'vbm_band_idx': vbm_band_idx, 'cbm_band_idx': cbm_band_idx,
        }

    def calculate_global_bandgap(self):
        global_vbm = -float('inf')
        global_cbm = float('inf')
        global_vbm_info = {}
        global_cbm_info = {}
        all_direct_gaps = []

        spins = [Spin.up, Spin.down] if self.is_spin_polarized else [Spin.up]

        for spin in spins:
            edges = self.find_band_edges(spin)
            label = "up" if spin == Spin.up else "down"
            if edges['vbm_energy'] > global_vbm:
                global_vbm = edges['vbm_energy']
                global_vbm_info = {**edges, 'spin': label}
            if edges['cbm_energy'] < global_cbm:
                global_cbm = edges['cbm_energy']
                global_cbm_info = {**edges, 'spin': label}
            direct_gaps = self.calculate_direct_gaps(spin)
            for g in direct_gaps:
                g['spin'] = label
            all_direct_gaps.extend(direct_gaps)

        fundamental_gap = global_cbm - global_vbm
        same_kpoint = np.allclose(global_vbm_info['vbm_kpoint'], global_cbm_info['cbm_kpoint'], atol=1e-6)
        fundamental_type = "direct" if same_kpoint else "indirect"

        min_direct_gap_info = min(all_direct_gaps, key=lambda x: x['gap']) if all_direct_gaps else None
        min_direct_gap = min_direct_gap_info['gap'] if min_direct_gap_info else None

        if fundamental_type == "direct":
            indirect_gap = None
            indirect_gap_info = {}
            for spin in spins:
                edges = self.find_band_edges(spin)
                label = "up" if spin == Spin.up else "down"
                if not np.allclose(edges['vbm_kpoint'], global_cbm_info['cbm_kpoint'], atol=1e-6):
                    gap = global_cbm - edges['vbm_energy']
                    if indirect_gap is None or gap < indirect_gap:
                        indirect_gap = gap
                        indirect_gap_info = {
                            'gap': gap, 'vbm_energy': edges['vbm_energy'], 'cbm_energy': global_cbm,
                            'vbm_kpoint': edges['vbm_kpoint'].tolist(), 'cbm_kpoint': global_cbm_info['cbm_kpoint'],
                            'vbm_spin': label, 'cbm_spin': global_cbm_info['spin'],
                        }
                if not np.allclose(edges['cbm_kpoint'], global_vbm_info['vbm_kpoint'], atol=1e-6):
                    gap = edges['cbm_energy'] - global_vbm
                    if indirect_gap is None or gap < indirect_gap:
                        indirect_gap = gap
                        indirect_gap_info = {
                            'gap': gap, 'vbm_energy': global_vbm, 'cbm_energy': edges['cbm_energy'],
                            'vbm_kpoint': global_vbm_info['vbm_kpoint'], 'cbm_kpoint': edges['cbm_kpoint'].tolist(),
                            'vbm_spin': global_vbm_info['spin'], 'cbm_spin': label,
                        }
        else:
            indirect_gap = fundamental_gap
            indirect_gap_info = {
                'gap': indirect_gap, 'vbm_energy': global_vbm, 'cbm_energy': global_cbm,
                'vbm_kpoint': global_vbm_info['vbm_kpoint'], 'cbm_kpoint': global_cbm_info['cbm_kpoint'],
                'vbm_spin': global_vbm_info['spin'], 'cbm_spin': global_cbm_info['spin'],
            }

        return {
            'fundamental_gap': fundamental_gap, 'fundamental_type': fundamental_type,
            'direct_gap': min_direct_gap,
            'direct_gap_kpoint': min_direct_gap_info['kpoint'] if min_direct_gap_info else None,
            'direct_gap_kpoint_idx': min_direct_gap_info['kpoint_idx'] if min_direct_gap_info else None,
            'direct_gap_spin': min_direct_gap_info['spin'] if min_direct_gap_info else None,
            'indirect_gap': indirect_gap,
            'indirect_gap_info': indirect_gap_info if indirect_gap else None,
            'is_metal': fundamental_gap <= 0.0,
            'global_vbm_energy': global_vbm, 'global_cbm_energy': global_cbm,
            'vbm_spin': global_vbm_info['spin'], 'cbm_spin': global_cbm_info['spin'],
            'vbm_kpoint_idx': global_vbm_info['vbm_kpoint_idx'],
            'cbm_kpoint_idx': global_cbm_info['cbm_kpoint_idx'],
            'vbm_kpoint': global_vbm_info['vbm_kpoint'], 'cbm_kpoint': global_cbm_info['cbm_kpoint'],
            'vbm_band_idx': global_vbm_info['vbm_band_idx'],
            'cbm_band_idx': global_cbm_info['cbm_band_idx'],
        }

    def calculate_direct_gaps(self, spin=Spin.up):
        eigenvals = self.eigenvalues[spin]  # type: ignore
        direct_gaps = []
        for k_idx, kpoint_data in enumerate(eigenvals):
            occupied = [e for e, o in kpoint_data if o > 0.5]
            unoccupied = [e for e, o in kpoint_data if o <= 0.5]
            if occupied and unoccupied:
                local_vbm = max(occupied)
                local_cbm = min(unoccupied)
                direct_gaps.append({
                    'gap': local_cbm - local_vbm, 'kpoint_idx': k_idx,
                    'kpoint': self.kpoints[k_idx], 'vbm': local_vbm, 'cbm': local_cbm,
                })
        return direct_gaps

    def calculate_spin_resolved_gaps(self):
        results = {}
        spins = [Spin.up, Spin.down] if self.is_spin_polarized else [Spin.up]
        for spin in spins:
            label = "up" if spin == Spin.up else "down"
            edges = self.find_band_edges(spin)
            direct_gaps = self.calculate_direct_gaps(spin)
            min_dg = min(direct_gaps, key=lambda x: x['gap']) if direct_gaps else None
            fg = edges['cbm_energy'] - edges['vbm_energy']
            same_k = np.allclose(edges['vbm_kpoint'], edges['cbm_kpoint'], atol=1e-6)
            results[label] = {
                'fundamental_gap': fg, 'gap_type': "direct" if same_k else "indirect",
                'is_metal': fg <= 0.0, 'vbm_energy': edges['vbm_energy'], 'cbm_energy': edges['cbm_energy'],
                'vbm_kpoint': edges['vbm_kpoint'], 'cbm_kpoint': edges['cbm_kpoint'],
                'vbm_kpoint_idx': edges['vbm_kpoint_idx'], 'cbm_kpoint_idx': edges['cbm_kpoint_idx'],
                'vbm_band_idx': edges['vbm_band_idx'], 'cbm_band_idx': edges['cbm_band_idx'],
                'min_direct_gap': min_dg['gap'] if min_dg else None,
                'min_direct_gap_kpoint': min_dg['kpoint'] if min_dg else None,
                'min_direct_gap_kpoint_idx': min_dg['kpoint_idx'] if min_dg else None,
                'all_direct_gaps': direct_gaps,
            }
        return results

    def analyze(self):
        global_info = self.calculate_global_bandgap()
        spin_resolved = self.calculate_spin_resolved_gaps()
        system_info = {
            'fermi_level': self.fermi_level,
            'is_spin_polarized': self.is_spin_polarized,
            'n_kpoints': len(self.kpoints),
            'n_bands': len(self.eigenvalues[Spin.up][0]),  # type: ignore
        }
        return {'system_info': system_info, 'global_bandgap': global_info, 'spin_resolved': spin_resolved}


def analyze_bandgap(vasprun_path="vasprun.xml"):
    """便捷函数：分析带隙并返回结果字典。"""
    analyzer = BandGapAnalyzer(vasprun_path)
    return analyzer.analyze()
