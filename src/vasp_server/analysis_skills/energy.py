"""能量相关分析技能"""
import re
from pathlib import Path
from typing import List, Optional


def extract_final_energy(work_dir: str) -> Optional[float]:
    """从 OUTCAR 提取最终总能量 (eV)"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    energy = None
    with open(outcar, errors="ignore") as f:
        for line in f:
            if "free energy    TOTEN" in line:
                try:
                    energy = float(line.split()[4])
                except (IndexError, ValueError):
                    pass
    return energy


def extract_energy_per_atom(work_dir: str) -> Optional[float]:
    """能量/原子数 (eV/atom)"""
    energy = extract_final_energy(work_dir)
    if energy is None:
        return None
    for fname in ("CONTCAR", "POSCAR"):
        p = Path(work_dir) / fname
        if p.exists() and p.stat().st_size > 0:
            try:
                lines = p.read_text(errors="ignore").splitlines()
                n_atoms = sum(int(x) for x in lines[6].split())
                return energy / n_atoms if n_atoms > 0 else None
            except Exception:
                continue
    return None


def get_ionic_steps_energy(work_dir: str) -> List[float]:
    """每个离子步结束时的总能量列表 (eV)"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return []
    energies: List[float] = []
    with open(outcar, errors="ignore") as f:
        for line in f:
            if "free energy    TOTEN" in line:
                try:
                    energies.append(float(line.split()[4]))
                except (IndexError, ValueError):
                    pass
    return energies


def get_max_force(work_dir: str) -> Optional[float]:
    """从 OUTCAR 获取最终离子步的最大原子力 (eV/Å)"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    content = outcar.read_text(errors="ignore")
    blocks = list(re.finditer(r"TOTAL-FORCE.*?\n(.*?)(?=\n-{10})", content, re.DOTALL))
    if not blocks:
        return None
    try:
        forces = []
        for line in blocks[-1].group(1).strip().splitlines():
            parts = line.split()
            if len(parts) == 6:
                fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                forces.append((fx**2 + fy**2 + fz**2) ** 0.5)
        return max(forces) if forces else None
    except Exception:
        return None


def get_computation_time(work_dir: str) -> Optional[float]:
    """从 OUTCAR 提取总计算时间 (秒)"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    content = outcar.read_text(errors="ignore")
    m = re.search(r"Total CPU time used \(sec\):\s*([\d.]+)", content)
    if m:
        return float(m.group(1))
    m = re.search(r"Elapsed time \(sec\):\s*([\d.]+)", content)
    return float(m.group(1)) if m else None
