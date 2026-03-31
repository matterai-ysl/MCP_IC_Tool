"""电子结构相关分析技能"""
import re
from pathlib import Path
from typing import Any, Dict, List, Optional


def is_converged(work_dir: str) -> bool:
    """检查 VASP 计算是否正常完成"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return False
    return "General timing and accounting informations for this job:" in outcar.read_text(errors="ignore")


def get_fermi_energy(work_dir: str) -> Optional[float]:
    """从 OUTCAR 提取最终费米能级 (eV)"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    content = outcar.read_text(errors="ignore")
    matches = re.findall(r"E-fermi\s*:\s*([-\d.]+)", content)
    return float(matches[-1]) if matches else None


def get_band_gap(work_dir: str) -> Optional[Dict[str, Any]]:
    """从 vasprun.xml 提取带隙信息 {energy_ev, direct, transition, vbm, cbm}"""
    vr_path = Path(work_dir) / "vasprun.xml"
    if not vr_path.exists():
        return None
    try:
        from pymatgen.io.vasp import Vasprun
        vr = Vasprun(str(vr_path), parse_dos=False, parse_eigen=True)
        bs = vr.get_band_structure()
        gap_info = bs.get_band_gap()
        vbm = bs.get_vbm()
        cbm = bs.get_cbm()
        return {
            "energy_ev": round(gap_info.get("energy", 0.0), 4),
            "direct": gap_info.get("direct", False),
            "transition": gap_info.get("transition"),
            "vbm_ev": round(vbm["energy"], 4) if vbm.get("energy") is not None else None,
            "cbm_ev": round(cbm["energy"], 4) if cbm.get("energy") is not None else None,
            "is_metal": bs.is_metal(),
        }
    except Exception as e:
        return {"error": str(e)}


def get_total_energy_from_vasprun(work_dir: str) -> Optional[float]:
    """从 vasprun.xml 提取最终总能量 (eV)"""
    vr_path = Path(work_dir) / "vasprun.xml"
    if not vr_path.exists():
        return None
    try:
        from pymatgen.io.vasp import Vasprun
        vr = Vasprun(str(vr_path), parse_dos=False, parse_eigen=False)
        return vr.final_energy
    except Exception:
        return None


def get_electronic_steps(work_dir: str) -> Optional[Dict[str, Any]]:
    """统计电子步信息：最后离子步的电子步数，是否收敛"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    content = outcar.read_text(errors="ignore")
    # 每个离子步的电子步数
    all_steps = [int(x) for x in re.findall(r"Iteration\s+\d+\(\s*(\d+)\)", content)]
    if not all_steps:
        return None
    return {
        "last_ionic_step_elec_steps": all_steps[-1],
        "total_elec_steps": len(all_steps),
        "converged": is_converged(work_dir),
    }


def get_magnetization(work_dir: str) -> Optional[Dict[str, Any]]:
    """提取总磁矩（自旋极化计算）"""
    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return None
    content = outcar.read_text(errors="ignore")
    matches = re.findall(r"number of electron\s+[\d.]+\s+magnetization\s+([-\d.]+)", content)
    if not matches:
        return None
    values = [float(v) for v in matches]
    return {"final_magnetization": values[-1], "all_ionic_steps": values}
