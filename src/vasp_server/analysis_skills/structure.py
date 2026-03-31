"""结构相关分析技能"""
from pathlib import Path
from typing import Any, Dict, Optional


def get_final_structure(work_dir: str):
    """获取最终结构 (pymatgen Structure)，优先 CONTCAR，其次 POSCAR"""
    from pymatgen.core import Structure
    for fname in ("CONTCAR", "POSCAR"):
        p = Path(work_dir) / fname
        if p.exists() and p.stat().st_size > 0:
            try:
                return Structure.from_file(str(p))
            except Exception:
                continue
    return None


def get_volume(work_dir: str) -> Optional[float]:
    """晶胞体积 (Å³)"""
    s = get_final_structure(work_dir)
    return float(s.volume) if s else None


def get_lattice_params(work_dir: str) -> Optional[Dict[str, float]]:
    """晶格参数 (a, b, c, alpha, beta, gamma, volume)"""
    s = get_final_structure(work_dir)
    if not s:
        return None
    lat = s.lattice
    return {
        "a": round(lat.a, 6), "b": round(lat.b, 6), "c": round(lat.c, 6),
        "alpha": round(lat.alpha, 4), "beta": round(lat.beta, 4), "gamma": round(lat.gamma, 4),
        "volume": round(lat.volume, 4),
    }


def get_formula(work_dir: str) -> Optional[str]:
    """化学式（约化）"""
    s = get_final_structure(work_dir)
    return s.composition.reduced_formula if s else None


def get_n_atoms(work_dir: str) -> Optional[int]:
    """结构中的原子总数"""
    s = get_final_structure(work_dir)
    return len(s) if s else None


def get_volume_change(work_dir: str) -> Optional[Dict[str, Any]]:
    """比较初始(POSCAR)与最终(CONTCAR)体积变化"""
    from pymatgen.core import Structure
    results: Dict[str, Any] = {}
    for label, fname in (("initial", "POSCAR"), ("final", "CONTCAR")):
        p = Path(work_dir) / fname
        if p.exists() and p.stat().st_size > 0:
            try:
                s = Structure.from_file(str(p))
                results[label] = round(s.volume, 4)
            except Exception:
                pass
    if "initial" in results and "final" in results:
        results["change_pct"] = round(
            (results["final"] - results["initial"]) / results["initial"] * 100, 3
        )
    return results or None
