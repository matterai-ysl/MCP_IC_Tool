"""
共享的 VASP 文件发现和结构解析工具。

从 optimization_analyzer.py 和 scf_analyzer.py 的重复代码中提取。
"""

import re
from pathlib import Path
from typing import Dict, List, Optional


def find_outcar_file(work_dir: Path) -> Optional[Path]:
    """智能寻找 OUTCAR 文件，兼容不同后缀和大小写。"""
    return find_vasp_file(work_dir, "OUTCAR")


def find_vasp_file(work_dir: Path, base_name: str) -> Optional[Path]:
    """
    智能寻找 VASP 文件，兼容不同后缀和大小写。

    尝试顺序: POSCAR → POSCAR.txt → poscar → poscar.txt → Poscar → Poscar.txt
    """
    possible_names = [
        base_name,
        f"{base_name}.txt",
        base_name.lower(),
        f"{base_name.lower()}.txt",
        base_name.capitalize(),
        f"{base_name.capitalize()}.txt",
    ]
    for name in possible_names:
        file_path = work_dir / name
        if file_path.exists() and file_path.is_file():
            return file_path
    return None


def resolve_work_dir_and_outcar(input_path: str) -> tuple[Path, Path]:
    """
    从 input_path (文件或目录) 解析 work_dir 和 outcar_path。

    Returns:
        (work_dir, outcar_path)

    Raises:
        FileNotFoundError
    """
    p = Path(input_path)
    if p.is_file():
        return p.parent, p
    elif p.is_dir():
        outcar = find_outcar_file(p)
        if not outcar:
            tried = ["OUTCAR", "OUTCAR.txt", "outcar", "outcar.txt", "Outcar", "Outcar.txt"]
            raise FileNotFoundError(
                f"文件夹中未找到 OUTCAR 文件，尝试了: {', '.join(tried)}\n路径: {p}"
            )
        return p, outcar
    else:
        raise FileNotFoundError(f"输入路径不存在: {input_path}")


def extract_composition_from_poscar(poscar_content: str) -> str:
    """
    从 POSCAR 内容提取化学组成字符串，如 'Li2O3'。

    POSCAR 格式:
        Line 1: 注释
        Lines 2-5: 标度因子 + 晶胞矢量
        Line 6: 元素名 (可选, VASP5+)
        Line 7: 原子数目
    """
    lines = poscar_content.strip().split('\n')
    if len(lines) < 7:
        return "Unknown"

    # 第 6 行通常是元素符号
    elements_line = lines[5].strip()
    counts_line = lines[6].strip()

    # 检查第 6 行是否全是字母（元素符号行）
    elements_parts = elements_line.split()
    counts_parts = counts_line.split()

    if not elements_parts:
        return "Unknown"

    # 如果第 6 行全是数字，说明没有元素符号行 (VASP4 格式)
    if all(_is_number(p) for p in elements_parts):
        return "Unknown"

    # 正常 VASP5 格式
    if len(elements_parts) != len(counts_parts):
        return "Unknown"

    try:
        parts = []
        for elem, cnt_str in zip(elements_parts, counts_parts):
            cnt = int(cnt_str)
            parts.append(f"{elem}{cnt}" if cnt > 1 else elem)
        return "".join(parts)
    except (ValueError, IndexError):
        return "Unknown"


def parse_poscar_elements_and_counts(work_dir: Path) -> tuple[List[str], List[int]]:
    """
    从 POSCAR 读取元素名和原子数。

    Returns:
        (["Li", "O", ...], [4, 6, ...])
    """
    poscar_path = find_vasp_file(work_dir, "POSCAR")
    if not poscar_path:
        return [], []

    with open(poscar_path, 'r', encoding='utf-8', errors='ignore') as f:
        lines = f.readlines()

    if len(lines) < 7:
        return [], []

    elements = lines[5].strip().split()
    counts_str = lines[6].strip().split()

    if all(_is_number(p) for p in elements):
        # VASP4 格式，没有元素符号行
        return [], []

    try:
        counts = [int(c) for c in counts_str]
    except ValueError:
        return [], []

    return elements, counts


def read_structure_files(work_dir: Path) -> Dict[str, str]:
    """读取 POSCAR 和 CONTCAR 的内容。"""
    result: Dict[str, str] = {}
    for name in ("POSCAR", "CONTCAR"):
        path = find_vasp_file(work_dir, name)
        if path and path.exists():
            try:
                result[name.lower()] = path.read_text(encoding='utf-8', errors='ignore')
            except Exception:
                pass
    return result


def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False
