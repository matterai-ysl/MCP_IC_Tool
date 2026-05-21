from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass(frozen=True)
class FailureGuidance:
    failure_type: Optional[str] = None
    reason: Optional[str] = None
    suggested_action: Optional[str] = None


ELECTRONIC_DIVERGENCE = "electronic_divergence"
IONIC_LINE_SEARCH_FAILURE = "ionic_line_search_failure"
LATTICE_INCONSISTENCY = "lattice_inconsistency"
INVALID_STRUCTURE_FILE = "invalid_structure_file"
MISSING_CHARGE_DENSITY = "missing_charge_density"
DISORDERED_STRUCTURE = "disordered_structure"
DISK_QUOTA_EXCEEDED = "disk_quota_exceeded"


_SCF_DIVERGENCE_ACTION = (
    "先对该结构做一次 SCF，确认电子收敛并生成稳定的 CHGCAR；如果这是原始上传结构、"
    "表面/slab、吸附体系或强极化体系，先做更稳健的预优化或结构优化，再继续 "
    "band_structure、DOS 或其他依赖 CHGCAR 的后续任务。"
)

_IONIC_LINE_SEARCH_ACTION = (
    "先将当前 CONTCAR 作为新的 POSCAR 续算；如果该结构还没有经过稳定 SCF 或预优化，"
    "先做一次 SCF 或轻量预优化，再重试结构优化。必要时减小 POTIM、收紧 EDIFF，"
    "或降低单步结构更新幅度。"
)

_MISSING_CHGCAR_ACTION = (
    "先对该结构执行一次稳定的 SCF，确保生成可复用的 CHGCAR；确认 SCF 已收敛后，"
    "再继续 band_structure、DOS 或其他依赖电荷密度的后续任务。"
)

_LATTICE_CONSISTENCY_ACTION = (
    "先规范化该结构的晶格参数和晶格向量，去掉应为 0 的极小噪声分量；"
    "如果结构来自外部导出文件，建议重新导出标准 CIF/POSCAR，必要时调大 SYMPREC，"
    "再重新提交结构优化。"
)

_INVALID_STRUCTURE_FILE_ACTION = (
    "请先确认上传的是有效的结构文件，而不是空文件、说明文本或缺少原子坐标的 CIF；"
    "优先重新导出标准 CIF/POSCAR 后再提交。"
)

_DISORDERED_STRUCTURE_ACTION = (
    "该 CIF 含有小数占位或混占晶位，不能直接转换成 VASP POSCAR。请先生成明确有序结构："
    "例如构建足够大的超胞后按化学计量选择具体原子、使用 pymatgen 的有序化/枚举占位流程，"
    "或从结构编辑软件导出已 ordered 的 POSCAR/CIF；确认所有站位占位为整数原子后再提交计算。"
)

_DISK_QUOTA_ACTION = (
    "HPC 账号或文件系统配额已经不足，VASP 写 WAVECAR/CHG/vasprun.xml 等大文件时失败。"
    "请先删除或归档旧的 vasp_calculations 数据，确认 quota 恢复后再重提；必要时关闭不需要的 "
    "LWAVE/LCHARG 输出或降低输出体积。"
)


def _read_tail_text(path: Path, max_bytes: int = 200_000) -> str:
    if not path.exists() or not path.is_file():
        return ""
    try:
        with path.open("rb") as fh:
            fh.seek(0, 2)
            size = fh.tell()
            fh.seek(max(size - max_bytes, 0))
            return fh.read().decode("utf-8", errors="ignore")
    except OSError:
        return ""


def derive_failure_guidance(error_text: str, task_type: Optional[str] = None) -> FailureGuidance:
    if not error_text:
        return FailureGuidance()

    normalized = error_text.upper()
    normalized_task_type = str(task_type or "").lower()

    if any(
        marker in normalized
        for marker in (
            "DISK QUOTA EXCEEDED",
            "NO SPACE LEFT ON DEVICE",
            "ERROR DURING WRITE",
        )
    ):
        return FailureGuidance(
            failure_type=DISK_QUOTA_EXCEEDED,
            reason="VASP 在写入大文件时遇到磁盘配额或剩余空间限制，计算已无法继续可靠产出结果。",
            suggested_action=_DISK_QUOTA_ACTION,
        )

    if any(
        marker in normalized
        for marker in (
            "BRMIX",
            "EDDDAV",
            "ZHEGV",
            "ORBITAL ORTHONORMALIZATION FAILED",
            "ZPOTRF",
            "ZTRTRI",
        )
    ):
        return FailureGuidance(
            failure_type=ELECTRONIC_DIVERGENCE,
            reason="电子自洽过程已经发散，当前结构很可能还没有经过稳定 SCF，或初始结构/参数导致电荷混合失稳。",
            suggested_action=_SCF_DIVERGENCE_ACTION,
        )

    if "ZBRENT" in normalized:
        return FailureGuidance(
            failure_type=IONIC_LINE_SEARCH_FAILURE,
            reason="离子优化的线搜索失败，当前结构沿优化方向没有找到稳定的能量下降区间。",
            suggested_action=_IONIC_LINE_SEARCH_ACTION,
        )

    if any(
        marker in normalized
        for marker in (
            "HNFORM",
            "INCONSISTENT BRAVAIS LATTICE TYPES",
            "RECIPROCAL LATTICE ARE INCOMMENSURATE",
        )
    ):
        return FailureGuidance(
            failure_type=LATTICE_INCONSISTENCY,
            reason="当前结构的晶格向量存在不一致或数值噪声，VASP 无法为该晶胞构造稳定的 k 点网格。",
            suggested_action=_LATTICE_CONSISTENCY_ACTION,
        )

    if any(
        marker in normalized
        for marker in (
            "DISORDERED STRUCTURE WITH PARTIAL OCCUPANCIES",
            "PARTIAL OCCUPANCIES CANNOT BE CONVERTED INTO POSCAR",
            "PARTIAL OCCUPANC",
            "DISORDERED STRUCTURE",
            "小数占位",
            "混占晶位",
        )
    ):
        return FailureGuidance(
            failure_type=DISORDERED_STRUCTURE,
            reason="当前 CIF 包含小数占位或混占晶位，无法直接写成 VASP 需要的确定性 POSCAR。",
            suggested_action=_DISORDERED_STRUCTURE_ACTION,
        )

    if any(
        marker in normalized
        for marker in (
            "INVALID CIF FILE WITH NO STRUCTURES",
            "_ATOM_SITE_LABEL",
            "无法从 CIF 解析结构".upper(),
        )
    ):
        return FailureGuidance(
            failure_type=INVALID_STRUCTURE_FILE,
            reason="当前上传的结构文件不完整或不符合常见 CIF 约定，解析阶段无法得到可计算的原子结构。",
            suggested_action=_INVALID_STRUCTURE_FILE_ACTION,
        )

    downstream_task_types = {
        "band_structure",
        "band_structure_analysis",
        "dos_calculation",
        "dos_analysis",
    }
    if (
        "未生成有效的 CHGCAR" in error_text
        or "FAILED TO GENERATE A VALID CHGCAR" in normalized
        or (
            normalized_task_type in downstream_task_types
            and "CHGCAR" in normalized
            and ("SCF" in normalized or "CONVERG" in normalized or "FAILED" in normalized)
        )
    ):
        return FailureGuidance(
            failure_type=MISSING_CHARGE_DENSITY,
            reason="当前任务缺少稳定的电荷密度输入，后续 NSCF/DOS 分析无法直接继续。",
            suggested_action=_MISSING_CHGCAR_ACTION,
        )

    return FailureGuidance()


def guidance_for_failure_type(failure_type: Optional[str]) -> FailureGuidance:
    if failure_type == ELECTRONIC_DIVERGENCE:
        return FailureGuidance(
            failure_type=ELECTRONIC_DIVERGENCE,
            reason="电子自洽过程已经发散，当前结构很可能还没有经过稳定 SCF，或初始结构/参数导致电荷混合失稳。",
            suggested_action=_SCF_DIVERGENCE_ACTION,
        )
    if failure_type == IONIC_LINE_SEARCH_FAILURE:
        return FailureGuidance(
            failure_type=IONIC_LINE_SEARCH_FAILURE,
            reason="离子优化的线搜索失败，当前结构沿优化方向没有找到稳定的能量下降区间。",
            suggested_action=_IONIC_LINE_SEARCH_ACTION,
        )
    if failure_type == LATTICE_INCONSISTENCY:
        return FailureGuidance(
            failure_type=LATTICE_INCONSISTENCY,
            reason="当前结构的晶格向量存在不一致或数值噪声，VASP 无法为该晶胞构造稳定的 k 点网格。",
            suggested_action=_LATTICE_CONSISTENCY_ACTION,
        )
    if failure_type == INVALID_STRUCTURE_FILE:
        return FailureGuidance(
            failure_type=INVALID_STRUCTURE_FILE,
            reason="当前上传的结构文件不完整或不符合常见 CIF 约定，解析阶段无法得到可计算的原子结构。",
            suggested_action=_INVALID_STRUCTURE_FILE_ACTION,
        )
    if failure_type == MISSING_CHARGE_DENSITY:
        return FailureGuidance(
            failure_type=MISSING_CHARGE_DENSITY,
            reason="当前任务缺少稳定的电荷密度输入，后续 NSCF/DOS 分析无法直接继续。",
            suggested_action=_MISSING_CHGCAR_ACTION,
        )
    if failure_type == DISORDERED_STRUCTURE:
        return FailureGuidance(
            failure_type=DISORDERED_STRUCTURE,
            reason="当前 CIF 包含小数占位或混占晶位，无法直接写成 VASP 需要的确定性 POSCAR。",
            suggested_action=_DISORDERED_STRUCTURE_ACTION,
        )
    if failure_type == DISK_QUOTA_EXCEEDED:
        return FailureGuidance(
            failure_type=DISK_QUOTA_EXCEEDED,
            reason="VASP 在写入大文件时遇到磁盘配额或剩余空间限制，计算已无法继续可靠产出结果。",
            suggested_action=_DISK_QUOTA_ACTION,
        )
    return FailureGuidance()


def derive_failure_guidance_from_workdir(work_dir: Path, task_type: Optional[str] = None) -> FailureGuidance:
    fragments: list[str] = []
    for candidate in ("result.log", "OUTCAR"):
        content = _read_tail_text(work_dir / candidate)
        if content:
            fragments.append(content)

    for error_file in sorted(work_dir.glob("*.err")):
        content = _read_tail_text(error_file)
        if content:
            fragments.append(content)
            break

    return derive_failure_guidance("\n".join(fragments), task_type=task_type)
