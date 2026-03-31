from enum import Enum
from pydantic import BaseModel, Field, HttpUrl
from typing import Optional, Dict, Any, List


class SelectionMode(str, Enum):
    auto = "auto"
    stable = "stable"
    most_stable = "most_stable"
    first = "first"


class StructOptInput(BaseModel):
    user_id: Optional[str] = Field(default=None)
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    kpoint_density: float = 30.0
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class SCFInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    optimized_task_id: Optional[str] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    kpoint_density: float = 30.0
    precision: str = "Accurate"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class DOSInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    kpoint_density: float = 30.0
    kpoint_multiplier: float = 2.0
    precision: str = "Accurate"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class BandStructureInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    kpoint_density: float = 30.0
    line_density: int = 20
    precision: str = "Accurate"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class MDInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    md_steps: int = 1000
    temperature: float = Field(default=300.0, description="目标温度(K)")
    time_step: float = 1.0
    ensemble: str = "NVT"
    precision: str = "Normal"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class NEBInput(BaseModel):
    user_id: str
    # 初始结构（三选一）
    initial_formula: Optional[str] = None
    initial_cif_url: Optional[HttpUrl] = None
    initial_task_id: Optional[str] = None
    # 终态结构（三选一）
    final_formula: Optional[str] = None
    final_cif_url: Optional[HttpUrl] = None
    final_task_id: Optional[str] = None
    # 结构搜索参数
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    # NEB 参数
    n_images: int = Field(default=5, description="中间图像数（不含端点），2-20")
    kpoint_density: float = 30.0
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class PhononInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    kpoint_density: float = 30.0
    displacement: float = Field(default=0.015, description="有限位移步长（Å）")
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class CustomCalcInput(BaseModel):
    """通用自定义计算输入 — 用户完全控制INCAR"""
    user_id: str
    # 输入源（三选一）
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    from_task_id: Optional[str] = Field(default=None, description="复用已完成任务的结构")
    # 材料搜索参数（formula时有效）
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    # INCAR — 用户完全控制
    incar: Dict[str, Any] = Field(..., description="完整的INCAR参数字典")
    # K点设置
    kpoint_density: float = 30.0
    kpoint_mode: str = Field(default="mesh", description="K点模式: mesh 或 gamma")


