from enum import Enum
from pydantic import BaseModel, Field, HttpUrl
from typing import Optional, Dict, Any, Union, List


class CalcType(str, Enum):
    OXC = "OXC"
    SSE = "SSE"
    ORC = "ORC"
    ECAT_OER = "ECAT_OER"
    ECAT_HER = "ECAT_HER"


class SelectionMode(str, Enum):
    auto = "auto"
    stable = "stable"
    most_stable = "most_stable"
    first = "first"


class StructOptInput(BaseModel):
    user_id: Optional[str] = Field(default=None)
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    calc_type: CalcType
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
    calc_type: CalcType
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
    calc_type: CalcType
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


class MDInput(BaseModel):
    user_id: str
    formula: Optional[str] = None
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    calc_type: CalcType
    spacegroup: Optional[str] = None
    max_energy_above_hull: Optional[float] = Field(default=0.1)
    min_band_gap: Optional[float] = None
    max_band_gap: Optional[float] = None
    max_nsites: Optional[int] = None
    min_nsites: Optional[int] = None
    stable_only: bool = True
    selection_mode: SelectionMode = SelectionMode.auto
    md_steps: int = 1000
    temperature: Union[float, List[float]] = Field(default=300.0, description="目标温度(K) - 支持单个温度或温度列表进行多温度扫描")
    time_step: float = 1.0
    ensemble: str = "NVT"
    precision: str = "Normal"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


