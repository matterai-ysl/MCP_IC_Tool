from pydantic import BaseModel, Field, HttpUrl
from typing import Optional, Dict, Any, List, Union


class StrictInputModel(BaseModel):
    model_config = {"extra": "forbid"}


class StructOptInput(StrictInputModel):
    user_id: Optional[str] = Field(default=None)
    cif_url: HttpUrl
    queue_name: Optional[str] = None
    kpoint_density: float = 15.0
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class SCFInput(StrictInputModel):
    user_id: str
    cif_url: Optional[HttpUrl] = None
    optimized_task_id: Optional[str] = None
    queue_name: Optional[str] = None
    kpoint_density: float = 30.0
    precision: str = "Accurate"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class DOSInput(StrictInputModel):
    user_id: str
    input_url: Optional[Union[HttpUrl, List[HttpUrl]]] = None
    scf_task_id: Optional[str] = None
    queue_name: Optional[str] = None
    kpoint_density: float = 20.0
    kpoint_multiplier: float = 1.5
    precision: str = "High"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class BandStructureInput(StrictInputModel):
    user_id: str
    input_url: Optional[Union[HttpUrl, List[HttpUrl]]] = None
    scf_task_id: Optional[str] = None
    queue_name: Optional[str] = None
    kpoint_density: float = 20.0
    line_density: int = 20
    precision: str = "High"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class MDInput(StrictInputModel):
    user_id: str
    input_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    queue_name: Optional[str] = None
    md_steps: int = 1000
    temperature: float = Field(default=300.0, description="目标温度(K)")
    time_step: float = 1.0
    ensemble: str = "NVT"
    precision: str = "Normal"
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class NEBInput(StrictInputModel):
    user_id: str
    # 初始结构（二选一）
    initial_cif_url: Optional[HttpUrl] = None
    initial_task_id: Optional[str] = None
    # 终态结构（二选一）
    final_cif_url: Optional[HttpUrl] = None
    final_task_id: Optional[str] = None
    # NEB 参数
    n_images: int = Field(default=5, description="中间图像数（不含端点），2-20")
    queue_name: Optional[str] = None
    kpoint_density: float = 15.0
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class PhononInput(StrictInputModel):
    user_id: str
    cif_url: Optional[HttpUrl] = None
    scf_task_id: Optional[str] = None
    queue_name: Optional[str] = None
    kpoint_density: float = 15.0
    displacement: float = Field(default=0.015, description="有限位移步长（Å）")
    custom_incar: Optional[Dict[str, Any]] = Field(default=None, description="自定义INCAR参数字典")


class CustomCalcInput(StrictInputModel):
    """通用自定义计算输入 — 用户完全控制INCAR"""
    user_id: str
    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = None
    from_task_id: Optional[str] = Field(default=None, description="复用已完成任务的结构")
    queue_name: Optional[str] = None
    # INCAR — 用户完全控制
    incar: Dict[str, Any] = Field(..., description="完整的INCAR参数字典")
    # K点设置
    kpoint_density: float = 30.0
    kpoint_mode: str = Field(default="mesh", description="K点模式: mesh 或 gamma")
