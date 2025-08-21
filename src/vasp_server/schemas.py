from pydantic import BaseModel, Field, HttpUrl
from typing import Optional, Literal, Union
from enum import Enum

class CalcType(str, Enum):
    """计算类型枚举"""
    OXC = "OXC"  # 氧化物/硫化物固体电解质
    SSE = "SSE"  # 固体电解质（等同于OXC）
    ORC = "ORC"  # 氧化物还原催化剂
    ECAT_OER = "ECAT_OER"  # 氧析出反应催化剂
    ECAT_HER = "ECAT_HER"  # 氢析出反应催化剂

class SelectionMode(str, Enum):
    """材料选择模式"""
    auto = "auto"  # 自动选择最稳定的
    stable = "stable"  # 只选择稳定材料
    most_stable = "most_stable"  # 选择最稳定的
    first = "first"  # 选择第一个匹配的

class StructOptRequest(BaseModel):
    """结构优化请求模型"""
    # 用户ID（必填）
    user_id: Optional[str] = Field(None, description="用户ID")
    
    # 输入源（二选一）
    formula: Optional[str] = Field(None, description="化学式，如 'Li2O', 'LiFePO4'")
    cif_url: Optional[HttpUrl] = Field(None, description="CIF文件的URL地址")
    
    # 计算类型（必填）
    calc_type: CalcType = Field(..., description="计算类型")
    
    # 材料搜索参数（仅当使用formula时有效）
    spacegroup: Optional[str] = Field(None, description="空间群符号")
    max_energy_above_hull: Optional[float] = Field(0.1, description="最大能量上凸包距离 (eV/atom)")
    min_band_gap: Optional[float] = Field(None, description="最小带隙 (eV)")
    max_band_gap: Optional[float] = Field(None, description="最大带隙 (eV)")
    max_nsites: Optional[int] = Field(None, description="最大原子数")
    min_nsites: Optional[int] = Field(None, description="最小原子数")
    stable_only: bool = Field(True, description="只选择稳定材料")
    selection_mode: SelectionMode = Field(SelectionMode.auto, description="选择模式")
    
    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    
    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        if not self.formula and not self.cif_url:
            raise ValueError("必须提供 formula 或 cif_url 中的一个")
        if self.formula and self.cif_url:
            raise ValueError("不能同时提供 formula 和 cif_url")

class TaskStatus(str, Enum):
    """任务状态"""
    queued = "queued"
    running = "running"
    completed = "completed"
    failed = "failed"
    canceled = "canceled"
    canceling = "canceling"

class StructOptResponse(BaseModel):
    """结构优化响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class TaskStatusResponse(BaseModel):
    """任务状态查询响应模型"""
    task_id: str
    user_id: str
    task_type: str
    status: TaskStatus
    progress: int = Field(..., description="进度百分比 (0-100)")
    params: Optional[dict] = None
    result_path: Optional[str] = None
    external_job_id: Optional[str] = None
    process_id: Optional[int] = Field(None, description="VASP进程ID")
    error_message: Optional[str] = None
    created_at: str
    updated_at: str

class StructOptResult(BaseModel):
    """结构优化结果模型"""
    optimized_structure: Optional[str] = Field(None, description="优化后的结构文件路径")
    energy: Optional[float] = Field(None, description="总能量 (eV)")
    final_forces: Optional[list] = Field(None, description="最终力矩阵")
    convergence: bool = Field(False, description="是否收敛")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")

class SCFRequest(BaseModel):
    """自洽场计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（三选一）
    formula: Optional[str] = Field(None, description="化学式，如 'Li2O', 'LiFePO4'")
    cif_url: Optional[HttpUrl] = Field(None, description="CIF文件的URL地址")
    optimized_task_id: Optional[str] = Field(None, description="已完成的结构优化任务ID")
    
    # 计算类型（必填）
    calc_type: CalcType = Field(..., description="计算类型")
    
    # 材料搜索参数（仅当使用formula时有效）
    spacegroup: Optional[str] = Field(None, description="空间群符号")
    max_energy_above_hull: Optional[float] = Field(0.1, description="最大能量上凸包距离 (eV/atom)")
    min_band_gap: Optional[float] = Field(None, description="最小带隙 (eV)")
    max_band_gap: Optional[float] = Field(None, description="最大带隙 (eV)")
    max_nsites: Optional[int] = Field(None, description="最大原子数")
    min_nsites: Optional[int] = Field(None, description="最小原子数")
    stable_only: bool = Field(True, description="只选择稳定材料")
    selection_mode: SelectionMode = Field(SelectionMode.auto, description="选择模式")
    
    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    precision: str = Field("Accurate", description="计算精度 (Normal, High, Accurate)")
    
    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([
            bool(self.formula),
            bool(self.cif_url), 
            bool(self.optimized_task_id)
        ])
        if input_count != 1:
            raise ValueError("必须提供 formula、cif_url 或 optimized_task_id 中的一个")

class SCFResponse(BaseModel):
    """自洽场计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class SCFResult(BaseModel):
    """自洽场计算结果模型"""
    scf_structure: Optional[str] = Field(None, description="自洽场计算的结构文件路径")
    total_energy: Optional[float] = Field(None, description="总能量 (eV)")
    fermi_energy: Optional[float] = Field(None, description="费米能级 (eV)")
    band_gap: Optional[float] = Field(None, description="带隙 (eV)")
    convergence: bool = Field(False, description="是否收敛")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    electronic_steps: Optional[int] = Field(None, description="电子步数")

class DOSRequest(BaseModel):
    """态密度计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（三选一）
    formula: Optional[str] = Field(None, description="化学式，如 'Li2O', 'LiFePO4'")
    cif_url: Optional[HttpUrl] = Field(None, description="CIF文件的URL地址")
    scf_task_id: Optional[str] = Field(None, description="已完成的自洽场计算任务ID")
    
    # 计算类型（必填）
    calc_type: CalcType = Field(..., description="计算类型")
    
    # 材料搜索参数（仅当使用formula时有效）
    spacegroup: Optional[str] = Field(None, description="空间群符号")
    max_energy_above_hull: Optional[float] = Field(0.1, description="最大能量上凸包距离 (eV/atom)")
    min_band_gap: Optional[float] = Field(None, description="最小带隙 (eV)")
    max_band_gap: Optional[float] = Field(None, description="最大带隙 (eV)")
    max_nsites: Optional[int] = Field(None, description="最大原子数")
    min_nsites: Optional[int] = Field(None, description="最小原子数")
    stable_only: bool = Field(True, description="只选择稳定材料")
    selection_mode: SelectionMode = Field(SelectionMode.auto, description="选择模式")
    
    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    kpoint_multiplier: float = Field(2.0, description="K点倍增因子 (相对于优化计算)")
    precision: str = Field("Accurate", description="计算精度 (Normal, High, Accurate)")
    
    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([
            bool(self.formula),
            bool(self.cif_url), 
            bool(self.scf_task_id)
        ])
        if input_count != 1:
            raise ValueError("必须提供 formula、cif_url 或 scf_task_id 中的一个")

class DOSResponse(BaseModel):
    """态密度计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class DOSResult(BaseModel):
    """态密度计算结果模型"""
    dos_structure: Optional[str] = Field(None, description="态密度计算的结构文件路径")
    doscar_path: Optional[str] = Field(None, description="DOSCAR文件路径")
    total_energy: Optional[float] = Field(None, description="总能量 (eV)")
    fermi_energy: Optional[float] = Field(None, description="费米能级 (eV)")
    band_gap: Optional[float] = Field(None, description="带隙 (eV)")
    dos_data: Optional[dict] = Field(None, description="态密度数据")
    convergence: bool = Field(False, description="是否收敛")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    kpoints_used: Optional[list] = Field(None, description="使用的K点网格")

class MDRequest(BaseModel):
    """分子动力学计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（三选一）
    formula: Optional[str] = Field(None, description="化学式，如 'Li2O', 'LiFePO4'")
    cif_url: Optional[HttpUrl] = Field(None, description="CIF文件的URL地址")
    scf_task_id: Optional[str] = Field(None, description="已完成的自洽场计算任务ID")
    
    # 计算类型（必填）
    calc_type: CalcType = Field(..., description="计算类型")
    
    # 材料搜索参数（仅当使用formula时有效）
    spacegroup: Optional[str] = Field(None, description="空间群符号")
    max_energy_above_hull: Optional[float] = Field(0.1, description="最大能量上凸包距离 (eV/atom)")
    min_band_gap: Optional[float] = Field(None, description="最小带隙 (eV)")
    max_band_gap: Optional[float] = Field(None, description="最大带隙 (eV)")
    max_nsites: Optional[int] = Field(None, description="最大原子数")
    min_nsites: Optional[int] = Field(None, description="最小原子数")
    stable_only: bool = Field(True, description="只选择稳定材料")
    selection_mode: SelectionMode = Field(SelectionMode.auto, description="选择模式")
    
    # MD计算参数
    md_steps: int = Field(1000, description="MD步数")
    temperature: float = Field(300.0, description="目标温度 (K)")
    time_step: float = Field(1.0, description="时间步长 (fs)")
    ensemble: str = Field("NVT", description="系综类型 (NVT, NVE, NPT)")
    precision: str = Field("Normal", description="计算精度 (Normal, High, Accurate)")
    
    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([
            bool(self.formula),
            bool(self.cif_url), 
            bool(self.scf_task_id)
        ])
        if input_count != 1:
            raise ValueError("必须提供 formula、cif_url 或 scf_task_id 中的一个")

class MDResponse(BaseModel):
    """分子动力学计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class MDResult(BaseModel):
    """分子动力学计算结果模型"""
    md_structure: Optional[str] = Field(None, description="MD计算的初始结构文件路径")
    xdatcar_path: Optional[str] = Field(None, description="XDATCAR轨迹文件路径")
    oszicar_path: Optional[str] = Field(None, description="OSZICAR能量文件路径")
    final_energy: Optional[float] = Field(None, description="最终能量 (eV)")
    average_temperature: Optional[float] = Field(None, description="平均温度 (K)")
    total_md_steps: Optional[int] = Field(None, description="完成的MD步数")
    convergence: bool = Field(False, description="是否正常完成")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    trajectory_data: Optional[dict] = Field(None, description="轨迹统计数据") 