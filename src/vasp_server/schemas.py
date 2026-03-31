from pydantic import BaseModel, Field, HttpUrl, field_validator
from typing import Optional, Literal, Union, Dict, Any, List
from enum import Enum
from datetime import datetime

TASK_STATES = {
    "queued",
    "leased",
    "running",
    "uploading",
    "analyzing",
    "completed",
    "failed",
    "cancel_requested",
    "canceled",
}

# custom_incar 中禁止覆盖的关键参数（会破坏工作流或文件路径逻辑）
_INCAR_KEY_BLACKLIST = frozenset({
    "SYSTEM",     # 内部用于标识计算类型
    "ISTART",     # 由工作流控制（续算 vs 新算）
    "ICHARG",     # 同上
    "IBRION",     # 由计算类型决定（优化/MD/SCF）
    "NSW",        # 由计算类型决定
    "LCHARG",     # 输出控制由分析流程依赖
    "LWAVE",      # 同上
})

def _validate_custom_incar(v: Dict[str, Any] | None) -> Dict[str, Any] | None:
    """校验 custom_incar 是否包含黑名单参数"""
    if v is None:
        return v
    blocked = {k for k in v if k.upper() in _INCAR_KEY_BLACKLIST}
    if blocked:
        raise ValueError(
            f"custom_incar 中不允许覆盖以下参数（由工作流自动管理）: {', '.join(sorted(blocked))}"
        )
    return v


class StrictRequestModel(BaseModel):
    model_config = {"extra": "forbid"}


class StructOptRequest(StrictRequestModel):
    """结构优化请求模型"""
    user_id: Optional[str] = Field(None, description="用户ID")
    cif_url: HttpUrl = Field(..., description="结构文件URL")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")
    kpoint_density: float = Field(30.0, description="K点密度参数")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典，用于覆盖默认参数")
    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

class TaskStatus(str, Enum):
    """任务状态"""
    queued = "queued"
    leased = "leased"
    running = "running"
    uploading = "uploading"
    analyzing = "analyzing"
    completed = "completed"
    failed = "failed"
    cancel_requested = "cancel_requested"
    canceled = "canceled"


class AnalysisStatus(str, Enum):
    """分析状态"""
    pending = "pending"
    running = "running"
    completed = "completed"
    failed = "failed"


class ArtifactInfo(BaseModel):
    """Artifact 摘要 (对外返回)"""
    id: str
    artifact_type: str
    mime_type: Optional[str] = None
    content_type: Optional[str] = None
    size_bytes: Optional[float] = None
    download_url: Optional[str] = None

class StructOptResponse(BaseModel):
    """结构优化响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class TaskStatusResponse(BaseModel):
    """任务状态查询响应模型 — v2-lite 增强版，保持向后兼容"""
    task_id: str
    user_id: str
    task_type: str
    status: TaskStatus
    progress: int = Field(..., description="进度百分比 (0-100)")

    # --- v2-lite 新增字段 ---
    analysis_status: Optional[AnalysisStatus] = Field(None, description="分析状态")
    result_summary: Optional[dict] = Field(None, description="结构化结果摘要")
    html_report_url: Optional[str] = Field(None, description="HTML分析报告URL")
    artifacts: Optional[List[ArtifactInfo]] = Field(None, description="产物列表")
    progress_message: Optional[str] = Field(None, description="进度消息")

    # --- 兼容旧字段 ---
    params: Optional[dict] = None
    result_path: Optional[str] = None
    external_job_id: Optional[str] = None
    process_id: Optional[int] = Field(None, description="VASP进程ID")
    error_message: Optional[str] = None
    result_data: Optional[dict] = Field(None, description="详细的计算结果数据")
    created_at: str
    updated_at: str


class WorkerRegisterRequest(BaseModel):
    worker_id: str
    queue_name: str = "default"


class WorkerRegisterResponse(BaseModel):
    worker_id: str
    queue_name: str
    poll_interval_seconds: int


class WorkerClaimRequest(BaseModel):
    worker_id: str
    queue_name: str = "default"


class WorkerClaimResponse(BaseModel):
    task_id: str
    status: TaskStatus
    worker_id: str
    lease_token: str
    lease_expires_at: datetime
    task_type: str
    queue_name: str
    params: Optional[dict] = None


class WorkerLeaseRequest(BaseModel):
    worker_id: str
    lease_token: str


class WorkerHeartbeatResponse(BaseModel):
    task_id: str
    status: TaskStatus
    lease_expires_at: datetime
    cancel_requested: bool = False


class WorkerTaskStatusResponse(BaseModel):
    task_id: str
    status: TaskStatus


class WorkerCompleteRequest(WorkerLeaseRequest):
    result_data: Optional[dict] = None
    artifact_manifest: List[dict] = Field(default_factory=list)


class WorkerFailRequest(WorkerLeaseRequest):
    error_message: str
    failure_code: Optional[str] = None
    artifact_manifest: List[dict] = Field(default_factory=list)

class StructOptResult(BaseModel):
    """结构优化结果模型"""
    optimized_structure: Optional[str] = Field(None, description="优化后的结构文件路径")
    energy: Optional[float] = Field(None, description="总能量 (eV)")
    final_forces: Optional[list] = Field(None, description="最终力矩阵")
    convergence: bool = Field(False, description="是否收敛")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    html_analysis_report: Optional[str] = Field(None, description="可视化分析报告HTML文件路径")

class SCFRequest(StrictRequestModel):
    """自洽场计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    optimized_task_id: Optional[str] = Field(None, description="已完成的结构优化任务ID")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    precision: str = Field("Accurate", description="计算精度 (Normal, High, Accurate)")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([bool(self.cif_url), bool(self.optimized_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 optimized_task_id 中的一个")

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

class DOSRequest(StrictRequestModel):
    """态密度计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    scf_task_id: Optional[str] = Field(None, description="已完成的自洽场计算任务ID")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    kpoint_multiplier: float = Field(2.0, description="K点倍增因子 (相对于优化计算)")
    precision: str = Field("Accurate", description="计算精度 (Normal, High, Accurate)")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([bool(self.cif_url), bool(self.scf_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 scf_task_id 中的一个")

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

class BandStructureRequest(StrictRequestModel):
    """能带结构计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")

    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    scf_task_id: Optional[str] = Field(None, description="已完成的自洽场计算任务ID")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # VASP计算参数
    kpoint_density: float = Field(30.0, description="K点密度参数 (用于SCF步骤)")
    line_density: int = Field(20, description="能带k路径上每段的k点数")
    precision: str = Field("Accurate", description="计算精度 (Normal, High, Accurate)")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([bool(self.cif_url), bool(self.scf_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 scf_task_id 中的一个")


class BandStructureResponse(BaseModel):
    """能带结构计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")


class BandStructureResult(BaseModel):
    """能带结构计算结果模型"""
    band_structure_path: Optional[str] = Field(None, description="能带结构数据文件路径")
    total_energy: Optional[float] = Field(None, description="总能量 (eV)")
    fermi_energy: Optional[float] = Field(None, description="费米能级 (eV)")
    band_gap: Optional[float] = Field(None, description="带隙 (eV)")
    is_direct: Optional[bool] = Field(None, description="是否为直接带隙")
    vbm: Optional[float] = Field(None, description="价带顶 (eV)")
    cbm: Optional[float] = Field(None, description="导带底 (eV)")
    convergence: bool = Field(False, description="是否收敛")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    kpath_info: Optional[dict] = Field(None, description="高对称k路径信息")


class MDRequest(StrictRequestModel):
    """分子动力学计算请求模型"""
    # 用户ID（必填）
    user_id: str = Field(..., description="用户ID")
    
    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    scf_task_id: Optional[str] = Field(None, description="已完成的自洽场计算任务ID")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # MD计算参数
    md_steps: int = Field(1000, description="MD步数")
    temperature: float = Field(300.0, description="目标温度 (K)")
    time_step: float = Field(1.0, description="时间步长 (fs)")
    ensemble: str = Field("NVT", description="系综类型 (NVT, NVE, NPT)")
    precision: str = Field("Normal", description="计算精度 (Normal, High, Accurate)")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        """验证输入参数"""
        input_count = sum([bool(self.cif_url), bool(self.scf_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 scf_task_id 中的一个")

class MDResponse(BaseModel):
    """分子动力学计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")

class MDResult(BaseModel):
    """分子动力学计算结果模型（单温度）"""
    md_structure: Optional[str] = Field(None, description="MD计算的初始结构文件路径")
    xdatcar_path: Optional[str] = Field(None, description="XDATCAR轨迹文件路径")
    oszicar_path: Optional[str] = Field(None, description="OSZICAR能量文件路径")
    final_energy: Optional[float] = Field(None, description="最终能量 (eV)")
    average_temperature: Optional[float] = Field(None, description="平均温度 (K)")
    total_md_steps: Optional[int] = Field(None, description="完成的MD步数")
    convergence: bool = Field(False, description="是否正常完成")
    computation_time: Optional[float] = Field(None, description="计算耗时 (秒)")
    trajectory_data: Optional[dict] = Field(None, description="轨迹统计数据")


class MDMultiAnalyzeRequest(BaseModel):
    """多任务MD聚合分析请求模型 — 接受多个已完成的MD task_id"""
    user_id: str = Field(..., description="用户ID")
    task_ids: List[str] = Field(..., description="多个已完成的MD任务ID列表", min_length=1)


# ====================================================================== #
#  独立分析请求/响应模型
# ====================================================================== #

class AnalyzeRequest(StrictRequestModel):
    """独立分析请求模型 — 支持 task_id 或 file_url 二选一"""
    user_id: Optional[str] = Field(None, description="用户ID")
    task_id: Optional[str] = Field(None, description="已有任务ID，直接分析其结果目录")
    file_url: Optional[HttpUrl] = Field(None, description="VASP输出文件压缩包URL（zip/tar.gz），包含OUTCAR等文件")

    def model_post_init(self, __context) -> None:
        input_count = sum([bool(self.task_id), bool(self.file_url)])
        if input_count != 1:
            raise ValueError("必须提供 task_id 或 file_url 中的一个")


class AnalyzeResponse(BaseModel):
    """独立分析响应模型"""
    success: bool = Field(..., description="分析是否成功")
    analysis_type: str = Field(..., description="分析类型")
    summary: Optional[dict] = Field(None, description="结构化分析摘要")
    html_report_url: Optional[str] = Field(None, description="HTML分析报告URL")
    error_message: Optional[str] = Field(None, description="错误信息")


class NEBRequest(StrictRequestModel):
    """NEB（过渡态）计算请求模型"""
    user_id: str = Field(..., description="用户ID")

    # 初始结构（二选一）
    initial_cif_url: Optional[HttpUrl] = Field(None, description="初始结构文件URL")
    initial_task_id: Optional[str] = Field(None, description="已完成的初始结构优化任务ID")

    # 终态结构（二选一）
    final_cif_url: Optional[HttpUrl] = Field(None, description="终态结构文件URL")
    final_task_id: Optional[str] = Field(None, description="已完成的终态结构优化任务ID")

    # NEB 参数
    n_images: int = Field(5, ge=2, le=20, description="中间图像数（不含端点）")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")
    kpoint_density: float = Field(30.0, description="K点密度参数")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        initial_count = sum([bool(self.initial_cif_url), bool(self.initial_task_id)])
        final_count = sum([bool(self.final_cif_url), bool(self.final_task_id)])
        if initial_count != 1:
            raise ValueError("必须提供 initial_cif_url 或 initial_task_id 中的一个")
        if final_count != 1:
            raise ValueError("必须提供 final_cif_url 或 final_task_id 中的一个")


class NEBResponse(BaseModel):
    """NEB 计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")


class PhononRequest(StrictRequestModel):
    """声子计算请求模型（DFPT/有限位移，Gamma 点）"""
    user_id: str = Field(..., description="用户ID")

    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    scf_task_id: Optional[str] = Field(None, description="已完成的SCF/优化任务ID（复用POSCAR+POTCAR）")
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # 声子参数
    kpoint_density: float = Field(30.0, description="K点密度参数")
    displacement: float = Field(0.015, description="有限位移步长（Å），影响精度")
    custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")

    _check_custom_incar = field_validator("custom_incar")(_validate_custom_incar)

    def model_post_init(self, __context) -> None:
        input_count = sum([bool(self.cif_url), bool(self.scf_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 scf_task_id 中的一个")


class PhononResponse(BaseModel):
    """声子计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")


class CustomCalcRequest(StrictRequestModel):
    """通用自定义计算请求模型 — 用户完全控制INCAR，覆盖 HSE06/DFPT/GW/Wannier90/SOC 等长尾需求
    注意：Bader 电荷分析和 ELF 已由 SCF 分析器原生支持，无需使用此接口。"""
    user_id: str = Field(..., description="用户ID")

    # 输入源（二选一）
    cif_url: Optional[HttpUrl] = Field(None, description="结构文件URL")
    from_task_id: Optional[str] = Field(
        None, description="复用已完成任务的结构（优先CONTCAR，其次POSCAR）"
    )
    queue_name: Optional[str] = Field(None, description="目标超算队列/共享存储标识")

    # INCAR — 由用户完全指定，不做黑名单过滤
    incar: Dict[str, Any] = Field(..., description="完整的INCAR参数字典，用户自行控制所有参数")

    # K点设置
    kpoint_density: float = Field(30.0, description="K点密度，用于 mesh 模式")
    kpoint_mode: str = Field(
        "mesh",
        description="K点模式: 'mesh'=自动Monkhorst-Pack网格, 'gamma'=仅Gamma点(1×1×1)",
    )

    def model_post_init(self, __context) -> None:
        input_count = sum([bool(self.cif_url), bool(self.from_task_id)])
        if input_count != 1:
            raise ValueError("必须提供 cif_url 或 from_task_id 中的一个")
        if not self.incar:
            raise ValueError("incar 参数不能为空，请提供至少一个INCAR参数")


class CustomCalcResponse(BaseModel):
    """通用自定义计算响应模型"""
    task_id: str = Field(..., description="任务ID")
    status: TaskStatus = Field(..., description="任务状态")
    message: str = Field(..., description="响应消息")


class AgentAnalyzeRequest(StrictRequestModel):
    """AI Agent 分析请求"""
    user_id: str = Field(..., description="用户ID")
    task_id: str = Field(..., description="已完成的任务ID")
    question: str = Field(..., description="分析问题，例如：'带隙是多少？绘制能量收敛曲线。'")
    model: str = Field("claude-haiku-4-5-20251001", description="内层分析 Agent 使用的模型")


class AgentAnalyzeResponse(BaseModel):
    """AI Agent 分析响应"""
    success: bool
    answer: str = Field("", description="Agent 的完整分析回答（含 JSON 摘要）")
    steps: int = Field(0, description="Agent 执行代码的次数")
    generated_plots: List[str] = Field(default_factory=list, description="生成图表的访问URL列表")
    error_message: Optional[str] = None
