from typing import Any, Dict, Optional

from fastmcp import FastMCP,Context

from .client import VaspAPIClient
from .config import mcp_config
from .models import StructOptInput, SCFInput, DOSInput, MDInput


mcp = FastMCP("VASP-MCP")
client = VaspAPIClient()


# 自定义INCAR参数使用说明
CUSTOM_INCAR_HELP = """
自定义INCAR参数功能使用说明:

custom_incar参数允许智能体直接指定VASP计算的INCAR参数，提供最大的灵活性。

常用INCAR参数示例:
- 电子结构参数:
  * EDIFF: 电子收敛精度 (如: 1e-6, 1e-7)
  * NELM: 最大电子步数 (如: 100, 200)
  * ALGO: 算法选择 (如: "Fast", "Normal", "Very_Fast")
  * ISMEAR: 展宽方法 (如: 0, -5)
  * SIGMA: 展宽参数 (如: 0.05, 0.1)

- DOS计算专用:
  * LORBIT: 轨道投影 (如: 11, 12)
  * NEDOS: 能量网格点数 (如: 2000, 3000)
  * EMIN/EMAX: 能量范围 (如: -20, 10)

- MD计算专用:
  * SMASS: 热浴质量 (如: 0, 1)
  * POTIM: 时间步长 (如: 0.5, 1.0)
  * ISYM: 对称性 (如: 0关闭, 1开启)
  * MDALGO: MD算法 (如: 1,2,3)

- 优化计算:
  * EDIFFG: 力收敛精度 (如: -0.01, -0.02)
  * IBRION: 离子运动方法 (如: 1,2,3)
  * ISIF: 应力张量 (如: 2,3,7)

使用方式:
custom_incar参数接受字典格式，键为INCAR参数名，值为参数值。
参数名不区分大小写，会自动转换为大写。
自定义参数会覆盖默认生成的同名参数。

示例:
{"custom_incar": {"EDIFF": 1e-7, "NELM": 100, "ALGO": "Fast"}}
{"custom_incar": {"LORBIT": 11, "NEDOS": 3000, "ISMEAR": -5}}
{"custom_incar": {"SMASS": 1, "POTIM": 0.5, "ISYM": 0}}
"""
def get_user_id(ctx: Context) -> str:
    if ctx is not None:
        user_id = ctx.request_context.request.headers.get("user_id", None)  # type: ignore
    else:
        user_id = "123"
    return user_id #type: ignore

@mcp.tool()
async def submit_structure_optimization(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    提交结构优化任务，返回任务信息。
    
    参数说明:
    - calc_type (必填): 计算类型，可选值:
      * "OXC" - 氧化物/硫化物固体电解质
      * "SSE" - 固体电解质(等同于OXC)
      * "ORC" - 氧化物还原催化剂
      * "ECAT_OER" - 氧析出反应催化剂
      * "ECAT_HER" - 氢析出反应催化剂
    - formula (可选): 化学式，如'Li2O', 'LiFePO4'，与cif_url二选一
    - cif_url (可选): CIF文件的URL地址，与formula二选一
    - spacegroup (可选): 空间群符号，如"P1", "Fm-3m"，仅当使用formula时有效，用于筛选结构
    - max_energy_above_hull (可选): 最大能量上凸包距离(eV/atom)，默认0.1，仅当使用formula时有效，用于筛选结构
    - min_band_gap (可选): 最小带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_band_gap (可选): 最大带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_nsites (可选): 最大原子数，仅当使用formula时有效，用于筛选结构
    - min_nsites (可选): 最小原子数，仅当使用formula时有效，用于筛选结构
    - stable_only (可选): 只选择稳定材料，默认true，仅当使用formula时有效，用于筛选结构
    - selection_mode (可选): 选择模式，可选"auto"/"stable"/"first"，仅当使用formula时有效，用于筛选结构
    - kpoint_density (可选): K点密度参数，默认30.0
    - custom_incar (可选): 自定义INCAR参数字典，用于覆盖或添加特定的VASP计算参数
    
    示例:
    submit_structure_optimization(calc_type="OXC", formula="Li2O", kpoint_density=25.0)
    submit_structure_optimization(calc_type="OXC", formula="Li2O", custom_incar={"EDIFF": 1e-7, "NELM": 100, "ALGO": "Fast"})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    print("submit_structure_optimization")
    print(params)
    payload = StructOptInput(**params).model_dump(mode="json", by_alias=True)
    if payload.get("user_id") is None:
        payload["user_id"] = "123"
    return await client.submit_structure_optimization(payload)


@mcp.tool()
async def submit_scf_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    optimized_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore

) -> Dict[str, Any]:
    """
    提交自洽场计算任务，返回任务信息。
    
    参数说明:
    - calc_type (必填): 计算类型，同结构优化
    - formula (可选): 化学式，如'Li2O', 'LiFePO4'，与cif_url和optimized_task_id三选一
    - cif_url (可选): CIF文件的URL地址，与formula和optimized_task_id三选一
    - optimized_task_id (可选): 已完成的结构优化任务ID，与formula和cif_url三选一
    - spacegroup (可选): 空间群符号，仅当使用formula时有效
    - max_energy_above_hull (可选): 最大能量上凸包距离(eV/atom)，默认0.1，仅当使用formula时有效，用于筛选结构
    - min_band_gap (可选): 最小带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_band_gap (可选): 最大带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_nsites (可选): 最大原子数，仅当使用formula时有效，用于筛选结构
    - min_nsites (可选): 最小原子数，仅当使用formula时有效，用于筛选结构
    - stable_only (可选): 只选择稳定材料，默认true，仅当使用formula时有效，用于筛选结构
    - selection_mode (可选): 选择模式，可选"auto"/"stable"/"first"，仅当使用formula时有效，用于筛选结构
    - kpoint_density (可选): K点密度参数，默认30.0
    - precision (可选): 计算精度，可选"Normal"/"High"/"Accurate"，默认"Accurate"
    - custom_incar (可选): 自定义INCAR参数字典，用于覆盖或添加特定的VASP计算参数
    
    示例:
    submit_scf_calculation(user_id="user123", calc_type="OXC", optimized_task_id="task_001", precision="High")
    submit_scf_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"ISMEAR": 0, "SIGMA": 0.05})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = SCFInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_scf(payload)


@mcp.tool()
async def submit_dos_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    kpoint_multiplier: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    提交态密度计算任务，返回任务信息。
    
    参数说明:
    - calc_type (必填): 计算类型，同结构优化
    - formula (可选): 化学式，如'Li2O', 'LiFePO4'，与cif_url和scf_task_id三选一
    - cif_url (可选): CIF文件的URL地址，与formula和scf_task_id三选一  
    - scf_task_id (可选): 已完成的自洽场计算任务ID，与formula和cif_url三选一
    - spacegroup (可选): 空间群符号，仅当使用formula时有效
    - max_energy_above_hull (可选): 最大能量上凸包距离(eV/atom)，默认0.1
    - min_band_gap (可选): 最小带隙(eV)
    - max_band_gap (可选): 最大带隙(eV)
    - max_nsites (可选): 最大原子数
    - min_nsites (可选): 最小原子数
    - stable_only (可选): 只选择稳定材料，默认true
    - selection_mode (可选): 选择模式，可选"auto"/"stable"/"first"
    - kpoint_density (可选): K点密度参数，默认30.0
    - kpoint_multiplier (可选): K点倍增因子，相对于优化计算，默认2.0
    - precision (可选): 计算精度，可选"Normal"/"High"/"Accurate"，默认"Accurate"
    - custom_incar (可选): 自定义INCAR参数字典，用于覆盖或添加特定的VASP计算参数
    
    示例:
    submit_dos_calculation(user_id="user123", calc_type="OXC", scf_task_id="scf_001", kpoint_multiplier=3.0)
    submit_dos_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"LORBIT": 11, "NEDOS": 3000})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = DOSInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_dos(payload)


@mcp.tool()
async def submit_md_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    md_steps: Optional[int] = None,
    temperature: Optional[Any] = None,  # Can be float or List[float]
    time_step: Optional[float] = None,
    ensemble: Optional[str] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    提交分子动力学计算任务，返回任务信息。
    
    参数说明:
    - calc_type (必填): 计算类型，同结构优化
    - formula (可选): 化学式，如'Li2O', 'LiFePO4'，与cif_url和scf_task_id三选一
    - cif_url (可选): CIF文件的URL地址，与formula和scf_task_id三选一
    - scf_task_id (可选): 已完成的自洽场计算任务ID，与formula和cif_url三选一
    - spacegroup (可选): 空间群符号，仅当使用formula时有效
    - max_energy_above_hull (可选): 最大能量上凸包距离(eV/atom)，默认0.1，仅当使用formula时有效，用于筛选结构
    - min_band_gap (可选): 最小带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_band_gap (可选): 最大带隙(eV)，仅当使用formula时有效，用于筛选结构
    - max_nsites (可选): 最大原子数，仅当使用formula时有效，用于筛选结构
    - min_nsites (可选): 最小原子数，仅当使用formula时有效，用于筛选结构
    - stable_only (可选): 只选择稳定材料，默认true，仅当使用formula时有效，用于筛选结构
    - selection_mode (可选): 选择模式，可选"auto"/"stable"/"first"，仅当使用formula时有效，用于筛选结构
    - md_steps (可选): MD步数，默认1000
    - temperature (可选): 目标温度(K)，支持单个温度或温度列表进行多温度扫描
      * 单温度: 300.0
      * 多温度: [200.0, 300.0, 400.0, 500.0] - 将为每个温度创建子任务
    - time_step (可选): 时间步长(fs)，默认1.0
    - ensemble (可选): 系综类型，可选"NVT"/"NVE"/"NPT"，默认"NVT"
    - precision (可选): 计算精度，可选"Normal"/"High"/"Accurate"，默认"Normal"
    - custom_incar (可选): 自定义INCAR参数字典，用于覆盖或添加特定的VASP计算参数
    
    🌡️ 多温度MD计算特性:
    - 每个温度点创建独立的子目录 (如: T_300K/, T_400K/)
    - 所有子任务共享相同的任务ID，便于管理
    - 支持子任务独立成功/失败状态
    - 自动生成多温度汇总分析报告
    
    示例:
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", md_steps=2000, temperature=350.0, ensemble="NPT")
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", temperature=[300.0, 400.0, 500.0], md_steps=1500)
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"SMASS": 1, "POTIM": 0.5, "ISYM": 0})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = MDInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_md(payload)


@mcp.tool()
async def get_task_status(task_id: str, 
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    查询任务状态。
    
    参数说明:
    - task_id (必填): 任务ID，由提交任务时返回
    
    返回信息包含:
    - status: 任务状态("queued"/"running"/"completed"/"failed"/"canceled")
    - progress: 进度百分比(0-100)
    - result_path: 结果文件路径(完成时)
    - error_message: 错误信息(失败时)
    - result_data: 详细结果数据
    
    🌡️ 对于多温度MD任务，result_data还包含:
    - multi_temperature_info: 多温度子任务状态汇总
      * is_multi_temperature: 是否为多温度计算
      * total_subtasks: 子任务总数
      * completed_subtasks: 完成的子任务数
      * failed_subtasks: 失败的子任务数
      * subtask_status: 各温度点详细状态列表
    
    示例:
    get_task_status("task_001")
    """
    return await client.get_task_status(task_id, get_user_id(ctx))


@mcp.tool()
async def cancel_task(task_id: str, 
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    取消正在运行或排队的任务。
    
    参数说明:
    - task_id (必填): 要取消的任务ID
    
    示例:
    cancel_task("task_001", "user123")
    """
    return await client.cancel_task(task_id, get_user_id(ctx))


@mcp.tool()
async def list_user_tasks(
ctx: Context = None #type: ignore
) -> Any:
    """
    列出用户的所有任务。
    
    参数说明:
    
    返回任务列表，每个任务包含:
    - task_id: 任务ID
    - task_type: 任务类型("structure_optimization"/"scf_calculation"/"dos_calculation"/"md_calculation")
    - status: 任务状态
    - created_at: 创建时间
    - updated_at: 更新时间
    
    示例:
    list_user_tasks()
    """
    return await client.list_tasks(get_user_id(ctx))


@mcp.tool()
async def get_task_result(task_id: str, 
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    获取已完成任务的结果路径信息。
    
    参数说明:
    - task_id (必填): 已完成的任务ID
    
    返回结果路径信息，包含:
    - result_path: 结果文件所在目录路径
    - message: 结果描述信息
    
    示例:
    get_task_result("task_001", "user123")
    """
    return await client.get_task_result(task_id, get_user_id(ctx))


@mcp.tool()
async def get_md_result(task_id: str, 
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    获取分子动力学任务的详细计算结果。
    
    参数说明:
    - task_id (必填): 已完成的MD任务ID
    
    返回详细MD结果，包含:
    
    📊 基本信息:
    - is_multi_temperature: 是否为多温度计算
    - total_subtasks: 子任务总数
    - completed_subtasks: 完成的子任务数
    - failed_subtasks: 失败的子任务数
    - convergence: 整体是否成功完成
    - computation_time: 总计算耗时(秒)
    
    🌡️ 多温度计算专用:
    - subtask_results: 各温度点详细结果列表，每项包含:
      * temperature: 温度(K)
      * subtask_dir: 子任务目录
      * md_structure: 初始结构文件路径
      * xdatcar_path: 轨迹文件路径(XDATCAR)
      * oszicar_path: 能量文件路径(OSZICAR)
      * final_energy: 最终能量(eV)
      * average_temperature: 平均温度(K)
      * total_md_steps: 完成的MD步数
      * convergence: 该温度点是否正常完成
      * status: 子任务状态
      * error_message: 错误信息(如有)
    
    📈 分析报告:
    - md_html_analysis_report: 多温度分析报告HTML文件路径
    - md_output_dir: 分析输出目录
    
    示例:
    get_md_result("md_task_001", "user123")
    """
    return await client.get_md_result(task_id, get_user_id(ctx))


@mcp.tool()
async def get_custom_incar_help() -> str:
    """
    获取自定义INCAR参数的详细使用说明和常用参数示例。
    
    这个工具提供了:
    - custom_incar参数的使用方法
    - 常用INCAR参数的详细说明
    - 不同计算类型的推荐参数
    - 实际使用示例
    
    调用此工具可以帮助智能体了解如何正确使用custom_incar参数来自定义VASP计算。
    """
    return CUSTOM_INCAR_HELP


# def run():
#     # 运行 FastMCP 服务器
#     mcp.run(host=mcp_config.host, port=mcp_config.port)


# if __name__ == "__main__":
#     run()


