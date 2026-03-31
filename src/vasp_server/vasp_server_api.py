from fastapi import FastAPI, HTTPException, Query
from fastapi.middleware.cors import CORSMiddleware
from typing import List
import os
import uvicorn
import logging
import sys
from pathlib import Path
# Config import moved to __main__ block
from .internal_worker_api import build_internal_worker_router
from .schemas import (
    StructOptRequest, StructOptResponse, TaskStatusResponse,
    TaskStatus, AnalysisStatus, ArtifactInfo,
    SCFRequest, SCFResponse,
    DOSRequest, DOSResponse,
    BandStructureRequest, BandStructureResponse,
    MDRequest, MDResponse, MDMultiAnalyzeRequest,
    AnalyzeRequest, AnalyzeResponse,
    NEBRequest, NEBResponse,
    PhononRequest, PhononResponse,
    CustomCalcRequest, CustomCalcResponse,
    AgentAnalyzeRequest, AgentAnalyzeResponse,
)
from .task_manager.manager import TaskManager
from .task_manager.database import check_and_init_db, engine

# 配置日志 - 包含线程信息，确保子线程日志也能正确输出
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - [%(threadName)s] - %(levelname)s - %(message)s',
    handlers=[
        logging.StreamHandler(sys.stdout),  # 输出到标准输出
    ],
    force=True  # 强制重新配置，覆盖已有配置
)
logger = logging.getLogger(__name__)

app = FastAPI(title="VASP计算服务API", version="1.0.0")

# --- Prometheus Metrics ---
try:
    from prometheus_fastapi_instrumentator import Instrumentator
    from prometheus_client import Gauge, Counter

    vasp_tasks_queued = Gauge(
        "vasp_tasks_queued", "Number of queued VASP tasks"
    )
    vasp_tasks_running = Gauge(
        "vasp_tasks_running", "Number of running VASP tasks"
    )
    vasp_tasks_total = Counter(
        "vasp_tasks_submitted_total",
        "Total VASP tasks submitted",
        ["task_type"],
    )
    vasp_tasks_completed = Counter(
        "vasp_tasks_completed_total",
        "Total VASP tasks completed",
        ["task_type", "status"],  # status: completed / failed
    )

    Instrumentator().instrument(app).expose(app, endpoint="/metrics")
    _METRICS_ENABLED = True
    logger.info("Prometheus metrics enabled at /metrics")
except ImportError:
    _METRICS_ENABLED = False
    logger.warning("prometheus-fastapi-instrumentator not installed, /metrics disabled")

# 配置CORS
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],  # 允许所有来源，生产环境中应该设置具体的域名
    allow_credentials=True,
    allow_methods=["*"],  # 允许所有方法
    allow_headers=["*"],  # 允许所有头部
)

# 挂载静态文件服务（在 startup 时根据配置决定）


@app.on_event("startup")
async def _startup_init():
    logger.info("🚀 启动VASP计算服务...")
    logger.info("🔧 检查并初始化数据库...")
    check_and_init_db()

# 创建全局任务管理器实例
task_manager = TaskManager()
app.include_router(build_internal_worker_router(task_manager))

@app.get("/")
async def root():
    return {"message": "VASP计算服务API", "version": "1.0.0"}


@app.get("/health")
async def health_check():
    """Health check endpoint: verifies database connectivity and reports task counts."""
    from sqlalchemy import text
    from fastapi.responses import JSONResponse

    version = "1.0.0"

    # 1. Check database connectivity
    try:
        with engine.connect() as conn:
            conn.execute(text("SELECT 1"))
        db_status = "connected"
    except Exception:
        return JSONResponse(
            status_code=503,
            content={
                "status": "unhealthy",
                "database": "disconnected",
                "tasks": None,
                "version": version,
            },
        )

    # 2. Count tasks by status
    try:
        with engine.connect() as conn:
            rows = conn.execute(
                text("SELECT status, COUNT(*) FROM tasks GROUP BY status")
            ).fetchall()
        counts = {row[0]: row[1] for row in rows}
        queued = counts.get("queued", 0)
        running = counts.get("running", 0)
        total = sum(counts.values())
    except Exception:
        queued = running = total = 0

    return {
        "status": "healthy",
        "database": db_status,
        "tasks": {"queued": queued, "running": running, "total": total},
        "version": version,
    }


def _record_task_submitted(task_type: str):
    """Record a task submission in Prometheus metrics."""
    if _METRICS_ENABLED:
        vasp_tasks_total.labels(task_type=task_type).inc()
        vasp_tasks_queued.inc()


@app.post("/vasp/structure-optimization", response_model=StructOptResponse)
async def submit_structure_optimization(request: StructOptRequest):
    """
    提交结构优化任务
    
    支持两种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件
    2. CIF URL：直接从指定URL下载CIF文件
    
    Returns:
        StructOptResponse: 包含任务ID和状态的响应
    """
    try:
        # 准备任务参数 (input validation is handled by StructOptRequest's model_post_init)
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,

            "kpoint_density": request.kpoint_density,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        if request.user_id is None:
            request.user_id = os.getenv("DEFAULT_USER_ID", "123")
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="structure_optimization",
            params=task_params
        )
        _record_task_submitted("structure_optimization")

        return StructOptResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"结构优化任务已提交，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交任务失败: {str(e)}")


@app.post("/vasp/scf-calculation", response_model=SCFResponse)
async def submit_scf_calculation(request: SCFRequest):
    """
    提交自洽场计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件
    2. CIF URL：直接从指定URL下载CIF文件
    3. 优化任务ID：基于已完成的结构优化任务的CONTCAR文件
    
    Returns:
        SCFResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.optimized_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 optimized_task_id 中的一个"
            )
        
        # 如果基于优化任务，验证任务存在性
        if request.optimized_task_id:
            opt_task = task_manager.get_task(request.optimized_task_id, request.user_id)
            if not opt_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"结构优化任务 {request.optimized_task_id} 未找到或无权限访问"
                )
            if str(opt_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"结构优化任务 {request.optimized_task_id} 尚未完成"
                )
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "optimized_task_id": request.optimized_task_id,

            "kpoint_density": request.kpoint_density,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="scf_calculation",
            params=task_params
        )
        _record_task_submitted("scf_calculation")

        input_source = "化学式" if request.formula else "CIF URL" if request.cif_url else "结构优化任务"
        
        return SCFResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"自洽场计算任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交自洽场计算任务失败: {str(e)}")


@app.post("/vasp/dos-calculation", response_model=DOSResponse)
async def submit_dos_calculation(request: DOSRequest):
    """
    提交态密度计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件（需要先完成自洽场计算）
    2. CIF URL：直接从指定URL下载CIF文件（需要先完成自洽场计算）
    3. 自洽场任务ID：基于已完成的自洽场计算任务结果
    
    Returns:
        DOSResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.scf_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 scf_task_id 中的一个"
            )
        
        # 如果基于自洽场任务，验证任务存在性
        if request.scf_task_id:
            scf_task = task_manager.get_task(request.scf_task_id, request.user_id)
            if not scf_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 未找到或无权限访问"
                )
            if str(scf_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 尚未完成"
                )
        
        # 准备任务参数
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,

            "kpoint_density": request.kpoint_density,
            "kpoint_multiplier": request.kpoint_multiplier,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="dos_calculation",
            params=task_params
        )
        _record_task_submitted("dos_calculation")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "态密度计算"
        elif request.formula:
            input_source = "化学式"
            calc_mode = "单点自洽+DOS计算"
        else:
            input_source = "CIF URL"
            calc_mode = "单点自洽+DOS计算"
        
        return DOSResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交态密度计算任务失败: {str(e)}")


@app.post("/vasp/band-structure", response_model=BandStructureResponse)
async def submit_band_structure_calculation(request: BandStructureRequest):
    """
    提交能带结构计算任务

    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件
    2. CIF URL：直接从指定URL下载CIF文件
    3. 自洽场任务ID：基于已完成的自洽场计算任务结果（推荐）

    Returns:
        BandStructureResponse: 包含任务ID和状态的响应
    """
    try:
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url),
            bool(request.scf_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400,
                detail="必须提供 formula、cif_url 或 scf_task_id 中的一个"
            )

        if request.scf_task_id:
            scf_task = task_manager.get_task(request.scf_task_id, request.user_id)
            if not scf_task:
                raise HTTPException(
                    status_code=404,
                    detail=f"自洽场计算任务 {request.scf_task_id} 未找到或无权限访问"
                )
            if str(scf_task.status) != "completed":
                raise HTTPException(
                    status_code=400,
                    detail=f"自洽场计算任务 {request.scf_task_id} 尚未完成"
                )

        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,
            "kpoint_density": request.kpoint_density,
            "line_density": request.line_density,
            "precision": request.precision,
        }

        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value

        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="band_structure",
            params=task_params
        )
        _record_task_submitted("band_structure")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
        elif request.formula:
            input_source = "化学式"
        else:
            input_source = "CIF URL"

        return BandStructureResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"能带结构计算任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )

    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交能带结构计算任务失败: {str(e)}")


@app.post("/vasp/md-calculation", response_model=MDResponse)
async def submit_md_calculation(request: MDRequest):
    """
    提交分子动力学计算任务
    
    支持三种输入方式：
    1. 化学式：从Materials Project数据库搜索和下载CIF文件（纯MD计算）
    2. CIF URL：直接从指定URL下载CIF文件（纯MD计算）
    3. 自洽场任务ID：基于已完成的自洽场计算任务结果
    
    Returns:
        MDResponse: 包含任务ID和状态的响应
    """
    try:
        # 验证输入参数
        input_count = sum([
            bool(request.formula),
            bool(request.cif_url), 
            bool(request.scf_task_id)
        ])
        if input_count != 1:
            raise HTTPException(
                status_code=400, 
                detail="必须提供 formula、cif_url 或 scf_task_id 中的一个"
            )
        
        # 如果基于自洽场任务，验证任务存在性
        if request.scf_task_id:
            scf_task = task_manager.get_task(request.scf_task_id, request.user_id)
            if not scf_task:
                raise HTTPException(
                    status_code=404, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 未找到或无权限访问"
                )
            if str(scf_task.status) != "completed":  # type: ignore
                raise HTTPException(
                    status_code=400, 
                    detail=f"自洽场计算任务 {request.scf_task_id} 尚未完成"
                )
        
        # 准备任务参数（temperature 始终为单个 float）
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,

            "md_steps": request.md_steps,
            "temperature": float(request.temperature),
            "time_step": request.time_step,
            "ensemble": request.ensemble,
            "precision": request.precision,
        }
        
        # 添加材料搜索参数（仅当使用formula时）
        if request.formula:
            search_params = {
                "spacegroup": request.spacegroup,
                "max_energy_above_hull": request.max_energy_above_hull,
                "min_band_gap": request.min_band_gap,
                "max_band_gap": request.max_band_gap,
                "max_nsites": request.max_nsites,
                "min_nsites": request.min_nsites,
                "stable_only": request.stable_only,
                "selection_mode": request.selection_mode.value,
            }
            # 只添加非None的参数
            for key, value in search_params.items():
                if value is not None:
                    task_params[key] = value
        
        # 提交任务
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="md_calculation",
            params=task_params
        )
        _record_task_submitted("md_calculation")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "分子动力学计算"
        elif request.formula:
            input_source = "化学式"
            calc_mode = "纯MD计算"
        else:
            input_source = "CIF URL"
            calc_mode = "纯MD计算"
        
        return MDResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交分子动力学计算任务失败: {str(e)}")

@app.post("/vasp/neb-calculation", response_model=NEBResponse)
async def submit_neb_calculation(request: NEBRequest):
    """提交 NEB（过渡态）计算任务"""
    try:
        task_params = {
            "initial_formula": request.initial_formula,
            "initial_cif_url": str(request.initial_cif_url) if request.initial_cif_url else None,
            "initial_task_id": request.initial_task_id,
            "final_formula": request.final_formula,
            "final_cif_url": str(request.final_cif_url) if request.final_cif_url else None,
            "final_task_id": request.final_task_id,
            "n_images": request.n_images,
            "kpoint_density": request.kpoint_density,
            "spacegroup": request.spacegroup,
            "max_energy_above_hull": request.max_energy_above_hull,
            "stable_only": request.stable_only,
            "selection_mode": request.selection_mode.value,
        }
        if request.custom_incar:
            task_params["custom_incar"] = request.custom_incar

        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="neb_calculation",
            params={k: v for k, v in task_params.items() if v is not None},
        )
        _record_task_submitted("neb_calculation")
        return NEBResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"NEB 计算任务已提交，任务ID: {task_id}",
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交 NEB 计算任务失败: {str(e)}")


@app.post("/vasp/phonon-calculation", response_model=PhononResponse)
async def submit_phonon_calculation(request: PhononRequest):
    """提交声子计算任务（IBRION=6，Gamma 点有限位移法）"""
    try:
        task_params = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,
            "kpoint_density": request.kpoint_density,
            "displacement": request.displacement,
        }
        if request.formula:
            for k in ("spacegroup", "max_energy_above_hull", "min_band_gap", "max_band_gap",
                      "max_nsites", "min_nsites", "stable_only", "selection_mode"):
                v = getattr(request, k)
                if v is not None:
                    task_params[k] = v.value if hasattr(v, 'value') else v
        if request.custom_incar:
            task_params["custom_incar"] = request.custom_incar

        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="phonon_calculation",
            params={k: v for k, v in task_params.items() if v is not None},
        )
        _record_task_submitted("phonon_calculation")
        return PhononResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"声子计算任务已提交，任务ID: {task_id}",
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交声子计算任务失败: {str(e)}")


@app.post("/vasp/custom-calculation", response_model=CustomCalcResponse)
async def submit_custom_calculation(request: CustomCalcRequest):
    """提交通用自定义计算任务 — 用户完全控制INCAR，适合HSE06、DFPT介电、ELF、Wannier等长尾需求"""
    try:
        task_params: Dict[str, Any] = {
            "formula": request.formula,
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "from_task_id": request.from_task_id,
            "incar": request.incar,
            "kpoint_density": request.kpoint_density,
            "kpoint_mode": request.kpoint_mode,
        }
        if request.formula:
            for k in ("spacegroup", "max_energy_above_hull", "min_band_gap", "max_band_gap",
                      "max_nsites", "min_nsites", "stable_only", "selection_mode"):
                v = getattr(request, k)
                if v is not None:
                    task_params[k] = v.value if hasattr(v, "value") else v

        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="custom_calculation",
            params={k: v for k, v in task_params.items() if v is not None},
        )
        _record_task_submitted("custom_calculation")
        return CustomCalcResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"自定义计算任务已提交，任务ID: {task_id}",
        )
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"提交自定义计算任务失败: {str(e)}")


@app.get("/vasp/task/{task_id}", response_model=TaskStatusResponse)
async def get_task_status(task_id: str, user_id: str = Query(..., description="用户ID")):
    """查询任务状态与任务结果"""
    try:
        task = task_manager.get_task(task_id, user_id)
        if not task:
            raise HTTPException(status_code=404, detail="任务未找到或无权限访问")

        return _build_task_response(task, task_id)

    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"查询任务状态失败: {e}", exc_info=True)
        raise HTTPException(status_code=500, detail=f"查询任务状态失败: {str(e)}")


def _safe_task_status(status_val) -> TaskStatus:
    """安全转换任务状态"""
    try:
        return TaskStatus(status_val) if status_val else TaskStatus.queued
    except (ValueError, AttributeError):
        return TaskStatus.queued


def _safe_analysis_status(status_val) -> AnalysisStatus | None:
    """安全转换分析状态"""
    try:
        return AnalysisStatus(status_val) if status_val else None
    except (ValueError, AttributeError):
        return None


def _is_public_url(value) -> bool:
    return isinstance(value, str) and value.startswith(("https://", "http://"))


def _sanitize_public_value(value):
    if isinstance(value, dict):
        return {key: _sanitize_public_value(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_sanitize_public_value(item) for item in value]
    if isinstance(value, str) and value.startswith("/"):
        return None
    return value


def _sanitize_report_urls(result_summary, result_data, html_report_url):
    sanitized_summary = _sanitize_public_value(result_summary)
    sanitized_data = _sanitize_public_value(result_data)

    if isinstance(sanitized_summary, dict) and "html_report_url" in sanitized_summary:
        sanitized_summary["html_report_url"] = html_report_url if _is_public_url(html_report_url) else None

    report_keys = [
        "analysis_report_html_path",
        "html_analysis_report",
        "md_analysis_report_html_path",
        "dos_analysis_report_html_path",
        "scf_analysis_report_html_path",
        "band_structure_report_html_path",
        "neb_report_html_path",
        "phonon_report_html_path",
    ]
    if isinstance(sanitized_data, dict):
        for key in report_keys:
            if key in sanitized_data:
                sanitized_data[key] = html_report_url if _is_public_url(html_report_url) else None

    return sanitized_summary, sanitized_data


def _build_task_response(task, task_id: str) -> TaskStatusResponse:
    """从 Task ORM 构建统一的 TaskStatusResponse"""
    # 从 result_summary 提取 html_report_url
    result_summary = getattr(task, 'result_summary', None)
    html_report_url = None
    if isinstance(result_summary, dict):
        html_report_url = result_summary.get('html_report_url')

    # 如果 result_summary 没有 html_report_url，从 result_data 兜底
    if not html_report_url:
        result_data = getattr(task, 'result_data', None)
        if isinstance(result_data, dict):
            for key in ['analysis_report_html_path', 'html_analysis_report',
                        'md_analysis_report_html_path', 'dos_analysis_report_html_path',
                        'scf_analysis_report_html_path', 'band_structure_report_html_path']:
                val = result_data.get(key)
                if val:
                    html_report_url = str(val)
                    break

    # 查询 artifacts
    artifacts_list = None
    try:
        raw_artifacts = task_manager.get_task_artifacts(task_id)
        if raw_artifacts:
            artifacts_list = []
            for a in raw_artifacts:
                download_url = task_manager.get_artifact_download_url(a)
                content_type = getattr(a, 'content_type', None) or getattr(a, 'mime_type', None)
                artifacts_list.append(
                    ArtifactInfo(
                        id=str(a.id),
                        artifact_type=str(a.artifact_type),
                        mime_type=getattr(a, 'mime_type', None),
                        content_type=content_type,
                        size_bytes=getattr(a, 'size_bytes', None),
                        download_url=download_url,
                    )
                )
                if str(getattr(a, 'artifact_type', '')) == "html_report" and download_url:
                    html_report_url = download_url
    except Exception:
        pass  # artifact 查询失败不影响主响应

    if not _is_public_url(html_report_url):
        html_report_url = None

    result_data = getattr(task, 'result_data', None)
    result_summary, result_data = _sanitize_report_urls(result_summary, result_data, html_report_url)
    result_path = getattr(task, 'result_path', None)
    if not _is_public_url(result_path):
        result_path = None

    response_data = {
        "task_id": getattr(task, 'id', None),
        "user_id": getattr(task, 'user_id', None),
        "task_type": getattr(task, 'task_type', None),
        "status": _safe_task_status(getattr(task, 'status', None)),
        "progress": getattr(task, 'progress', 0) or 0,
        "analysis_status": _safe_analysis_status(getattr(task, 'analysis_status', None)),
        "result_summary": result_summary,
        "html_report_url": html_report_url,
        "artifacts": artifacts_list,
        "progress_message": getattr(task, 'progress_message', None),
        "params": getattr(task, 'params', None),
        "result_path": result_path,
        "external_job_id": getattr(task, 'external_job_id', None),
        "process_id": getattr(task, 'process_id', None),
        "error_message": getattr(task, 'error_message', None),
        "result_data": result_data,
        "created_at": getattr(task, 'created_at', None),
        "updated_at": getattr(task, 'updated_at', None),
    }

    # 时间字段转字符串
    for tf in ("created_at", "updated_at"):
        val = response_data[tf]
        response_data[tf] = val.isoformat() if val else ""

    return TaskStatusResponse(**response_data)


@app.post("/vasp/task/{task_id}/cancel")
async def cancel_task(task_id: str, user_id: str = Query(..., description="用户ID")):
    """
    取消任务
    
    Args:
        task_id: 任务ID
        user_id: 用户ID
        
    Returns:
        dict: 取消结果
    """
    try:
        success = task_manager.cancel_task(task_id, user_id)
        
        if not success:
            raise HTTPException(status_code=404, detail="任务未找到或无法取消")
        
        return {"message": f"任务 {task_id} 已请求取消"}
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"取消任务失败: {str(e)}")


@app.get("/vasp/tasks", response_model=List[TaskStatusResponse])
async def list_user_tasks(user_id: str = Query(..., description="用户ID")):
    """列出用户的所有任务"""
    try:
        tasks = task_manager.list_tasks(user_id)
        return [_build_task_response(task, str(task.id)) for task in tasks]
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"获取任务列表失败: {str(e)}")


@app.get("/download/file/{file_path:path}")
async def download_file(file_path: str):
    raise HTTPException(
        status_code=410,
        detail="本地文件下载已下线，请使用任务状态中的对象存储签名 URL",
    )


# ====================================================================== #
#  独立分析端点 — 支持 task_id 或 file_url
# ====================================================================== #

import tempfile, zipfile, tarfile, shutil, requests as _requests


def _resolve_analysis_input_dir(req: AnalyzeRequest) -> tuple[str, str | None]:
    """解析 VASP 输出目录。返回 (work_dir, tmp_root)，tmp_root 非 None 时调用方需清理。"""
    if req.task_id:
        user_id = req.user_id or os.getenv("DEFAULT_USER_ID", "123")
        task = task_manager.get_task(req.task_id, user_id)
        if task is None:
            raise HTTPException(status_code=404, detail=f"任务 {req.task_id} 未找到或无权限访问")
        work_dir = getattr(task, "result_path", None)
        if not work_dir:
            result_data = getattr(task, "result_data", None) or {}
            work_dir = result_data.get("work_directory")
        if not work_dir or not Path(work_dir).exists():
            raise HTTPException(status_code=400, detail=f"任务 {req.task_id} 的工作目录不存在或任务尚未完成")
        return str(work_dir), None

    # --- file_url: 下载压缩包并解压 ---
    url = str(req.file_url)
    tmp_dir = tempfile.mkdtemp(prefix="vasp_analyze_")
    archive_path = Path(tmp_dir) / "archive"
    try:
        resp = _requests.get(url, timeout=120, stream=True)
        resp.raise_for_status()
        with open(archive_path, "wb") as f:
            for chunk in resp.iter_content(chunk_size=8192):
                f.write(chunk)
    except Exception as e:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise HTTPException(status_code=400, detail=f"下载文件失败: {e}")

    extract_dir = Path(tmp_dir) / "extracted"
    extract_dir.mkdir()
    try:
        if zipfile.is_zipfile(archive_path):
            with zipfile.ZipFile(archive_path, "r") as zf:
                for member in zf.namelist():
                    member_path = (extract_dir / member).resolve()
                    if not str(member_path).startswith(str(extract_dir.resolve())):
                        raise HTTPException(status_code=400, detail=f"压缩包含非法路径: {member}")
                zf.extractall(extract_dir)
        elif tarfile.is_tarfile(str(archive_path)):
            with tarfile.open(str(archive_path), "r:*") as tf:
                for member in tf.getmembers():
                    member_path = (extract_dir / member.name).resolve()
                    if not str(member_path).startswith(str(extract_dir.resolve())):
                        raise HTTPException(status_code=400, detail=f"压缩包含非法路径: {member.name}")
                tf.extractall(extract_dir)
        else:
            raise HTTPException(status_code=400, detail="文件格式不支持，请上传 zip 或 tar.gz 压缩包")
    except (zipfile.BadZipFile, tarfile.TarError) as e:
        shutil.rmtree(tmp_dir, ignore_errors=True)
        raise HTTPException(status_code=400, detail=f"解压文件失败: {e}")

    # 如果解压后只有一个子目录，进入该子目录
    children = list(extract_dir.iterdir())
    if len(children) == 1 and children[0].is_dir():
        return str(children[0]), tmp_dir
    return str(extract_dir), tmp_dir


def _make_analyze_response(analysis_type: str, summary: dict, html_path: str | None) -> AnalyzeResponse:
    html_url = None
    if html_path:
        from .Config import get_static_url
        html_url = get_static_url(html_path)
    return AnalyzeResponse(
        success=True,
        analysis_type=analysis_type,
        summary=summary,
        html_report_url=html_url,
    )


@app.post("/vasp/analyze/optimization", response_model=AnalyzeResponse)
async def analyze_structure_optimization(request: AnalyzeRequest):
    """独立分析结构优化结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.optimization import OUTCARAnalyzer, generate_optimization_report
        analyzer = OUTCARAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_optimization_report(work_dir)
        summary = {
            "force_convergence": data.get("convergence_analysis", {}).get("force_converged"),
            "final_max_force": data.get("convergence_analysis", {}).get("final_max_force"),
            "energy_convergence": data.get("convergence_analysis", {}).get("energy_converged"),
            "final_energy": data.get("final_results", {}).get("total_energy"),
        }
        return _make_analyze_response("optimization", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析结构优化失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="optimization", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/scf", response_model=AnalyzeResponse)
async def analyze_scf(request: AnalyzeRequest):
    """独立分析自洽场计算结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.scf import SCFAnalyzer, generate_scf_report
        analyzer = SCFAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_scf_report(work_dir)
        summary = {
            "convergence": data.get("electronic_convergence", {}).get("converged"),
            "total_energy": data.get("final_results", {}).get("total_energy"),
            "fermi_energy": data.get("electronic_structure", {}).get("fermi_energy"),
            "band_gap": data.get("electronic_structure", {}).get("band_gap"),
        }
        return _make_analyze_response("scf", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析SCF失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="scf", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/dos", response_model=AnalyzeResponse)
async def analyze_dos(request: AnalyzeRequest):
    """独立分析态密度计算结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.dos import PyMatGenDOSAnalyzer, generate_pymatgen_dos_report
        analyzer = PyMatGenDOSAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_pymatgen_dos_report(work_dir)
        summary = {
            "band_gap": data.get("final_results", {}).get("band_gap"),
            "fermi_energy": data.get("final_results", {}).get("fermi_energy"),
            "is_metal": data.get("final_results", {}).get("is_metal"),
        }
        return _make_analyze_response("dos", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析DOS失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="dos", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/band-structure", response_model=AnalyzeResponse)
async def analyze_band_structure(request: AnalyzeRequest):
    """独立分析能带结构计算结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.band_structure import BandStructureAnalyzer, generate_band_structure_report
        analyzer = BandStructureAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_band_structure_report(work_dir)
        bs_data = data.get("band_structure", {})
        results = data.get("final_results", {})
        summary = {
            "band_gap": bs_data.get("band_gap") or results.get("band_gap"),
            "is_direct": bs_data.get("is_direct"),
            "is_metal": bs_data.get("is_metal"),
            "vbm": results.get("vbm"),
            "cbm": results.get("cbm"),
            "fermi_energy": results.get("fermi_energy"),
        }
        return _make_analyze_response("band_structure", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析能带结构失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="band_structure", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/md", response_model=AnalyzeResponse)
async def analyze_md(request: AnalyzeRequest):
    """独立分析单个MD任务结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.md import generate_md_analysis_report
        html_path = generate_md_analysis_report(work_dir)
        from .analyzers.md import VASP_MDAnalyzer
        analyzer = VASP_MDAnalyzer(work_dir)
        data = analyzer.analyze()
        summary = {
            "convergence": data.get("final_results", {}).get("convergence"),
            "final_energy": data.get("final_results", {}).get("final_energy"),
            "average_temperature": data.get("final_results", {}).get("average_temperature"),
            "total_md_steps": data.get("final_results", {}).get("total_md_steps"),
        }
        return _make_analyze_response("md", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析MD失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="md", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/md-multi", response_model=AnalyzeResponse)
async def analyze_md_multi(request: MDMultiAnalyzeRequest):
    """多任务MD聚合分析 — 接受多个已完成的MD task_id，汇聚到临时目录做Arrhenius等聚合分析"""
    import tempfile
    tmp_root = None
    try:
        # 收集各任务的工作目录
        task_dirs = []
        for tid in request.task_ids:
            task = task_manager.get_task(tid, request.user_id)
            if not task:
                raise HTTPException(status_code=404, detail=f"任务 {tid} 未找到或无权限访问")
            if str(task.status) != "completed":  # type: ignore
                raise HTTPException(status_code=400, detail=f"任务 {tid} 尚未完成 (status={task.status})")  # type: ignore

            # 获取工作目录
            from .task_manager.database import SessionLocal
            from .task_manager.models import ExecutionAttempt
            db = SessionLocal()
            try:
                attempt = (
                    db.query(ExecutionAttempt)
                    .filter(ExecutionAttempt.task_id == tid, ExecutionAttempt.status == "succeeded")
                    .order_by(ExecutionAttempt.attempt_no.desc())
                    .first()
                )
                if attempt and attempt.work_directory:
                    task_dirs.append((tid, str(attempt.work_directory)))
                else:
                    raise HTTPException(status_code=400, detail=f"任务 {tid} 无可用的工作目录")
            finally:
                db.close()

        if not task_dirs:
            raise HTTPException(status_code=400, detail="没有可分析的任务")

        # 创建临时聚合目录：将各任务目录软链接为 T_{temperature}K 子目录
        tmp_root = tempfile.mkdtemp(prefix="md_multi_")
        aggregate_dir = Path(tmp_root) / "aggregate"
        aggregate_dir.mkdir()

        for tid, work_dir_str in task_dirs:
            work_dir_path = Path(work_dir_str)
            # 从 INCAR 读取温度来命名子目录
            incar_path = work_dir_path / "INCAR"
            temp_label = tid[:8]  # fallback
            if incar_path.exists():
                try:
                    from pymatgen.io.vasp.inputs import Incar
                    incar = Incar.from_file(str(incar_path))
                    tebeg = incar.get("TEBEG", None)
                    if tebeg is not None:
                        temp_label = f"T_{float(tebeg)}K"
                except Exception:
                    pass
            link_path = aggregate_dir / temp_label
            if link_path.exists():
                # 避免重名，追加 task_id 前缀
                link_path = aggregate_dir / f"{temp_label}_{tid[:8]}"
            link_path.symlink_to(work_dir_path)

        # 执行聚合分析
        from .analyzers.md import generate_md_analysis_report
        html_path = generate_md_analysis_report(str(aggregate_dir), task_id="md_multi")

        from .analyzers.md import VASP_MDAnalyzer
        analyzer = VASP_MDAnalyzer(str(aggregate_dir))
        data = analyzer.analyze()

        summary = {
            "task_count": len(task_dirs),
            "task_ids": request.task_ids,
            "arrhenius": data.get("arrhenius"),
        }

        # 复制 HTML 报告到持久路径（不放在 tmp 里）
        from .settings import settings as _settings
        persistent_dir = Path(_settings.vasp_work_root) / request.user_id / "md_multi_analysis"
        persistent_dir.mkdir(parents=True, exist_ok=True)
        import shutil as _shutil
        persistent_html = persistent_dir / "md_multi_analysis_report.html"
        _shutil.copy2(html_path, str(persistent_html))
        from .Config import get_static_url
        report_url = get_static_url(str(persistent_html))

        return _make_analyze_response("md_multi", summary, report_url)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"多任务MD聚合分析失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="md_multi", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/neb", response_model=AnalyzeResponse)
async def analyze_neb(request: AnalyzeRequest):
    """独立分析 NEB 计算结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.neb import NEBAnalyzer, generate_neb_report
        analyzer = NEBAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_neb_report(work_dir)
        neb = data.get("neb", {})
        summary = {
            "forward_barrier_eV": neb.get("forward_barrier_eV"),
            "backward_barrier_eV": neb.get("backward_barrier_eV"),
            "reaction_energy_eV": neb.get("reaction_energy_eV"),
            "n_images": neb.get("n_images"),
            "ts_image_index": neb.get("ts_image_index"),
        }
        return _make_analyze_response("neb", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析 NEB 失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="neb", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/analyze/phonon", response_model=AnalyzeResponse)
async def analyze_phonon(request: AnalyzeRequest):
    """独立分析声子计算结果 — 传入 task_id 或 file_url"""
    work_dir, tmp_root = _resolve_analysis_input_dir(request)
    try:
        from .analyzers.phonon import PhononAnalyzer, generate_phonon_report
        analyzer = PhononAnalyzer(work_dir)
        data = analyzer.analyze()
        html_path = generate_phonon_report(work_dir)
        phonon = data.get("phonon", {})
        summary = {
            "n_modes": phonon.get("n_modes"),
            "n_imaginary": phonon.get("n_imaginary"),
            "dynamically_stable": phonon.get("dynamically_stable"),
            "max_imaginary_freq_cm1": phonon.get("max_imaginary_freq_cm1"),
            "max_real_freq_cm1": phonon.get("max_real_freq_cm1"),
        }
        return _make_analyze_response("phonon", summary, html_path)
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"分析声子失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="phonon", error_message=str(e))
    finally:
        if tmp_root:
            shutil.rmtree(tmp_root, ignore_errors=True)


@app.post("/vasp/agent/analyze", response_model=AgentAnalyzeResponse)
async def agent_analyze(request: AgentAnalyzeRequest):
    """
    使用 AI Agent（LangGraph + Python REPL + skills 库）分析任意 VASP 计算结果。
    支持自然语言提问，Agent 自主编写和执行 Python 代码完成分析。
    """
    try:
        task = task_manager.get_task(request.task_id, request.user_id)
        if not task:
            raise HTTPException(status_code=404, detail=f"任务 {request.task_id} 未找到")

        # 获取工作目录
        work_dir = getattr(task, "result_path", None)
        if not work_dir:
            result_data = getattr(task, "result_data", None) or {}
            work_dir = result_data.get("work_directory")
        if not work_dir or not Path(work_dir).exists():
            raise HTTPException(
                status_code=400,
                detail="任务工作目录不存在或任务尚未完成，请确认任务状态为 completed",
            )

        from .analysis_agent import run_analysis
        from .Config import get_static_url

        result = await run_analysis(work_dir, request.question, request.model)

        # 将生成的 PNG 文件名转换为可访问 URL
        plot_urls = [
            get_static_url(str(Path(work_dir) / fname))
            for fname in result.get("generated_plots", [])
        ]

        return AgentAnalyzeResponse(
            success=True,
            answer=result["answer"],
            steps=result["steps"],
            generated_plots=plot_urls,
        )

    except HTTPException:
        raise
    except Exception as e:
        logger.error("Agent 分析失败: %s", e, exc_info=True)
        return AgentAnalyzeResponse(success=False, error_message=str(e))


if __name__ == "__main__":
    from .Config import VASP_remote_run_port
    uvicorn.run(app, host="0.0.0.0", port=VASP_remote_run_port)
