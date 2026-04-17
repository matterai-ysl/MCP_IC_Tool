from fastapi import FastAPI, HTTPException, Query
from fastapi.responses import FileResponse, JSONResponse
from fastapi.middleware.cors import CORSMiddleware
from typing import Any, Dict, List, Optional
from urllib.parse import urlsplit, urlunsplit
import csv
from html import escape
import json
import math
import os
import uvicorn
import logging
import sys
import uuid
from pathlib import Path
# Config import moved to __main__ block
from .failure_guidance import derive_failure_guidance, derive_failure_guidance_from_workdir
from .internal_worker_api import build_internal_worker_router
from .errors import APIError, build_api_error
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
from .settings import settings
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


@app.exception_handler(APIError)
async def handle_api_error(_request, exc: APIError):
    return JSONResponse(status_code=exc.status_code, content={"error": exc.to_dict()})


def _wrap_unexpected_api_error(exc: Exception, *, code: str, message: str) -> APIError:
    if isinstance(exc, APIError):
        return exc
    return build_api_error(
        status_code=500,
        code=code,
        message=message,
        retryable=False,
        details={"reason": str(exc)},
        suggested_action="请稍后重试；如果问题持续存在，请联系管理员排查服务日志",
    )


def _raise_validation_error(message: str, *, code: str, details: Dict[str, Any] | None = None) -> None:
    raise build_api_error(
        status_code=400,
        code=code,
        message=message,
        retryable=False,
        details=details,
    )


def _raise_not_found_error(message: str, *, code: str, details: Dict[str, Any] | None = None) -> None:
    raise build_api_error(
        status_code=404,
        code=code,
        message=message,
        retryable=False,
        details=details,
    )

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


@app.get("/dft/{file_path:path}")
async def serve_public_artifact(file_path: str):
    root = Path(settings.local_public_artifact_root).resolve()
    target = (root / file_path).resolve()
    if root not in target.parents and target != root:
        raise HTTPException(status_code=404, detail="文件不存在")
    if not target.exists() or not target.is_file():
        raise HTTPException(status_code=404, detail="文件不存在")
    return FileResponse(target)


def _record_task_submitted(task_type: str):
    """Record a task submission in Prometheus metrics."""
    if _METRICS_ENABLED:
        vasp_tasks_total.labels(task_type=task_type).inc()
        vasp_tasks_queued.inc()


def _require_completed_upstream_task(task_id: str, user_id: str, task_label: str):
    task = task_manager.get_task(task_id, user_id)
    if not task:
        raise HTTPException(
            status_code=404,
            detail=f"{task_label} {task_id} 未找到或无权限访问",
        )
    if str(task.status) != "completed":
        raise HTTPException(
            status_code=400,
            detail=f"{task_label} {task_id} 尚未完成",
        )
    return task


def _resolve_queue_name(requested_queue_name: str | None, upstream_tasks: List[Any] | None = None) -> str:
    requested = str(requested_queue_name).strip() if requested_queue_name else None
    upstream_queue_names = {
        str(getattr(task, "queue_name", "default") or "default")
        for task in (upstream_tasks or [])
        if task is not None
    }
    if not upstream_queue_names:
        return requested or "default"

    if len(upstream_queue_names) > 1:
        raise HTTPException(
            status_code=400,
            detail="上游任务位于不同 queue_name，无法在不同超算之间直接续算",
        )

    upstream_queue_name = next(iter(upstream_queue_names))
    if requested and requested != upstream_queue_name:
        raise HTTPException(
            status_code=400,
            detail=f"queue_name 必须与上游任务一致: {upstream_queue_name}",
        )
    return upstream_queue_name


@app.post("/vasp/structure-optimization", response_model=StructOptResponse)
async def submit_structure_optimization(request: StructOptRequest):
    """提交结构优化任务。"""
    try:
        task_params = {
            "cif_url": str(request.cif_url),
            "queue_name": _resolve_queue_name(request.queue_name),
            "kpoint_density": request.kpoint_density,
        }
        if request.custom_incar:
            task_params["custom_incar"] = request.custom_incar
        if request.user_id is None:
            request.user_id = os.getenv("DEFAULT_USER_ID", "123")
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
        
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="STRUCTURE_OPTIMIZATION_SUBMIT_FAILED",
            message="提交结构优化任务失败",
        )


@app.post("/vasp/scf-calculation", response_model=SCFResponse)
async def submit_scf_calculation(request: SCFRequest):
    """提交自洽场计算任务。"""
    try:
        input_count = sum([bool(request.cif_url), bool(request.optimized_task_id)])
        if input_count != 1:
            _raise_validation_error(
                "必须提供 cif_url 或 optimized_task_id 中的一个",
                code="INVALID_SCF_INPUT_SOURCE",
                details={"fields": ["cif_url", "optimized_task_id"]},
            )
        
        # 如果基于优化任务，验证任务存在性
        upstream_task = None
        if request.optimized_task_id:
            upstream_task = _require_completed_upstream_task(
                request.optimized_task_id,
                request.user_id,
                "结构优化任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])
        
        task_params = {
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "optimized_task_id": request.optimized_task_id,
            "queue_name": queue_name,
            "kpoint_density": request.kpoint_density,
            "precision": request.precision,
        }
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="scf_calculation",
            params=task_params
        )
        _record_task_submitted("scf_calculation")

        input_source = "结构 URL" if request.cif_url else "结构优化任务"
        
        return SCFResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"自洽场计算任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="SCF_SUBMIT_FAILED",
            message="提交自洽场计算任务失败",
        )


@app.post("/vasp/dos-calculation", response_model=DOSResponse)
async def submit_dos_calculation(request: DOSRequest):
    """提交态密度计算任务。"""
    try:
        input_count = sum([bool(request.input_url), bool(request.scf_task_id)])
        if input_count != 1:
            _raise_validation_error(
                "必须提供 input_url 或 scf_task_id 中的一个",
                code="INVALID_DOS_INPUT_SOURCE",
                details={"fields": ["input_url", "scf_task_id"]},
            )
        
        # 如果基于自洽场任务，验证任务存在性
        upstream_task = None
        if request.scf_task_id:
            upstream_task = _require_completed_upstream_task(
                request.scf_task_id,
                request.user_id,
                "自洽场计算任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])
        
        task_params = {
            "input_url": (
                [str(item) for item in request.input_url]
                if isinstance(request.input_url, list)
                else (str(request.input_url) if request.input_url else None)
            ),
            "scf_task_id": request.scf_task_id,
            "queue_name": queue_name,
            "kpoint_density": request.kpoint_density,
            "kpoint_multiplier": request.kpoint_multiplier,
            "precision": request.precision,
        }
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="dos_calculation",
            params=task_params
        )
        _record_task_submitted("dos_calculation")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "态密度计算"
        else:
            input_source = "输入 URL"
            calc_mode = "单点自洽+DOS计算"
        
        return DOSResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="DOS_SUBMIT_FAILED",
            message="提交态密度计算任务失败",
        )


@app.post("/vasp/band-structure", response_model=BandStructureResponse)
async def submit_band_structure_calculation(request: BandStructureRequest):
    """提交能带结构计算任务。"""
    try:
        input_url = request.input_url
        input_count = sum([bool(input_url), bool(request.scf_task_id)])
        if input_count != 1:
            _raise_validation_error(
                "必须提供 input_url 或 scf_task_id 中的一个",
                code="INVALID_BAND_INPUT_SOURCE",
                details={"fields": ["input_url", "scf_task_id"]},
            )

        upstream_task = None
        if request.scf_task_id:
            upstream_task = _require_completed_upstream_task(
                request.scf_task_id,
                request.user_id,
                "自洽场计算任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])

        task_params = {
            "input_url": (
                [str(item) for item in input_url]
                if isinstance(input_url, list)
                else (str(input_url) if input_url else None)
            ),
            "scf_task_id": request.scf_task_id,
            "queue_name": queue_name,
            "kpoint_density": request.kpoint_density,
            "line_density": request.line_density,
            "precision": request.precision,
        }

        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="band_structure",
            params=task_params
        )
        _record_task_submitted("band_structure")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
        else:
            input_source = "输入 URL"

        return BandStructureResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"能带结构计算任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )

    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="BAND_STRUCTURE_SUBMIT_FAILED",
            message="提交能带结构计算任务失败",
        )


@app.post("/vasp/md-calculation", response_model=MDResponse)
async def submit_md_calculation(request: MDRequest):
    """提交分子动力学计算任务。"""
    try:
        input_count = sum([bool(request.input_url), bool(request.scf_task_id)])
        if input_count != 1:
            _raise_validation_error(
                "必须提供 input_url 或 scf_task_id 中的一个",
                code="INVALID_MD_INPUT_SOURCE",
                details={"fields": ["input_url", "scf_task_id"]},
            )
        
        # 如果基于自洽场任务，验证任务存在性
        upstream_task = None
        if request.scf_task_id:
            upstream_task = _require_completed_upstream_task(
                request.scf_task_id,
                request.user_id,
                "自洽场计算任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])
        
        task_params = {
            "input_url": str(request.input_url) if request.input_url else None,
            "scf_task_id": request.scf_task_id,
            "queue_name": queue_name,
            "md_steps": request.md_steps,
            "temperature": float(request.temperature),
            "time_step": request.time_step,
            "ensemble": request.ensemble,
            "precision": request.precision,
        }
        task_id = task_manager.submit_task(
            user_id=request.user_id,
            task_type="md_calculation",
            params=task_params
        )
        _record_task_submitted("md_calculation")

        if request.scf_task_id:
            input_source = "自洽场计算任务"
            calc_mode = "分子动力学计算"
        else:
            input_source = "输入 URL"
            calc_mode = "纯MD计算"
        
        return MDResponse(
            task_id=task_id,
            status=TaskStatus.queued,
            message=f"{calc_mode}任务已提交，输入源：{input_source}，任务ID: {task_id}"
        )
        
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="MD_SUBMIT_FAILED",
            message="提交分子动力学计算任务失败",
        )

@app.post("/vasp/neb-calculation", response_model=NEBResponse)
async def submit_neb_calculation(request: NEBRequest):
    """提交 NEB（过渡态）计算任务"""
    try:
        upstream_tasks = []
        if request.initial_task_id:
            upstream_tasks.append(_require_completed_upstream_task(
                request.initial_task_id,
                request.user_id,
                "初始态任务",
            ))
        if request.final_task_id:
            upstream_tasks.append(_require_completed_upstream_task(
                request.final_task_id,
                request.user_id,
                "终态任务",
            ))
        queue_name = _resolve_queue_name(request.queue_name, upstream_tasks)

        task_params = {
            "initial_cif_url": str(request.initial_cif_url) if request.initial_cif_url else None,
            "initial_task_id": request.initial_task_id,
            "final_cif_url": str(request.final_cif_url) if request.final_cif_url else None,
            "final_task_id": request.final_task_id,
            "queue_name": queue_name,
            "n_images": request.n_images,
            "kpoint_density": request.kpoint_density,
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
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="NEB_SUBMIT_FAILED",
            message="提交 NEB 计算任务失败",
        )


@app.post("/vasp/phonon-calculation", response_model=PhononResponse)
async def submit_phonon_calculation(request: PhononRequest):
    """提交声子计算任务（IBRION=6，Gamma 点有限位移法）"""
    try:
        upstream_task = None
        if request.scf_task_id:
            upstream_task = _require_completed_upstream_task(
                request.scf_task_id,
                request.user_id,
                "上游任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])
        task_params = {
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "scf_task_id": request.scf_task_id,
            "queue_name": queue_name,
            "kpoint_density": request.kpoint_density,
            "displacement": request.displacement,
        }
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
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="PHONON_SUBMIT_FAILED",
            message="提交声子计算任务失败",
        )


@app.post("/vasp/custom-calculation", response_model=CustomCalcResponse)
async def submit_custom_calculation(request: CustomCalcRequest):
    """提交通用自定义计算任务 — 用户完全控制INCAR，适合HSE06、DFPT介电、ELF、Wannier等长尾需求"""
    try:
        upstream_task = None
        if request.from_task_id:
            upstream_task = _require_completed_upstream_task(
                request.from_task_id,
                request.user_id,
                "上游任务",
            )
        queue_name = _resolve_queue_name(request.queue_name, [upstream_task] if upstream_task else [])
        task_params: Dict[str, Any] = {
            "cif_url": str(request.cif_url) if request.cif_url else None,
            "from_task_id": request.from_task_id,
            "queue_name": queue_name,
            "incar": request.incar,
            "kpoint_density": request.kpoint_density,
            "kpoint_mode": request.kpoint_mode,
        }

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
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="CUSTOM_CALC_SUBMIT_FAILED",
            message="提交自定义计算任务失败",
        )


@app.get("/vasp/task/{task_id}", response_model=TaskStatusResponse)
async def get_task_status(task_id: str, user_id: str = Query(..., description="用户ID")):
    """查询任务状态与任务结果"""
    try:
        task = task_manager.get_task(task_id, user_id)
        if not task:
            _raise_not_found_error(
                "任务未找到或无权限访问",
                code="TASK_NOT_FOUND",
                details={"task_id": task_id},
            )

        return _build_task_response(task, task_id)

    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        logger.error(f"查询任务状态失败: {e}", exc_info=True)
        raise _wrap_unexpected_api_error(
            e,
            code="TASK_STATUS_FETCH_FAILED",
            message="查询任务状态失败",
        )


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
    if not isinstance(value, str):
        return False
    if value.startswith(("https://", "http://")):
        parsed = urlsplit(value)
        if parsed.path.startswith(("/static/", "/download/file/")):
            return False
        return True
    if settings.legacy_local_artifact_urls_enabled and value.startswith(("/static/", "/download/file/")):
        return True
    return False


_LEGACY_PUBLIC_HOSTS = {"api.matterai.tech"}


def _normalize_public_url(value: str):
    if not isinstance(value, str):
        return value
    if not value.startswith(("https://", "http://")):
        return value

    parsed = urlsplit(value)
    if parsed.netloc not in _LEGACY_PUBLIC_HOSTS:
        return value
    if not parsed.path.startswith(("/download/file/", "/static/")):
        return value

    target = urlsplit(settings.vasp_public_base_url)
    scheme = target.scheme or parsed.scheme
    netloc = target.netloc or parsed.netloc
    return urlunsplit((scheme, netloc, parsed.path, parsed.query, parsed.fragment))


def _sanitize_public_value(value):
    if isinstance(value, dict):
        return {key: _sanitize_public_value(val) for key, val in value.items()}
    if isinstance(value, list):
        return [_sanitize_public_value(item) for item in value]
    if isinstance(value, str):
        value = _normalize_public_url(value)
        if value.startswith("/"):
            return None
        return value
    return value


def _sanitize_report_urls(result_summary, result_data, html_report_url):
    sanitized_summary = _sanitize_public_value(result_summary)
    sanitized_data = _sanitize_public_value(result_data)

    if isinstance(sanitized_summary, dict):
        sanitized_summary.pop("html_report_url", None)

    report_keys = [
        "analysis_report_html_path",
        "html_analysis_report",
        "md_analysis_report_html_path",
        "dos_analysis_report_html_path",
        "scf_analysis_report_html_path",
        "band_structure_report_html_path",
        "neb_report_html_path",
        "phonon_report_html_path",
        "agent_analysis_report_html_path",
    ]
    if isinstance(sanitized_data, dict):
        sanitized_data.pop("html_report_url", None)
        for key in report_keys:
            sanitized_data.pop(key, None)
        for key in ("work_directory", "source_work_directory"):
            sanitized_data.pop(key, None)
        # Large analyzer payloads belong in the HTML report and downloadable CSV/JSON
        # artifacts, not in the task status JSON returned to agents/frontends.
        for key in ("analysis_data", "dos_data", "raw_dos_data", "trajectory_data"):
            sanitized_data.pop(key, None)

    return sanitized_summary, sanitized_data


def _apply_artifact_download_overrides(result_data, artifacts):
    if not isinstance(result_data, dict) or not artifacts:
        return result_data

    by_filename = {}
    for artifact in artifacts:
        download_url = getattr(artifact, "download_url", None)
        filename = getattr(artifact, "filename", None)
        if not download_url or not filename:
            continue
        lowered = str(filename).lower()
        by_filename.setdefault(lowered, download_url)
        by_filename.setdefault(os.path.basename(lowered), download_url)

    if not by_filename:
        return result_data

    def _rewrite(value):
        if isinstance(value, dict):
            return {key: _rewrite(val) for key, val in value.items()}
        if isinstance(value, list):
            return [_rewrite(item) for item in value]
        if isinstance(value, str):
            basename = os.path.basename(urlsplit(value).path).lower()
            if basename in by_filename:
                return by_filename[basename]
        return value

    return _rewrite(result_data)


def _public_url_for_local_path(path: str | None) -> Optional[str]:
    if not path:
        return None
    candidate = Path(path).resolve()
    public_root = Path(settings.local_public_artifact_root).resolve()
    try:
        rel = candidate.relative_to(public_root).as_posix()
    except ValueError:
        from .Config import get_static_url
        return get_static_url(str(candidate))
    return f"{settings.public_artifact_base_url.rstrip('/')}/{rel}"


_IMAGE_EXTENSIONS = (".png", ".jpg", ".jpeg", ".svg", ".webp")
_DATA_EXTENSIONS = (".csv", ".json")


def _artifact_display_filename(artifact) -> Optional[str]:
    raw_path = getattr(artifact, "object_key", None) or getattr(artifact, "storage_key", None)
    if not raw_path:
        return None
    text = str(raw_path).strip("/")
    marker = "/attempts/"
    if marker in text:
        suffix = text.split(marker, 1)[1]
        parts = suffix.split("/", 1)
        if len(parts) == 2:
            return parts[1]
    return os.path.basename(text) or text


def _is_html_artifact(artifact_type: str, filename: Optional[str], content_type: Optional[str]) -> bool:
    if artifact_type == "html_report":
        return True
    if content_type and content_type.startswith("text/html"):
        return True
    return bool(filename and filename.lower().endswith(".html"))


def _is_preview_image_artifact(filename: Optional[str], content_type: Optional[str]) -> bool:
    if content_type and content_type.startswith("image/"):
        return True
    return bool(filename and filename.lower().endswith(_IMAGE_EXTENSIONS))


def _is_data_download_artifact(filename: Optional[str], content_type: Optional[str]) -> bool:
    if content_type in {"text/csv", "application/json"}:
        return True
    return bool(filename and filename.lower().endswith(_DATA_EXTENSIONS))


def _derive_task_suggested_action(task, task_id: str) -> Optional[str]:
    guidance = _derive_task_failure_guidance(task, task_id)
    return guidance.suggested_action


def _derive_task_failure_type(task, task_id: str) -> Optional[str]:
    guidance = _derive_task_failure_guidance(task, task_id)
    return guidance.failure_type


def _derive_task_failure_guidance(task, task_id: str):
    fragments: list[str] = []
    task_type = getattr(task, "task_type", None)
    error_message = getattr(task, "error_message", None)
    if error_message:
        fragments.append(str(error_message))

    try:
        attempt = task_manager.get_execution_attempt(task_id)
    except Exception:
        attempt = None
    if attempt is not None:
        failure_detail = getattr(attempt, "failure_detail", None)
        if failure_detail:
            fragments.append(str(failure_detail))

    guidance = derive_failure_guidance(
        "\n".join(fragments),
        task_type=task_type,
    )
    if guidance.failure_type:
        return guidance

    if attempt is not None:
        work_directory = getattr(attempt, "work_directory", None)
        if work_directory:
            try:
                workdir_guidance = derive_failure_guidance_from_workdir(
                    Path(str(work_directory)),
                    task_type=task_type,
                )
            except OSError:
                workdir_guidance = None
            if workdir_guidance and workdir_guidance.failure_type:
                return workdir_guidance

    return guidance


def _build_task_response(task, task_id: str) -> TaskStatusResponse:
    """从 Task ORM 构建统一的 TaskStatusResponse"""
    guidance = _derive_task_failure_guidance(task, task_id)

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
                        'scf_analysis_report_html_path', 'band_structure_report_html_path',
                        'agent_analysis_report_html_path']:
                val = result_data.get(key)
                if val:
                    html_report_url = str(val)
                    break

    # 查询 artifacts
    artifacts_list = None
    preview_images = None
    data_downloads = None
    try:
        raw_artifacts = task_manager.get_task_artifacts(task_id)
        if raw_artifacts:
            artifacts_list = []
            preview_images = []
            data_downloads = []
            for a in raw_artifacts:
                download_url = task_manager.get_artifact_download_url(a)
                content_type = getattr(a, 'content_type', None) or getattr(a, 'mime_type', None)
                filename = _artifact_display_filename(a)
                artifact_info = ArtifactInfo(
                    id=str(a.id),
                    artifact_type=str(a.artifact_type),
                    filename=filename,
                    mime_type=getattr(a, 'mime_type', None),
                    content_type=content_type,
                    size_bytes=getattr(a, 'size_bytes', None),
                    download_url=download_url,
                )
                artifact_type = str(getattr(a, 'artifact_type', ''))
                if _is_html_artifact(artifact_type, filename, content_type) and download_url:
                    html_report_url = download_url
                is_report_asset = (
                    _is_html_artifact(artifact_type, filename, content_type)
                    or _is_preview_image_artifact(filename, content_type)
                    or _is_data_download_artifact(filename, content_type)
                )
                if download_url and not is_report_asset:
                    artifacts_list.append(artifact_info)
                if _is_preview_image_artifact(filename, content_type) and download_url:
                    preview_images.append(artifact_info)
                if _is_data_download_artifact(filename, content_type) and download_url:
                    data_downloads.append(artifact_info)
            if not artifacts_list:
                artifacts_list = None
            if not preview_images:
                preview_images = None
            if not data_downloads:
                data_downloads = None
    except Exception:
        pass  # artifact 查询失败不影响主响应

    if not _is_public_url(html_report_url):
        html_report_url = None

    result_data = getattr(task, 'result_data', None)
    result_summary, result_data = _sanitize_report_urls(result_summary, result_data, html_report_url)
    result_data = _apply_artifact_download_overrides(result_data, artifacts_list)
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
        "preview_images": preview_images,
        "data_downloads": data_downloads,
        "progress_message": getattr(task, 'progress_message', None),
        "params": getattr(task, 'params', None),
        "result_path": result_path,
        "external_job_id": getattr(task, 'external_job_id', None),
        "process_id": getattr(task, 'process_id', None),
        "error_message": getattr(task, 'error_message', None),
        "failure_type": guidance.failure_type,
        "suggested_action": guidance.suggested_action,
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
            _raise_not_found_error(
                "任务未找到或无法取消",
                code="TASK_CANCEL_NOT_FOUND",
                details={"task_id": task_id},
            )
        
        return {"message": f"任务 {task_id} 已请求取消"}
        
    except APIError:
        raise
    except HTTPException:
        raise
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="TASK_CANCEL_FAILED",
            message="取消任务失败",
        )


@app.get("/vasp/tasks", response_model=List[TaskStatusResponse])
async def list_user_tasks(user_id: str = Query(..., description="用户ID")):
    """列出用户的所有任务"""
    try:
        tasks = task_manager.list_tasks(user_id)
        return [_build_task_response(task, str(task.id)) for task in tasks]
    except Exception as e:
        raise _wrap_unexpected_api_error(
            e,
            code="TASK_LIST_FAILED",
            message="获取任务列表失败",
        )


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


ANALYSIS_TYPE_TO_TASK_TYPE = {
    "dos": "dos_analysis",
    "band_structure": "band_structure_analysis",
    "md": "md_analysis",
    "agent": "agent_analysis",
}


def _get_completed_source_task(task_id: str, user_id: str):
    task = task_manager.get_task(task_id, user_id)
    if task is None:
        raise HTTPException(status_code=404, detail=f"任务 {task_id} 未找到或无权限访问")
    if str(getattr(task, "status", "")) != "completed":
        raise HTTPException(status_code=400, detail=f"任务 {task_id} 尚未完成 (status={getattr(task, 'status', '')})")
    return task


def _submit_remote_analysis_task(
    *,
    analysis_type: str,
    user_id: str,
    source_task_id: str,
    extra_params: Optional[dict[str, Any]] = None,
):
    source_task = _get_completed_source_task(source_task_id, user_id)
    task_type = ANALYSIS_TYPE_TO_TASK_TYPE[analysis_type]
    source_result_data = getattr(source_task, "result_data", None) or {}
    source_work_directory = getattr(source_task, "result_path", None) or source_result_data.get("work_directory")
    params = {
        "source_task_id": source_task_id,
        "source_task_type": getattr(source_task, "task_type", None),
        "source_work_directory": source_work_directory,
        "queue_name": getattr(source_task, "queue_name", None) or "default",
    }
    if extra_params:
        params.update(extra_params)

    analysis_task_id = task_manager.submit_task(user_id, task_type, params)
    analysis_task = task_manager.get_task(analysis_task_id, user_id)
    if analysis_task is None:
        raise HTTPException(status_code=500, detail=f"分析任务 {analysis_task_id} 创建后无法读取")
    return analysis_task


def _analysis_task_message(status: str, analysis_type: str) -> str:
    if status == "completed":
        return f"{analysis_type} 分析已完成"
    if status == "failed":
        return f"{analysis_type} 分析失败"
    if status == "running":
        return f"{analysis_type} 分析任务正在执行"
    return f"{analysis_type} 分析任务已提交"


def _build_remote_analysis_response(task, analysis_type: str) -> AnalyzeResponse:
    task_id = str(getattr(task, "id", ""))
    status = _safe_task_status(getattr(task, "status", None))
    task_response = _build_task_response(task, task_id)
    success = status != "failed"
    return AnalyzeResponse(
        success=success,
        analysis_type=analysis_type,
        analysis_task_id=task_id,
        status=status,
        message=_analysis_task_message(str(status), analysis_type),
        summary=task_response.result_summary,
        html_report_url=task_response.html_report_url,
        error_message=getattr(task, "error_message", None) if status == "failed" else None,
        failure_type=task_response.failure_type if status == "failed" else None,
        suggested_action=task_response.suggested_action if status == "failed" else None,
    )


def _build_remote_agent_response(task) -> AgentAnalyzeResponse:
    task_id = str(getattr(task, "id", ""))
    status = _safe_task_status(getattr(task, "status", None))
    task_response = _build_task_response(task, task_id)
    result_data = getattr(task, "result_data", None) or {}
    generated_plots = [
        artifact.download_url
        for artifact in (task_response.preview_images or [])
        if artifact.download_url
    ]
    success = status != "failed"
    return AgentAnalyzeResponse(
        success=success,
        analysis_task_id=task_id,
        status=status,
        message=_analysis_task_message(str(status), "agent"),
        answer=str(result_data.get("answer") or ""),
        steps=int(result_data.get("steps") or 0),
        generated_plots=generated_plots,
        html_report_url=task_response.html_report_url,
        error_message=getattr(task, "error_message", None) if status == "failed" else None,
        failure_type=task_response.failure_type if status == "failed" else None,
        suggested_action=task_response.suggested_action if status == "failed" else None,
    )


def _make_analyze_response(analysis_type: str, summary: dict, html_path: str | None) -> AnalyzeResponse:
    html_url = _public_url_for_local_path(html_path)
    return AnalyzeResponse(
        success=True,
        analysis_type=analysis_type,
        summary=summary,
        html_report_url=html_url,
    )


def _analysis_error_message(exc: Exception) -> str:
    detail = getattr(exc, "detail", None)
    if detail:
        return str(detail)
    return str(exc)


def _load_local_or_remote_json(local_path: str | None, download_url: str | None) -> Dict[str, Any]:
    if local_path and Path(local_path).exists():
        return json.loads(Path(local_path).read_text(encoding="utf-8"))
    if download_url:
        response = _requests.get(download_url, timeout=30)
        response.raise_for_status()
        return response.json()
    raise FileNotFoundError("未找到可访问的 JSON 产物")


def _load_local_or_remote_text(local_path: str | None, download_url: str | None) -> str:
    if local_path and Path(local_path).exists():
        return Path(local_path).read_text(encoding="utf-8")
    if download_url:
        response = _requests.get(download_url, timeout=30)
        response.raise_for_status()
        response.encoding = response.encoding or "utf-8"
        return response.text
    raise FileNotFoundError("未找到可访问的文本产物")


def _find_task_artifact(task_id: str, *suffixes: str):
    normalized_suffixes = [suffix.lower() for suffix in suffixes]
    for artifact in reversed(task_manager.get_task_artifacts(task_id)):
        display_name = (_artifact_display_filename(artifact) or "").lower()
        object_key = str(getattr(artifact, "object_key", "") or "").lower()
        storage_key = str(getattr(artifact, "storage_key", "") or "").lower()
        if any(
            display_name.endswith(suffix)
            or object_key.endswith(suffix)
            or storage_key.endswith(suffix)
            for suffix in normalized_suffixes
        ):
            return artifact
    return None


def _load_md_diffusion_results(task_id: str) -> Dict[str, Any]:
    artifact = _find_task_artifact(task_id, "md_output/data/diffusion_results.json", "diffusion_results.json")
    if artifact is not None:
        try:
            local_path = getattr(artifact, "storage_key", None)
            download_url = task_manager.get_artifact_download_url(artifact)
            return _load_local_or_remote_json(local_path, download_url)
        except Exception:
            pass

    raise FileNotFoundError(f"任务 {task_id} 缺少可用的 diffusion_results.json")


def _extract_md_diffusion_point(task_id: str, payload: Dict[str, Any]) -> Dict[str, Any]:
    temperature = payload.get("temperature_K")
    try:
        temperature = float(temperature)
    except (TypeError, ValueError):
        raise ValueError(f"任务 {task_id} 的 diffusion_results.json 缺少 temperature_K")

    diffusion_by_element = payload.get("diffusion_by_element") or {}
    conductivity_by_element = payload.get("ionic_conductivity_S_per_m") or {}
    mobile_species = payload.get("mobile_species") or list(diffusion_by_element.keys())

    diffusivities: List[float] = []
    conductivities: List[float] = []
    used_species: List[str] = []
    for specie in mobile_species:
        specie_key = str(specie)
        specie_data = diffusion_by_element.get(specie_key) or {}
        try:
            diffusivity = float(specie_data.get("D_m2_per_s"))
        except (TypeError, ValueError):
            continue
        if diffusivity <= 0:
            continue
        diffusivities.append(diffusivity)
        used_species.append(specie_key)
        try:
            conductivity = float(conductivity_by_element.get(specie_key))
        except (TypeError, ValueError):
            conductivity = None
        if conductivity is not None and conductivity > 0:
            conductivities.append(conductivity)

    if not diffusivities:
        raise ValueError(f"任务 {task_id} 的 diffusion_results.json 缺少有效扩散系数")

    avg_diffusivity = sum(diffusivities) / len(diffusivities)
    avg_conductivity = (sum(conductivities) / len(conductivities)) if conductivities else None
    diffusion_detail = {
        specie: float((diffusion_by_element.get(specie) or {}).get("D_m2_per_s"))
        for specie in used_species
        if (diffusion_by_element.get(specie) or {}).get("D_m2_per_s") is not None
    }
    conductivity_detail = {
        specie: float(conductivity_by_element[specie])
        for specie in used_species
        if conductivity_by_element.get(specie) is not None
    }
    return {
        "task_id": task_id,
        "temperature_K": temperature,
        "mobile_species": used_species or [str(s) for s in mobile_species],
        "diffusivity_m2_per_s": avg_diffusivity,
        "ionic_conductivity_S_per_m": avg_conductivity,
        "diffusion_by_element": diffusion_detail,
        "conductivity_by_element": conductivity_detail,
    }


def _safe_float(value: Any) -> Optional[float]:
    try:
        if value is None or value == "":
            return None
        return float(value)
    except (TypeError, ValueError):
        return None


def _load_md_time_series_rows(task_id: str) -> List[Dict[str, str]]:
    artifact = _find_task_artifact(task_id, "md_output/data/time_series.csv", "time_series.csv")
    if artifact is None:
        return []
    try:
        text = _load_local_or_remote_text(
            getattr(artifact, "storage_key", None),
            task_manager.get_artifact_download_url(artifact),
        )
    except Exception:
        return []
    return list(csv.DictReader(text.splitlines()))


def _load_md_rdf_summary(task_id: str) -> Dict[str, Any]:
    artifact = _find_task_artifact(task_id, "md_output/data/rdf_summary.json", "rdf_summary.json")
    if artifact is None:
        return {}
    try:
        return _load_local_or_remote_json(
            getattr(artifact, "storage_key", None),
            task_manager.get_artifact_download_url(artifact),
        )
    except Exception:
        return {}


def _extract_md_stability_point(task_id: str, rows: List[Dict[str, str]]) -> Optional[Dict[str, Any]]:
    if not rows:
        return None

    temperatures = [_safe_float(row.get("temperature_K")) for row in rows]
    energies = [_safe_float(row.get("energy_eV")) for row in rows]
    pressures = [_safe_float(row.get("pressure_kB")) for row in rows]
    densities = [_safe_float(row.get("density_kg_m3")) for row in rows]
    volumes = [_safe_float(row.get("volume_A3")) for row in rows]

    def _valid(values: List[Optional[float]]) -> List[float]:
        return [float(v) for v in values if v is not None]

    valid_temperatures = _valid(temperatures)
    valid_energies = _valid(energies)
    valid_pressures = _valid(pressures)
    valid_densities = _valid(densities)
    valid_volumes = _valid(volumes)

    def _mean(values: List[float]) -> Optional[float]:
        return (sum(values) / len(values)) if values else None

    def _std(values: List[float]) -> Optional[float]:
        if len(values) < 2:
            return 0.0 if values else None
        mean = sum(values) / len(values)
        return math.sqrt(sum((value - mean) ** 2 for value in values) / len(values))

    energy_drift = None
    if len(valid_energies) >= 2:
        energy_drift = valid_energies[-1] - valid_energies[0]

    return {
        "task_id": task_id,
        "average_temperature_K": _mean(valid_temperatures),
        "temperature_std_K": _std(valid_temperatures),
        "average_pressure_kB": _mean(valid_pressures),
        "average_density_kg_m3": _mean(valid_densities),
        "average_volume_A3": _mean(valid_volumes),
        "energy_drift_eV": energy_drift,
    }


def _extract_md_rdf_point(task_id: str, payload: Dict[str, Any]) -> Optional[Dict[str, Any]]:
    coordination = payload.get("coordination_numbers") or {}
    peak_positions = payload.get("peak_positions") or {}
    if not coordination and not peak_positions:
        return None

    rows = []
    all_pairs = sorted(set(coordination.keys()) | set(peak_positions.keys()))
    for pair_name in all_pairs:
        first_peak_values = peak_positions.get(pair_name) or []
        first_peak = None
        if first_peak_values:
            first_peak = _safe_float(first_peak_values[0])
        rows.append(
            {
                "pair": pair_name,
                "coordination_number": _safe_float(coordination.get(pair_name)),
                "first_peak_A": first_peak,
            }
        )
    return {"task_id": task_id, "pairs": rows}


def _linear_fit(points: List[Dict[str, Any]], value_key: str, prefactor_key: str) -> Optional[Dict[str, float]]:
    valid_points = [
        point for point in points
        if point.get("temperature_K") and point.get(value_key) and point[value_key] > 0
    ]
    if len(valid_points) < 2:
        return None

    x_values = [1.0 / float(point["temperature_K"]) for point in valid_points]
    y_values = [math.log(float(point[value_key])) for point in valid_points]
    n = float(len(valid_points))
    x_mean = sum(x_values) / n
    y_mean = sum(y_values) / n
    numerator = sum((x - x_mean) * (y - y_mean) for x, y in zip(x_values, y_values))
    denominator = sum((x - x_mean) ** 2 for x in x_values)
    if denominator == 0:
        return None
    slope = numerator / denominator
    intercept = y_mean - slope * x_mean
    y_pred = [slope * x + intercept for x in x_values]
    ss_res = sum((y - pred) ** 2 for y, pred in zip(y_values, y_pred))
    ss_tot = sum((y - y_mean) ** 2 for y in y_values)
    r_squared = 1.0 if ss_tot == 0 else max(0.0, 1.0 - (ss_res / ss_tot))

    kb_ev_per_k = 8.617333262145e-5
    fit = {
        "slope": float(slope),
        "intercept": float(intercept),
        "r2": float(r_squared),
        prefactor_key: float(math.exp(intercept)),
    }
    if value_key == "diffusivity_m2_per_s":
        fit["Ea_eV"] = float(-slope * kb_ev_per_k)
    return fit


def _build_md_multi_summary(points: List[Dict[str, Any]]) -> Dict[str, Any]:
    ordered_points = sorted(points, key=lambda item: float(item["temperature_K"]))
    return {
        "task_count": len(ordered_points),
        "points": ordered_points,
        "fit": _linear_fit(ordered_points, "diffusivity_m2_per_s", "D0_m2_per_s"),
        "conductivity_fit": _linear_fit(ordered_points, "ionic_conductivity_S_per_m", "sigma0_S_per_m"),
    }


def _build_md_multi_stability_summary(stability_points: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not stability_points:
        return {"task_count": 0}

    density_values = [point["average_density_kg_m3"] for point in stability_points if point.get("average_density_kg_m3") is not None]
    energy_drifts = [point["energy_drift_eV"] for point in stability_points if point.get("energy_drift_eV") is not None]
    temperature_std_values = [point["temperature_std_K"] for point in stability_points if point.get("temperature_std_K") is not None]

    summary: Dict[str, Any] = {
        "task_count": len(stability_points),
        "average_temperature_std_K": (sum(temperature_std_values) / len(temperature_std_values)) if temperature_std_values else None,
        "average_energy_drift_eV": (sum(energy_drifts) / len(energy_drifts)) if energy_drifts else None,
        "rows": stability_points,
    }
    if density_values:
        summary["density_range_kg_m3"] = [min(density_values), max(density_values)]
    return summary


def _build_md_multi_rdf_summary(rdf_points: List[Dict[str, Any]]) -> Dict[str, Any]:
    if not rdf_points:
        return {"task_count": 0, "pair_names": [], "rows": []}

    pair_names = sorted({row["pair"] for point in rdf_points for row in point.get("pairs", [])})
    flat_rows = []
    for point in rdf_points:
        for row in point.get("pairs", []):
            flat_rows.append({"task_id": point["task_id"], **row})
    return {
        "task_count": len(rdf_points),
        "pair_names": pair_names,
        "rows": flat_rows,
    }


def _build_md_multi_public_summary(summary: Dict[str, Any], task_ids: List[str]) -> Dict[str, Any]:
    arrhenius = summary.get("arrhenius") or {}
    stability = summary.get("stability") or {}
    rdf = summary.get("rdf") or {}
    points = arrhenius.get("points") or []
    temperatures = [
        float(point["temperature_K"])
        for point in points
        if point.get("temperature_K") is not None
    ]
    mobile_species = sorted(
        {
            str(specie)
            for point in points
            for specie in (point.get("mobile_species") or [])
            if specie
        }
    )
    public_arrhenius = {
        "task_count": arrhenius.get("task_count", len(points)),
        "fit": arrhenius.get("fit"),
        "conductivity_fit": arrhenius.get("conductivity_fit"),
        "mobile_species": mobile_species or None,
    }
    if temperatures:
        public_arrhenius["temperature_range_K"] = [min(temperatures), max(temperatures)]

    public_stability = {
        "task_count": stability.get("task_count", 0),
        "average_temperature_std_K": stability.get("average_temperature_std_K"),
        "average_energy_drift_eV": stability.get("average_energy_drift_eV"),
    }
    if stability.get("density_range_kg_m3") is not None:
        public_stability["density_range_kg_m3"] = stability.get("density_range_kg_m3")

    public_rdf = {
        "task_count": rdf.get("task_count", 0),
        "pair_names": rdf.get("pair_names") or [],
    }

    return {
        "task_count": len(task_ids),
        "task_ids": task_ids,
        "arrhenius": public_arrhenius,
        "stability": public_stability,
        "rdf": public_rdf,
    }


def _svg_line_chart(
    title: str,
    x_label: str,
    y_label: str,
    series: List[Dict[str, Any]],
    *,
    width: int = 760,
    height: int = 320,
) -> str:
    valid_series = []
    for item in series:
        points = [(float(x), float(y)) for x, y in item.get("points", []) if x is not None and y is not None]
        if points:
            valid_series.append({**item, "points": points})

    if not valid_series:
        return '<div class="muted">暂无可绘制数据</div>'

    all_x = [x for item in valid_series for x, _ in item["points"]]
    all_y = [y for item in valid_series for _, y in item["points"]]
    x_min, x_max = min(all_x), max(all_x)
    y_min, y_max = min(all_y), max(all_y)
    if x_min == x_max:
        x_min -= 1.0
        x_max += 1.0
    if y_min == y_max:
        y_min -= 1.0
        y_max += 1.0

    left, right, top, bottom = 68, 18, 26, 56
    plot_w = width - left - right
    plot_h = height - top - bottom

    def sx(value: float) -> float:
        return left + (value - x_min) / (x_max - x_min) * plot_w

    def sy(value: float) -> float:
        return top + (y_max - value) / (y_max - y_min) * plot_h

    def fmt(value: float) -> str:
        magnitude = abs(value)
        if magnitude >= 1000 or (magnitude and magnitude < 0.01):
            return f"{value:.2e}"
        if magnitude >= 10:
            return f"{value:.2f}"
        return f"{value:.3f}"

    x_ticks = [x_min + (x_max - x_min) * i / 4 for i in range(5)]
    y_ticks = [y_min + (y_max - y_min) * i / 4 for i in range(5)]
    palette = ["#2563eb", "#ef4444", "#10b981", "#f59e0b", "#8b5cf6"]

    grid = []
    for value in x_ticks:
        x_pos = sx(value)
        grid.append(f'<line x1="{x_pos:.1f}" y1="{top}" x2="{x_pos:.1f}" y2="{top + plot_h}" stroke="#e2e8f0" stroke-width="1" />')
        grid.append(f'<text x="{x_pos:.1f}" y="{height - 20}" font-size="11" text-anchor="middle" fill="#475569">{escape(fmt(value))}</text>')
    for value in y_ticks:
        y_pos = sy(value)
        grid.append(f'<line x1="{left}" y1="{y_pos:.1f}" x2="{left + plot_w}" y2="{y_pos:.1f}" stroke="#e2e8f0" stroke-width="1" />')
        grid.append(f'<text x="{left - 8}" y="{y_pos + 4:.1f}" font-size="11" text-anchor="end" fill="#475569">{escape(fmt(value))}</text>')

    series_markup = []
    legends = []
    for index, item in enumerate(valid_series):
        color = item.get("color") or palette[index % len(palette)]
        points_attr = " ".join(f"{sx(x):.1f},{sy(y):.1f}" for x, y in item["points"])
        series_markup.append(
            f'<polyline fill="none" stroke="{color}" stroke-width="2.5" points="{points_attr}" />'
        )
        for x, y in item["points"]:
            series_markup.append(f'<circle cx="{sx(x):.1f}" cy="{sy(y):.1f}" r="3.2" fill="{color}" />')
        legends.append(
            f'<g transform="translate({left + index * 150},{height - 8})"><rect width="12" height="12" fill="{color}" rx="2" />'
            f'<text x="18" y="10" font-size="11" fill="#334155">{escape(str(item.get("name") or f"series-{index + 1}"))}</text></g>'
        )

    return f"""
    <div class="chart-card">
      <h3>{escape(title)}</h3>
      <svg viewBox="0 0 {width} {height}" width="100%" height="auto" role="img" aria-label="{escape(title)}">
        <rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff" rx="10" />
        {''.join(grid)}
        <line x1="{left}" y1="{top}" x2="{left}" y2="{top + plot_h}" stroke="#94a3b8" stroke-width="1.4" />
        <line x1="{left}" y1="{top + plot_h}" x2="{left + plot_w}" y2="{top + plot_h}" stroke="#94a3b8" stroke-width="1.4" />
        {''.join(series_markup)}
        <text x="{left + plot_w / 2:.1f}" y="{height - 34}" text-anchor="middle" font-size="12" fill="#334155">{escape(x_label)}</text>
        <text x="18" y="{top + plot_h / 2:.1f}" text-anchor="middle" font-size="12" fill="#334155" transform="rotate(-90 18 {top + plot_h / 2:.1f})">{escape(y_label)}</text>
        {''.join(legends)}
      </svg>
    </div>
    """


def _format_optional_float(value: Any, fmt: str) -> str:
    if value is None:
        return "--"
    return format(float(value), fmt)


def _write_md_multi_report(user_id: str, summary: Dict[str, Any]) -> str:
    report_dir = (
        Path(settings.local_public_artifact_root).resolve()
        / "md_multi_analysis"
        / str(user_id)
        / uuid.uuid4().hex
    )
    report_dir.mkdir(parents=True, exist_ok=True)

    points = summary.get("points") or []
    points_csv = report_dir / "md_multi_points.csv"
    with points_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "task_id",
                "temperature_K",
                "diffusivity_m2_per_s",
                "ionic_conductivity_S_per_m",
                "mobile_species",
            ],
        )
        writer.writeheader()
        for point in points:
            writer.writerow(
                {
                    "task_id": point["task_id"],
                    "temperature_K": point["temperature_K"],
                    "diffusivity_m2_per_s": point["diffusivity_m2_per_s"],
                    "ionic_conductivity_S_per_m": point.get("ionic_conductivity_S_per_m"),
                    "mobile_species": ",".join(point.get("mobile_species") or []),
                }
            )

    stability_summary = summary.get("stability") or {}
    rdf_summary = summary.get("rdf") or {}

    stability_rows = stability_summary.get("rows") or []
    stability_csv = report_dir / "md_multi_stability.csv"
    with stability_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "task_id",
                "average_temperature_K",
                "temperature_std_K",
                "average_pressure_kB",
                "average_density_kg_m3",
                "average_volume_A3",
                "energy_drift_eV",
            ],
        )
        writer.writeheader()
        for row in stability_rows:
            writer.writerow(row)

    rdf_rows = rdf_summary.get("rows") or []
    rdf_csv = report_dir / "md_multi_rdf_summary.csv"
    with rdf_csv.open("w", encoding="utf-8", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["task_id", "pair", "coordination_number", "first_peak_A"],
        )
        writer.writeheader()
        for row in rdf_rows:
            writer.writerow(row)

    fit = summary.get("fit") or {}
    conductivity_fit = summary.get("conductivity_fit") or {}
    row_html: List[str] = []
    for point in points:
        conductivity = point.get("ionic_conductivity_S_per_m")
        conductivity_text = "" if conductivity is None else f"{float(conductivity):.6e}"
        row_html.append(
            "<tr>"
            f"<td>{point['task_id']}</td>"
            f"<td>{float(point['temperature_K']):.1f}</td>"
            f"<td>{float(point['diffusivity_m2_per_s']):.6e}</td>"
            f"<td>{conductivity_text}</td>"
            f"<td>{', '.join(point.get('mobile_species') or [])}</td>"
            "</tr>"
        )
    rows = "\n".join(row_html)

    arrhenius_series = [
        {
            "name": "任务点",
            "color": "#2563eb",
            "points": [
                (1000.0 / float(point["temperature_K"]), math.log(float(point["diffusivity_m2_per_s"])))
                for point in points
            ],
        }
    ]
    if fit and len(points) >= 2:
        min_x = min(1000.0 / float(point["temperature_K"]) for point in points)
        max_x = max(1000.0 / float(point["temperature_K"]) for point in points)
        arrhenius_series.append(
            {
                "name": "线性拟合",
                "color": "#ef4444",
                "points": [
                    (min_x, float(fit["slope"]) * (min_x / 1000.0) + float(fit["intercept"])),
                    (max_x, float(fit["slope"]) * (max_x / 1000.0) + float(fit["intercept"])),
                ],
            }
        )

    diffusivity_series = [
        {
            "name": "平均扩散系数",
            "color": "#0f766e",
            "points": [
                (float(point["temperature_K"]), math.log10(float(point["diffusivity_m2_per_s"])))
                for point in points
            ],
        }
    ]
    conductivity_points = [
        (float(point["temperature_K"]), math.log10(float(point["ionic_conductivity_S_per_m"])))
        for point in points
        if point.get("ionic_conductivity_S_per_m") and float(point["ionic_conductivity_S_per_m"]) > 0
    ]
    conductivity_series = [{"name": "平均离子电导率", "color": "#7c3aed", "points": conductivity_points}] if conductivity_points else []

    density_points = [
        (float(row["average_temperature_K"]), float(row["average_density_kg_m3"]))
        for row in stability_rows
        if row.get("average_temperature_K") is not None and row.get("average_density_kg_m3") is not None
    ]
    energy_drift_points = [
        (float(row["average_temperature_K"]), float(row["energy_drift_eV"]))
        for row in stability_rows
        if row.get("average_temperature_K") is not None and row.get("energy_drift_eV") is not None
    ]

    stability_table_rows = "\n".join(
        "<tr>"
        f"<td>{escape(str(row['task_id']))}</td>"
        f"<td>{_format_optional_float(row.get('average_temperature_K'), '.2f')}</td>"
        f"<td>{_format_optional_float(row.get('temperature_std_K'), '.3f')}</td>"
        f"<td>{_format_optional_float(row.get('average_pressure_kB'), '.3f')}</td>"
        f"<td>{_format_optional_float(row.get('energy_drift_eV'), '.4f')}</td>"
        f"<td>{_format_optional_float(row.get('average_density_kg_m3'), '.2f')}</td>"
        f"<td>{_format_optional_float(row.get('average_volume_A3'), '.3f')}</td>"
        "</tr>"
        for row in stability_rows
    )

    rdf_table_rows = "\n".join(
        "<tr>"
        f"<td>{escape(str(row['task_id']))}</td>"
        f"<td>{escape(str(row['pair']))}</td>"
        f"<td>{_format_optional_float(row.get('coordination_number'), '.3f')}</td>"
        f"<td>{_format_optional_float(row.get('first_peak_A'), '.3f')}</td>"
        "</tr>"
        for row in rdf_rows
    )
    html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="utf-8" />
  <title>MD 多任务聚合分析</title>
  <style>
    body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', sans-serif; margin: 32px; color: #0f172a; }}
    h1, h2 {{ margin-bottom: 12px; }}
    .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 16px; margin: 20px 0; }}
    .card {{ border: 1px solid #dbeafe; border-left: 4px solid #2563eb; border-radius: 12px; padding: 16px; background: #f8fbff; }}
    .chart-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(340px, 1fr)); gap: 20px; margin: 20px 0; }}
    .chart-card {{ border: 1px solid #e2e8f0; border-radius: 12px; padding: 12px; background: #fff; }}
    table {{ width: 100%; border-collapse: collapse; margin-top: 16px; }}
    th, td {{ border-bottom: 1px solid #e2e8f0; padding: 10px; text-align: left; font-size: 14px; }}
    th {{ background: #f8fafc; }}
    .muted {{ color: #475569; }}
    img {{ max-width: 760px; width: 100%; border: 1px solid #e2e8f0; border-radius: 12px; margin-top: 16px; }}
  </style>
</head>
<body>
  <h1>MD 多任务聚合分析</h1>
  <p class="muted">基于已上传的 diffusion_results.json 聚合，不依赖 HPC 工作目录。</p>
  <div class="cards">
    <div class="card"><strong>任务数</strong><div>{len(points)}</div></div>
    <div class="card"><strong>扩散激活能</strong><div>{'--' if not fit else f"{fit.get('Ea_eV', float('nan')):.4f} eV"}</div></div>
    <div class="card"><strong>扩散拟合 R²</strong><div>{'--' if not fit else f"{fit.get('r2', float('nan')):.4f}"}</div></div>
    <div class="card"><strong>电导拟合 R²</strong><div>{'--' if not conductivity_fit else f"{conductivity_fit.get('r2', float('nan')):.4f}"}</div></div>
  </div>
  <p><a href="{points_csv.name}">下载聚合数据 CSV</a> | <a href="{stability_csv.name}">下载稳定性摘要 CSV</a> | <a href="{rdf_csv.name}">下载 RDF 摘要 CSV</a></p>
  <h2>扩散与电导图</h2>
  <div class="chart-grid">
    {_svg_line_chart("Arrhenius 图", "1000 / T (K^-1)", "ln D", arrhenius_series)}
    {_svg_line_chart("扩散系数-温度", "温度 (K)", "log10 D (m²/s)", diffusivity_series)}
    {_svg_line_chart("离子电导率-温度", "温度 (K)", "log10 σ (S/m)", conductivity_series)}
  </div>
  <h2>任务点</h2>
  <table>
    <thead><tr><th>任务ID</th><th>温度 (K)</th><th>扩散系数 (m²/s)</th><th>离子电导率 (S/m)</th><th>迁移离子</th></tr></thead>
    <tbody>{rows}</tbody>
  </table>
  <h2>稳定性摘要</h2>
  <div class="chart-grid">
    {_svg_line_chart("平均密度-温度", "温度 (K)", "密度 (kg/m³)", [{"name": "平均密度", "color": "#0891b2", "points": density_points}])}
    {_svg_line_chart("能量漂移-温度", "温度 (K)", "能量漂移 (eV)", [{"name": "能量漂移", "color": "#dc2626", "points": energy_drift_points}])}
  </div>
  <table>
    <thead><tr><th>任务ID</th><th>平均温度 (K)</th><th>温度波动 (K)</th><th>平均压力 (kB)</th><th>能量漂移 (eV)</th><th>平均密度 (kg/m³)</th><th>平均体积 (Å³)</th></tr></thead>
    <tbody>{stability_table_rows}</tbody>
  </table>
  <h2>RDF 摘要</h2>
  <p class="muted">多任务聚合只展示每个任务已生成的 RDF 摘要指标，不在此处重复展开整条 RDF 曲线；完整曲线仍在单任务 HTML 中查看。</p>
  <table>
    <thead><tr><th>任务ID</th><th>元素对</th><th>配位数</th><th>首峰位置 (Å)</th></tr></thead>
    <tbody>{rdf_table_rows}</tbody>
  </table>
</body>
</html>"""
    html_path = report_dir / "md_multi_analysis_report.html"
    html_path.write_text(html, encoding="utf-8")
    return str(html_path)


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
    if request.task_id:
        analysis_task = _submit_remote_analysis_task(
            analysis_type="dos",
            user_id=request.user_id or os.getenv("DEFAULT_USER_ID", "123"),
            source_task_id=request.task_id,
        )
        return _build_remote_analysis_response(analysis_task, "dos")

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
    if request.task_id:
        analysis_task = _submit_remote_analysis_task(
            analysis_type="band_structure",
            user_id=request.user_id or os.getenv("DEFAULT_USER_ID", "123"),
            source_task_id=request.task_id,
        )
        return _build_remote_analysis_response(analysis_task, "band_structure")

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
    if request.task_id:
        analysis_task = _submit_remote_analysis_task(
            analysis_type="md",
            user_id=request.user_id or os.getenv("DEFAULT_USER_ID", "123"),
            source_task_id=request.task_id,
        )
        return _build_remote_analysis_response(analysis_task, "md")

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
    try:
        points: List[Dict[str, Any]] = []
        stability_points: List[Dict[str, Any]] = []
        rdf_points: List[Dict[str, Any]] = []
        for tid in request.task_ids:
            task = task_manager.get_task(tid, request.user_id)
            if not task:
                raise HTTPException(status_code=404, detail=f"任务 {tid} 未找到或无权限访问")
            if str(task.status) != "completed":  # type: ignore
                raise HTTPException(status_code=400, detail=f"任务 {tid} 尚未完成 (status={task.status})")  # type: ignore

            diffusion_payload = _load_md_diffusion_results(tid)
            points.append(_extract_md_diffusion_point(tid, diffusion_payload))
            stability_point = _extract_md_stability_point(tid, _load_md_time_series_rows(tid))
            if stability_point:
                stability_points.append(stability_point)
            rdf_point = _extract_md_rdf_point(tid, _load_md_rdf_summary(tid))
            if rdf_point:
                rdf_points.append(rdf_point)

        if not points:
            raise FileNotFoundError("没有可用于聚合分析的 diffusion_results.json")

        arrhenius = _build_md_multi_summary(points)
        full_summary = {
            "task_count": len(points),
            "task_ids": request.task_ids,
            "points": arrhenius.get("points") or [],
            "fit": arrhenius.get("fit"),
            "conductivity_fit": arrhenius.get("conductivity_fit"),
            "arrhenius": arrhenius,
            "stability": _build_md_multi_stability_summary(stability_points),
            "rdf": _build_md_multi_rdf_summary(rdf_points),
        }
        html_path = _write_md_multi_report(request.user_id, full_summary)
        public_summary = _build_md_multi_public_summary(full_summary, request.task_ids)
        return AnalyzeResponse(
            success=True,
            analysis_type="md_multi",
            summary=public_summary,
            html_report_url=_public_url_for_local_path(html_path),
        )
    except HTTPException as exc:
        return AnalyzeResponse(
            success=False,
            analysis_type="md_multi",
            error_message=_analysis_error_message(exc),
        )
    except Exception as e:
        logger.error(f"多任务MD聚合分析失败: {e}", exc_info=True)
        return AnalyzeResponse(success=False, analysis_type="md_multi", error_message=str(e))


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
        analysis_task = _submit_remote_analysis_task(
            analysis_type="agent",
            user_id=request.user_id,
            source_task_id=request.task_id,
            extra_params={
                "question": request.question,
                "model": request.model,
            },
        )
        return _build_remote_agent_response(analysis_task)

    except HTTPException:
        raise
    except Exception as e:
        logger.error("Agent 分析失败: %s", e, exc_info=True)
        return AgentAnalyzeResponse(success=False, error_message=str(e))


if __name__ == "__main__":
    from .Config import VASP_remote_run_port
    uvicorn.run(app, host="0.0.0.0", port=VASP_remote_run_port)
