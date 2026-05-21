import hashlib
import json
import mimetypes
import os
import threading
import uuid
import logging
from datetime import datetime, timedelta, timezone
from types import SimpleNamespace
from typing import Any, Dict, List, Optional

from sqlalchemy.orm import Session
from sqlalchemy import inspect as sa_inspect

from sqlalchemy import and_, func as sa_func, or_

from .database import SessionLocal
from ..errors import APIError
from .models import Task, ExecutionAttempt, AnalysisRun, Artifact
from ..notification_client import (
    TASK_TYPE_TO_DISPLAY_NAME,
    build_task_notification_payload,
    send_notification_async,
)
from ..storage import ObjectStorageService

logger = logging.getLogger(__name__)

# 可通过环境变量调整
MAX_CONCURRENT_TASKS = int(os.getenv("MAX_CONCURRENT_TASKS", "4"))
MAX_TASK_RETRIES = int(os.getenv("MAX_TASK_RETRIES", "1"))

# 可重试的错误关键词（HPC/基础设施故障）
_RETRYABLE_KEYWORDS = [
    "slurmstepd", "node failure", "oom", "out of memory",
    "connection", "timeout", "resource temporarily unavailable",
    "job killed", "segfault", "bus error",
]


def _orm_to_ns(obj):
    """Convert a SQLAlchemy ORM object to a SimpleNamespace with all column values.

    This must be called while the session is still open so that lazy-loaded
    attributes are accessible.  The returned SimpleNamespace supports
    attribute access (obj.status) and getattr() just like the original ORM
    object, but is fully detached from the session.
    """
    if obj is None:
        return None
    mapper = sa_inspect(type(obj))
    data = {col.key: getattr(obj, col.key) for col in mapper.column_attrs}
    return SimpleNamespace(**data)

# 支持的任务类型 → 分析类型映射
TASK_TYPE_TO_ANALYSIS = {
    "structure_optimization": "optimization",
    "scf_calculation": "scf",
    "dos_calculation": "dos",
    "dos_analysis": "dos",
    "md_calculation": "md",
    "md_analysis": "md",
    "band_structure": "band_structure",
    "band_structure_analysis": "band_structure",
    "neb_calculation": "neb",
    "phonon_calculation": "phonon",
    "custom_calculation": "custom",
    "agent_analysis": "agent",
}

# 支持的任务类型 → VaspWorker 方法名映射
TASK_TYPE_TO_METHOD = {
    "structure_optimization": "run_structure_optimization",
    "scf_calculation": "run_scf_calculation",
    "dos_calculation": "run_dos_calculation",
    "dos_analysis": "run_dos_analysis",
    "md_calculation": "run_md_calculation",
    "md_analysis": "run_md_analysis",
    "band_structure": "run_band_structure_calculation",
    "band_structure_analysis": "run_band_structure_analysis",
    "neb_calculation": "run_neb_calculation",
    "phonon_calculation": "run_phonon_calculation",
    "custom_calculation": "run_custom_calculation",
    "agent_analysis": "run_agent_analysis",
}

TASK_TYPE_TO_NOTIFICATION_METRICS = {
    "structure_optimization": (
        ("final_energy", "最终能量"),
        ("final_max_force", "最终最大受力"),
    ),
    "scf_calculation": (
        ("total_energy", "总能量"),
        ("band_gap", "带隙"),
    ),
    "dos_calculation": (
        ("band_gap", "带隙"),
        ("total_energy", "总能量"),
    ),
    "dos_analysis": (
        ("band_gap", "带隙"),
        ("total_energy", "总能量"),
    ),
    "band_structure": (
        ("band_gap", "带隙"),
        ("total_energy", "总能量"),
    ),
    "band_structure_analysis": (
        ("band_gap", "带隙"),
        ("total_energy", "总能量"),
    ),
    "md_calculation": (
        ("final_energy", "最终能量"),
        ("average_temperature", "平均温度"),
    ),
    "md_analysis": (
        ("final_energy", "最终能量"),
        ("average_temperature", "平均温度"),
    ),
    "neb_calculation": (
        ("forward_barrier_eV", "正向势垒"),
        ("reaction_energy_eV", "反应能"),
    ),
    "phonon_calculation": (
        ("n_imaginary", "虚频数量"),
        ("max_real_freq_cm1", "最大实频"),
    ),
    "custom_calculation": (
        ("total_energy", "总能量"),
        ("band_gap", "带隙"),
    ),
    "agent_analysis": (
        ("steps", "分析步数"),
    ),
}


def _is_retryable_error(error_msg: str) -> bool:
    """判断错误是否为可重试的基础设施故障"""
    lower = error_msg.lower()
    return any(kw in lower for kw in _RETRYABLE_KEYWORDS)


def _compute_input_hash(task_type: str, params: Dict[str, Any]) -> str:
    """基于任务类型和关键输入参数计算 SHA-256 hash，用于去重。"""
    # 只取影响计算结果的关键参数，忽略 user_id 等
    key_fields = {
        "task_type": task_type,
        "queue_name": params.get("queue_name"),
        "cif_url": params.get("cif_url"),
        "input_url": params.get("input_url"),
        "source_task_id": params.get("source_task_id"),
        "source_work_directory": params.get("source_work_directory"),
        "question": params.get("question"),
        "model": params.get("model"),
        "optimized_task_id": params.get("optimized_task_id"),
        "scf_task_id": params.get("scf_task_id"),
        "custom_incar": params.get("custom_incar"),
        "kpoint_density": params.get("kpoint_density"),
        "precision": params.get("precision"),
        "md_steps": params.get("md_steps"),
        "temperature": params.get("temperature"),
        "time_step": params.get("time_step"),
        "ensemble": params.get("ensemble"),
        "line_density": params.get("line_density"),
        # NEB
        "initial_task_id": params.get("initial_task_id"),
        "initial_cif_url": params.get("initial_cif_url"),
        "final_task_id": params.get("final_task_id"),
        "final_cif_url": params.get("final_cif_url"),
        "n_images": params.get("n_images"),
        # Phonon
        "displacement": params.get("displacement"),
        # Custom
        "from_task_id": params.get("from_task_id"),
        "incar": params.get("incar"),
        "kpoint_mode": params.get("kpoint_mode"),
    }
    # 去掉 None 值，保证 hash 稳定
    filtered = {k: v for k, v in key_fields.items() if v is not None}
    canonical = json.dumps(filtered, sort_keys=True, ensure_ascii=False)
    return hashlib.sha256(canonical.encode()).hexdigest()


def _format_notification_value(value: Any) -> str:
    if isinstance(value, float):
        return f"{value:.4f}"
    return str(value)


class TaskManager:
    def __init__(self, max_workers: int = MAX_CONCURRENT_TASKS) -> None:
        self._max_workers = max_workers
        self._cancel_flags: Dict[str, threading.Event] = {}
        self._lock = threading.Lock()
        self.storage_service = ObjectStorageService.from_settings()

    @staticmethod
    def _has_active_cancel_handshake(task: Task) -> bool:
        return any(
            [
                getattr(task, "worker_id", None),
                getattr(task, "lease_token", None),
                getattr(task, "heartbeat_at", None),
                getattr(task, "external_job_id", None),
            ]
        )

    @staticmethod
    def _active_task_filter():
        return or_(
            Task.status.in_([
                "queued",
                "leased",
                "running",
                "uploading",
                "analyzing",
            ]),
            and_(
                Task.status == "cancel_requested",
                or_(
                    Task.worker_id.is_not(None),
                    Task.lease_token.is_not(None),
                    Task.heartbeat_at.is_not(None),
                    Task.external_job_id.is_not(None),
                ),
            ),
        )

    @classmethod
    def _normalize_unclaimed_cancel_requests(cls, db: Session, user_id: Optional[str] = None) -> int:
        query = db.query(Task).filter(Task.status == "cancel_requested")
        if user_id is not None:
            query = query.filter(Task.user_id == user_id)

        stuck_tasks = list(
            query.filter(
                Task.worker_id.is_(None),
                Task.lease_token.is_(None),
                Task.heartbeat_at.is_(None),
                Task.external_job_id.is_(None),
            ).all()
        )
        if not stuck_tasks:
            return 0

        now = datetime.now(timezone.utc)
        for task in stuck_tasks:
            task.status = "canceled"  # type: ignore
            task.finalized_at = task.finalized_at or now  # type: ignore
            db.add(task)
        db.flush()
        return len(stuck_tasks)

    @staticmethod
    def _notification_language(task: Task) -> str:
        params = task.params or {}
        language = str(params.get("notification_language") or params.get("language") or "zh")
        return language if language in {"zh", "en"} else "zh"

    @staticmethod
    def _notification_email(task: Task) -> Optional[str]:
        params = task.params or {}
        email = params.get("to_email") or params.get("notification_email")
        return str(email).strip() if email else None

    @staticmethod
    def _notification_execution_time(*sources: Optional[Dict[str, Any]]) -> Optional[float]:
        for source in sources:
            if not source:
                continue
            value = source.get("computation_time")
            if value is None:
                value = source.get("execution_time_seconds")
            if value is None:
                continue
            try:
                return float(value)
            except (TypeError, ValueError):
                continue
        return None

    @staticmethod
    def _notification_result_summary(task_type: str, result: Optional[Dict[str, Any]]) -> str:
        display_name = TASK_TYPE_TO_DISPLAY_NAME.get(task_type, task_type)
        data = result or {}
        for key, label in TASK_TYPE_TO_NOTIFICATION_METRICS.get(task_type, ()):
            value = data.get(key)
            if value is not None:
                return f"{display_name}完成，{label} = {_format_notification_value(value)}"
        return f"{display_name}完成"

    @staticmethod
    def _notification_error_message(error_message: str) -> str:
        return " ".join(str(error_message).strip().splitlines())[:300]

    @staticmethod
    def _is_public_notification_html_url(value: Any) -> bool:
        if not isinstance(value, str):
            return False
        return value.startswith(("https://", "http://"))

    def _notification_html_report_url(
        self,
        task: Task,
        *sources: Optional[Dict[str, Any]],
    ) -> Optional[str]:
        db: Session = SessionLocal()
        try:
            public_artifact_url = self._extract_html_report_url_from_artifacts(db, str(task.id))
        finally:
            db.close()
        if public_artifact_url:
            return public_artifact_url

        for source in (*sources, task.result_summary, task.result_data):
            if not source or not isinstance(source, dict):
                continue
            html_report_url = self._extract_html_report_url(source)
            if html_report_url and self._is_public_notification_html_url(html_report_url):
                return html_report_url
        return None

    def _notify_task_success(
        self,
        task: Task,
        result_data: Optional[Dict[str, Any]] = None,
        result_summary: Optional[Dict[str, Any]] = None,
    ) -> None:
        html_report_url = self._notification_html_report_url(task, result_summary, result_data)
        summary_text = self._notification_result_summary(
            str(task.task_type),
            result_summary or result_data,
        )
        if html_report_url:
            summary_text = f"{summary_text}，报告链接: {html_report_url}"

        payload = build_task_notification_payload(
            user_id=str(task.user_id),
            task_id=str(task.id),
            task_type=str(task.task_type),
            notification_type="tool_complete",
            status="success",
            language=self._notification_language(task),
            execution_time_seconds=self._notification_execution_time(result_summary, result_data),
            result_summary=summary_text,
            html_report_url=html_report_url,
            to_email=self._notification_email(task),
        )
        send_notification_async(payload)

    def _notify_task_error(
        self,
        task: Task,
        error_message: str,
        result_data: Optional[Dict[str, Any]] = None,
    ) -> None:
        html_report_url = self._notification_html_report_url(task, result_data)
        payload = build_task_notification_payload(
            user_id=str(task.user_id),
            task_id=str(task.id),
            task_type=str(task.task_type),
            notification_type="tool_error",
            status="error",
            language=self._notification_language(task),
            execution_time_seconds=self._notification_execution_time(result_data),
            error_message=self._notification_error_message(error_message),
            html_report_url=html_report_url,
            to_email=self._notification_email(task),
        )
        send_notification_async(payload)

    # ------------------------------------------------------------------ #
    #  提交任务
    # ------------------------------------------------------------------ #
    def submit_task(self, user_id: str, task_type: str, params: Optional[Dict[str, Any]]) -> str:
        params = params or {}
        input_hash = _compute_input_hash(task_type, params)

        # --- 去重：查找同用户下已完成且 input_hash 相同的任务 ---
        db: Session = SessionLocal()
        try:
            existing = (
                db.query(Task)
                .filter(
                    Task.user_id == user_id,
                    Task.input_hash == input_hash,
                    Task.status == "completed",
                )
                .order_by(Task.created_at.desc())
                .first()
            )
            if existing is not None:
                logger.info(
                    "任务去重命中: 用户 %s 已有相同计算 %s (hash=%s)",
                    user_id, existing.id, input_hash[:12],
                )
                return str(existing.id)

            # --- 配额检查 ---
            from ..settings import settings as _settings

            normalized = self._normalize_unclaimed_cancel_requests(db, user_id=user_id)
            if normalized:
                logger.info("提交前自动清理了 %s 个未领取的 cancel_requested 任务", normalized)

            # 1) 用户并发任务数限制
            active_count = (
                db.query(sa_func.count(Task.id))
                .filter(
                    Task.user_id == user_id,
                    self._active_task_filter(),
                )
                .scalar()
            ) or 0
            if active_count >= _settings.max_user_concurrent_tasks:
                active_task_ids = [
                    str(row[0])
                    for row in (
                        db.query(Task.id)
                        .filter(
                            Task.user_id == user_id,
                            self._active_task_filter(),
                        )
                        .order_by(Task.created_at.desc())
                        .all()
                    )
                ]
                raise APIError(
                    status_code=429,
                    code="USER_CONCURRENT_TASK_LIMIT",
                    message="已达到并发任务上限",
                    retryable=False,
                    details={
                        "limit": _settings.max_user_concurrent_tasks,
                        "active_count": active_count,
                        "active_task_ids": active_task_ids,
                    },
                    suggested_action="请等待现有任务完成，或先清理取消中的任务后再提交",
                )

            # 2) 每日任务总数限制
            today_start = datetime.now(timezone.utc).replace(
                hour=0, minute=0, second=0, microsecond=0
            )
            daily_count = (
                db.query(sa_func.count(Task.id))
                .filter(
                    Task.user_id == user_id,
                    Task.created_at >= today_start,
                )
                .scalar()
            ) or 0
            if daily_count >= _settings.max_user_total_tasks_per_day:
                raise APIError(
                    status_code=429,
                    code="USER_DAILY_TASK_LIMIT",
                    message="已达到每日任务上限",
                    retryable=False,
                    details={
                        "limit": _settings.max_user_total_tasks_per_day,
                        "daily_count": daily_count,
                    },
                    suggested_action="请明天再试，或减少重复提交",
                )

            task_id = uuid.uuid4().hex

            task = Task(
                id=task_id,
                user_id=user_id,
                task_type=task_type,
                tenant_id=str(params.get("tenant_id", "default")),
                queue_name=str(params.get("queue_name", "default")),
                priority=int(params.get("priority", 0)),
                status="queued",
                analysis_status="pending",
                progress=0,
                params=params,
                input_snapshot=params,
                input_hash=input_hash,
                retry_count=0,
                max_retries=int(params.get("max_retries", MAX_TASK_RETRIES)),
            )
            db.add(task)
            db.commit()
        finally:
            db.close()

        return task_id

    def _lease_expiration(self, now: datetime) -> datetime:
        from ..settings import settings as _settings

        return now + timedelta(seconds=_settings.task_lease_seconds)

    def register_artifact(
        self,
        task_id: str,
        artifact_type: str,
        owner_type: str,
        owner_id: Optional[str],
        filename: str,
        content_type: Optional[str] = None,
        size_bytes: Optional[float] = None,
        etag: Optional[str] = None,
        sha256: Optional[str] = None,
        attempt_no: int = 1,
    ) -> SimpleNamespace:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None:
                raise ValueError(f"任务不存在: {task_id}")

            location = self.storage_service.build_location(
                tenant_id=str(task.tenant_id),
                task_id=task_id,
                attempt_no=attempt_no,
                filename=filename,
            )
            artifact = Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type=owner_type,
                owner_id=owner_id,
                artifact_type=artifact_type,
                storage_backend=location.storage_backend,
                storage_key=location.object_key,
                bucket=location.bucket,
                object_key=location.object_key,
                mime_type=content_type,
                content_type=content_type,
                size_bytes=size_bytes,
                checksum=sha256,
                etag=etag,
                sha256=sha256,
            )
            db.add(artifact)
            db.commit()
            db.refresh(artifact)
            return _orm_to_ns(artifact)
        finally:
            db.close()

    def get_artifact_download_url(self, artifact: SimpleNamespace) -> Optional[str]:
        object_key = getattr(artifact, "object_key", None) or getattr(artifact, "storage_key", None)
        storage_backend = getattr(artifact, "storage_backend", None)
        if storage_backend in {"oss", "s3"} and object_key:
            return self.storage_service.create_download_url(str(object_key))
        if storage_backend == "local_public" and object_key:
            return self.storage_service.create_public_download_url(str(object_key))
        return None

    @staticmethod
    def _infer_artifact_type(filename: str) -> str:
        lower = filename.lower()
        if lower.endswith(".html"):
            return "html_report"
        if lower.endswith((".png", ".jpg", ".jpeg", ".svg", ".webp")):
            return "image"
        if lower.endswith(".csv"):
            return "csv"
        if lower.endswith(".json"):
            return "json"
        mapping = {
            "outcar": "outcar",
            "contcar": "contcar",
            "poscar": "poscar",
            "potcar": "potcar",
            "doscar": "doscar",
            "chgcar": "chgcar",
            "chg": "chg",
            "wavecar": "wavecar",
            "kpoints": "kpoints",
            "incar": "incar",
            "result.log": "result_log",
        }
        return mapping.get(lower, mapping.get(os.path.basename(lower), "file"))

    def _persist_artifact_manifest(
        self,
        db: Session,
        task: Task,
        artifact_manifest: List[Dict[str, Any]] | None,
    ) -> None:
        if not artifact_manifest:
            return

        existing_filenames = {
            str(a.storage_key)
            for a in db.query(Artifact).filter(Artifact.task_id == task.id).all()
        }

        for item in artifact_manifest:
            filename = str(item.get("filename") or os.path.basename(str(item.get("local_path") or "")) or "artifact.bin")
            if not filename:
                continue
            item_storage_key = str(item.get("storage_key") or item.get("local_path") or "")
            if any(existing.endswith(f"/{filename}") or existing == filename or (item_storage_key and existing == item_storage_key) for existing in existing_filenames):
                continue

            content_type = item.get("content_type")
            if content_type is None:
                content_type = mimetypes.guess_type(filename)[0]

            storage_backend = str(item.get("storage_backend") or "").strip().lower()
            storage_key = None
            object_key = None
            bucket = None
            if storage_backend == "local_public":
                object_key = str(item.get("object_key") or filename)
                storage_key = str(item.get("storage_key") or object_key)
            elif storage_backend == "local" or item.get("local_path"):
                storage_backend = "local"
                storage_key = str(item.get("storage_key") or item.get("local_path") or filename)
                object_key = item.get("object_key")
            else:
                location = self.storage_service.build_location(
                    tenant_id=str(task.tenant_id),
                    task_id=str(task.id),
                    attempt_no=int(item.get("attempt_no", 1)),
                    filename=filename,
                )
                storage_backend = location.storage_backend
                storage_key = location.object_key
                bucket = location.bucket
                object_key = location.object_key

            artifact = Artifact(
                id=uuid.uuid4().hex,
                task_id=task.id,
                owner_type=str(item.get("owner_type", "execution")),
                owner_id=item.get("owner_id"),
                artifact_type=str(item.get("artifact_type") or self._infer_artifact_type(filename)),
                storage_backend=storage_backend,
                storage_key=str(storage_key),
                bucket=bucket,
                object_key=str(object_key) if object_key is not None else None,
                mime_type=content_type,
                content_type=content_type,
                size_bytes=item.get("size_bytes"),
                checksum=item.get("sha256"),
                etag=item.get("etag"),
                sha256=item.get("sha256"),
            )
            db.add(artifact)
            existing_filenames.add(str(storage_key))

    @staticmethod
    def _extract_attempt_scheduler_job_id(payload: Optional[Dict[str, Any]]) -> Optional[str]:
        if not payload:
            return None
        for key in ("scheduler_job_id", "slurm_job_id", "process_id", "external_job_id"):
            value = payload.get(key)
            if value is not None and value != "":
                return str(value)
        return None

    @staticmethod
    def _extract_attempt_work_directory(payload: Optional[Dict[str, Any]]) -> Optional[str]:
        if not payload:
            return None
        work_directory = payload.get("work_directory")
        if work_directory:
            return str(work_directory)
        return None

    def _upsert_distributed_execution_attempt(
        self,
        db: Session,
        task: Task,
        *,
        worker_id: str,
        lease_token: str,
        status: str,
        artifact_manifest: Optional[List[Dict[str, Any]]] = None,
        result_data: Optional[Dict[str, Any]] = None,
        failure_context: Optional[Dict[str, Any]] = None,
        failure_code: Optional[str] = None,
        failure_detail: Optional[str] = None,
    ) -> ExecutionAttempt:
        now = datetime.now(timezone.utc)
        current_attempt_id = getattr(task, "current_execution_attempt_id", None)
        attempt = db.get(ExecutionAttempt, current_attempt_id) if current_attempt_id else None

        if attempt is not None and str(getattr(attempt, "status", "")) in {"succeeded", "runtime_failed", "canceled", "submit_failed"}:
            attempt = None

        if attempt is None:
            attempt = ExecutionAttempt(
                id=uuid.uuid4().hex,
                task_id=str(task.id),
                attempt_no=int(getattr(task, "retry_count", 0) or 0) + 1,
                executor_type="slurm",
                status="submitted",
                worker_id=worker_id,
                lease_token=lease_token,
                submitted_at=now,
                started_at=now,
            )
            db.add(attempt)
            db.flush()
            task.current_execution_attempt_id = attempt.id  # type: ignore
        else:
            if not attempt.submitted_at:
                attempt.submitted_at = now  # type: ignore
            if not attempt.started_at:
                attempt.started_at = now  # type: ignore

        attempt.worker_id = worker_id  # type: ignore
        attempt.lease_token = lease_token  # type: ignore
        attempt.status = status  # type: ignore

        scheduler_job_id = (
            self._extract_attempt_scheduler_job_id(result_data)
            or self._extract_attempt_scheduler_job_id(failure_context)
        )
        if scheduler_job_id:
            attempt.scheduler_job_id = scheduler_job_id  # type: ignore
            attempt.external_job_id = scheduler_job_id  # type: ignore

        work_directory = (
            self._extract_attempt_work_directory(result_data)
            or self._extract_attempt_work_directory(failure_context)
        )
        if work_directory:
            attempt.work_directory = work_directory  # type: ignore

        if artifact_manifest is not None:
            attempt.artifact_manifest = artifact_manifest  # type: ignore
        if failure_code is not None:
            attempt.failure_code = failure_code  # type: ignore
        if failure_detail is not None:
            attempt.failure_detail = failure_detail  # type: ignore
        if status in {"succeeded", "runtime_failed", "canceled", "submit_failed"}:
            attempt.finished_at = now  # type: ignore

        db.add(attempt)
        db.add(task)
        return attempt

    def _validate_lease(
        self,
        db: Session,
        task_id: str,
        lease_token: str,
        worker_id: str,
    ) -> Task:
        task: Optional[Task] = db.get(Task, task_id)
        if task is None:
            raise ValueError(f"任务不存在: {task_id}")
        if str(task.worker_id) != worker_id or str(task.lease_token) != lease_token:
            raise ValueError("lease token 或 worker_id 无效")
        return task

    def claim_next_task(self, worker_id: str, queue_name: str) -> Optional[SimpleNamespace]:
        self.mark_orphaned_running_tasks()
        self.requeue_expired_leases(queue_name=queue_name)

        db: Session = SessionLocal()
        try:
            now = datetime.now(timezone.utc)
            query = (
                db.query(Task)
                .filter(
                    Task.status == "queued",
                    Task.queue_name == queue_name,
                )
                .order_by(Task.priority.desc(), Task.created_at.asc())
            )
            if db.bind is not None and db.bind.dialect.name == "postgresql":
                query = query.with_for_update(skip_locked=True)

            task = query.first()
            if task is None:
                return None

            task.status = "leased"  # type: ignore
            task.worker_id = worker_id  # type: ignore
            task.lease_token = uuid.uuid4().hex  # type: ignore
            task.lease_expires_at = self._lease_expiration(now)  # type: ignore
            task.heartbeat_at = now  # type: ignore
            task.error_message = None  # type: ignore
            db.add(task)
            db.commit()
            db.refresh(task)
            return _orm_to_ns(task)
        finally:
            db.close()

    def resume_task(self, task_id: str, worker_id: str, queue_name: str) -> SimpleNamespace:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None:
                raise ValueError(f"任务不存在: {task_id}")
            if str(task.queue_name) != queue_name:
                raise ValueError(f"任务不属于队列 {queue_name}")
            if str(task.status) in {"completed", "failed", "canceled"}:
                raise ValueError(f"任务已结束，无法恢复: {task.status}")

            now = datetime.now(timezone.utc)
            if str(task.status) in {"queued", "leased"}:
                task.status = "leased"  # type: ignore

            task.worker_id = worker_id  # type: ignore
            task.lease_token = uuid.uuid4().hex  # type: ignore
            task.lease_expires_at = self._lease_expiration(now)  # type: ignore
            task.heartbeat_at = now  # type: ignore
            task.error_message = None  # type: ignore
            db.add(task)
            db.commit()
            db.refresh(task)
            return _orm_to_ns(task)
        finally:
            db.close()

    def heartbeat_task(self, task_id: str, lease_token: str, worker_id: str) -> None:
        db: Session = SessionLocal()
        try:
            task = self._validate_lease(db, task_id, lease_token, worker_id)
            if str(task.status) not in {"leased", "running", "uploading", "analyzing", "cancel_requested"}:
                raise ValueError(f"任务状态不允许续租: {task.status}")

            now = datetime.now(timezone.utc)
            task.heartbeat_at = now  # type: ignore
            task.lease_expires_at = self._lease_expiration(now)  # type: ignore
            db.add(task)
            db.commit()
        finally:
            db.close()

    def mark_task_running(self, task_id: str, lease_token: str, worker_id: str) -> None:
        db: Session = SessionLocal()
        try:
            task = self._validate_lease(db, task_id, lease_token, worker_id)
            task.status = "running"  # type: ignore
            task.error_message = None  # type: ignore
            db.add(task)
            db.commit()
        finally:
            db.close()

    def request_cancel(self, task_id: str, user_id: str) -> bool:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None or str(task.user_id) != user_id:
                return False
            if str(task.status) in {"completed", "failed", "canceled"}:
                return True

            now = datetime.now(timezone.utc)
            task.cancel_requested_at = now  # type: ignore
            if str(task.status) == "queued" and not self._has_active_cancel_handshake(task):
                task.status = "canceled"  # type: ignore
                task.finalized_at = now  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                task.heartbeat_at = None  # type: ignore
            else:
                task.status = "cancel_requested"  # type: ignore
            db.add(task)
            db.commit()
            return True
        finally:
            db.close()

    def ack_cancel(self, task_id: str, lease_token: str, worker_id: str) -> None:
        db: Session = SessionLocal()
        try:
            task = self._validate_lease(db, task_id, lease_token, worker_id)
            task.status = "canceled"  # type: ignore
            task.finalized_at = datetime.now(timezone.utc)  # type: ignore
            task.worker_id = None  # type: ignore
            task.lease_token = None  # type: ignore
            task.lease_expires_at = None  # type: ignore
            db.add(task)
            db.commit()
        finally:
            db.close()

    def complete_execution(
        self,
        task_id: str,
        lease_token: str,
        worker_id: str,
        result_data: Optional[Dict[str, Any]] = None,
        artifact_manifest: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        db: Session = SessionLocal()
        try:
            task = self._validate_lease(db, task_id, lease_token, worker_id)
            self._upsert_distributed_execution_attempt(
                db,
                task,
                worker_id=worker_id,
                lease_token=lease_token,
                status="succeeded",
                artifact_manifest=artifact_manifest,
                result_data=result_data,
            )
            self._persist_artifact_manifest(db, task, artifact_manifest)
            task.status = "completed"  # type: ignore
            task.result_data = result_data  # type: ignore
            task.finalized_at = datetime.now(timezone.utc)  # type: ignore
            task.worker_id = None  # type: ignore
            task.lease_token = None  # type: ignore
            task.lease_expires_at = None  # type: ignore
            summary = None
            self._finalize_inline_analysis(db, task, result_data or {})
            db.add(task)
            db.commit()
            summary = dict(task.result_summary or {}) if task.result_summary else None
            self._notify_task_success(task, result_data=result_data, result_summary=summary)
        finally:
            db.close()

    def fail_execution(
        self,
        task_id: str,
        lease_token: str,
        worker_id: str,
        error_message: str,
        failure_code: Optional[str] = None,
        failure_context: Optional[Dict[str, Any]] = None,
        artifact_manifest: Optional[List[Dict[str, Any]]] = None,
    ) -> None:
        db: Session = SessionLocal()
        try:
            task = self._validate_lease(db, task_id, lease_token, worker_id)
            attempt = self._upsert_distributed_execution_attempt(
                db,
                task,
                worker_id=worker_id,
                lease_token=lease_token,
                status="runtime_failed",
                artifact_manifest=artifact_manifest,
                failure_context=failure_context,
                failure_code=failure_code,
                failure_detail=error_message,
            )
            self._persist_artifact_manifest(db, task, artifact_manifest)
            if str(task.status) == "cancel_requested":
                task.status = "canceled"  # type: ignore
                task.finalized_at = datetime.now(timezone.utc)  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                attempt.status = "canceled"  # type: ignore
                attempt.finished_at = datetime.now(timezone.utc)  # type: ignore
                db.add(attempt)
                db.add(task)
                db.commit()
                return

            should_requeue = _is_retryable_error(error_message) and task.retry_count < task.max_retries
            if should_requeue:
                task.retry_count = (task.retry_count or 0) + 1  # type: ignore
                task.status = "queued"  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                task.heartbeat_at = None  # type: ignore
                task.current_execution_attempt_id = None  # type: ignore
                task.error_message = None  # type: ignore
            else:
                task.status = "failed"  # type: ignore
                task.finalized_at = datetime.now(timezone.utc)  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                task.error_message = error_message  # type: ignore

            db.add(task)
            db.commit()
            if not should_requeue and str(task.status) == "failed":
                self._notify_task_error(task, error_message)
        finally:
            db.close()

    def requeue_expired_leases(self, queue_name: Optional[str] = None) -> int:
        db: Session = SessionLocal()
        try:
            now = datetime.now(timezone.utc)
            query = db.query(Task).filter(
                Task.status == "leased",
                Task.lease_expires_at.is_not(None),
                Task.lease_expires_at <= now,
            )
            if queue_name is not None:
                query = query.filter(Task.queue_name == queue_name)

            expired_tasks = list(query.all())

            for task in expired_tasks:
                task.status = "queued"  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                task.heartbeat_at = None  # type: ignore
                db.add(task)

            db.commit()
            return len(expired_tasks)
        finally:
            db.close()

    def mark_orphaned_running_tasks(self) -> int:
        from ..settings import settings as _settings

        db: Session = SessionLocal()
        try:
            now = datetime.now(timezone.utc)
            heartbeat_deadline = now - timedelta(seconds=_settings.worker_heartbeat_timeout_seconds)
            orphaned_tasks = list(
                db.query(Task)
                .filter(
                    Task.status.in_(["running", "uploading", "analyzing"]),
                    Task.heartbeat_at.is_not(None),
                    Task.heartbeat_at <= heartbeat_deadline,
                )
                .all()
            )

            for task in orphaned_tasks:
                if task.retry_count < task.max_retries:
                    task.retry_count = (task.retry_count or 0) + 1  # type: ignore
                    task.status = "queued"  # type: ignore
                    task.error_message = None  # type: ignore
                    task.current_execution_attempt_id = None  # type: ignore
                else:
                    task.status = "failed"  # type: ignore
                    task.finalized_at = now  # type: ignore
                    task.error_message = "Worker 心跳超时，任务失去执行进度，已终止"  # type: ignore
                task.worker_id = None  # type: ignore
                task.lease_token = None  # type: ignore
                task.lease_expires_at = None  # type: ignore
                task.heartbeat_at = None  # type: ignore
                db.add(task)

            db.commit()
            return len(orphaned_tasks)
        finally:
            db.close()

    # ------------------------------------------------------------------ #
    #  任务执行主循环 (daemon thread)
    # ------------------------------------------------------------------ #
    def _run_task_worker(self, task_id: str, cancel_event: threading.Event) -> None:
        import asyncio

        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None:
                return

            # --- 标记 running ---
            task.status = "running"  # type: ignore
            task.progress = 1  # type: ignore
            task.progress_message = "初始化..."  # type: ignore
            db.add(task)
            db.commit()

            # --- 创建 VaspWorker ---
            from ..vasp_worker import VaspWorker
            from ..Config import DOWNLOAD_URL
            user_vasp_worker = VaspWorker(
                user_id=str(task.user_id),
                base_work_dir=DOWNLOAD_URL
            )

            task_type = str(task.task_type)
            method_name = TASK_TYPE_TO_METHOD.get(task_type)
            if not method_name:
                raise Exception(f"不支持的任务类型: {task_type}")

            # --- 带重试的执行循环 ---
            max_attempts = MAX_TASK_RETRIES
            last_error = ""

            for attempt_no in range(1, max_attempts + 1):
                if cancel_event.is_set():
                    raise Exception("任务已被取消")

                attempt_id = uuid.uuid4().hex
                attempt = ExecutionAttempt(
                    id=attempt_id,
                    task_id=task_id,
                    attempt_no=attempt_no,
                    executor_type="slurm",
                    status="submitting",
                )
                db.add(attempt)
                task.current_execution_attempt_id = attempt_id  # type: ignore
                if attempt_no > 1:
                    task.progress = 1  # type: ignore
                    task.progress_message = f"第 {attempt_no} 次尝试..."  # type: ignore
                    logger.info(f"任务 {task_id[:8]}... 第 {attempt_no}/{max_attempts} 次重试")
                db.add(task)
                db.commit()

                # --- 进度回调 ---
                async def progress_callback(progress: int, message: str, pid: Optional[int] = None):
                    if cancel_event.is_set():
                        raise Exception("任务已被取消")
                    task.progress = progress  # type: ignore
                    task.progress_message = message  # type: ignore

                    if pid is not None:
                        task.process_id = pid  # type: ignore
                        attempt.external_job_id = str(pid)  # type: ignore
                        attempt.status = "running"  # type: ignore
                        if not attempt.started_at:
                            attempt.started_at = datetime.now(timezone.utc)  # type: ignore
                        db.add(attempt)
                        logger.info(f"任务 {task_id[:8]}... 作业ID: {pid}")

                    db.add(task)
                    db.commit()

                # --- 执行计算 ---
                exec_success = False
                result: Dict[str, Any] = {}
                try:
                    loop = asyncio.new_event_loop()
                    try:
                        attempt.status = "submitted"  # type: ignore
                        attempt.submitted_at = datetime.now(timezone.utc)  # type: ignore
                        db.add(attempt)
                        db.commit()

                        worker_method = getattr(user_vasp_worker, method_name)
                        result = loop.run_until_complete(
                            worker_method(task_id, task.params or {}, progress_callback)
                        )
                    finally:
                        loop.close()

                    logger.info(f"任务 {task_id[:8]}... 计算返回 (attempt {attempt_no})")

                    # --- 判断执行结果 ---
                    exec_success = result.get('success', False)
                    if not exec_success and task_type == "md_calculation":
                        if result.get('convergence', False):
                            exec_success = True

                    now = datetime.now(timezone.utc)
                    attempt.finished_at = now  # type: ignore
                    attempt.work_directory = result.get(
                        'work_directory', str(user_vasp_worker.base_work_dir / task_id)
                    )  # type: ignore

                    if exec_success:
                        attempt.status = "succeeded"  # type: ignore
                        db.add(attempt)
                        db.commit()
                        break  # 成功，退出重试循环

                    # 执行失败
                    last_error = result.get('error', result.get('error_message', '计算失败'))
                    attempt.status = "runtime_failed"  # type: ignore
                    attempt.failure_detail = last_error  # type: ignore
                    db.add(attempt)
                    db.commit()

                except Exception as attempt_exc:
                    last_error = str(attempt_exc)
                    attempt.status = "runtime_failed"  # type: ignore
                    attempt.failure_detail = last_error  # type: ignore
                    attempt.finished_at = datetime.now(timezone.utc)  # type: ignore
                    db.add(attempt)
                    db.commit()

                # 检查是否可重试
                if attempt_no < max_attempts and _is_retryable_error(last_error):
                    logger.warning(f"任务 {task_id[:8]}... 可重试错误: {last_error}")
                    continue
                else:
                    break  # 不可重试或已到最大次数，退出

            # --- 重试循环结束后处理最终结果 ---
            if exec_success:
                task.result_path = result.get('work_directory')  # type: ignore
                task.result_data = result  # type: ignore
                if result.get('process_id') and not task.process_id:
                    task.process_id = result.get('process_id')  # type: ignore

                self._run_analysis(db, task, attempt, result, user_vasp_worker)
                self._emit_completion_metrics(task_type, "completed")
            else:
                task.status = "failed"  # type: ignore
                task.error_message = last_error  # type: ignore
                db.add(task)
                db.commit()
                logger.info(f"任务 {task_id[:8]}... 最终失败 (尝试 {attempt_no} 次)")
                self._emit_completion_metrics(task_type, "failed")
                self._notify_task_error(task, last_error, result_data=result)

        except Exception as exc:
            logger.error(f"任务 {task_id[:8]}... 异常: {exc}", exc_info=True)
            task = db.get(Task, task_id)
            if task is not None:
                task.status = "failed"  # type: ignore
                task.error_message = str(exc)  # type: ignore
                db.add(task)
                db.commit()
                self._notify_task_error(task, str(exc))
        finally:
            db.close()
            with self._lock:
                self._cancel_flags.pop(task_id, None)

    # ------------------------------------------------------------------ #
    #  Prometheus 指标
    # ------------------------------------------------------------------ #
    @staticmethod
    def _emit_completion_metrics(task_type: str, status: str) -> None:
        """Update Prometheus gauges/counters when a task finishes."""
        try:
            from prometheus_client import REGISTRY
            # Check if metrics are registered
            vasp_tasks_queued = REGISTRY._names_to_collectors.get("vasp_tasks_queued")
            vasp_tasks_running = REGISTRY._names_to_collectors.get("vasp_tasks_running")
            vasp_tasks_completed = REGISTRY._names_to_collectors.get("vasp_tasks_completed_total")
            if vasp_tasks_queued:
                vasp_tasks_queued.dec()
            if vasp_tasks_completed:
                vasp_tasks_completed.labels(task_type=task_type, status=status).inc()
        except Exception:
            pass  # metrics are optional

    # ------------------------------------------------------------------ #
    #  分析阶段
    # ------------------------------------------------------------------ #
    def _run_analysis(
        self,
        db: Session,
        task: Task,
        attempt: ExecutionAttempt,
        exec_result: Dict[str, Any],
        worker,
    ) -> None:
        """执行成功后自动触发分析。分析失败不影响执行记录。"""
        task_id = str(task.id)
        task_type = str(task.task_type)
        analysis_type = TASK_TYPE_TO_ANALYSIS.get(task_type, task_type)

        # --- 创建 AnalysisRun ---
        run_id = uuid.uuid4().hex
        analysis_run = AnalysisRun(
            id=run_id,
            task_id=task_id,
            analysis_type=analysis_type,
            status="running",
            started_at=datetime.now(timezone.utc),
        )
        db.add(analysis_run)

        task.status = "analyzing"  # type: ignore
        task.analysis_status = "running"  # type: ignore
        task.current_analysis_run_id = run_id  # type: ignore
        db.add(task)
        db.commit()

        logger.info(f"任务 {task_id[:8]}... 开始分析 (run={run_id[:8]})")

        summary: Optional[Dict[str, Any]] = None
        try:
            # 分析逻辑：从 exec_result 中提取摘要
            summary = self._extract_analysis_summary(task_type, exec_result)

            # 从 exec_result 推导 html_report_url
            html_report_url = self._extract_html_report_url(exec_result)

            # 注册 artifacts (将在 Phase 4 细化)
            work_dir = exec_result.get('work_directory', str(worker.base_work_dir / task_id))
            self._register_artifacts(db, task_id, attempt.id, run_id, work_dir, exec_result)

            # --- 分析成功 ---
            analysis_run.status = "completed"  # type: ignore
            analysis_run.summary = summary  # type: ignore
            analysis_run.finished_at = datetime.now(timezone.utc)  # type: ignore
            db.add(analysis_run)

            task.analysis_status = "completed"  # type: ignore
            task.result_summary = summary  # type: ignore
            if html_report_url:
                # 在 summary 中保存 html_report_url
                if task.result_summary is None:
                    task.result_summary = {}  # type: ignore
                # JSON column — 需要重新赋值触发 dirty 检测
                updated_summary = dict(task.result_summary or {})
                updated_summary['html_report_url'] = html_report_url
                task.result_summary = updated_summary  # type: ignore

            task.status = "completed"  # type: ignore
            task.progress = 100  # type: ignore
            task.progress_message = "完成"  # type: ignore
            task.error_message = None  # type: ignore
            db.add(task)
            db.commit()

            logger.info(f"任务 {task_id[:8]}... 分析完成")
            self._notify_task_success(task, result_data=exec_result, result_summary=summary)

        except Exception as exc:
            logger.error(f"任务 {task_id[:8]}... 分析失败: {exc}", exc_info=True)

            analysis_run.status = "failed"  # type: ignore
            analysis_run.failure_detail = str(exc)  # type: ignore
            analysis_run.finished_at = datetime.now(timezone.utc)  # type: ignore
            db.add(analysis_run)

            # 分析失败 — 不覆盖成功的执行记录
            task.analysis_status = "failed"  # type: ignore
            task.status = "completed"  # type: ignore  # 执行成功，分析失败仍标记 completed
            task.progress = 100  # type: ignore
            task.error_message = f"分析失败: {str(exc)}"  # type: ignore
            # 保留 result_data (原始执行结果)
            db.add(task)
            db.commit()
            self._notify_task_success(task, result_data=exec_result, result_summary=summary)

    # ------------------------------------------------------------------ #
    #  提取分析摘要
    # ------------------------------------------------------------------ #
    @staticmethod
    def _extract_analysis_summary(task_type: str, result: Dict[str, Any]) -> Dict[str, Any]:
        """从执行结果中提取结构化摘要"""
        summary: Dict[str, Any] = {}

        if task_type == "structure_optimization":
            summary = {
                "success": result.get("success", False),
                "force_convergence": result.get("force_convergence"),
                "final_max_force": result.get("final_max_force"),
                "energy_convergence": result.get("energy_convergence"),
                "final_energy": result.get("final_energy") or result.get("energy"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "scf_calculation":
            summary = {
                "success": result.get("success", False),
                "convergence": result.get("convergence", False),
                "total_energy": result.get("total_energy") or result.get("energy"),
                "fermi_energy": result.get("fermi_energy"),
                "band_gap": result.get("band_gap"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "dos_calculation":
            summary = {
                "success": result.get("success", False),
                "convergence": result.get("convergence", False),
                "total_energy": result.get("total_energy") or result.get("energy"),
                "fermi_energy": result.get("fermi_energy"),
                "band_gap": result.get("band_gap"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "dos_analysis":
            summary = {
                "success": result.get("success", False),
                "band_gap": result.get("band_gap"),
                "fermi_energy": result.get("fermi_energy"),
                "is_metal": result.get("is_metal"),
                "total_energy": result.get("total_energy") or result.get("energy"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "band_structure":
            summary = {
                "success": result.get("success", False),
                "convergence": result.get("convergence", False),
                "total_energy": result.get("total_energy") or result.get("energy"),
                "fermi_energy": result.get("fermi_energy"),
                "band_gap": result.get("band_gap"),
                "is_direct": result.get("is_direct"),
                "vbm": result.get("vbm"),
                "cbm": result.get("cbm"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "band_structure_analysis":
            summary = {
                "success": result.get("success", False),
                "total_energy": result.get("total_energy") or result.get("energy"),
                "fermi_energy": result.get("fermi_energy"),
                "band_gap": result.get("band_gap"),
                "is_direct": result.get("is_direct"),
                "vbm": result.get("vbm"),
                "cbm": result.get("cbm"),
                "is_metal": result.get("is_metal"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "md_calculation":
            summary = {
                "convergence": result.get("convergence", False),
                "final_energy": result.get("final_energy"),
                "average_temperature": result.get("average_temperature"),
                "total_md_steps": result.get("total_md_steps"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "md_analysis":
            summary = {
                "success": result.get("success", False),
                "convergence": result.get("convergence", False),
                "final_energy": result.get("final_energy"),
                "average_temperature": result.get("average_temperature"),
                "total_md_steps": result.get("total_md_steps"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "neb_calculation":
            summary = {
                "success": result.get("success", False),
                "forward_barrier_eV": result.get("forward_barrier_eV"),
                "backward_barrier_eV": result.get("backward_barrier_eV"),
                "reaction_energy_eV": result.get("reaction_energy_eV"),
                "n_images": result.get("n_images"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "phonon_calculation":
            summary = {
                "success": result.get("success", False),
                "n_modes": result.get("n_modes"),
                "n_imaginary": result.get("n_imaginary"),
                "dynamically_stable": result.get("dynamically_stable"),
                "max_imaginary_freq_cm1": result.get("max_imaginary_freq_cm1"),
                "max_real_freq_cm1": result.get("max_real_freq_cm1"),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "custom_calculation":
            summary = {
                "success": result.get("success", False),
                "convergence": result.get("convergence", False),
                "final_energy": result.get("final_energy"),
                "fermi_energy": result.get("fermi_energy"),
                "n_ionic_steps": result.get("n_ionic_steps"),
                "n_electronic_steps": result.get("n_electronic_steps"),
                "output_files": result.get("output_files", []),
                "computation_time": result.get("computation_time"),
            }
        elif task_type == "agent_analysis":
            summary = {
                "success": result.get("success", False),
                "steps": result.get("steps"),
                "computation_time": result.get("computation_time"),
                "html_report_url": result.get("html_report_url"),
            }
        else:
            summary = {"success": result.get("success", False)}

        return summary

    # ------------------------------------------------------------------ #
    #  提取 HTML 报告 URL
    # ------------------------------------------------------------------ #
    @staticmethod
    def _extract_html_report_url(result: Dict[str, Any]) -> Optional[str]:
        """从结果中提取 html_report_url"""
        # 各种计算类型把报告路径放在不同的 key 里
        for key in [
            'html_report_url',
            'analysis_report_html_path',
            'html_analysis_report',
            'md_analysis_report_html_path',
            'dos_analysis_report_html_path',
            'scf_analysis_report_html_path',
            'band_structure_report_html_path',
            'neb_report_html_path',
            'phonon_report_html_path',
            'agent_analysis_report_html_path',
        ]:
            val = result.get(key)
            if val:
                return str(val)
        return None

    def _extract_html_report_url_from_artifacts(self, db: Session, task_id: str) -> Optional[str]:
        artifacts = (
            db.query(Artifact)
            .filter(Artifact.task_id == task_id)
            .order_by(Artifact.created_at.asc())
            .all()
        )
        for artifact in artifacts:
            storage_key = str(getattr(artifact, "storage_key", "") or "")
            object_key = str(getattr(artifact, "object_key", "") or "")
            artifact_type = str(getattr(artifact, "artifact_type", "") or "")
            if artifact_type != "html_report" and not storage_key.lower().endswith(".html") and not object_key.lower().endswith(".html"):
                continue
            download_url = self.get_artifact_download_url(_orm_to_ns(artifact))
            if download_url:
                return download_url
        return None

    def _finalize_inline_analysis(
        self,
        db: Session,
        task: Task,
        exec_result: Dict[str, Any],
    ) -> None:
        task_id = str(task.id)
        task_type = str(task.task_type)
        analysis_type = TASK_TYPE_TO_ANALYSIS.get(task_type, task_type)

        run_id = uuid.uuid4().hex
        analysis_run = AnalysisRun(
            id=run_id,
            task_id=task_id,
            analysis_type=analysis_type,
            status="running",
            started_at=datetime.now(timezone.utc),
        )
        db.add(analysis_run)

        task.analysis_status = "running"  # type: ignore
        task.current_analysis_run_id = run_id  # type: ignore
        db.add(task)
        db.flush()

        try:
            summary = self._extract_analysis_summary(task_type, exec_result)
            html_report_url = self._extract_html_report_url(exec_result) or self._extract_html_report_url_from_artifacts(db, task_id)

            analysis_run.status = "completed"  # type: ignore
            analysis_run.summary = summary  # type: ignore
            analysis_run.finished_at = datetime.now(timezone.utc)  # type: ignore
            db.add(analysis_run)

            task.analysis_status = "completed"  # type: ignore
            task.result_summary = summary  # type: ignore
            if html_report_url:
                updated_summary = dict(task.result_summary or {})
                updated_summary["html_report_url"] = html_report_url
                task.result_summary = updated_summary  # type: ignore
            db.add(task)
        except Exception as exc:
            analysis_run.status = "failed"  # type: ignore
            analysis_run.failure_detail = str(exc)  # type: ignore
            analysis_run.finished_at = datetime.now(timezone.utc)  # type: ignore
            db.add(analysis_run)

            task.analysis_status = "failed"  # type: ignore
            task.error_message = f"分析失败: {str(exc)}"  # type: ignore
            db.add(task)

    # ------------------------------------------------------------------ #
    #  注册 artifacts (最小版本)
    # ------------------------------------------------------------------ #
    def _register_artifacts(
        self,
        db: Session,
        task_id: str,
        attempt_id: str,
        analysis_run_id: str,
        work_dir: str,
        exec_result: Dict[str, Any],
    ) -> None:
        """注册产物到 artifacts 表"""
        import os

        # 执行产物
        exec_files = {
            "outcar": "OUTCAR",
            "contcar": "CONTCAR",
            "poscar": "POSCAR",
            "potcar": "POTCAR",
            "doscar": "DOSCAR",
            "xdatcar": "XDATCAR",
            "eigenval": "EIGENVAL",
            "procar": "PROCAR",
            "result_log": "result.log",
            "incar": "INCAR",
            "kpoints": "KPOINTS",
        }

        for artifact_type, filename in exec_files.items():
            filepath = os.path.join(work_dir, filename)
            if os.path.exists(filepath):
                try:
                    size = os.path.getsize(filepath)
                except OSError:
                    size = None
                artifact = Artifact(
                    id=uuid.uuid4().hex,
                    task_id=task_id,
                    owner_type="execution",
                    owner_id=attempt_id,
                    artifact_type=artifact_type,
                    storage_backend="local",
                    storage_key=filepath,
                    object_key=None,
                    size_bytes=size,
                )
                db.add(artifact)

        # 分析产物: HTML 报告
        html_url = self._extract_html_report_url(exec_result)
        if html_url:
            artifact = Artifact(
                id=uuid.uuid4().hex,
                task_id=task_id,
                owner_type="analysis",
                owner_id=analysis_run_id,
                artifact_type="html_report",
                storage_backend="local",
                storage_key=html_url,
                mime_type="text/html",
                content_type="text/html",
            )
            db.add(artifact)

        db.commit()

    # ------------------------------------------------------------------ #
    #  取消任务
    # ------------------------------------------------------------------ #
    def cancel_task(self, task_id: str, user_id: str) -> bool:
        return self.request_cancel(task_id, user_id)

    # ------------------------------------------------------------------ #
    #  查询
    # ------------------------------------------------------------------ #
    def get_task(self, task_id: str, user_id: str) -> Optional[SimpleNamespace]:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None or str(task.user_id) != user_id:
                return None
            return _orm_to_ns(task)
        finally:
            db.close()

    def list_tasks(self, user_id: str) -> List[SimpleNamespace]:
        db: Session = SessionLocal()
        try:
            tasks = list(
                db.query(Task)
                .filter(Task.user_id == user_id)
                .order_by(Task.created_at.desc())
                .all()
            )
            return [_orm_to_ns(t) for t in tasks]
        finally:
            db.close()

    def get_task_artifacts(self, task_id: str) -> List[SimpleNamespace]:
        """查询任务的所有 artifacts"""
        db: Session = SessionLocal()
        try:
            rows = list(
                db.query(Artifact)
                .filter(Artifact.task_id == task_id)
                .order_by(Artifact.created_at)
                .all()
            )
            return [_orm_to_ns(a) for a in rows]
        finally:
            db.close()

    def get_task_artifact_manifest(self, task_id: str) -> List[Dict[str, Any]]:
        """构造可供下游任务复用的产物清单。"""
        manifest: List[Dict[str, Any]] = []
        for artifact in self.get_task_artifacts(task_id):
            storage_key = getattr(artifact, "storage_key", None)
            object_key = getattr(artifact, "object_key", None)
            filename = os.path.basename(str(object_key or storage_key or getattr(artifact, "artifact_type", "artifact")))
            content_type = getattr(artifact, "content_type", None) or getattr(artifact, "mime_type", None)
            item: Dict[str, Any] = {
                "id": str(getattr(artifact, "id")),
                "artifact_type": str(getattr(artifact, "artifact_type")),
                "filename": filename,
                "storage_backend": getattr(artifact, "storage_backend", None),
                "storage_key": storage_key,
                "bucket": getattr(artifact, "bucket", None),
                "object_key": object_key,
                "content_type": content_type,
                "size_bytes": getattr(artifact, "size_bytes", None),
                "etag": getattr(artifact, "etag", None),
                "sha256": getattr(artifact, "sha256", None),
            }

            if getattr(artifact, "storage_backend", None) == "local" and storage_key:
                item["local_path"] = str(storage_key)

            download_url = self.get_artifact_download_url(artifact)
            if download_url:
                item["download_url"] = download_url

            manifest.append({key: value for key, value in item.items() if value is not None})

        return manifest

    def get_analysis_run(self, task_id: str) -> Optional[SimpleNamespace]:
        """获取任务最新的分析记录"""
        db: Session = SessionLocal()
        try:
            row = (
                db.query(AnalysisRun)
                .filter(AnalysisRun.task_id == task_id)
                .order_by(AnalysisRun.created_at.desc())
                .first()
            )
            return _orm_to_ns(row)
        finally:
            db.close()

    def get_execution_attempt(self, task_id: str) -> Optional[SimpleNamespace]:
        """获取任务最新的执行记录"""
        db: Session = SessionLocal()
        try:
            row = (
                db.query(ExecutionAttempt)
                .filter(ExecutionAttempt.task_id == task_id)
                .order_by(ExecutionAttempt.created_at.desc())
                .first()
            )
            return _orm_to_ns(row)
        finally:
            db.close()
