import uuid
from sqlalchemy import Column, String, Integer, DateTime, Text, Float, ForeignKey
from sqlalchemy.sql import func
from sqlalchemy.types import JSON

from .database import Base


def _gen_id():
    return uuid.uuid4().hex


class Task(Base):
    """主任务表 — 对外主对象，保持与现有 API 兼容"""
    __tablename__ = "tasks"

    id = Column(String, primary_key=True, index=True)
    user_id = Column(String, index=True, nullable=False)
    task_type = Column(String, nullable=False)

    # --- 状态 ---
    # queued / running / analyzing / completed / failed / cancel_requested / canceled
    status = Column(String, index=True, nullable=False, default="queued")
    # pending / running / completed / failed
    analysis_status = Column(String, nullable=False, default="pending")

    progress = Column(Integer, nullable=False, default=0)
    # 进度消息 (不再写到 error_message)
    progress_message = Column(Text, nullable=True)

    # --- 输入快照 ---
    params = Column(JSON, nullable=True)
    input_snapshot = Column(JSON, nullable=True)

    # --- 结果 ---
    result_summary = Column(JSON, nullable=True)  # 来自 AnalysisRun.summary 回写
    result_data = Column(JSON, nullable=True)      # 兼容：保留原始完整结果
    result_path = Column(String, nullable=True)    # 兼容：保留

    # --- 当前关联 ---
    current_execution_attempt_id = Column(String, nullable=True)
    current_analysis_run_id = Column(String, nullable=True)

    # --- 兼容字段 (逐步废弃) ---
    external_job_id = Column(String, nullable=True)
    process_id = Column(Integer, nullable=True)
    error_message = Column(Text, nullable=True)

    # --- 时间戳 ---
    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False
    )


class ExecutionAttempt(Base):
    """执行尝试记录 — 一个 Task 可以有多次执行尝试"""
    __tablename__ = "execution_attempts"

    id = Column(String, primary_key=True, default=_gen_id)
    task_id = Column(String, ForeignKey("tasks.id"), index=True, nullable=False)
    attempt_no = Column(Integer, nullable=False, default=1)
    executor_type = Column(String, nullable=False, default="slurm")  # slurm / local / ...

    # submitting / submitted / running / succeeded / submit_failed / runtime_failed / canceled
    status = Column(String, nullable=False, default="submitting")

    external_job_id = Column(String, nullable=True)  # e.g. SLURM job id
    work_directory = Column(String, nullable=True)

    failure_code = Column(String, nullable=True)
    failure_detail = Column(Text, nullable=True)

    submitted_at = Column(DateTime(timezone=True), nullable=True)
    started_at = Column(DateTime(timezone=True), nullable=True)
    finished_at = Column(DateTime(timezone=True), nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)


class AnalysisRun(Base):
    """分析运行记录 — 独立于执行，支持重跑"""
    __tablename__ = "analysis_runs"

    id = Column(String, primary_key=True, default=_gen_id)
    task_id = Column(String, ForeignKey("tasks.id"), index=True, nullable=False)

    analysis_type = Column(String, nullable=False)  # optimization / scf / dos / md
    # pending / running / completed / failed
    status = Column(String, nullable=False, default="pending")

    parser_version = Column(String, nullable=True)
    analyzer_version = Column(String, nullable=True)

    summary = Column(JSON, nullable=True)  # 结构化摘要，回写 Task.result_summary

    failure_code = Column(String, nullable=True)
    failure_detail = Column(Text, nullable=True)

    started_at = Column(DateTime(timezone=True), nullable=True)
    finished_at = Column(DateTime(timezone=True), nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)


class Artifact(Base):
    """产物索引 — 统一管理所有文件产物"""
    __tablename__ = "artifacts"

    id = Column(String, primary_key=True, default=_gen_id)
    task_id = Column(String, ForeignKey("tasks.id"), index=True, nullable=False)

    # owner_type: execution / analysis / task
    owner_type = Column(String, nullable=False, default="task")
    owner_id = Column(String, nullable=True)  # ExecutionAttempt.id 或 AnalysisRun.id

    # artifact_type: outcar / contcar / doscar / xdatcar / result_log /
    #                html_report / png / json / poscar / potcar / incar / ...
    artifact_type = Column(String, nullable=False)

    storage_backend = Column(String, nullable=False, default="local")
    storage_key = Column(String, nullable=False)  # 本地路径 (内部)，对外不直接暴露

    mime_type = Column(String, nullable=True)
    size_bytes = Column(Float, nullable=True)
    checksum = Column(String, nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
