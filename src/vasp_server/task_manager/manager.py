import threading
import time
import uuid
import logging
from datetime import datetime, timezone
from typing import Any, Dict, List, Optional

from sqlalchemy.orm import Session

from .database import SessionLocal
from .models import Task, ExecutionAttempt, AnalysisRun, Artifact

logger = logging.getLogger(__name__)

# 支持的任务类型 → 分析类型映射
TASK_TYPE_TO_ANALYSIS = {
    "structure_optimization": "optimization",
    "scf_calculation": "scf",
    "dos_calculation": "dos",
    "md_calculation": "md",
}

# 支持的任务类型 → VaspWorker 方法名映射
TASK_TYPE_TO_METHOD = {
    "structure_optimization": "run_structure_optimization",
    "scf_calculation": "run_scf_calculation",
    "dos_calculation": "run_dos_calculation",
    "md_calculation": "run_md_calculation",
}


class TaskManager:
    def __init__(self) -> None:
        self._cancel_flags: Dict[str, threading.Event] = {}
        self._lock = threading.Lock()

    # ------------------------------------------------------------------ #
    #  提交任务
    # ------------------------------------------------------------------ #
    def submit_task(self, user_id: str, task_type: str, params: Optional[Dict[str, Any]]) -> str:
        task_id = uuid.uuid4().hex
        cancel_event = threading.Event()
        with self._lock:
            self._cancel_flags[task_id] = cancel_event

        db: Session = SessionLocal()
        try:
            task = Task(
                id=task_id,
                user_id=user_id,
                task_type=task_type,
                status="queued",
                analysis_status="pending",
                progress=0,
                params=params or {},
                input_snapshot=params or {},
            )
            db.add(task)
            db.commit()
        finally:
            db.close()

        worker_thread = threading.Thread(
            target=self._run_task_worker, args=(task_id, cancel_event), daemon=True
        )
        worker_thread.start()

        return task_id

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
            user_vasp_worker = VaspWorker(
                user_id=str(task.user_id),
                base_work_dir="/data/home/ysl9527/vasp_calculations"
            )

            task_type = str(task.task_type)
            method_name = TASK_TYPE_TO_METHOD.get(task_type)
            if not method_name:
                raise Exception(f"不支持的任务类型: {task_type}")

            # --- 创建 ExecutionAttempt ---
            attempt_id = uuid.uuid4().hex
            attempt = ExecutionAttempt(
                id=attempt_id,
                task_id=task_id,
                attempt_no=1,
                executor_type="slurm",
                status="submitting",
            )
            db.add(attempt)
            task.current_execution_attempt_id = attempt_id  # type: ignore
            db.add(task)
            db.commit()

            # --- 进度回调 ---
            async def progress_callback(progress: int, message: str, pid: Optional[int] = None):
                if cancel_event.is_set():
                    raise Exception("任务已被取消")
                task.progress = progress  # type: ignore
                task.progress_message = message  # type: ignore
                # 不再写到 error_message

                if pid is not None:
                    task.process_id = pid  # type: ignore
                    # 也更新 ExecutionAttempt 的 external_job_id
                    attempt.external_job_id = str(pid)  # type: ignore
                    attempt.status = "running"  # type: ignore
                    if not attempt.started_at:
                        attempt.started_at = datetime.now(timezone.utc)  # type: ignore
                    db.add(attempt)
                    logger.info(f"任务 {task_id[:8]}... 作业ID: {pid}")

                db.add(task)
                db.commit()

            # --- 执行计算 ---
            loop = asyncio.new_event_loop()
            asyncio.set_event_loop(loop)
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

            logger.info(f"任务 {task_id[:8]}... 计算返回")

            # --- 判断执行结果 ---
            # 对 MD 任务做兼容：multi-temperature 没有顶层 success 字段
            exec_success = result.get('success', False)
            if not exec_success and task_type == "md_calculation":
                # MD multi-temp: 有 completed_subtasks > 0 就算成功
                if result.get('is_multi_temperature') and result.get('completed_subtasks', 0) > 0:
                    exec_success = True
                # 单温度 MD: convergence=True 也算成功
                if result.get('convergence', False):
                    exec_success = True

            now = datetime.now(timezone.utc)
            attempt.finished_at = now  # type: ignore
            attempt.work_directory = result.get('work_directory', str(user_vasp_worker.base_work_dir / task_id))  # type: ignore

            if exec_success:
                attempt.status = "succeeded"  # type: ignore
                db.add(attempt)
                db.commit()

                # 保存兼容字段
                task.result_path = result.get('work_directory')  # type: ignore
                task.result_data = result  # type: ignore
                if result.get('process_id') and not task.process_id:
                    task.process_id = result.get('process_id')  # type: ignore

                # --- 自动触发分析 ---
                self._run_analysis(db, task, attempt, result, user_vasp_worker)
            else:
                attempt.status = "runtime_failed"  # type: ignore
                attempt.failure_detail = result.get('error', result.get('error_message', '计算失败'))  # type: ignore
                db.add(attempt)

                task.status = "failed"  # type: ignore
                task.error_message = result.get('error', result.get('error_message', '计算失败'))  # type: ignore
                db.add(task)
                db.commit()
                logger.info(f"任务 {task_id[:8]}... 执行失败")

        except Exception as exc:
            logger.error(f"任务 {task_id[:8]}... 异常: {exc}", exc_info=True)
            task = db.get(Task, task_id)
            if task is not None:
                task.status = "failed"  # type: ignore
                task.error_message = str(exc)  # type: ignore
                db.add(task)
                db.commit()
        finally:
            db.close()
            with self._lock:
                self._cancel_flags.pop(task_id, None)

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
        elif task_type == "md_calculation":
            if result.get("is_multi_temperature"):
                summary = {
                    "is_multi_temperature": True,
                    "total_subtasks": result.get("total_subtasks"),
                    "completed_subtasks": result.get("completed_subtasks"),
                    "failed_subtasks": result.get("failed_subtasks"),
                    "computation_time": result.get("computation_time"),
                }
            else:
                summary = {
                    "convergence": result.get("convergence", False),
                    "final_energy": result.get("final_energy"),
                    "average_temperature": result.get("average_temperature"),
                    "total_md_steps": result.get("total_md_steps"),
                    "computation_time": result.get("computation_time"),
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
            'analysis_report_html_path',
            'html_analysis_report',
            'md_analysis_report_html_path',
            'dos_analysis_report_html_path',
            'scf_analysis_report_html_path',
        ]:
            val = result.get(key)
            if val:
                return str(val)
        return None

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
            )
            db.add(artifact)

        db.commit()

    # ------------------------------------------------------------------ #
    #  取消任务
    # ------------------------------------------------------------------ #
    def cancel_task(self, task_id: str, user_id: str) -> bool:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None or str(task.user_id) != user_id:
                return False
            if str(task.status) in {"completed", "failed", "canceled"}:
                return True

            with self._lock:
                cancel_event = self._cancel_flags.get(task_id)
                if cancel_event:
                    cancel_event.set()

            task.status = "cancel_requested"  # type: ignore
            db.add(task)
            db.commit()
            return True
        finally:
            db.close()

    # ------------------------------------------------------------------ #
    #  查询
    # ------------------------------------------------------------------ #
    def get_task(self, task_id: str, user_id: str) -> Optional[Task]:
        db: Session = SessionLocal()
        try:
            task: Optional[Task] = db.get(Task, task_id)
            if task is None or str(task.user_id) != user_id:
                return None
            # Eagerly load all attributes while session is open
            _ = task.id, task.status, task.result_data, task.result_summary
            _ = task.analysis_status, task.progress_message
            return task
        finally:
            db.close()

    def list_tasks(self, user_id: str) -> List[Task]:
        db: Session = SessionLocal()
        try:
            tasks = list(
                db.query(Task)
                .filter(Task.user_id == user_id)
                .order_by(Task.created_at.desc())
                .all()
            )
            # Eagerly load attributes
            for t in tasks:
                _ = t.id, t.status, t.result_summary, t.analysis_status
            return tasks
        finally:
            db.close()

    def get_task_artifacts(self, task_id: str) -> List[Artifact]:
        """查询任务的所有 artifacts"""
        db: Session = SessionLocal()
        try:
            return list(
                db.query(Artifact)
                .filter(Artifact.task_id == task_id)
                .order_by(Artifact.created_at)
                .all()
            )
        finally:
            db.close()

    def get_analysis_run(self, task_id: str) -> Optional[AnalysisRun]:
        """获取任务最新的分析记录"""
        db: Session = SessionLocal()
        try:
            return (
                db.query(AnalysisRun)
                .filter(AnalysisRun.task_id == task_id)
                .order_by(AnalysisRun.created_at.desc())
                .first()
            )
        finally:
            db.close()

    def get_execution_attempt(self, task_id: str) -> Optional[ExecutionAttempt]:
        """获取任务最新的执行记录"""
        db: Session = SessionLocal()
        try:
            return (
                db.query(ExecutionAttempt)
                .filter(ExecutionAttempt.task_id == task_id)
                .order_by(ExecutionAttempt.created_at.desc())
                .first()
            )
        finally:
            db.close()
