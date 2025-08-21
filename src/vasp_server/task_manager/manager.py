import threading
import time
import uuid
from typing import Any, Dict, Optional

from sqlalchemy.orm import Session

from .database import SessionLocal
from .models import Task
from ..vasp_worker import VaspWorker


class TaskManager:
    def __init__(self) -> None:
        self._cancel_flags: Dict[str, threading.Event] = {}
        self._lock = threading.Lock()

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
                progress=0,
                params=params or {},
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

    def _run_task_worker(self, task_id: str, cancel_event: threading.Event) -> None:
        db: Session = SessionLocal()
        try:
            task: Task = db.get(Task, task_id)  # type: ignore
            if task is None:
                return
            task.status = "running"  # type: ignore
            task.progress = 1  # type: ignore
            db.add(task)
            db.commit()

            # 创建进度回调函数
            async def progress_callback(progress: int, message: str, pid: Optional[int] = None):
                if cancel_event.is_set():
                    raise Exception("任务已被取消")
                task.progress = progress  # type: ignore
                task.error_message = message  # type: ignore  # 临时存储进度信息
                
                # 如果提供了PID，则设置到任务中
                if pid is not None:
                    task.process_id = pid  # type: ignore
                    print(f"📝 任务 {task_id[:8]}... 的进程ID已设置: {pid}")
                
                db.add(task)
                db.commit()
            
            # 运行相应的任务类型
            if str(task.task_type) == "structure_optimization":  # type: ignore
                import asyncio
                from ..vasp_worker import VaspWorker
                
                # 为当前用户创建专用的 VaspWorker
                user_vasp_worker = VaspWorker(
                    user_id=str(task.user_id),  # type: ignore
                    base_work_dir="/data/home/ysl9527/vasp_calculations"
                )
                
                # 创建新的事件循环
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                try:
                    result = loop.run_until_complete(
                        user_vasp_worker.run_structure_optimization(
                            task_id, task.params or {}, progress_callback  # type: ignore
                        )
                    )
                    
                    if result.get('success'):
                        task.status = "completed"  # type: ignore
                        task.progress = 100  # type: ignore
                        task.result_path = result.get('work_directory')  # type: ignore
                        task.error_message = None  # type: ignore
                        # 确保PID被设置（如果还没有设置的话）
                        if result.get('process_id') and not task.process_id:  # type: ignore
                            task.process_id = result.get('process_id')  # type: ignore
                    else:
                        task.status = "failed"  # type: ignore
                        task.error_message = result.get('error', '计算失败')  # type: ignore
                finally:
                    loop.close()
            elif str(task.task_type) == "scf_calculation":  # type: ignore
                import asyncio
                from ..vasp_worker import VaspWorker
                
                # 为当前用户创建专用的 VaspWorker
                user_vasp_worker = VaspWorker(
                    user_id=str(task.user_id),  # type: ignore
                    base_work_dir="/data/home/ysl9527/vasp_calculations"
                )
                
                # 创建新的事件循环
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                try:
                    result = loop.run_until_complete(
                        user_vasp_worker.run_scf_calculation(
                            task_id, task.params or {}, progress_callback  # type: ignore
                        )
                    )
                    
                    if result.get('success'):
                        task.status = "completed"  # type: ignore
                        task.progress = 100  # type: ignore
                        task.result_path = result.get('work_directory')  # type: ignore
                        task.error_message = None  # type: ignore
                        # 确保PID被设置（如果还没有设置的话）
                        if result.get('process_id') and not task.process_id:  # type: ignore
                            task.process_id = result.get('process_id')  # type: ignore
                    else:
                        task.status = "failed"  # type: ignore
                        task.error_message = result.get('error', '自洽场计算失败')  # type: ignore
                finally:
                    loop.close()
            elif str(task.task_type) == "dos_calculation":  # type: ignore
                import asyncio
                from ..vasp_worker import VaspWorker
                
                # 为当前用户创建专用的 VaspWorker
                user_vasp_worker = VaspWorker(
                    user_id=str(task.user_id),  # type: ignore
                    base_work_dir="/data/home/ysl9527/vasp_calculations"
                )
                
                # 创建新的事件循环
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
                try:
                    result = loop.run_until_complete(
                        user_vasp_worker.run_dos_calculation(
                            task_id, task.params or {}, progress_callback  # type: ignore
                        )
                    )
                    
                    if result.get('success'):
                        task.status = "completed"  # type: ignore
                        task.progress = 100  # type: ignore
                        task.result_path = result.get('work_directory')  # type: ignore
                        task.error_message = None  # type: ignore
                        # 确保PID被设置（如果还没有设置的话）
                        if result.get('process_id') and not task.process_id:  # type: ignore
                            task.process_id = result.get('process_id')  # type: ignore
                    else:
                        task.status = "failed"  # type: ignore
                        task.error_message = result.get('error', '态密度计算失败')  # type: ignore
                finally:
                    loop.close()
            else:
                # 模拟其他类型的任务
                for step in range(1, 21):
                    if cancel_event.is_set():
                        task.status = "canceled"  # type: ignore
                        db.add(task)
                        db.commit()
                        return
                    time.sleep(1)
                    task.progress = min(100, int(step * 5))  # type: ignore
                    db.add(task)
                    db.commit()
                
                task.status = "completed"  # type: ignore
                task.progress = 100  # type: ignore
                task.result_path = None  # type: ignore
                
            db.add(task)
            db.commit()
        except Exception as exc:  # noqa: BLE001
            task = db.get(Task, task_id)  # type: ignore
            if task is not None:
                task.status = "failed"  # type: ignore
                task.error_message = str(exc)  # type: ignore
                db.add(task)
                db.commit()
        finally:
            db.close()
            with self._lock:
                self._cancel_flags.pop(task_id, None)

    def cancel_task(self, task_id: str, user_id: str) -> bool:
        db: Session = SessionLocal()
        try:
            task: Task = db.get(Task, task_id)  # type: ignore
            if task is None or str(task.user_id) != user_id:  # type: ignore
                return False
            if str(task.status) in {"completed", "failed", "canceled"}:  # type: ignore
                return True

            with self._lock:
                cancel_event = self._cancel_flags.get(task_id)
                if cancel_event:
                    cancel_event.set()

            task.status = "canceling"  # type: ignore
            db.add(task)
            db.commit()
            return True
        finally:
            db.close()

    def get_task(self, task_id: str, user_id: str) -> Optional[Task]:
        db: Session = SessionLocal()
        try:
            task: Task = db.get(Task, task_id)  # type: ignore
            if task is None or str(task.user_id) != user_id:  # type: ignore
                return None
            return task
        finally:
            db.close()

    def list_tasks(self, user_id: str) -> list[Task]:
        db: Session = SessionLocal()
        try:
            return list(db.query(Task).filter(Task.user_id == user_id).order_by(Task.created_at.desc()).all())
        finally:
            db.close() 