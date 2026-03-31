from fastapi import APIRouter, Depends, Header, HTTPException, Response, status

from .schemas import (
    TaskStatus,
    WorkerClaimRequest,
    WorkerClaimResponse,
    WorkerCompleteRequest,
    WorkerFailRequest,
    WorkerHeartbeatResponse,
    WorkerLeaseRequest,
    WorkerRegisterRequest,
    WorkerRegisterResponse,
    WorkerTaskStatusResponse,
)
from .settings import settings
from .task_manager.database import SessionLocal
from .task_manager.models import Task


def build_internal_worker_router(task_manager) -> APIRouter:
    router = APIRouter(tags=["internal-workers"])

    def require_internal_worker_auth(
        authorization: str | None = Header(default=None),
    ) -> None:
        expected = f"Bearer {settings.internal_worker_token}"
        if authorization != expected:
            raise HTTPException(
                status_code=status.HTTP_401_UNAUTHORIZED,
                detail="invalid worker token",
            )

    def load_task_snapshot(task_id: str) -> dict:
        db = SessionLocal()
        try:
            task = db.get(Task, task_id)
            if task is None:
                raise HTTPException(status_code=404, detail=f"task not found: {task_id}")
            return {
                "id": str(task.id),
                "status": str(task.status),
                "worker_id": str(task.worker_id) if task.worker_id else None,
                "lease_token": str(task.lease_token) if task.lease_token else None,
                "lease_expires_at": task.lease_expires_at,
                "task_type": str(task.task_type),
                "queue_name": str(task.queue_name),
                "params": task.params or {},
            }
        finally:
            db.close()

    @router.post(
        "/internal/workers/register",
        response_model=WorkerRegisterResponse,
    )
    def register_worker(
        request: WorkerRegisterRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerRegisterResponse:
        return WorkerRegisterResponse(
            worker_id=request.worker_id,
            queue_name=request.queue_name,
            poll_interval_seconds=settings.worker_poll_interval_seconds,
        )

    @router.post(
        "/internal/workers/claim",
        response_model=WorkerClaimResponse,
    )
    def claim_task(
        request: WorkerClaimRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ):
        claimed = task_manager.claim_next_task(
            worker_id=request.worker_id,
            queue_name=request.queue_name,
        )
        if claimed is None:
            return Response(status_code=status.HTTP_204_NO_CONTENT)

        params = getattr(claimed, "params", {}) or {}
        return WorkerClaimResponse(
            task_id=str(claimed.id),
            status=TaskStatus(str(claimed.status)),
            worker_id=str(claimed.worker_id),
            lease_token=str(claimed.lease_token),
            lease_expires_at=claimed.lease_expires_at,
            task_type=str(claimed.task_type),
            queue_name=str(claimed.queue_name),
            params=params,
        )

    @router.post(
        "/internal/tasks/{task_id}/heartbeat",
        response_model=WorkerHeartbeatResponse,
    )
    def heartbeat_task(
        task_id: str,
        request: WorkerLeaseRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerHeartbeatResponse:
        try:
            task_manager.heartbeat_task(task_id, request.lease_token, request.worker_id)
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        snapshot = load_task_snapshot(task_id)
        return WorkerHeartbeatResponse(
            task_id=task_id,
            status=TaskStatus(snapshot["status"]),
            lease_expires_at=snapshot["lease_expires_at"],
            cancel_requested=snapshot["status"] == "cancel_requested",
        )

    @router.post(
        "/internal/tasks/{task_id}/running",
        response_model=WorkerTaskStatusResponse,
    )
    def mark_task_running(
        task_id: str,
        request: WorkerLeaseRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerTaskStatusResponse:
        try:
            task_manager.mark_task_running(task_id, request.lease_token, request.worker_id)
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        return WorkerTaskStatusResponse(task_id=task_id, status=TaskStatus.running)

    @router.post(
        "/internal/tasks/{task_id}/complete",
        response_model=WorkerTaskStatusResponse,
    )
    def complete_task(
        task_id: str,
        request: WorkerCompleteRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerTaskStatusResponse:
        try:
            task_manager.complete_execution(
                task_id=task_id,
                lease_token=request.lease_token,
                worker_id=request.worker_id,
                result_data=request.result_data,
                artifact_manifest=request.artifact_manifest,
            )
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        return WorkerTaskStatusResponse(task_id=task_id, status=TaskStatus.completed)

    @router.post(
        "/internal/tasks/{task_id}/fail",
        response_model=WorkerTaskStatusResponse,
    )
    def fail_task(
        task_id: str,
        request: WorkerFailRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerTaskStatusResponse:
        try:
            task_manager.fail_execution(
                task_id=task_id,
                lease_token=request.lease_token,
                worker_id=request.worker_id,
                error_message=request.error_message,
                artifact_manifest=request.artifact_manifest,
            )
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        snapshot = load_task_snapshot(task_id)
        return WorkerTaskStatusResponse(
            task_id=task_id,
            status=TaskStatus(snapshot["status"]),
        )

    @router.post(
        "/internal/tasks/{task_id}/cancel-ack",
        response_model=WorkerTaskStatusResponse,
    )
    def cancel_ack(
        task_id: str,
        request: WorkerLeaseRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerTaskStatusResponse:
        try:
            task_manager.ack_cancel(task_id, request.lease_token, request.worker_id)
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        return WorkerTaskStatusResponse(task_id=task_id, status=TaskStatus.canceled)

    return router
