from pathlib import Path, PurePosixPath

from fastapi import APIRouter, Depends, Header, HTTPException, Request, Response, status

from .schemas import (
    TaskStatus,
    WorkerClaimRequest,
    WorkerClaimResponse,
    WorkerArtifactUploadResponse,
    WorkerCompleteRequest,
    WorkerFailRequest,
    WorkerHeartbeatResponse,
    WorkerLeaseRequest,
    WorkerRegisterRequest,
    WorkerRegisterResponse,
    WorkerResumeRequest,
    WorkerTaskStatusResponse,
)
from .settings import settings
from .task_manager.database import SessionLocal
from .task_manager.models import ExecutionAttempt, Task


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

    def current_attempt_no(task_id: str) -> tuple[Task, int]:
        db = SessionLocal()
        try:
            task = db.get(Task, task_id)
            if task is None:
                raise HTTPException(status_code=404, detail=f"task not found: {task_id}")
            attempt_no = int(getattr(task, "retry_count", 0) or 0) + 1
            current_attempt_id = getattr(task, "current_execution_attempt_id", None)
            if current_attempt_id:
                attempt = db.get(ExecutionAttempt, current_attempt_id)
                if attempt is not None and getattr(attempt, "attempt_no", None):
                    attempt_no = int(attempt.attempt_no)
            return task, attempt_no
        finally:
            db.close()

    def sanitize_artifact_path(artifact_path: str) -> str:
        pure = PurePosixPath(artifact_path)
        if pure.is_absolute() or ".." in pure.parts:
            raise HTTPException(status_code=400, detail="invalid artifact path")
        normalized = pure.as_posix().lstrip("./")
        if not normalized:
            raise HTTPException(status_code=400, detail="empty artifact path")
        return normalized

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
        "/internal/tasks/{task_id}/resume",
        response_model=WorkerClaimResponse,
    )
    def resume_task(
        task_id: str,
        request: WorkerResumeRequest,
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerClaimResponse:
        try:
            resumed = task_manager.resume_task(
                task_id=task_id,
                worker_id=request.worker_id,
                queue_name=request.queue_name,
            )
        except ValueError as exc:
            raise HTTPException(status_code=409, detail=str(exc)) from exc

        params = getattr(resumed, "params", {}) or {}
        return WorkerClaimResponse(
            task_id=str(resumed.id),
            status=TaskStatus(str(resumed.status)),
            worker_id=str(resumed.worker_id),
            lease_token=str(resumed.lease_token),
            lease_expires_at=resumed.lease_expires_at,
            task_type=str(resumed.task_type),
            queue_name=str(resumed.queue_name),
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
            snapshot = load_task_snapshot(task_id)
            if snapshot["status"] == "completed":
                return WorkerTaskStatusResponse(task_id=task_id, status=TaskStatus.completed)
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
                failure_code=request.failure_code,
                failure_context=request.failure_context,
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

    @router.put(
        "/internal/tasks/{task_id}/artifacts/{artifact_path:path}",
        response_model=WorkerArtifactUploadResponse,
    )
    async def upload_public_artifact(
        task_id: str,
        artifact_path: str,
        request: Request,
        content_type: str | None = Header(default=None, alias="Content-Type"),
        _auth: None = Depends(require_internal_worker_auth),
    ) -> WorkerArtifactUploadResponse:
        task, attempt_no = current_attempt_no(task_id)
        artifact_rel_path = sanitize_artifact_path(artifact_path)
        body = await request.body()
        if not body:
            raise HTTPException(status_code=400, detail="empty artifact payload")

        location = task_manager.storage_service.build_public_location(
            tenant_id=str(task.tenant_id),
            task_id=str(task.id),
            attempt_no=attempt_no,
            filename=artifact_rel_path,
        )
        target_path = Path(location.storage_key)
        target_path.parent.mkdir(parents=True, exist_ok=True)
        target_path.write_bytes(body)

        download_url = task_manager.storage_service.create_public_download_url(location.object_key)
        return WorkerArtifactUploadResponse(
            storage_backend=location.storage_backend,
            storage_key=location.storage_key,
            object_key=location.object_key,
            download_url=download_url,
            content_type=content_type,
            size_bytes=float(len(body)),
        )

    return router
