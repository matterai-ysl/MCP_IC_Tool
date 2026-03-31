from pathlib import Path

from src.vasp_server.storage import ObjectStorageService
from src.vasp_server.task_manager.database import SessionLocal
from src.vasp_server.task_manager.models import Artifact, ExecutionAttempt, Task


def backfill_local_artifacts(session, storage_service: ObjectStorageService | None = None) -> int:
    storage_service = storage_service or ObjectStorageService.from_settings()
    migrated = 0

    rows = list(
        session.query(Artifact)
        .filter(Artifact.storage_backend == "local")
        .all()
    )

    for artifact in rows:
        if artifact.object_key:
            continue

        task = session.get(Task, artifact.task_id)
        if task is None:
            continue

        attempt_no = 1
        if artifact.owner_type == "execution" and artifact.owner_id:
            attempt = session.get(ExecutionAttempt, artifact.owner_id)
            if attempt is not None and attempt.attempt_no:
                attempt_no = attempt.attempt_no

        filename = Path(str(artifact.storage_key)).name or f"{artifact.artifact_type}.bin"
        location = storage_service.build_location(
            tenant_id=str(task.tenant_id),
            task_id=str(task.id),
            attempt_no=attempt_no,
            filename=filename,
        )

        artifact.storage_backend = storage_service.backend
        artifact.bucket = location.bucket
        artifact.object_key = location.object_key
        artifact.storage_key = location.object_key
        artifact.content_type = artifact.content_type or artifact.mime_type
        session.add(artifact)
        migrated += 1

    session.commit()
    return migrated


def main() -> int:
    session = SessionLocal()
    try:
        migrated = backfill_local_artifacts(session)
        print(f"migrated_artifacts={migrated}")
        return 0
    finally:
        session.close()


if __name__ == "__main__":
    raise SystemExit(main())
