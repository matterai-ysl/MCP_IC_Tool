from dataclasses import dataclass
from typing import Optional


@dataclass(frozen=True)
class ObjectStorageLocation:
    bucket: str
    object_key: str
    storage_backend: str = "oss"


@dataclass(frozen=True)
class RegisteredArtifact:
    task_id: str
    artifact_type: str
    owner_type: str
    owner_id: Optional[str]
    location: ObjectStorageLocation
    content_type: Optional[str] = None
    size_bytes: Optional[float] = None
    etag: Optional[str] = None
    sha256: Optional[str] = None
