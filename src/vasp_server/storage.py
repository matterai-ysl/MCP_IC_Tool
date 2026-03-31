import hashlib
import hmac
from urllib.parse import quote, urlencode

from .settings import settings
from .storage_models import ObjectStorageLocation


class ObjectStorageService:
    def __init__(
        self,
        endpoint: str,
        bucket: str,
        region: str,
        access_key_id: str,
        access_key_secret: str,
        backend: str = "oss",
    ) -> None:
        self.endpoint = endpoint.rstrip("/")
        self.bucket = bucket
        self.region = region
        self.access_key_id = access_key_id
        self.access_key_secret = access_key_secret
        self.backend = backend

    @classmethod
    def from_settings(cls) -> "ObjectStorageService":
        endpoint = settings.s3_endpoint or settings.oss_endpoint
        return cls(
            endpoint=endpoint,
            bucket=settings.oss_bucket,
            region=settings.oss_region,
            access_key_id=settings.oss_access_key_id,
            access_key_secret=settings.oss_access_key_secret,
            backend=settings.artifact_storage_backend,
        )

    def build_task_prefix(self, tenant_id: str, task_id: str, attempt_no: int) -> str:
        return f"tenant/{tenant_id}/tasks/{task_id}/attempts/{attempt_no}"

    def build_location(self, tenant_id: str, task_id: str, attempt_no: int, filename: str) -> ObjectStorageLocation:
        prefix = self.build_task_prefix(tenant_id=tenant_id, task_id=task_id, attempt_no=attempt_no)
        clean_filename = filename.lstrip("/")
        return ObjectStorageLocation(
            bucket=self.bucket,
            object_key=f"{prefix}/{clean_filename}",
            storage_backend=self.backend,
        )

    def create_upload_url(self, object_key: str, content_type: str) -> str:
        return self._signed_url(
            method="PUT",
            object_key=object_key,
            content_type=content_type,
            expires_in=settings.artifact_url_expire_seconds,
        )

    def create_download_url(self, object_key: str, expires_in: int = 3600) -> str:
        return self._signed_url(
            method="GET",
            object_key=object_key,
            content_type=None,
            expires_in=expires_in,
        )

    def _signed_url(
        self,
        method: str,
        object_key: str,
        content_type: str | None,
        expires_in: int,
    ) -> str:
        payload = "\n".join(
            [
                method,
                self.bucket,
                object_key,
                str(expires_in),
                content_type or "",
            ]
        )
        signature = hmac.new(
            self.access_key_secret.encode("utf-8"),
            payload.encode("utf-8"),
            hashlib.sha256,
        ).hexdigest()
        query = urlencode(
            {
                "access_key": self.access_key_id,
                "expires": expires_in,
                "signature": signature,
                "method": method,
                "content_type": content_type,
                "region": self.region,
            }
        )
        quoted_key = quote(object_key)
        return f"{self.endpoint}/{self.bucket}/{quoted_key}?{query}"
