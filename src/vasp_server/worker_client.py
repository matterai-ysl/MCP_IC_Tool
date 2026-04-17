from typing import Any, Optional

import requests


class ControlPlaneClient:
    def __init__(
        self,
        base_url: str,
        worker_token: str,
        worker_id: str,
        queue_name: str = "default",
        timeout: int = 30,
    ) -> None:
        self.base_url = base_url.rstrip("/")
        self.worker_token = worker_token
        self.worker_id = worker_id
        self.queue_name = queue_name
        self.timeout = timeout
        self.session = requests.Session()

    @property
    def _headers(self) -> dict[str, str]:
        return {"Authorization": f"Bearer {self.worker_token}"}

    def claim(self) -> Optional[dict[str, Any]]:
        response = self.session.post(
            f"{self.base_url}/internal/workers/claim",
            json={"worker_id": self.worker_id, "queue_name": self.queue_name},
            headers=self._headers,
            timeout=self.timeout,
        )
        if response.status_code == 204:
            return None
        response.raise_for_status()
        return response.json()

    def resume(self, task_id: str) -> dict[str, Any]:
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/resume",
            json={"worker_id": self.worker_id, "queue_name": self.queue_name},
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def mark_running(self, task_id: str, lease_token: str) -> dict[str, Any]:
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/running",
            json={"worker_id": self.worker_id, "lease_token": lease_token},
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def heartbeat(self, task_id: str, lease_token: str) -> dict[str, Any]:
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/heartbeat",
            json={"worker_id": self.worker_id, "lease_token": lease_token},
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def complete(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> dict[str, Any]:
        body = {"worker_id": self.worker_id, "lease_token": lease_token, **payload}
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/complete",
            json=body,
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def fail(self, task_id: str, lease_token: str, payload: dict[str, Any]) -> dict[str, Any]:
        body = {"worker_id": self.worker_id, "lease_token": lease_token, **payload}
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/fail",
            json=body,
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def cancel_ack(self, task_id: str, lease_token: str) -> dict[str, Any]:
        response = self.session.post(
            f"{self.base_url}/internal/tasks/{task_id}/cancel-ack",
            json={"worker_id": self.worker_id, "lease_token": lease_token},
            headers=self._headers,
            timeout=self.timeout,
        )
        response.raise_for_status()
        return response.json()

    def upload_public_artifact(
        self,
        task_id: str,
        filename: str,
        local_path: str,
        content_type: str | None = None,
    ) -> dict[str, Any]:
        with open(local_path, "rb") as fh:
            headers = dict(self._headers)
            if content_type:
                headers["Content-Type"] = content_type
            response = self.session.put(
                f"{self.base_url}/internal/tasks/{task_id}/artifacts/{filename.lstrip('/')}",
                data=fh.read(),
                headers=headers,
                timeout=self.timeout,
            )
        response.raise_for_status()
        return response.json()
