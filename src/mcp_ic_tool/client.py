from typing import Any, Dict, List, Optional
import logging
import httpx

from .config import vasp_config

logger = logging.getLogger(__name__)


class VaspAPIError(Exception):
    def __init__(
        self,
        *,
        status_code: int,
        code: str,
        message: str,
        retryable: bool = False,
        details: Optional[Dict[str, Any]] = None,
        suggested_action: Optional[str] = None,
    ) -> None:
        super().__init__(message)
        self.status_code = status_code
        self.code = code
        self.message = message
        self.retryable = retryable
        self.details = details or {}
        self.suggested_action = suggested_action

    def __str__(self) -> str:
        parts = [f"[{self.code}] {self.message}"]
        if self.details:
            detail_parts = []
            for key, value in self.details.items():
                if value in (None, "", [], {}):
                    continue
                detail_parts.append(f"{key}={value}")
            if detail_parts:
                parts.append("details: " + ", ".join(detail_parts))
        parts.append(f"status={self.status_code}")
        parts.append(f"retryable={self.retryable}")
        if self.suggested_action:
            parts.append(f"建议: {self.suggested_action}")
        return " | ".join(parts)


class VaspAPIClient:
    """轻量级 HTTP 客户端，封装 vasp_server_api.py 的端点。"""

    def __init__(self, base_url: Optional[str] = None) -> None:
        self.base_url = base_url or vasp_config.base_url
        logger.debug("VaspAPIClient base_url: %s", self.base_url)

    def _handle_response(self, resp: httpx.Response) -> Dict[str, Any]:
        if resp.is_success:
            if not resp.content:
                return {}
            return resp.json()

        payload: Dict[str, Any] = {}
        try:
            payload = resp.json()
        except Exception:
            payload = {}

        error_payload = payload.get("error") if isinstance(payload, dict) else None
        if isinstance(error_payload, dict):
            raise VaspAPIError(
                status_code=resp.status_code,
                code=str(error_payload.get("code") or f"HTTP_{resp.status_code}"),
                message=str(error_payload.get("message") or "请求失败"),
                retryable=bool(error_payload.get("retryable", False)),
                details=error_payload.get("details") if isinstance(error_payload.get("details"), dict) else {},
                suggested_action=error_payload.get("suggested_action"),
            )

        detail = None
        if isinstance(payload, dict):
            detail = payload.get("detail")
        if isinstance(detail, dict):
            raise VaspAPIError(
                status_code=resp.status_code,
                code=str(detail.get("code") or f"HTTP_{resp.status_code}"),
                message=str(detail.get("message") or detail.get("detail") or "请求失败"),
                retryable=bool(detail.get("retryable", False)),
                details=detail.get("details") if isinstance(detail.get("details"), dict) else {},
                suggested_action=detail.get("suggested_action"),
            )

        raise VaspAPIError(
            status_code=resp.status_code,
            code=f"HTTP_{resp.status_code}",
            message=str(detail or resp.text or "请求失败"),
            retryable=resp.status_code in {408, 429, 500, 502, 503, 504},
        )

    async def _apost(self, path: str, json: Dict[str, Any]) -> Dict[str, Any]:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.post(f"{self.base_url}{path}", json=json)
            return self._handle_response(resp)

    async def _aget(self, path: str, params: Dict[str, Any]) -> Dict[str, Any]:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.get(f"{self.base_url}{path}", params=params)
            return self._handle_response(resp)

    # --- 提交任务 ---
    async def submit_structure_optimization(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/structure-optimization", payload)

    async def submit_scf(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/scf-calculation", payload)

    async def submit_dos(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/dos-calculation", payload)

    async def submit_band_structure(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/band-structure", payload)

    async def submit_md(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/md-calculation", payload)

    # --- 查询/控制 ---
    async def get_task_status(self, task_id: str, user_id: str) -> Dict[str, Any]:
        return await self._aget(f"/vasp/task/{task_id}", {"user_id": user_id})

    async def cancel_task(self, task_id: str, user_id: str) -> Dict[str, Any]:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.post(
                f"{self.base_url}/vasp/task/{task_id}/cancel",
                params={"user_id": user_id},
            )
            return self._handle_response(resp)

    async def list_tasks(self, user_id: str) -> List[Dict[str, Any]]:
        return await self._aget("/vasp/tasks", {"user_id": user_id})

    # --- 独立分析 ---
    async def analyze_optimization(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/optimization", payload)

    async def analyze_scf(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/scf", payload)

    async def analyze_dos(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/dos", payload)

    async def analyze_band_structure(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/band-structure", payload)

    async def analyze_md(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/md", payload)

    async def analyze_md_multi(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/md-multi", payload)

    async def submit_neb(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/neb-calculation", payload)

    async def submit_phonon(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/phonon-calculation", payload)

    async def analyze_neb(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/neb", payload)

    async def analyze_phonon(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/analyze/phonon", payload)

    async def submit_custom_calculation(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/custom-calculation", payload)

    async def agent_analyze(self, payload: Dict[str, Any]) -> Dict[str, Any]:
        return await self._apost("/vasp/agent/analyze", payload)
