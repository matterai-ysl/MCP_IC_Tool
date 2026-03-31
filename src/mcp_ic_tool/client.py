from typing import Any, Dict, List, Optional
import logging
import httpx

from .config import vasp_config

logger = logging.getLogger(__name__)


class VaspAPIClient:
    """轻量级 HTTP 客户端，封装 vasp_server_api.py 的端点。"""

    def __init__(self, base_url: Optional[str] = None) -> None:
        self.base_url = base_url or vasp_config.base_url
        logger.debug("VaspAPIClient base_url: %s", self.base_url)

    async def _apost(self, path: str, json: Dict[str, Any]) -> Dict[str, Any]:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.post(f"{self.base_url}{path}", json=json)
            resp.raise_for_status()
            return resp.json()

    async def _aget(self, path: str, params: Dict[str, Any]) -> Dict[str, Any]:
        async with httpx.AsyncClient(timeout=60.0) as client:
            resp = await client.get(f"{self.base_url}{path}", params=params)
            resp.raise_for_status()
            return resp.json()

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
            resp.raise_for_status()
            return resp.json()

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
