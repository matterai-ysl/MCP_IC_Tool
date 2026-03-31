import logging
import threading
from typing import Any, Dict, Optional

import requests

from .settings import VaspServerSettings

logger = logging.getLogger(__name__)

TASK_TYPE_TO_DISPLAY_NAME = {
    "structure_optimization": "结构优化计算",
    "scf_calculation": "SCF 计算",
    "dos_calculation": "DOS 计算",
    "band_structure": "能带计算",
    "md_calculation": "分子动力学计算",
    "neb_calculation": "NEB 计算",
    "phonon_calculation": "声子计算",
    "custom_calculation": "自定义计算",
}


def build_task_notification_payload(
    *,
    user_id: str,
    task_id: str,
    task_type: str,
    notification_type: str,
    status: str,
    language: str = "zh",
    execution_time_seconds: Optional[float] = None,
    result_summary: Optional[str] = None,
    error_message: Optional[str] = None,
    to_email: Optional[str] = None,
) -> Dict[str, Any]:
    context: Dict[str, Any] = {
        "tool_name": "preset-vasp-calculation",
        "tool_display_name": TASK_TYPE_TO_DISPLAY_NAME.get(task_type, task_type),
        "status": status,
        "task_id": task_id,
        "session_id": task_id,
    }
    if execution_time_seconds is not None:
        context["execution_time_seconds"] = execution_time_seconds
    if result_summary:
        context["result_summary"] = result_summary
    if error_message:
        context["error_message"] = error_message

    payload: Dict[str, Any] = {
        "notification_type": notification_type,
        "channel": "email",
        "language": language,
        "context": context,
    }
    if to_email:
        payload["to_email"] = to_email
    else:
        payload["user_id"] = user_id
    return payload


def send_notification(payload: Dict[str, Any]) -> None:
    settings = VaspServerSettings()
    if not settings.notification_service_base_url or not settings.notification_service_api_key:
        logger.debug("Notification service is not configured; skipping notification.")
        return

    url = settings.notification_service_base_url.rstrip("/") + "/api/v1/notify"
    headers = {"X-API-Key": settings.notification_service_api_key}

    try:
        response = requests.post(
            url,
            json=payload,
            headers=headers,
            timeout=settings.notification_timeout_seconds,
        )
        response.raise_for_status()
    except Exception as exc:
        logger.warning("Notification delivery failed: %s", exc)


def send_notification_async(payload: Dict[str, Any]) -> None:
    thread = threading.Thread(target=send_notification, args=(payload,), daemon=True)
    thread.start()
