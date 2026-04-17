from __future__ import annotations

from typing import Any, Dict, Optional


class APIError(Exception):
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

    def to_dict(self) -> Dict[str, Any]:
        payload: Dict[str, Any] = {
            "code": self.code,
            "message": self.message,
            "retryable": self.retryable,
            "details": self.details,
        }
        if self.suggested_action:
            payload["suggested_action"] = self.suggested_action
        return payload


def build_api_error(
    *,
    status_code: int,
    code: str,
    message: str,
    retryable: bool = False,
    details: Optional[Dict[str, Any]] = None,
    suggested_action: Optional[str] = None,
) -> APIError:
    return APIError(
        status_code=status_code,
        code=code,
        message=message,
        retryable=retryable,
        details=details,
        suggested_action=suggested_action,
    )
