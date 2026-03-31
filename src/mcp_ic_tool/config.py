import os
from pathlib import Path

BASE_URL = os.getenv("MCP_BASE_URL", "http://localhost:8130")
DOWNLOAD_URL = os.getenv("VASP_WORK_ROOT", "/data/home/ysl9527/vasp_calculations")
MCP_PORT = int(os.getenv("MCP_PORT", "8130"))
VASP_SERVER_BASE_URL = os.getenv("VASP_SERVER_BASE_URL", "https://api.matterai.tech")
VASP_SERVER_BASE_File_URL = os.getenv("VASP_PUBLIC_BASE_URL", "https://api.matterai.tech")


def get_download_url(path: str) -> str:
    p = Path(path)
    try:
        rel = p.relative_to(DOWNLOAD_URL).as_posix()
    except ValueError:
        rel = p.as_posix().lstrip("/")
    return f"{VASP_SERVER_BASE_File_URL}/vasp/download/{rel}"


def get_static_url(path: str) -> str:
    p = Path(path)
    try:
        rel = p.relative_to(DOWNLOAD_URL).as_posix()
    except ValueError:
        rel = p.as_posix().lstrip("/")
    return f"{VASP_SERVER_BASE_File_URL}/vasp/static/{rel}"


class VaspServerConfig:
    """VASP 服务端 API 配置。"""

    def __init__(self) -> None:
        self.base_url: str = VASP_SERVER_BASE_URL


class MCPServerConfig:
    """MCP 服务器配置。"""

    def __init__(self) -> None:
        self.host: str = os.getenv("MCP_HOST", "0.0.0.0") or "0.0.0.0"
        self.port: int = MCP_PORT


vasp_config = VaspServerConfig()
mcp_config = MCPServerConfig()
