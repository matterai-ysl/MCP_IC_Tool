import os
from typing import Optional
BASE_URL = "http://localhost:8130"
DOWNLOAD_URL = "/data/home/ysl9527/vasp_calculations"
MCP_PORT = 8130
VASP_SERVER_BASE_URL = "http://localhost:8135"
VASP_SERVER_BASE_File_URL = "https://47.99.180.80/vasp_server"
VASP_SERVER_PORT = 8135
from pathlib import Path


def get_download_url(path:str):
    return f"{VASP_SERVER_BASE_File_URL}/vasp/download/{Path(path).relative_to(DOWNLOAD_URL).as_posix()}"

def get_static_url(path:str):
    return f"{VASP_SERVER_BASE_File_URL}/vasp/static/{Path(path).relative_to(DOWNLOAD_URL).as_posix()}"

def get_env(name: str, default: Optional[str] = None) -> Optional[str]:
    return os.environ.get(name, default)


class VaspServerConfig:
    """VASP 服务端 API 配置。"""

    def __init__(self) -> None:
        # 允许通过环境变量覆盖，默认指向本机或内网可达地址
        self.base_url: str = VASP_SERVER_BASE_URL


class MCPServerConfig:
    """MCP 服务器配置。"""

    def __init__(self) -> None:
        self.host: str = get_env("MCP_HOST", "0.0.0.0") or "0.0.0.0"
        self.port: int = MCP_PORT


vasp_config = VaspServerConfig()
mcp_config = MCPServerConfig()



