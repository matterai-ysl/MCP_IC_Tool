"""
Centralized configuration for the VASP server using Pydantic BaseSettings.

All environment variables and defaults are managed here. Other modules
(Config.py, database.py) import from this module instead of calling
os.getenv() directly.
"""

from pydantic_settings import BaseSettings
from typing import Optional


class VaspServerSettings(BaseSettings):
    # ---- Database ----
    database_url: Optional[str] = None
    db_host: str = "localhost"
    db_port: int = 5432
    db_name: str = "postgres"
    db_user: str = "postgres"
    db_password: Optional[str] = None
    sqlite_fallback_url: str = "sqlite:////tmp/mcp_ic_tool_tasks.db"

    # ---- VASP paths ----
    vasp_path: str = "/data/app/vasp/6.3.2-intel/bin/vasp_std"
    pseudo_path: str = "/data/home/ysl9527/software/psudopotential"
    u_values_json: str = "/data/home/ysl9527/software/u_values.json"
    bader_path: str = "/data/home/ysl9527/software/bader"
    chgsum_pl_path: str = "/data/home/ysl9527/software/chgsum.pl"
    vasp_work_root: str = "/data/home/ysl9527/vasp_calculations"

    # ---- Server ----
    vasp_server_base_url: str = "http://localhost:8130"
    vasp_remote_run_url: str = "http://localhost:8140"
    vasp_remote_run_port: int = 8140
    vasp_public_base_url: str = "https://api.matterai.tech"
    mcp_port: int = 8130
    internal_worker_token: str = "worker-token"

    # ---- Object storage ----
    artifact_storage_backend: str = "oss"
    oss_endpoint: str = "https://oss-cn-hangzhou.aliyuncs.com"
    s3_endpoint: Optional[str] = None
    oss_bucket: str = "vasp-artifacts"
    oss_region: str = "cn-hangzhou"
    oss_access_key_id: str = "test-access-key"
    oss_access_key_secret: str = "test-secret"
    artifact_url_expire_seconds: int = 3600
    legacy_local_artifact_urls_enabled: bool = False

    # ---- Task management ----
    max_concurrent_tasks: int = 4
    default_user_id: str = "123"
    task_lease_seconds: int = 300
    worker_poll_interval_seconds: int = 10
    worker_heartbeat_timeout_seconds: int = 120

    # ---- Notification service ----
    notification_service_base_url: Optional[str] = None
    notification_service_api_key: Optional[str] = None
    notification_timeout_seconds: float = 5.0
    notification_language: str = "zh"

    # ---- Quota ----
    max_user_concurrent_tasks: int = 3  # 单用户最大并发任务数
    max_user_total_tasks_per_day: int = 50  # 单用户每日最大任务提交数

    model_config = {
        "env_file": ".env",
        "env_file_encoding": "utf-8",
        "extra": "ignore",
    }

    @property
    def effective_database_url(self) -> str:
        """Resolve the database URL with the same priority as the old logic:
        1. DATABASE_URL env / field
        2. Build from DB_HOST/DB_PORT/DB_NAME/DB_USER/DB_PASSWORD
        3. SQLite fallback
        """
        if self.database_url:
            return self.database_url
        if self.db_password:
            return (
                f"postgresql://{self.db_user}:{self.db_password}"
                f"@{self.db_host}:{self.db_port}/{self.db_name}"
            )
        return self.sqlite_fallback_url


settings = VaspServerSettings()
