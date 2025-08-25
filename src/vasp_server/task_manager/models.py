from sqlalchemy import Column, String, Integer, DateTime, Text
from sqlalchemy.sql import func
from sqlalchemy.types import JSON

from .database import Base


class Task(Base):
    __tablename__ = "tasks"

    id = Column(String, primary_key=True, index=True)
    user_id = Column(String, index=True, nullable=False)
    task_type = Column(String, nullable=False)
    status = Column(String, index=True, nullable=False, default="queued")
    progress = Column(Integer, nullable=False, default=0)
    params = Column(JSON, nullable=True)
    result_path = Column(String, nullable=True)
    result_data = Column(JSON, nullable=True)  # 存储分析结果数据
    external_job_id = Column(String, nullable=True)
    process_id = Column(Integer, nullable=True)  # 进程ID
    error_message = Column(Text, nullable=True)

    created_at = Column(DateTime(timezone=True), server_default=func.now(), nullable=False)
    updated_at = Column(
        DateTime(timezone=True), server_default=func.now(), onupdate=func.now(), nullable=False
    ) 