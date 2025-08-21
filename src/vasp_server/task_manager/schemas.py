from typing import Optional, Any, Dict, List
from pydantic import BaseModel, Field


class TaskCreate(BaseModel):
    task_type: str = Field(..., description="计算类型，例如 'vasp_scf'")
    params: Optional[Dict[str, Any]] = None


class TaskRead(BaseModel):
    id: str
    user_id: str
    task_type: str
    status: str
    progress: int
    result_path: Optional[str] = None
    external_job_id: Optional[str] = None
    error_message: Optional[str] = None

    class Config:
        from_attributes = True


class TaskStatus(BaseModel):
    id: str
    status: str
    progress: int
    error_message: Optional[str] = None


class TaskList(BaseModel):
    tasks: List[TaskRead] 