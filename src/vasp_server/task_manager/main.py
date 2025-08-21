from typing import Annotated, Optional

from fastapi import Depends, FastAPI, Header, HTTPException, status
from fastapi.middleware.cors import CORSMiddleware

from .database import init_db
from .manager import TaskManager
from .schemas import TaskCreate, TaskRead, TaskStatus, TaskList

app = FastAPI(title="MCP Task Manager", version="0.1.0")

app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

init_db()
manager = TaskManager()


async def get_current_user_id(x_user_id: Annotated[Optional[str], Header(alias="X-User-Id")]):
    if not x_user_id:
        raise HTTPException(status_code=status.HTTP_401_UNAUTHORIZED, detail="Missing X-User-Id header")
    return x_user_id


@app.post("/tasks", response_model=TaskStatus, status_code=status.HTTP_202_ACCEPTED)
async def create_task(payload: TaskCreate, user_id: str = Depends(get_current_user_id)):
    task_id = manager.submit_task(user_id=user_id, task_type=payload.task_type, params=payload.params)
    task = manager.get_task(task_id, user_id)
    if task is None:
        raise HTTPException(status_code=status.HTTP_500_INTERNAL_SERVER_ERROR, detail="Task creation failed")
    return TaskStatus(id=task.id, status=task.status, progress=task.progress, error_message=task.error_message)


@app.get("/tasks/{task_id}", response_model=TaskRead)
async def get_task(task_id: str, user_id: str = Depends(get_current_user_id)):
    task = manager.get_task(task_id, user_id)
    if task is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Task not found")
    return TaskRead.model_validate(task)


@app.get("/me/tasks", response_model=TaskList)
async def list_my_tasks(user_id: str = Depends(get_current_user_id)):
    tasks = manager.list_tasks(user_id)
    return TaskList(tasks=[TaskRead.model_validate(t) for t in tasks])


@app.post("/tasks/{task_id}/cancel", response_model=TaskStatus)
async def cancel_task(task_id: str, user_id: str = Depends(get_current_user_id)):
    ok = manager.cancel_task(task_id, user_id)
    if not ok:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Task not found or not owned by user")
    task = manager.get_task(task_id, user_id)
    if task is None:
        raise HTTPException(status_code=status.HTTP_404_NOT_FOUND, detail="Task not found")
    return TaskStatus(id=task.id, status=task.status, progress=task.progress, error_message=task.error_message) 