# MCP Task Manager

一个简单的基于 FastAPI + SQLite 的任务管理服务。提交任务会返回 `task_id`，用户可通过 `task_id` 查询进度并进行取消。

## 快速开始

1. 安装依赖：
```bash
pip install -r requirements.txt
```

2. 启动服务：
```bash
uvicorn task_manager.main:app --reload
```

3. 使用时需要在请求头中携带 `X-User-Id` 标识当前用户。

## API

- 创建任务
```bash
curl -X POST "http://127.0.0.1:8000/tasks" \
  -H "Content-Type: application/json" \
  -H "X-User-Id: user_a" \
  -d '{"task_type":"vasp_scf","params":{"ENCUT":520}}'
```

- 查询任务
```bash
curl -X GET "http://127.0.0.1:8000/tasks/<task_id>" -H "X-User-Id: user_a"
```

- 列出我的任务
```bash
curl -X GET "http://127.0.0.1:8000/me/tasks" -H "X-User-Id: user_a"
```

- 取消任务
```bash
curl -X POST "http://127.0.0.1:8000/tasks/<task_id>/cancel" -H "X-User-Id: user_a"
```

## 说明
- 当前实现为本地线程模拟长任务（约20秒），并持续更新进度。
- 生产环境可替换为调度系统提交（如 LSF/Slurm），并将 `external_job_id` 写入数据库，以轮询或回调更新 `status/progress`。 