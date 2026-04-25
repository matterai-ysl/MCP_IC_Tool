# MCP IC Tool / VASP Control Plane

这个仓库提供一套面向 VASP 计算任务的控制面与 MCP 工具层，当前主链路是：

```text
MCP 工具层 -> FastAPI 控制面 -> Task Manager -> Pull Worker -> VASP / 分析任务
```

仓库里的核心角色：

- `single_port_server.py`：单端口 MCP 服务入口，对外暴露 `/mcp`
- `src/vasp_server/vasp_server_api.py`：FastAPI 控制面，负责提交任务、查询状态、取消任务、worker 内部接口
- `src/vasp_server/task_manager/`：任务状态机、租约、重试、产物登记
- `src/vasp_server/hpc_pull_worker.py`：拉取式 worker，主动 claim 任务并回传结果
- `src/vasp_server/vasp_worker.py`：实际准备输入、执行 VASP、分析输出

## Repository Status

这个仓库当前不是 `uv sync` 风格项目：

- 没有 `pyproject.toml`
- 没有 `uv.lock`
- 依赖入口仍然是根目录的 `requirements.txt`
- `.venv` 已经存在，并且是 `uv` 创建的虚拟环境

我在仓库本地核对到的真实环境信息如下：

- `uv`: `0.10.10`
- `.venv` Python: `3.12.13`
- `.venv/pyvenv.cfg` 显示解释器来自 `uv` 管理的 CPython 3.12

## Requirements

如果你要在这台机器上继续使用现有 `.venv`，当前环境里已经装有一部分基础包，例如：

- `fastapi`
- `httpx`
- `pydantic`
- `pydantic-settings`
- `sqlalchemy`
- `pytest`
- `pytest-asyncio`
- `requests`

但按仓库当前代码实际情况，现有 `.venv` 还不能直接完整启动服务，也不能通过全部相关测试。缺失项里至少包括：

- `uvicorn`
- `fastmcp`
- `psycopg2-binary`
- `prometheus-fastapi-instrumentator`
- `mp-api`
- `langchain-openai`
- `langchain-anthropic`
- `langgraph`
- `pymatgen`

说明：

- `single_port_server.py` 依赖 `uvicorn` 和 `fastmcp`
- `src/vasp_server/vasp_server_api.py` 在 PostgreSQL 场景下依赖 `psycopg2-binary`
- `src/vasp_server/base.py` 的 CIF -> POSCAR 转换依赖 `pymatgen`
- 当前 `requirements.txt` 包含大部分运行依赖，但没有显式列出 `pymatgen`

如果要把当前 `.venv` 补到可启动状态，建议使用 `uv` 继续安装：

```bash
uv venv .venv --python 3.12
source .venv/bin/activate
uv pip install -r requirements.txt
uv pip install pymatgen
```

如果你只是想确认当前虚拟环境里到底装了什么，可以运行：

```bash
uv pip list --python .venv/bin/python
```

## Quick Start

### 1. 准备环境

```bash
cd /Users/songlin/Desktop/Code/MCP_IC_Tool
uv venv .venv --python 3.12
source .venv/bin/activate
uv pip install -r requirements.txt
uv pip install pymatgen
cp .env.example .env
```

本仓库默认会从 `.env` 读取配置。请只填写你本地或部署环境所需的值，不要把真实密钥、路径或内部连接方式提交到仓库。

### 2. 启动 FastAPI 控制面

```bash
source .venv/bin/activate
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140
```

健康检查：

```bash
curl http://127.0.0.1:8140/health
```

### 3. 启动 MCP 服务

```bash
source .venv/bin/activate
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 \
python single_port_server.py --host 0.0.0.0 --port 8130
```

默认入口：

- `http://127.0.0.1:8130/mcp`

### 4. 启动 Pull Worker

本仓库当前没有单独的 worker CLI 包装脚本，开发态可以直接运行：

```bash
source .venv/bin/activate
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 \
INTERNAL_WORKER_TOKEN=worker-token \
python -c "from src.vasp_server.hpc_pull_worker import PullWorker; PullWorker(worker_id='local-dev', queue_name='default').run_loop()"
```

如果你只是做 API 或 MCP 联调，不一定需要真的跑完整 VASP 作业；很多改动可以先通过单元测试验证。

## Development Notes

推荐先看这些文件理解代码结构：

- `CLAUDE.md`
- `src/vasp_server/task_manager/README.md`
- `VASP_API_README.md`

当前常见开发流程：

```bash
source .venv/bin/activate
python -m pytest tests -q
```

如果出现与结构解析相关的失败，先确认 `pymatgen` 是否已经安装。

## Safety

以下内容不要提交到 GitHub：

- `.env`
- `.venv/`
- `env/`
- `venv/`
- `tasks.db`
- `public_artifacts/`
- 任何真实密钥、令牌、证书或数据库连接信息
- 任何超算、集群、VPS 的登录方式、跳板方式或运维凭据

可以提交：

- `.env.example`，但只能保留示例值
- 代码、测试、公开文档
- 不包含敏感信息的启动说明

## Known Gap

仓库当前存在一个“文档与环境不完全一致”的现实情况：

- `.venv` 是 `uv` 创建的，但依赖并未完全补齐
- `requirements.txt` 不是锁文件，不能代表“当前机器已安装集合”
- CIF 相关能力实际依赖 `pymatgen`，这一点需要在使用时额外注意

如果后续要把环境管理收敛得更稳定，建议补齐 `pyproject.toml` 与 `uv.lock`，再把运行依赖和可选依赖拆开管理。
