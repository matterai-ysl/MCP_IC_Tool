# CLAUDE.md

This file provides guidance to Claude Code when working in this repository.

## Current Architecture

这个仓库当前的主链路已经从“API 进程内直接起线程执行 VASP”切到三层结构：

1. `src/mcp_ic_tool/`
   - FastMCP 工具层
   - 负责把 LLM 的工具调用转换成对 VASP 控制面的 HTTP 请求
   - 默认通过 [single_port_server.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/single_port_server.py) 暴露 `/mcp`

2. `src/vasp_server/vasp_server_api.py`
   - VPS 控制面 / 公网 API
   - 负责校验请求、写入 `Task` 记录、查询状态、取消任务
   - 同时挂载内部 worker API：`/internal/workers/*`、`/internal/tasks/*`

3. `src/vasp_server/hpc_pull_worker.py`
   - HPC 侧拉取式 worker
   - 主动向控制面 claim 任务、续租、执行 VASP、回传结果和产物清单

4. `src/vasp_server/task_manager/`
   - 数据库编排层
   - 负责任务状态机、租约、取消、重试、回收过期租约、产物索引

5. `src/vasp_server/storage.py`
   - 对象存储 URL / object key 约定层
   - 负责生成对象 key 和签名样式 URL

## Real End-To-End Flow

当前真实链路如下：

1. LLM 调用 MCP 工具，例如 `submit_structure_optimization`
2. [src/mcp_ic_tool/mcp_server.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/mcp_server.py) 组装 payload
3. [src/mcp_ic_tool/client.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/client.py) 把请求发到 `VASP_SERVER_BASE_URL`
4. [src/vasp_server/vasp_server_api.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py) 校验参数后调用 `TaskManager.submit_task()`
5. [src/vasp_server/task_manager/manager.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py) 只写数据库，把任务置为 `queued`
6. HPC 侧 [src/vasp_server/hpc_pull_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py) 通过 [src/vasp_server/worker_client.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/worker_client.py) 调 `/internal/workers/claim`
7. 控制面给 worker 分配 `lease_token`，任务进入 `leased`
8. worker 标记 `running`，调用 [src/vasp_server/vasp_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py) 执行实际计算
9. 执行过程中 worker 通过 heartbeat 续租；如果用户取消，控制面会把状态改成 `cancel_requested`，worker 在 heartbeat 响应里感知并调用 `scancel`
10. worker 完成后把 `result_data` 和 `artifact_manifest` 回传给控制面
11. 控制面把产物登记到 `Artifact` 表，并在状态接口里把产物转成签名 URL
12. 用户或 MCP 再查询 `/vasp/task/{task_id}` 时，看到的是数据库状态和产物 URL，而不是本地文件路径

当前输入契约里已经不再支持 `formula` 或 Materials Project 自动找结构。
结构来源统一是：

- 外部提供的 `cif_url`
- 或已完成上游任务的 `*_task_id`

## Canonical Task States

当前任务主状态机定义在 [src/vasp_server/task_manager/models.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py)：

- `queued`: 已入队，等待 worker claim
- `leased`: 已被某个 worker 领取，租约有效
- `running`: HPC 作业已提交并运行
- `uploading`: 预留给上传阶段
- `analyzing`: 预留给分析阶段
- `completed`: 成功完成
- `failed`: 失败且不再重试
- `cancel_requested`: 已请求取消，等待 worker 确认
- `canceled`: worker 已确认取消

当前分布式主路径是：

`queued -> leased -> running -> completed`

取消路径是：

`queued|leased|running -> cancel_requested -> canceled`

恢复路径是：

- `leased -> queued`：租约过期回收
- `running -> queued`：仅在可重试基础设施错误下允许

## Upstream Artifact Reuse

SCF / DOS 不再优先依赖共享目录，而是优先走 `upstream_artifact_manifest`：

- 解析逻辑在 [src/vasp_server/input_resolver.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/input_resolver.py)
- worker claim 响应会带上 `upstream_artifact_manifest`
- `VaspWorker` 会把上游 `CONTCAR` / `POSCAR` / `CHGCAR` / `WAVECAR` 等物化到当前工作目录

这意味着控制面是任务编排真源，上游输入不应再靠人工拼本地路径。

## Object Storage Contract

当前对外契约已经切换到“对象存储元数据 + 签名 URL”：

- `Artifact.storage_backend` / `bucket` / `object_key` 由 [src/vasp_server/storage.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/storage.py) 生成
- `/vasp/task/{task_id}` 会优先把产物转成下载 URL 返回
- 默认禁用旧的本地 `/static/` / `/download/file/` URL 暴露

需要注意的现实情况：

- 目前代码已经完成“对象 key 建模、URL 生成、状态接口出参切换”
- 但 `PullWorker` 现在回传的是本地产物 `artifact_manifest`
- 控制面会根据 manifest 生成对象存储元数据并写库
- 仓库里还没有看到真正把二进制文件 `PUT` 到 OSS/S3 的通用上传实现

也就是说，当前代码把“对外接口契约”切到了对象存储模型，但如果要在生产中真正下载这些对象，仍需要保证文件已经被上传，或者补上实际上传逻辑。

## Key Files

- MCP 入口： [single_port_server.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/single_port_server.py)
- MCP 工具定义： [src/mcp_ic_tool/mcp_server.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/mcp_server.py)
- MCP HTTP client： [src/mcp_ic_tool/client.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/client.py)
- 控制面 API： [src/vasp_server/vasp_server_api.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py)
- 内部 worker API： [src/vasp_server/internal_worker_api.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/internal_worker_api.py)
- 任务编排： [src/vasp_server/task_manager/manager.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py)
- 状态机模型： [src/vasp_server/task_manager/models.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py)
- HPC pull worker： [src/vasp_server/hpc_pull_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py)
- 实际执行器： [src/vasp_server/vasp_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py)
- 上游输入解析： [src/vasp_server/input_resolver.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/input_resolver.py)
- 对象存储 URL： [src/vasp_server/storage.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/storage.py)
- 切换 runbook： [docs/runbooks/distributed-cutover.md](/Users/ysl/Desktop/Code/MCP_IC_Tool/docs/runbooks/distributed-cutover.md)

## How To Run

### 1. Prepare Environment

```bash
cd /Users/ysl/Desktop/Code/MCP_IC_Tool
cp .env.example .env
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

最重要的环境变量：

- `DATABASE_URL` 或 `DB_HOST/DB_PORT/DB_NAME/DB_USER/DB_PASSWORD`
- `VASP_SERVER_BASE_URL`
- `INTERNAL_WORKER_TOKEN`
- `OSS_ENDPOINT` / `OSS_BUCKET` / `OSS_ACCESS_KEY_ID` / `OSS_ACCESS_KEY_SECRET`
- `VASP_WORK_ROOT`
- `VASP_PATH` / `PSEUDO_PATH`
- `DEFAULT_USER_ID`

本地联调时，建议至少改成：

```bash
VASP_SERVER_BASE_URL=http://127.0.0.1:8140
INTERNAL_WORKER_TOKEN=worker-token
LEGACY_LOCAL_ARTIFACT_URLS_ENABLED=false
```

如果你不想连 PostgreSQL，可以直接用 `.env.example` 里的 `SQLITE_FALLBACK_URL`。

### 2. Start The VASP Control Plane

```bash
cd /Users/ysl/Desktop/Code/MCP_IC_Tool
source .venv/bin/activate
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140
```

作用：

- 提供公网任务提交接口
- 提供任务状态 / 取消接口
- 提供内部 worker claim / heartbeat / complete / fail 接口

健康检查：

```bash
curl http://127.0.0.1:8140/health
```

### 3. Start The MCP Server

```bash
cd /Users/ysl/Desktop/Code/MCP_IC_Tool
source .venv/bin/activate
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 \
python single_port_server.py --host 0.0.0.0 --port 8130
```

MCP 入口地址是：

- `http://127.0.0.1:8130/mcp`

说明：

- 这个进程本身不执行 VASP
- 它只是把 MCP 工具调用转发到控制面
- `single_port_server.py` 里仍保留了一些旧的文件服务路由，但它们不是当前分布式 artifact 主链路的一部分

### 4. Start The HPC Pull Worker

当前仓库里还没有独立的 worker CLI 入口，但可以直接这样跑：

```bash
cd /Users/ysl/Desktop/Code/MCP_IC_Tool
source .venv/bin/activate
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 \
INTERNAL_WORKER_TOKEN=worker-token \
python -c "from src.vasp_server.hpc_pull_worker import PullWorker; PullWorker(worker_id='hpc-dev', queue_name='default').run_loop()"
```

如果是在真实 HPC 节点上运行，还需要把以下变量配置正确：

- `VASP_WORK_ROOT`
- `VASP_PATH`
- `PSEUDO_PATH`
- `SLURM_*`
- `ONEAPI_ENV_SCRIPT`
- `VASP_SLURM_RUN_LINE`

### 5. Submit And Query A Task

不经过 MCP，直接调用控制面可以这样验证：

```bash
curl -X POST http://127.0.0.1:8140/vasp/structure-optimization \
  -H "Content-Type: application/json" \
  -d '{"user_id":"test_user","cif_url":"https://structures.example.com/Li2O.cif"}'
```

查询状态：

```bash
curl "http://127.0.0.1:8140/vasp/task/<task_id>?user_id=test_user"
```

取消任务：

```bash
curl -X POST "http://127.0.0.1:8140/vasp/task/<task_id>/cancel?user_id=test_user"
```

### 6. Backfill And Recovery Scripts

历史产物迁移：

```bash
python scripts/backfill_artifacts.py
```

回收过期租约 / 心跳超时任务：

```bash
python scripts/requeue_stuck_tasks.py
```

### 7. Run Tests

```bash
cd /Users/ysl/Desktop/Code/MCP_IC_Tool
source .venv/bin/activate
./.venv/bin/python -m pytest -q tests
```

当前这套升级链路的重点测试包括：

- [tests/test_public_api_upgrade.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_public_api_upgrade.py)
- [tests/test_backfill_scripts.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_backfill_scripts.py)
- [tests/e2e/test_distributed_flow.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/e2e/test_distributed_flow.py)

## Development Notes

- 公网 API 现在应该只负责“校验 + 入库”，不要再把新功能接回进程内直接执行
- 结构输入只接受外部 URL 或上游任务 ID，不要重新引入 `formula` / Materials Project 下载路径
- HPC worker 不应该直接读写 PostgreSQL，只能走内部 API
- 新产物接口应返回签名 URL，不要新增本地文件路径出参
- SCF / DOS 上游输入应优先走 `upstream_artifact_manifest`
- `agent_analyze` 这类旧逻辑仍依赖控制面本地工作目录，和新的对象存储主链路不是一套模型，改这部分时要格外小心

## Adding A New Calculation Type

新增计算类型时，至少同步修改这些位置：

1. 在 [src/vasp_server/schemas.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py) 增加请求/响应模型
2. 在 [src/vasp_server/vasp_server_api.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py) 增加公网提交接口
3. 在 [src/vasp_server/vasp_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py) 增加执行方法
4. 在 [src/vasp_server/hpc_pull_worker.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py) 的 `TASK_TYPE_TO_METHOD` 注册映射
5. 在 [src/mcp_ic_tool/models.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/models.py) 和 [src/mcp_ic_tool/mcp_server.py](/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/mcp_server.py) 增加 MCP 工具
6. 补充对应单测和必要的端到端测试
