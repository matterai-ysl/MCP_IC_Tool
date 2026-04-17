# MCP IC Tool / VASP Control Plane

这个仓库是 MatterAI 的 VASP 计算与分析服务，当前主链路是：

```text
MCP 工具层 -> VPS 控制面/API -> HPC pull worker -> VASP/分析任务 -> VPS 产物与状态接口
```

如果你是接手修 bug 的开发者，建议先读：

- [CLAUDE.md](CLAUDE.md)：代码架构、任务状态机、主要文件和本地运行方式
- [docs/runbooks/session-handoff.md](docs/runbooks/session-handoff.md)：VPS/HPC 登录、tmux、健康检查、近期 bugfix 交接
- [docs/architecture/vps-hpc-control-plane.md](docs/architecture/vps-hpc-control-plane.md)：VPS 控制面和 HPC worker 架构说明
- [docs/runbooks/distributed-cutover.md](docs/runbooks/distributed-cutover.md)：分布式执行切换 runbook

## 当前架构摘要

- `src/mcp_ic_tool/`：MCP server 和 HTTP client，将智能体工具调用转发到 VASP 控制面。
- `src/vasp_server/vasp_server_api.py`：FastAPI 控制面，负责提交任务、查询状态、取消任务、worker 内部协议。
- `src/vasp_server/hpc_pull_worker.py`：HPC 侧 worker，主动领取任务、提交 Slurm/VASP、回传结果。
- `src/vasp_server/vasp_worker.py`：实际准备 VASP 输入、执行计算、分析输出。
- `src/vasp_server/task_manager/`：任务状态机、租约、重试、产物索引。
- `src/vasp_server/analyzers/`：DOS、Band、SCF、MD、NEB、Phonon 等分析器。

## 本地快速启动

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
cp .env.example .env
```

启动控制面：

```bash
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140
curl http://127.0.0.1:8140/health
```

启动 MCP 服务：

```bash
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 \
python single_port_server.py --host 0.0.0.0 --port 8130
```

本地开发一般不直接跑真实 VASP，可以优先用测试锁定行为。

## 近期已修复的重点

- MD 分析兼容 `numpy >= 2.0`：不再依赖已移除的 `np.trapz`，改为优先 `np.trapezoid`。
- CIF 解析增强：兼容缺少 `_atom_site_label` 的 CIF，并清理应为 0 的微小晶格/坐标噪声。
- VASP 失败诊断增强：任务状态、远程分析和 Agent 分析返回 `failure_type` 与 `suggested_action`。
- 数值类失败只反馈给智能体，不由后端自动重提：智能体根据 `failure_type` 决定是否重新提交和如何改参数。
- 失败证据链增强：如果任务错误信息过泛，但 execution attempt 有 `work_directory`，API 会回读 `result.log` / `OUTCAR` 做兜底诊断。

当前 `failure_type` 主要包括：

- `electronic_divergence`：如 `BRMIX`、`EDDDAV`、`ZHEGV`，建议先做稳定 SCF 或预优化。
- `ionic_line_search_failure`：如 `ZBRENT`，建议用 `CONTCAR` 续算并减小 `POTIM`。
- `lattice_inconsistency`：如 `HNFORM` 或 reciprocal lattice 不一致，建议规范化晶格或重新导出结构。
- `invalid_structure_file`：结构文件无效或 CIF 缺关键结构信息。
- `missing_charge_density`：后续 DOS/Band 类任务缺少稳定 `CHGCAR`。

## 推荐测试

接手 bugfix 时，优先跑和当前修改相关的定向测试：

```bash
./.venv/bin/python -m pytest \
  tests/test_base.py \
  tests/test_md_analyzer.py \
  tests/test_task_status_api.py \
  tests/test_hpc_pull_worker.py::test_scheduler_failure_message_includes_agent_guidance_for_unstable_scf \
  tests/test_hpc_pull_worker.py::test_scheduler_failure_message_includes_agent_guidance_for_lattice_inconsistency \
  -q
```

如果修改分布式 worker 或任务状态机，再补跑：

```bash
./.venv/bin/python -m pytest \
  tests/test_hpc_pull_worker.py \
  tests/test_internal_worker_api.py \
  tests/test_task_manager_claims.py \
  tests/e2e/test_distributed_flow.py \
  -q
```

## GitHub 上传前注意

不要提交真实密钥、数据库、计算产物或用户上传文件。尤其检查：

- `.env`
- `tasks.db`
- `public_artifacts/`
- `*.zip` / `*.tar.gz`
- HPC/VPS 登录凭据

`.env.example` 可以提交，但里面只能保留示例值。

## 旧文档提示

[VASP_API_README.md](VASP_API_README.md) 是早期 API 文档，其中部分 `formula` / Materials Project 相关描述已经不是当前主链路。当前以 `README.md`、`CLAUDE.md` 和 `docs/runbooks/session-handoff.md` 为准。
