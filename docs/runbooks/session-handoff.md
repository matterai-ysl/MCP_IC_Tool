# VPS / HPC 会话交接手册

## 目的

这份文档用于在新窗口里快速接手当前这套 `VPS 控制面 + HPC pull worker` 环境，重点包含：

- 如何登录 VPS / HPC
- 各个 `tmux` 会话分别负责什么
- 常用重启与健康检查命令
- 当前这套环境里最容易踩的坑

## 当前架构

- `VPS`
  - `8130`: MCP 入口服务
  - `8140`: 控制面 / API / 内部 worker 协议
  - `nginx`:
    - `/hcp/` 反代到 `8140`
    - `/dft/` 对外暴露静态报告与轻量产物
- `HPC`
  - pull worker 从 VPS 领取任务
  - 真正的 VASP 计算与分析在 HPC 本地目录执行
  - 轻量产物回传到 VPS，对外统一走 `/dft/...`

## 登录方式

### VPS

```bash
ssh matterai@47.99.180.80
```

> GitHub 交接注意：不要在仓库中提交真实密码或私钥。VPS 密码/密钥请通过安全渠道单独交接。



登录后建议：

```bash
export TERM=xterm
cd ~/MCP_IC_Tool
```

### HPC

先登录登录节点：

```bash
ssh -i ~/.ssh/chinahpc_ysl9527 -p 1014 ysl9527@hpc3.chinahpc.com
```

再切到计算节点：

```bash
ssh c2ln6
cd ~/mcp_ic_tool_api
```

## 关键 tmux 会话

### VPS

查看：

```bash
tmux ls
```

常用会话：

- `ic_mcp`
  - `8130` MCP 服务
- `hpc_worker`
  - `8140` 控制面服务

进入：

```bash
tmux attach -t ic_mcp
tmux attach -t hpc_worker
```

### HPC

查看：

```bash
tmux ls
```

常用会话：

- `mcp_ic_api`
  - 主 worker

进入：

```bash
tmux attach -t mcp_ic_api
```

说明：

- 正常情况下只保留一个主 worker。
- 只有在排队较多、需要临时并发验证时，才短时额外拉第二个 worker。

## 常用启动 / 重启命令

### VPS: 8140 控制面

```bash
cd ~/MCP_IC_Tool
source .venv/bin/activate
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140
```

如果放进 `tmux`：

```bash
tmux kill-session -t hpc_worker 2>/dev/null || true
tmux new-session -d -s hpc_worker 'cd ~/MCP_IC_Tool && source .venv/bin/activate && python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140'
```

### VPS: 8130 MCP 服务

```bash
cd ~/MCP_IC_Tool
source .venv/bin/activate
VASP_SERVER_BASE_URL=http://127.0.0.1:8140 python single_port_server.py
```

如果放进 `tmux`：

```bash
tmux kill-session -t ic_mcp 2>/dev/null || true
tmux new-session -d -s ic_mcp 'cd ~/MCP_IC_Tool && source .venv/bin/activate && VASP_SERVER_BASE_URL=http://127.0.0.1:8140 python single_port_server.py'
```

### HPC: 主 worker

当前环境统一使用 `uv` 管理的虚拟环境：

```bash
cd ~/mcp_ic_tool_api
source .venv/bin/activate
python run_hcp_worker.py
```

如果放进 `tmux`：

```bash
tmux kill-session -t mcp_ic_api 2>/dev/null || true
tmux new-session -d -s mcp_ic_api 'cd ~/mcp_ic_tool_api && source .venv/bin/activate && python run_hcp_worker.py'
```

## 常用健康检查

### VPS 本机

```bash
curl http://127.0.0.1:8140/health
curl http://127.0.0.1:8130/api/health
```

### 公网静态报告

报告与轻量产物统一走：

```text
https://www.matterai.cn/dft/...
```

### HPC 到 VPS 控制面

当前实际走的是 `nginx` 的 `/hcp/` 转发，不依赖域名时可用 IP：

```bash
curl http://47.99.180.80/hcp/health
```

## 当前重要配置约定

### VPS

- `8130` 必须指向：

```text
VASP_SERVER_BASE_URL=http://127.0.0.1:8140
```

- 公网 MCP / nginx 必须透传 `user_id` header，否则 VPS 内部 MCP 正常但公网工具会看不到用户任务。
- 如果怀疑公网 MCP header 透传异常，可以在 VPS 上先用本地控制面验证，再测公网路径。

本地控制面验证示例：

```bash
curl http://127.0.0.1:8130/mcp \
  -X POST \
  -H "Content-Type: application/json" \
  -H "user_id: 10004" \
  -d '{"jsonrpc":"2.0","method":"tools/call","params":{"name":"list_user_tasks","arguments":{}},"id":1}'
```

### HPC

建议确认 `.env` 中至少有：

```text
VASP_SERVER_BASE_URL=http://47.99.180.80/hcp
INTERNAL_WORKER_TOKEN=...
VASP_PATH=/data/app/vasp-intel2022.3/6.3.2/bin/vasp_std
ONEAPI_ENV_SCRIPT=/data/app/intel/oneapi-2023.2/setvars.sh
SLURM_SRUN_MPI=pmi2
I_MPI_PMI_LIBRARY=/usr/lib64/libpmi2.so
```

当前 VASP 作业默认走：

```text
srun --mpi=pmi2 <VASP_PATH> >result.log 2>&1
```

## 常见注意事项

### 1. 不要再依赖 HPC 本地路径给前端

对外可访问地址必须是：

- `/dft/...`

不要再把这类路径直接返回给前端或邮件：

- `/data/home/...`
- `https://www.matterai.cn/static/...`
- `/download/file/...`

### 2. 分析任务并不都一样

目前已经适配成“VPS 创建分析任务，HPC 执行分析，再回传 VPS”的主要是：

- `DOS`
- `Band`
- `MD`
- `Agent Analyze`

如果未来继续扩展其它分析器，要优先沿用这条模式，不要让 VPS 直接依赖 HPC 本地 `work_directory`。

### 3. 更新时最好先看有没有活跃任务

虽然现在已经补了 runtime state 和 resume / reconcile 能力，但更稳妥的操作仍然是：

1. 先看有没有 `running / leased`
2. 尽量等当前任务结束
3. 再重启服务或 worker

### 4. HPC 到 VPS 有偶发网络抖动

现象通常是：

- `Read timed out`
- `ConnectTimeout`

当前 worker 已经有重试与补发机制，但排障时不要把偶发 warning 和真正失败混为一谈。

### 5. 使用 `uv`，不要随手混用系统 Python

统一：

```bash
source .venv/bin/activate
```

需要安装包时优先：

```bash
uv pip install -i https://pypi.tuna.tsinghua.edu.cn/simple <package>
```

### 6. 声子工具当前不是完整声子谱

当前 `phonon` 只适合：

- `Gamma` 点稳定性检查

不等于完整 phonon dispersion。  
如果后续要做完整声子谱，需要新工作流：

- 超胞
- 位移
- 力常数
- `phonopy`

### 7. Band / DOS / MD 报告都应该以用户体验为先

近期已经统一做过这些约束：

- 不显示 HPC 本地路径
- 不把超长原始数据直接塞给前端
- 尽量提供 HTML 报告、预览图、数据下载
- 图下增加“怎么看这张图”的通俗说明

### 8. 数值类 VASP 失败只反馈给智能体，不由后端自动重提

当前设计边界是：

- 后端负责识别失败类型并返回明确修改方向。
- 智能体负责判断是否重新提交、提交什么任务、如何调整参数。
- 后端不要因为 `BRMIX`、`ZBRENT` 等数值错误自动重提任务，避免隐式多烧机时。

接口返回里重点看：

- `failure_type`
- `suggested_action`

当前主要分类：

- `electronic_divergence`
  - 常见日志：`BRMIX`、`EDDDAV`、`ZHEGV`
  - 给智能体的方向：先做稳定 SCF 或轻量预优化，确认电子收敛后再继续。
- `ionic_line_search_failure`
  - 常见日志：`ZBRENT`
  - 给智能体的方向：优先用当前 `CONTCAR` 作为新 `POSCAR` 续算，并减小 `POTIM`。
- `lattice_inconsistency`
  - 常见日志：`HNFORM`、`Inconsistent Bravais lattice`
  - 给智能体的方向：规范化晶格、去掉近零噪声分量，必要时重新导出结构。
- `invalid_structure_file`
  - 常见日志：`Invalid CIF file with no structures`、`_atom_site_label`
  - 给智能体的方向：不要盲目重提计算，先换成有效 CIF/POSCAR。
- `missing_charge_density`
  - 常见场景：DOS/Band 任务缺少稳定 `CHGCAR`
  - 给智能体的方向：先补稳定 SCF，再做后续 NSCF/DOS/Band。

相关文件：

- `src/vasp_server/failure_guidance.py`
- `src/vasp_server/vasp_server_api.py`
- `src/vasp_server/schemas.py`
- `src/vasp_server/hpc_pull_worker.py`

### 9. 近期 bugfix 交接

最近已经修过并部署过的点：

- `module 'numpy' has no attribute 'trapz'`
  - 根因：HPC 上的 NumPy 版本已经移除 `np.trapz`。
  - 修复：MD 分析使用兼容 helper，优先 `np.trapezoid`，再回退 `np.trapz`。
  - 相关文件：`src/vasp_server/analyzers/md.py`、`tests/test_md_analyzer.py`
- CIF 缺 `_atom_site_label`
  - 根因：部分外部 CIF 不符合 pymatgen 默认期待。
  - 修复：读取 CIF 文本后补齐 atom site label，再交给 `CifParser`。
  - 相关文件：`src/vasp_server/base.py`、`tests/test_base.py`
- `HNFORM` / reciprocal lattice 不一致
  - 根因：POSCAR 里出现应为 0 的极小晶格噪声，VASP 判定晶格/k 点不一致。
  - 修复：CIF 转 POSCAR 后清理极小 lattice / fractional coordinate 噪声。
  - 相关文件：`src/vasp_server/base.py`、`tests/test_base.py`
- 失败诊断字段缺失
  - 修复：任务状态、远程分析、Agent 分析都透传 `failure_type` 和 `suggested_action`。
  - 相关文件：`src/vasp_server/failure_guidance.py`、`src/vasp_server/schemas.py`、`src/vasp_server/vasp_server_api.py`、`tests/test_task_status_api.py`

建议下一位开发者不要优先处理“历史老任务缺少 workdir / 原始日志导致无法分类”的问题，除非明确要做历史数据回填。对新任务和保留了 `work_directory` 的任务，失败诊断链路已经能工作。

### 10. GitHub 上传前检查

上传前不要提交：

- `.env`
- `tasks.db`
- `public_artifacts/`
- `*.zip`、`*.tar.gz`
- 真实 VPS/HPC 密码、私钥、token
- 用户上传文件、计算中间结果、报告产物

建议上传前执行：

```bash
git status --short
git diff -- .gitignore README.md docs/runbooks/session-handoff.md VASP_API_README.md
```

## 典型排障清单

### 任务一直 `queued`

先查：

1. HPC worker 是否真在跑

```bash
pgrep -af 'python run_hcp_worker.py'
tmux ls
```

2. VPS 控制面是否健康

```bash
curl http://127.0.0.1:8140/health
```

3. HPC 到 VPS 链路是否通

```bash
curl http://47.99.180.80/hcp/health
```

如果 worker 存在但不领取任务，再查：

- `queue_name` 是否一致。
- worker 是否启动在正确代码目录。
- worker 的 `VASP_SERVER_BASE_URL` 是否指向 `http://47.99.180.80/hcp`。
- `INTERNAL_WORKER_TOKEN` 是否与 VPS 控制面一致。
- 控制面 `/health` 中 `queued/running` 数量是否符合预期。

### 任务卡在 `leased`

常见原因：

- claim 成功了，但 worker 在 `mark_running` 前掉线

优先检查：

- 是否有 runtime state 文件
- 是否有对应 Slurm 作业
- 是否需要回收过期 lease

### 报告有了，但前端打不开

优先检查：

1. 是否真上传到了 VPS `public_artifacts`
2. 返回的是不是 `/dft/...`
3. `nginx` 的 `/dft/` 配置是否还在

## 新窗口接手建议

建议新窗口先按这个顺序来：

1. 先读根目录 `README.md` 和 `CLAUDE.md`。
2. 登录 VPS，确认 `ic_mcp` 和 `hpc_worker`。
3. 登录 HPC，确认 `mcp_ic_api`；如果看到 `mcp_ic_api_fix*` 之类临时 worker，先确认是否仍需要。
4. 跑健康检查。
5. 再开始提交任务或排障。

最小检查集合：

```bash
# VPS
curl http://127.0.0.1:8140/health
curl http://127.0.0.1:8130/api/health

# HPC
tmux ls
pgrep -af 'python run_hcp_worker.py'
curl http://47.99.180.80/hcp/health
```
