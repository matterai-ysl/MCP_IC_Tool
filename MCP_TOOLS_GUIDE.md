# MCP 工具功能说明文档

本文档梳理当前 MCP 服务（`src/mcp_ic_tool/`）已实现的工具功能、已知局限，以及后续改进方向。

---

## 一、已实现功能

### 1. 计算提交工具

| 工具名 | 功能 | 结构来源 |
|--------|------|----------|
| `submit_structure_optimization` | 结构弛豫（ionic + cell relaxation） | `formula` 或 `cif_url` |
| `submit_scf_calculation` | 自洽场计算（SCF） | `formula` / `cif_url` / `optimized_task_id` |
| `submit_dos_calculation` | 态密度计算（DOS） | `formula` / `cif_url` / `scf_task_id` |
| `submit_band_structure_calculation` | 能带色散关系 E(k) | `formula` / `cif_url` / `scf_task_id` |
| `submit_md_calculation` | 分子动力学（MD） | `formula` / `cif_url` / `scf_task_id` |

**共同能力：**
- 所有提交工具均支持 `custom_incar` 字典，可覆盖任意 VASP INCAR 参数
- 结构选择支持丰富过滤器：`spacegroup`、`max_energy_above_hull`、`min/max_band_gap`、`min/max_nsites`、`stable_only`、`selection_mode`
- `submit_md_calculation` 支持多温度批量提交：`temperature=[300, 400, 500]` 一次调用提交多个独立任务

**INCAR 默认值：**

| 计算类型 | 关键默认参数 |
|----------|-------------|
| 结构优化 | ENCUT=520, EDIFF=1E-5, EDIFFG=-0.01, NSW=500, ISIF=3, IBRION=2 |
| SCF | 同优化，但 NSW=0, IBRION=-1, LWAVE=.TRUE., LCHARG=.TRUE. |
| DOS | 同SCF，另加 ISMEAR=-5, LORBIT=11, NEDOS=2001, ICHARG=11 |
| 能带 | 同SCF，另加 ISMEAR=0, LORBIT=11, ICHARG=11，k点为线模式 |
| MD | IBRION=0, SMASS=0（Nosé-Hoover）, ISYM=0, 1×1×1 Gamma |

---

### 2. 任务管理工具

| 工具名 | 功能 |
|--------|------|
| `get_task_status` | 查询单个任务状态、进度百分比、结果摘要 |
| `list_user_tasks` | 列出当前用户全部任务 |
| `cancel_task` | 取消排队中或运行中的任务 |
| `get_custom_incar_help` | 返回 INCAR 参数使用指南 |

**状态字段：**
- `status`: `queued` / `running` / `completed` / `failed` / `canceled`
- `analysis_status`: `pending` / `running` / `completed` / `failed`
- `progress`: 0–100 进度百分比
- `result_summary`: 计算完成后的结果摘要 JSON
- `html_report_url`: HTML 可视化报告下载链接

---

### 3. 独立分析工具

无需重新计算，可对已完成任务或用户上传的输出文件做分析：

| 工具名 | 依赖文件 | 输出 |
|--------|----------|------|
| `analyze_structure_optimization_results` | OUTCAR（必须）, POSCAR, CONTCAR | 力收敛、最终能量、HTML报告 |
| `analyze_scf_results` | OUTCAR（必须）, POSCAR+POTCAR+ACF.dat（可选Bader） | 收敛性、总能量、费米能级、带隙 |
| `analyze_dos_results` | vasprun.xml（必须）, OUTCAR（可选） | 带隙、费米能级、是否金属 |
| `analyze_band_structure_results` | vasprun.xml（必须）, OUTCAR, EIGENVAL（可选） | 带隙、直/间接、VBM/CBM |
| `analyze_md_results` | XDATCAR（必须）, OUTCAR, vasprun.xml, INCAR（可选） | 扩散系数、均方位移 |
| `analyze_md_multi_results` | 多个已完成 MD task_ids | Arrhenius 拟合、活化能、各温度离子电导率 |

分析输入可以是：
- `task_id`：直接引用服务器上已完成的任务
- `file_url`：指向 `.zip` 或 `.tar.gz` 压缩包的 URL（自动下载解压后分析）

---

### 4. 后端基础设施（已实现）

| 功能 | 实现方式 |
|------|----------|
| 任务持久化 | PostgreSQL（阿里云RDS），SQLAlchemy ORM |
| 任务重试 | 最多 `MAX_TASK_RETRIES` 次，区分可重试错误（节点故障）和不可重试错误（收敛失败） |
| 并发控制 | `ThreadPoolExecutor`，`MAX_CONCURRENT_TASKS` 环境变量（默认4） |
| 用户配额 | 单用户并发上限 `MAX_USER_CONCURRENT_TASKS=3`，每日总量 `MAX_USER_TOTAL_TASKS_PER_DAY=50` |
| 计算去重 | 基于 task_type + 关键参数的 SHA-256 哈希，命中已完成任务直接返回 |
| 健康检查 | `GET /health` 检查数据库连通性，异常返回 503 |
| 监控指标 | Prometheus `/metrics` 端点，含 queued/running/submitted_total/completed_total |
| 结构化日志 | 全局 `logging` 模块，含线程名，覆盖所有子模块 |
| 临时文件清理 | 分析结束后 `finally` 块自动 `shutil.rmtree` |
| 输入安全校验 | formula 通过 pymatgen Composition 校验；custom_incar 关键字黑名单；cif_url 域名白名单 |

---

## 二、已知局限 / 暂不支持

### 2.1 任务进度只能轮询
`get_task_status` 需要 Agent 反复调用才能知道任务完成。对于长达数小时的 VASP 计算，这会消耗大量 token，且时机难以把握。

**没有**：SSE 推送、WebSocket、Webhook 回调机制。

---

### 2.2 多温度 MD 工作流需要手动串联
`submit_md_calculation(temperature=[300,400,500])` 会返回多个 task_id，但：
- Agent 必须自己等所有任务完成后，再手动调用 `analyze_md_multi_results`
- 没有"等待一组任务全部完成"的工具，只能逐个轮询

---

### 2.3 没有资源预估
提交计算前无法预知：
- 大概需要多少小时
- 所需内存/CPU核心数
- 是否会因原子数过多导致计算超时

---

### 2.4 没有工作流编排工具
常见的计算流程（Opt → SCF → DOS + Band）需要 Agent 自己串联多个工具，并管理中间 task_id。没有一键式"完整电子结构分析"工具。

---

### 2.5 结构可视化缺失
MCP 工具只返回文本摘要和 HTML 报告 URL，没有直接在对话中展示结构图或能带图的能力。

---

### 2.6 不支持的计算类型
以下常见 DFT 任务当前没有对应工具：
- NEB（过渡态/反应路径）
- 声子谱（Phonon）
- 弹性常数
- 电荷密度差分图
- 表面/界面模型（slab）计算

---

### 2.7 结果文件下载需借助外部 URL
HTML 报告通过 `html_report_url` 返回，但：
- 该 URL 指向 VASP Server 的静态文件服务，Agent 无法直接"读取"报告内容
- 没有 MCP 工具能把报告内容作为文本返回给 Agent

---

## 三、后续改进方向

### 优先级 P0（影响可用性）

**任务进度推送**
- 在 VASP Server 添加 SSE 端点 `GET /vasp/task/{task_id}/stream`
- 或支持 Webhook：任务完成时主动 POST 回调到指定 URL
- MCP 侧添加 `wait_for_task(task_id, timeout_minutes)` 工具，内部轮询但减少 Agent 侧调用次数

**批量状态查询**
- 添加 `get_tasks_status(task_ids: List[str])` 工具，一次调用查询多个任务，避免 MD 多温度场景下反复调用

---

### 优先级 P1（提升体验）

**工作流编排工具**
```
submit_full_electronic_structure(formula, run_dos=True, run_band=True)
```
内部自动串联 Opt → SCF → DOS + Band，返回统一结果，减少 Agent 需要管理的中间状态。

**资源预估工具**
```
estimate_calculation_cost(formula, calc_type, custom_incar)
```
基于原子数、k点密度、计算类型估算运行时间和内存需求，让 Agent 在提交前给用户预警。

**报告内容工具**
```
get_task_report_text(task_id)  # 返回 HTML 报告的文本摘要
```
将报告关键数据（能量、带隙、扩散系数等）直接以结构化 JSON 返回，不依赖用户打开 URL。

---

### 优先级 P2（功能扩展）

**新计算类型**

| 工具名（规划） | 计算类型 | 关键 INCAR |
|----------------|----------|-----------|
| `submit_neb_calculation` | 过渡态 NEB | IMAGES, SPRING |
| `submit_phonon_calculation` | 声子谱（DFPT） | IBRION=8, LEPSILON=.TRUE. |
| `submit_elastic_calculation` | 弹性常数 | IBRION=6, ISIF=3 |

**结构操作工具**
- `build_slab_structure(formula, miller_index, vacuum)` — 基于 pymatgen 构建表面模型
- `get_structure_info(task_id)` — 返回优化后结构的晶格参数、空间群、键长等

**cif_url 来源扩展**
当前域名白名单较严格，建议增加：
- `materialsproject.org` 直接 API 获取
- `crystallography.net`（COD 数据库）
- 用户自建对象存储（OSS/S3）

---

### 优先级 P3（工程改进）

| 问题 | 建议 |
|------|------|
| CORS 允许全部来源 | 生产环境配置具体域名白名单 |
| 数据库凭证硬编码兜底 | 强制要求环境变量，启动时无配置则拒绝启动 |
| 没有 API 版本控制 | 所有端点加 `/v1/` 前缀，为后续兼容性预留空间 |
| 任务日志不可查询 | 添加 `get_task_logs(task_id)` 接口返回 VASP 运行日志片段 |

---

## 四、典型工作流示例

### 完整电子结构分析（当前需 Agent 手动串联）

```
1. submit_structure_optimization(formula="LiFePO4")
   → 等待 → task_id_opt

2. submit_scf_calculation(optimized_task_id=task_id_opt)
   → 等待 → task_id_scf

3. 并行提交：
   submit_dos_calculation(scf_task_id=task_id_scf)
   submit_band_structure_calculation(scf_task_id=task_id_scf)

4. 等待两者完成，读取 html_report_url
```

### 多温度离子扩散分析

```
1. submit_md_calculation(
       formula="Li2O",
       temperature=[600, 800, 1000, 1200],
       md_steps=5000
   )
   → 返回 4 个 task_id

2. 逐个轮询 get_task_status(task_id) 直到全部 completed

3. analyze_md_multi_results(task_ids=[...])
   → Arrhenius 拟合 + 活化能 + 离子电导率
```
