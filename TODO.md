# 生产级改进 TODO

## P0 — 已完成
- [x] DetachedInstanceError: manager.py 返回 ORM 对象改为 SimpleNamespace
- [x] tarfile/zipfile 路径穿越漏洞: extractall 前校验每个成员路径
- [x] SQLite 并发: check_same_thread=False + WAL mode + scoped_session

## P1 — 已完成
- [x] 任务重试机制: manager.py 支持 MAX_TASK_RETRIES，自动区分可重试错误（HPC节点故障/超时）和不可重试错误（VASP收敛失败），每次重试创建新 ExecutionAttempt
- [x] 任务并发控制: daemon thread 改为 ThreadPoolExecutor，MAX_CONCURRENT_TASKS 环境变量控制（默认 4）
- [x] 健康检查: 添加 /health 端点，检查数据库连接、报告任务计数，异常返回 503
- [x] 日志规范化: 全局替换 print() 为 logger（vasp_worker.py, database.py, base.py, client.py, mp.py 等），删除 emoji
- [x] 配置管理: 新增 settings.py (Pydantic BaseSettings)，Config.py 和 database.py 统一从 settings 读取，requirements.txt 添加 pydantic-settings

## P2 — 已完成
- [x] MCP get_user_id 可返回 None: 缺失时抛出 ValueError 而非静默传 None
- [x] 临时文件清理: _resolve_analysis_input_dir 返回 (work_dir, tmp_root)，分析端点在 finally 中 shutil.rmtree 清理
- [x] 输入校验增强: formula 通过 pymatgen Composition 校验，custom_incar key 黑名单（SYSTEM/ISTART/ICHARG/IBRION/NSW/LCHARG/LWAVE），cif_url 域名白名单
- [x] 计算结果缓存/去重: 基于 task_type + 关键参数计算 SHA-256 input_hash，提交前查询已完成任务，命中则直接返回已有 task_id
- [x] 费用/配额管理: 单用户并发任务数限制 (MAX_USER_CONCURRENT_TASKS=3)，每日任务总数限制 (MAX_USER_TOTAL_TASKS_PER_DAY=50)
- [x] Prometheus Metrics: 集成 prometheus_fastapi_instrumentator，/metrics 端点，自定义指标 vasp_tasks_queued/running/submitted_total/completed_total

## P3 — 未来改进

### 任务进度推送
当前只能轮询 get_task_status，智能体需反复调用浪费 token。方案：
- SSE 端点 `/vasp/task/{task_id}/stream`
- 或 webhook 回调机制，任务完成时主动通知

### Band Structure 计算 — 已完成
- [x] 添加 submit_band_structure_calculation MCP 工具
- [x] 添加 BandStructureRequest/Response/Result 模型
- [x] 添加 band_structure 任务类型到 manager.py 映射
- [x] VaspWorker: run_band_structure_calculation（SCF+能带 / 纯能带）
- [x] pymatgen HighSymmKpath 自动生成高对称k路径
- [x] BandStructureAnalyzer + HTML 报告（能带图绘制）
- [x] 独立分析端点 POST /vasp/analyze/band-structure
- [x] MCP 工具: submit_band_structure_calculation + analyze_band_structure_results

### 任务资源估算
- 原子数 → 预估计算时间/内存
- 提交前给用户预警大规模计算
