# 代码审查报告 - MCP_IC_Tool

> 审查日期：2026-03-16
> 审查范围：全部 Python 源文件

---

## 概览

| 严重级别 | 数量 |
|---------|------|
| 严重 (Critical) | 5 |
| 高危 (High) | 6 |
| 中危 (Medium) | 9 |
| 低危 (Low) | 6 |
| 架构问题 | 4 |

---

## 严重问题 (Critical)

### C1. 数据库密码明文硬编码在源代码中

**文件：** `src/vasp_server/task_manager/database.py`, 第 6-10 行

```python
DB_HOST = os.getenv("DB_HOST", "pgm-uf69uij17z9vh123jo.pg.rds.aliyuncs.com")
DB_USER = os.getenv("DB_USER", "a2252222223")
DB_PASSWORD = os.getenv("DB_PASSWORD", "Jixiaobei123")  # 明文密码
```

**风险：** 代码提交到 Git 仓库后，密码永久存在于历史记录中，任何能看到仓库的人均可直接获取数据库访问权限。

**修复：**
- 立即轮换数据库密码
- 从代码中删除所有默认值，改为强制要求环境变量
- 使用 `.env` 文件（加入 `.gitignore`）或密钥管理系统

```python
DB_PASSWORD = os.environ["DB_PASSWORD"]  # 不设默认值，缺失时直接报错
```

---

### C2. Materials Project API Key 明文硬编码

**文件：** `src/vasp_server/mp.py`, 第 6 行

```python
API_KEY = 'g8j3tc9BUugnPSzgJ2ppCaxPgEo5W8H7'
```

**风险：** API Key 已暴露于版本控制，可被滥用，导致配额耗尽或账号封禁。

**修复：**
- 立即到 Materials Project 官网撤销并重新生成此 Key
- 改为从环境变量读取：`API_KEY = os.environ["MP_API_KEY"]`

---

### C3. 文件下载接口存在路径穿越漏洞

**文件：** `src/vasp_server/vasp_server_api.py`, 第 570-578 行

```python
full_path = Path(DOWNLOAD_URL) / file_path
if not full_path.exists():
    raise HTTPException(status_code=404, ...)
# 缺少：验证 full_path 确实在 DOWNLOAD_URL 目录内
return FileResponse(full_path, ...)
```

**风险：** 攻击者可通过 `../../../etc/passwd` 等路径读取服务器上的任意文件。

**修复：**
```python
resolved = full_path.resolve()
allowed = Path(DOWNLOAD_URL).resolve()
if not str(resolved).startswith(str(allowed)):
    raise HTTPException(status_code=403, detail="禁止访问")
```

---

### C4. CORS 配置允许所有来源且同时开启 credentials

**文件：** `src/vasp_server/vasp_server_api.py`, 第 38-44 行

```python
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_credentials=True,  # 与 "*" 同时使用会导致 CSRF 漏洞
    allow_methods=["*"],
    allow_headers=["*"],
)
```

**风险：** `allow_origins=["*"]` 与 `allow_credentials=True` 同时设置，实际上 FastAPI/Starlette 会拒绝这种组合（浏览器不接受），但这暴露了没有访问控制设计的问题。任意第三方网站可发起跨域请求。

**修复：** 将 `allow_origins` 限制为具体的可信域名列表。

---

### C5. 无认证机制，user_id 默认值为 "123"

**文件：** `src/vasp_server/vasp_server_api.py` 第 103 行；`src/mcp_ic_tool/mcp_server.py` 第 54-59 行

```python
# vasp_server_api.py
if request.user_id is None:
    request.user_id = "123"

# mcp_server.py
def get_user_id(ctx: Context) -> str:
    user_id = ctx.request_context.request.headers.get("user_id", None)
    else:
        user_id = "123"
    return user_id
```

**风险：** 所有未携带 `user_id` 的请求都归属同一用户，导致任意用户可以查看、取消他人的任务。系统完全没有认证和授权机制。

**修复：** 实现 JWT 或 API Key 认证，未认证请求返回 401。

---

## 高危问题 (High)

### H1. 静态文件目录无访问控制，所有用户数据公开可访问

**文件：** `src/vasp_server/vasp_server_api.py`, 第 48 行

```python
app.mount("/static", StaticFiles(directory="/data/home/ysl9527/vasp_calculations"), name="static")
```

**风险：** `/static/` 路径下可直接浏览和下载所有用户的全部计算文件，无需任何身份验证。

**修复：** 删除此静态挂载，仅通过有访问控制的 `/download/file/` 端点提供文件下载。

---

### H2. TaskManager 中大量重复代码，且工作目录路径硬编码

**文件：** `src/vasp_server/task_manager/manager.py`, 第 77-244 行

四种任务类型（`structure_optimization`、`scf_calculation`、`dos_calculation`、`md_calculation`）的处理逻辑几乎完全相同，唯一的区别是调用不同的 `VaspWorker` 方法。每段都重复以下操作：
- 导入 `asyncio`、`VaspWorker`
- 实例化 `VaspWorker(user_id=..., base_work_dir="/data/home/ysl9527/vasp_calculations")`
- 创建事件循环、运行、关闭
- 检查 result 并更新 task 状态

工作目录 `/data/home/ysl9527/vasp_calculations` 在 `manager.py` 中硬编码了 4 次，在 `vasp_worker.py` 构造函数默认值中又出现 1 次。

**修复：** 重构为单个通用方法：

```python
TASK_METHOD_MAP = {
    "structure_optimization": "run_structure_optimization",
    "scf_calculation": "run_scf_calculation",
    "dos_calculation": "run_dos_calculation",
    "md_calculation": "run_md_calculation",
}

def _run_vasp_task(self, task, cancel_event, db):
    method_name = TASK_METHOD_MAP.get(str(task.task_type))
    if not method_name:
        raise ValueError(f"未知任务类型: {task.task_type}")
    worker = VaspWorker(user_id=str(task.user_id), base_work_dir=BASE_WORK_DIR)
    loop = asyncio.new_event_loop()
    try:
        return loop.run_until_complete(getattr(worker, method_name)(task.id, task.params or {}, ...))
    finally:
        loop.close()
```

将 `base_work_dir` 移至 `Config.py` 统一管理。

---

### H3. 进度回调中的数据库异常会导致任务静默失败

**文件：** `src/vasp_server/task_manager/manager.py`, 第 62-74 行

```python
async def progress_callback(progress: int, message: str, pid=None):
    if cancel_event.is_set():
        raise Exception("任务已被取消")
    task.progress = progress
    db.add(task)
    db.commit()  # 若此处抛出异常，整个计算任务会被当作失败
```

**风险：** 进度更新时数据库临时不可用，会导致正在运行的 VASP 计算被中断并标记为失败，造成算力浪费。

**修复：** 进度回调中的数据库错误应只记录日志，不应向上抛出，计算本身应继续运行。

---

### H4. 取消任务存在竞态条件

**文件：** `src/vasp_server/task_manager/manager.py`, 第 277-296 行

```python
task = db.get(Task, task_id)          # 读取状态
if str(task.status) in {"completed", ...}:
    return True
# 此处，worker 线程可能已将状态改为 completed
with self._lock:
    cancel_event.set()
task.status = "canceling"              # 覆盖了实际已完成的状态
db.commit()
```

**风险：** 已完成的任务可能被错误地标记为 `canceling`。

**修复：** 使用数据库级别的条件更新：
```sql
UPDATE tasks SET status='canceling' WHERE id=? AND status IN ('queued','running')
```

---

### H5. get_task 接口中 task 为 None 时的判断顺序错误

**文件：** `src/vasp_server/vasp_server_api.py`, 第 418-433 行

```python
task = task_manager.get_task(task_id, user_id)
logger.debug(f"Task object: {task}")
for attr in dir(task):          # 若 task 是 None，此处直接 AttributeError
    ...
if not task:                    # 这个判断在 dir(task) 之后，顺序错误
    raise HTTPException(status_code=404, ...)
```

**风险：** 当任务不存在时，会抛出 `AttributeError: 'NoneType' object has no attribute ...`，而不是返回 404，且会暴露内部堆栈信息。

**修复：** 将 `if not task` 检查移到所有属性访问之前。

---

### H6. 子进程执行时 AWK 命令使用字符串插值

**文件：** `src/vasp_server/base.py`（POTCAR 生成相关代码）

```python
cmd = ["awk", f'$1 == "{element}" {{print $2}}', paw_pick_list]
```

**风险：** 若 `element` 包含引号或特殊字符，可导致命令注入。虽然元素名来自 CIF 文件而非直接用户输入，但仍存在风险。

**修复：** 改用 Python 原生读取文件代替 AWK，或对 `element` 值做严格的白名单校验（只允许字母）。

---

## 中危问题 (Medium)

### M1. `vasp_server_api.py` 中保留了大量调试代码

**文件：** `src/vasp_server/vasp_server_api.py`, 第 420-431、475-477 行

`get_task_status` 接口中包含大量 `logger.debug` 调用，遍历并打印 task 对象的所有属性，包括潜在的敏感结果数据。这些代码在生产环境中不应存在。

**修复：** 删除 `dir(task)` 遍历日志，保留必要的业务日志。

---

### M2. 数据库连接池未配置，高并发下可能耗尽连接

**文件：** `src/vasp_server/task_manager/database.py`, 第 14-18 行

```python
engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    pool_pre_ping=True,
    pool_recycle=300,
    # 缺少 pool_size 和 max_overflow
)
```

每个任务会启动一个线程，每个线程持有一个数据库 Session，高并发时可能耗尽默认连接池（默认 5 个连接）。

**修复：**
```python
engine = create_engine(
    SQLALCHEMY_DATABASE_URL,
    pool_pre_ping=True,
    pool_recycle=300,
    pool_size=20,
    max_overflow=10,
)
```

---

### M3. 无数据库迁移机制

**文件：** `src/vasp_server/task_manager/database.py`, 第 23-28 行

使用 `Base.metadata.create_all()` 管理表结构，对已有数据库执行时不会修改已存在的表结构，导致字段变更无法应用。

**修复：** 引入 [Alembic](https://alembic.sqlalchemy.org/) 管理数据库迁移版本。

---

### M4. `custom_incar` 参数未做白名单校验

**文件：** `src/vasp_server/schemas.py`

```python
custom_incar: Optional[Dict[str, Any]] = Field(None, description="自定义INCAR参数字典")
```

任意键值均可传入，会被写入 INCAR 文件并提交给 VASP 执行。恶意或错误的参数可能导致计算结果错误或 VASP 崩溃。

**修复：** 建立 VASP 合法参数白名单，拒绝不在白名单内的参数键名。

---

### M5. `task_manager/manager.py` 中 `# type: ignore` 泛滥

**文件：** `src/vasp_server/task_manager/manager.py`（全文约 30+ 处）

```python
task: Task = db.get(Task, task_id)  # type: ignore
task.status = "running"  # type: ignore
task.progress = 1  # type: ignore
```

根本原因是 SQLAlchemy ORM 模型的字段在运行时是 Column 类型，但赋值时是普通 Python 值，类型检查器无法正确推断。可通过正确配置 SQLAlchemy 的 Mapped 注解解决。

**修复：** 升级到 SQLAlchemy 2.x 的 `Mapped` 类型注解：
```python
from sqlalchemy.orm import Mapped, mapped_column
class Task(Base):
    id: Mapped[str] = mapped_column(String, primary_key=True)
    status: Mapped[str] = mapped_column(String, default="queued")
```

---

### M6. `TaskManager._run_task_worker` 中异常处理可能导致二次异常

**文件：** `src/vasp_server/task_manager/manager.py`, 第 265-273 行

```python
except Exception as exc:
    task = db.get(Task, task_id)  # 若 db Session 已因异常损坏，此处会再次抛出
    if task is not None:
        task.status = "failed"
        db.add(task)
        db.commit()
```

若异常来自数据库操作，Session 可能已处于不可用状态，`db.get()` 会抛出新异常，导致任务状态永远停留在 `running`。

**修复：** 在 except 块中新建一个独立的 Session 来更新任务状态。

---

### M7. MCP 工具未对返回数据做错误结构校验

**文件：** `src/mcp_ic_tool/client.py`

HTTP 请求失败时（非 2xx 状态码），直接返回 `resp.json()` 或 `{}` 给 LLM 代理，没有统一的错误结构，代理难以区分成功和失败。

**修复：** 统一错误返回格式：
```python
if resp.status_code != 200:
    return {"success": False, "error": resp.text, "status_code": resp.status_code}
```

---

### M8. Materials Project API 请求无超时控制

**文件：** `src/vasp_server/mp.py`

`MPRester` 调用没有设置超时，若 MP 服务不可用，整个计算线程会无限阻塞，导致该 task 永远处于 `running` 状态。

**修复：** 设置超时并在超时后将任务标记为失败。

---

### M9. `single_port_server.py` 路径安全检查使用字符串前缀比较，存在绕过风险

**文件：** `single_port_server.py`, 第 122-125 行

```python
if not str(resolved_path).startswith(str(models_resolved)):
```

在某些平台（如 macOS 的大小写不敏感文件系统）或存在同名路径前缀时（如 `/data/models` 和 `/data/models_backup`），`startswith` 可能产生误判。

**修复：** 使用 Python 3.9+ 的 `Path.is_relative_to()`：
```python
if not resolved_path.is_relative_to(models_resolved):
    raise HTTPException(status_code=403)
```

---

## 低危问题 (Low)

### L1. `Config.py` 中 `import os` 重复两次

**文件：** `src/vasp_server/Config.py`, 第 1-2 行

```python
import os
import os  # 重复导入
```

**修复：** 删除重复的 `import os`。

---

### L2. `base.py` 中遗留硬编码测试路径

**文件：** `src/vasp_server/base.py`（`if __name__ == "__main__"` 块）

```python
cif_to_poscar("/Users/ysl/Desktop/Code/...")
```

遗留了开发者本地路径，说明测试代码未清理。

**修复：** 删除或移至单独的测试文件。

---

### L3. 错误信息中英文混用

全项目中，API 返回的错误信息、日志消息混用中文和英文（例如 `"任务未找到或无权限访问"` 与 `"Task not found"` 并存）。

**修复：** 统一使用一种语言，或引入 i18n 机制。

---

### L4. `vasp_server_api.py` 中注释掉的代码块过多

**文件：** `src/vasp_server/vasp_server_api.py`, 第 598-769 行

超过 170 行被注释掉的代码（`get_md_result` 接口及 `_analyze_single_temp_result` 函数）仍保留在文件中，严重影响可读性。

**修复：** 若确认不再需要，直接删除；若将来可能恢复，使用 Git 保存历史记录即可。

---

### L5. POST 接口返回状态码应为 201 而非默认 200

**文件：** `src/vasp_server/vasp_server_api.py`（所有 `@app.post` 装饰器）

提交任务的接口语义上是"创建资源"，应返回 HTTP 201，但当前默认返回 200。

**修复：** `@app.post("/vasp/structure-optimization", response_model=StructOptResponse, status_code=201)`

---

### L6. `mcp_server.py` 中 `print` 调试语句未清理

**文件：** `src/mcp_ic_tool/mcp_server.py`, 第 107-108 行

```python
print("submit_structure_optimization")
print(params)
```

生产代码中应使用 logger，不应直接 print。

**修复：** 改为 `logger.debug(...)` 或删除。

---

## 架构问题

### A1. 无任何认证授权体系

整个系统（VASP Server API 和 MCP Server）均无认证机制。`user_id` 通过 HTTP Header 或请求体明文传递，任何人可伪造任意 `user_id` 访问他人数据。

**建议：** 在 MCP Server 层实现 API Key 认证，VASP Server API 仅监听内网地址（不对外暴露）。

---

### A2. 任务执行模型不适合生产环境

每个任务独占一个 Python daemon 线程 + 一个 asyncio 事件循环。VASP 计算可能持续数小时，高并发时会导致：
- 线程数量失控
- 进程重启后所有运行中任务状态丢失（数据库中停留在 `running`）
- 无法跨进程或跨机器分发任务

**建议：** 引入 Celery + Redis/RabbitMQ 作为任务队列，或使用专门的 HPC 作业调度（已有 SLURM/LSF 相关代码，可进一步整合）。增加启动时将 `running` 状态的孤儿任务恢复为 `failed` 的逻辑。

---

### A3. 服务重启后 running 状态任务成为孤儿

**文件：** `src/vasp_server/task_manager/database.py` / `vasp_server_api.py`

服务重启后，数据库中处于 `running` 状态的任务永远不会再被更新，但实际计算进程已经终止。没有任何启动时的恢复逻辑。

**修复：** 在服务启动时，将所有 `running` 状态的任务标记为 `failed`，或实现断点续算机制。

---

### A4. 没有 API 版本管理

所有接口路径均为 `/vasp/xxx`，没有版本号（如 `/v1/vasp/xxx`）。未来 API 变更将导致所有客户端（包括 MCP 工具）同步失效。

**建议：** 在路由中引入版本前缀，或通过 FastAPI 的 `APIRouter` 管理多版本。

---

## 立即需要执行的操作

以下两项需要**立刻**处理，因为敏感信息已存在于 Git 历史中：

1. **立即轮换 Materials Project API Key**（当前值已通过本报告可见）
2. **立即修改数据库密码**（当前密码已通过本报告可见）
3. 在 `.gitignore` 中加入 `.env`，并将所有凭证迁移至环境变量
4. 修复 `vasp_server_api.py` 中文件下载接口的路径穿越漏洞（C3）
5. 将 `if not task` 判断移至属性访问之前（H5）
