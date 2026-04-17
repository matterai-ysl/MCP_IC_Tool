# Task Failure Fixes Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 修复 band_structure 内部 seed-SCF 的误判失败，并补齐分布式 worker 在成功/失败路径上的 `ExecutionAttempt` 元数据落库。

**Architecture:** 保持现有控制面与 pull-worker 协议不大改。`band_structure` 只调整 seed-SCF 成功判定，避免依赖不存在的返回字段；分布式执行记录则通过 fail/complete 路径补建或更新 `ExecutionAttempt`，并把失败上下文从 HPC worker 带回控制面。

**Tech Stack:** Python, pytest, FastAPI schemas, SQLAlchemy task manager

---

### Task 1: Lock band_structure regression with a failing test

**Files:**
- Modify: `tests/test_vasp_worker_slurm.py`
- Modify: `src/vasp_server/vasp_worker.py`

**Step 1: Write the failing test**

添加一个 `run_band_structure_calculation()` 用例，模拟 `_run_vasp_calculation()` 成功返回但不包含 `convergence` 字段，同时 seed SCF 已经产出非空 `CHGCAR`，期望 band 任务继续执行而不是抛错。

**Step 2: Run test to verify it fails**

Run: `pytest tests/test_vasp_worker_slurm.py -k seed_scf -q`

**Step 3: Write minimal implementation**

让 seed-SCF 判定基于实际产物是否存在，而不是依赖 `_run_vasp_calculation()` 未返回的 `convergence` 字段。

**Step 4: Run test to verify it passes**

Run: `pytest tests/test_vasp_worker_slurm.py -k seed_scf -q`

### Task 2: Lock distributed attempt persistence with failing tests

**Files:**
- Modify: `tests/test_hpc_pull_worker.py`
- Modify: `tests/test_task_manager_claims.py`
- Modify: `src/vasp_server/hpc_pull_worker.py`
- Modify: `src/vasp_server/schemas.py`
- Modify: `src/vasp_server/internal_worker_api.py`
- Modify: `src/vasp_server/task_manager/manager.py`

**Step 1: Write the failing tests**

添加两类用例：
- pull-worker 失败上报时会把 `scheduler_job_id` 与 `work_directory` 一起带回控制面；
- `complete_execution()` / `fail_execution()` 在分布式 worker 路径下会创建或更新 `ExecutionAttempt`。

**Step 2: Run tests to verify they fail**

Run: `pytest tests/test_hpc_pull_worker.py tests/test_task_manager_claims.py -k 'failure_context or execution_attempt' -q`

**Step 3: Write minimal implementation**

给失败协议增加 `failure_context`，并在控制面侧统一补建/回填 `ExecutionAttempt` 的状态、目录、调度作业号和失败详情。

**Step 4: Run tests to verify they pass**

Run: `pytest tests/test_hpc_pull_worker.py tests/test_task_manager_claims.py -k 'failure_context or execution_attempt' -q`

### Task 3: Run targeted regression

**Files:**
- Test: `tests/test_vasp_worker_slurm.py`
- Test: `tests/test_hpc_pull_worker.py`
- Test: `tests/test_task_manager_claims.py`

**Step 1: Run focused regression suite**

Run: `pytest tests/test_vasp_worker_slurm.py tests/test_hpc_pull_worker.py tests/test_task_manager_claims.py -q`

**Step 2: Check for unrelated failures**

如果有失败，确认是否与本次修改直接相关；不扩 scope 修 unrelated 问题。
