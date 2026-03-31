# Distributed Execution Cutover Runbook

## Goal

把当前 VASP 服务从“单机本地执行 + 本地文件暴露”平滑切换到“VPS 控制面 + HPC pull worker + 对象存储”。

## Rollout Phases

### Phase A: Database Queue

- VPS 公网 API 只负责校验请求和写入任务记录。
- 老的本地执行链路仍可保留给开发环境。
- 验证项：
  - 新提交任务状态为 `queued`
  - claim / heartbeat / cancel / retry 状态流可在数据库中正确反映

### Phase B: Internal Worker API Canary

- 启用 `/internal/workers/*` 端点和内部 token。
- 只放一个 HPC worker canary。
- 先放 `structure_optimization`。
- 验证项：
  - worker 能 claim、续租、上报 running、complete、fail
  - cancel 请求能经过 heartbeat 传播到 worker

### Phase C: Object Storage Required

- 新任务产物必须写对象存储元数据。
- 任务状态接口只返回签名 URL。
- 验证项：
  - 公共 API 中不再出现 `/static/` 或 `/download/file/`
  - HTML 报告和 artifacts 均返回 `https://` 签名地址

### Phase D: Legacy Path Removal

- 禁用本地文件下载接口。
- 停止依赖共享目录复用上游任务产物。
- 仅保留 pull-worker 执行模式。

## Backfill

- 先运行 `scripts/backfill_artifacts.py`，把历史 `storage_backend=local` 的产物记录迁移到对象存储 key。
- 再运行 `scripts/requeue_stuck_tasks.py`，回收过期租约和心跳超时任务。

## Production Acceptance

- 没有 `queued` 任务超过 2 个 claim interval
- 没有 `leased` 任务超过 lease timeout 且未回收
- 公共 API 响应中没有本地文件系统路径
- canary task type 在至少一个真实 HPC worker 上稳定成功
