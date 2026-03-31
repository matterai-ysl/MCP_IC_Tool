# VPS Control Plane + HPC Pull Worker Contract

## Goal

把当前单进程 VASP 服务升级成“VPS 控制面 + HPC 拉取式 worker + 对象存储”的分布式架构，并先锁定任务生命周期契约，避免后续 API、数据库和 worker 行为各自漂移。

## Canonical Task States

- `queued`: 任务已入队，等待 worker claim
- `leased`: 某个 worker 已领取任务，并持有有效租约
- `running`: HPC 侧作业已提交并开始执行
- `uploading`: worker 正在上传执行产物
- `analyzing`: VPS 正在做结果分析和摘要生成
- `completed`: 整个流程成功结束
- `failed`: 流程失败且不再重试
- `cancel_requested`: 用户或系统请求取消，等待 worker 确认
- `canceled`: worker 已确认取消并完成收尾

## Required Metadata

`Task` 需要携带租约和路由字段，以支持数据库驱动的分布式编排：

- `tenant_id`
- `queue_name`
- `priority`
- `worker_id`
- `lease_token`
- `lease_expires_at`
- `heartbeat_at`
- `cancel_requested_at`
- `finalized_at`
- `retry_count`
- `max_retries`

`ExecutionAttempt` 需要记录具体执行实例上下文：

- `worker_id`
- `lease_token`
- `scheduler_job_id`
- `artifact_manifest`

## Transition Rules

- 正常主路径：`queued -> leased -> running -> uploading -> analyzing -> completed`
- 取消路径：`queued|leased|running -> cancel_requested -> canceled`
- 租约过期恢复：`leased -> queued`
- 运行失败重试：`running -> queued`，仅在 `retry_count < max_retries` 时允许
- 运行失败终态：`running -> failed`，当重试预算耗尽时生效

## Design Constraints

- 只有 VPS 控制面可以读写任务编排状态。
- HPC worker 只能通过内部 API 认领任务、续租、上报状态。
- 对象存储是产物唯一可信来源，公共 API 不再暴露本地路径。
