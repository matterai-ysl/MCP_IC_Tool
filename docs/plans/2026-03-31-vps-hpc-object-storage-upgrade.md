# VPS Control Plane + HPC Pull Worker + Object Storage Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Refactor the current VASP service into a VPS-hosted control plane with database-backed orchestration, HPC workers that actively claim tasks from the VPS, and object storage as the canonical artifact store.

**Architecture:** The FastAPI service on the VPS becomes the only public control plane and the only component that reads and writes orchestration state in PostgreSQL. HPC workers no longer share a filesystem or rely on in-process `ThreadPoolExecutor`; instead they poll internal worker endpoints, acquire leases for queued tasks, run Slurm jobs locally on the cluster, upload artifacts to object storage, and report completion back to the VPS. Frontend and MCP clients only access VPS APIs plus signed object-storage URLs.

**Tech Stack:** FastAPI, SQLAlchemy, PostgreSQL, Pydantic settings, existing `TaskManager` models, Alibaba OSS or S3-compatible object storage, Slurm, pytest.

---

### Task 1: Lock The New Control-Plane Contract

**Files:**
- Create: `docs/architecture/vps-hpc-control-plane.md`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_state_machine.py`

**Step 1: Write the failing test**

```python
def test_task_lifecycle_for_pull_workers():
    assert [
        "queued",
        "leased",
        "running",
        "uploading",
        "analyzing",
        "completed",
    ]
```

Add tests that assert:
- `queued -> leased -> running -> uploading -> analyzing -> completed`
- `queued -> leased -> cancel_requested -> canceled`
- `leased -> lease_expired -> queued`
- `running -> runtime_failed -> queued` only when retry budget remains

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_state_machine.py`
Expected: FAIL because the new states and model fields do not exist yet.

**Step 3: Write minimal implementation**

Update the task schema contract first.

- In `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py`, extend `Task` with:
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
- Extend `ExecutionAttempt` with:
  - `worker_id`
  - `lease_token`
  - `scheduler_job_id`
  - `artifact_manifest`
- Add or document allowed statuses in `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py`.

Use these canonical task states:

```python
TASK_STATES = {
    "queued",
    "leased",
    "running",
    "uploading",
    "analyzing",
    "completed",
    "failed",
    "cancel_requested",
    "canceled",
}
```

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_state_machine.py`
Expected: PASS

**Step 5: Commit**

```bash
git add docs/architecture/vps-hpc-control-plane.md src/vasp_server/schemas.py src/vasp_server/task_manager/models.py tests/test_task_state_machine.py
git commit -m "refactor: define pull-worker task state machine"
```

### Task 2: Replace In-Process Execution With Database-Backed Orchestration

**Files:**
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/database.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/settings.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_manager_claims.py`

**Step 1: Write the failing test**

```python
def test_claim_task_assigns_worker_and_lease_token(task_manager):
    claimed = task_manager.claim_next_task(worker_id="hpc-a", queue_name="default")
    assert claimed.worker_id == "hpc-a"
    assert claimed.lease_token
    assert claimed.status == "leased"
```

Add tests for:
- only one worker can claim the same queued task
- expired lease returns task to `queued`
- `cancel_requested` task is never newly claimed

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_manager_claims.py`
Expected: FAIL because `TaskManager` still uses `_executor` and `_cancel_flags`.

**Step 3: Write minimal implementation**

Refactor `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py`:

- Remove in-process execution ownership from:
  - `self._executor`
  - `self._cancel_flags`
  - `_run_task_worker(...)`
- Keep `submit_task(...)`, but make it only insert a row with `status="queued"`.
- Add methods:
  - `claim_next_task(worker_id: str, queue_name: str) -> Optional[ClaimedTask]`
  - `heartbeat_task(task_id: str, lease_token: str, worker_id: str) -> None`
  - `mark_task_running(...)`
  - `request_cancel(task_id: str, user_id: str) -> bool`
  - `ack_cancel(...)`
  - `complete_execution(...)`
  - `fail_execution(...)`
  - `requeue_expired_leases(...)`
- Use a single transactional claim query with row locking when PostgreSQL is active.
- In `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/settings.py`, add:
  - `task_lease_seconds`
  - `worker_poll_interval_seconds`
  - `worker_heartbeat_timeout_seconds`

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_manager_claims.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/task_manager/manager.py src/vasp_server/task_manager/database.py src/vasp_server/settings.py tests/test_task_manager_claims.py
git commit -m "refactor: move task orchestration fully into database"
```

### Task 3: Add Internal Worker API On The VPS

**Files:**
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/internal_worker_api.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_internal_worker_api.py`

**Step 1: Write the failing test**

```python
def test_worker_claim_endpoint_returns_leased_task(client):
    resp = client.post("/internal/workers/claim", headers={"Authorization": "Bearer worker-token"})
    assert resp.status_code == 200
    assert resp.json()["status"] == "leased"
```

Add tests for:
- worker registration or token validation
- claim response includes `lease_token`, task params, and upstream artifact manifest
- heartbeat refreshes `lease_expires_at`
- complete/fail endpoints reject invalid lease tokens

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_internal_worker_api.py`
Expected: FAIL because internal worker endpoints do not exist.

**Step 3: Write minimal implementation**

Create `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/internal_worker_api.py` with endpoints:

- `POST /internal/workers/register`
- `POST /internal/workers/claim`
- `POST /internal/tasks/{task_id}/heartbeat`
- `POST /internal/tasks/{task_id}/running`
- `POST /internal/tasks/{task_id}/complete`
- `POST /internal/tasks/{task_id}/fail`
- `POST /internal/tasks/{task_id}/cancel-ack`

Use an internal worker token from settings, not public end-user auth.

Mount the router from `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py`.

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_internal_worker_api.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/internal_worker_api.py src/vasp_server/vasp_server_api.py src/vasp_server/schemas.py tests/test_internal_worker_api.py
git commit -m "feat: add internal pull-worker control plane endpoints"
```

### Task 4: Introduce Object Storage As The Canonical Artifact Store

**Files:**
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/storage.py`
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/storage_models.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/Config.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/settings.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_artifact_storage.py`

**Step 1: Write the failing test**

```python
def test_register_artifact_uses_object_storage_key(task_manager):
    artifact = task_manager.register_artifact(...)
    assert artifact.storage_backend == "oss"
    assert artifact.storage_key.startswith("tenant/")
```

Add tests for:
- presigned upload URL generation
- presigned download URL generation
- artifact metadata persisted in DB
- public task status never returns raw local filesystem paths

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_artifact_storage.py`
Expected: FAIL because artifacts still use `storage_backend="local"` and raw `storage_key` paths.

**Step 3: Write minimal implementation**

Create `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/storage.py` with:

```python
class ObjectStorageService:
    def create_upload_url(self, object_key: str, content_type: str) -> str: ...
    def create_download_url(self, object_key: str, expires_in: int = 3600) -> str: ...
    def build_task_prefix(self, tenant_id: str, task_id: str, attempt_no: int) -> str: ...
```

Update artifact metadata in `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/models.py`:
- `bucket`
- `object_key`
- `etag`
- `sha256`
- `content_type`

Update `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/Config.py` and `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/settings.py` with object-storage settings:
- `OSS_ENDPOINT` or `S3_ENDPOINT`
- `OSS_BUCKET`
- `OSS_REGION`
- `OSS_ACCESS_KEY_ID`
- `OSS_ACCESS_KEY_SECRET`
- `ARTIFACT_URL_EXPIRE_SECONDS`

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_artifact_storage.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/storage.py src/vasp_server/storage_models.py src/vasp_server/task_manager/models.py src/vasp_server/Config.py src/vasp_server/settings.py tests/test_artifact_storage.py
git commit -m "feat: add object storage artifact backend"
```

### Task 5: Remove Public Filesystem Serving From The API

**Files:**
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_status_api.py`

**Step 1: Write the failing test**

```python
def test_task_status_returns_signed_artifact_urls_not_local_paths(client):
    data = client.get("/vasp/task/task123", params={"user_id": "u1"}).json()
    assert data["html_report_url"].startswith("https://")
    assert "/static/" not in data["html_report_url"]
```

Add tests for:
- `/download/file/*` is no longer used in task payloads
- `/static` is no longer mounted for task artifacts
- artifact list returns typed metadata and signed URLs only

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_status_api.py`
Expected: FAIL because the API still returns `/static/...` and `/download/file/...`.

**Step 3: Write minimal implementation**

Update `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py`:
- remove startup `app.mount("/static", ...)`
- deprecate `/download/file/{file_path:path}`
- build task response URLs from object storage service
- include `artifacts[].download_url` only when authorized

Update `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/schemas.py` so `ArtifactInfo` includes:
- `artifact_type`
- `content_type`
- `size_bytes`
- `download_url`

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_task_status_api.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/vasp_server_api.py src/vasp_server/schemas.py tests/test_task_status_api.py
git commit -m "refactor: serve artifacts through signed object storage urls"
```

### Task 6: Decouple VaspWorker From Shared Local Directories

**Files:**
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py`
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/input_resolver.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_input_resolver.py`

**Step 1: Write the failing test**

```python
def test_scf_uses_upstream_artifact_manifest_instead_of_local_task_dir():
    resolver = UpstreamInputResolver(...)
    paths = resolver.materialize_inputs(task_type="scf_calculation", upstream_artifacts=[...])
    assert "POSCAR" in paths
```

Add tests for:
- SCF resolves `CONTCAR` from artifact metadata instead of `self.base_work_dir / optimized_task_id`
- DOS resolves `CHGCAR`, `WAVECAR`, `POTCAR` from object storage
- fallback for formula/cif_url input still works

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_input_resolver.py`
Expected: FAIL because upstream reuse still copies from local task folders.

**Step 3: Write minimal implementation**

Create `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/input_resolver.py`:

```python
class UpstreamInputResolver:
    def resolve_for_scf(self, artifacts: list[dict], work_dir: Path) -> dict: ...
    def resolve_for_dos(self, artifacts: list[dict], work_dir: Path) -> dict: ...
```

Refactor `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py`:
- replace direct reads from:
  - `self.base_work_dir / optimized_task_id`
  - `self.base_work_dir / scf_task_id`
- accept materialized upstream artifacts in task payload or worker claim payload
- keep local scratch space only for the active attempt

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_input_resolver.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/vasp_worker.py src/vasp_server/input_resolver.py tests/test_input_resolver.py
git commit -m "refactor: resolve upstream calculation inputs from artifact manifests"
```

### Task 7: Build The HPC Pull Worker Runtime

**Files:**
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py`
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/worker_client.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_hpc_pull_worker.py`

**Step 1: Write the failing test**

```python
def test_pull_worker_claims_runs_uploads_and_completes(mock_server):
    worker = PullWorker(...)
    worker.run_once()
    assert mock_server.completed_tasks == ["task-1"]
```

Add tests for:
- periodic claim loop
- heartbeat while Slurm job is running
- fail path uploads logs then reports `failed`
- cancel path calls `scancel` and reports `canceled`

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_hpc_pull_worker.py`
Expected: FAIL because the standalone pull worker does not exist.

**Step 3: Write minimal implementation**

Create `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/worker_client.py`:

```python
class ControlPlaneClient:
    def claim(self) -> dict | None: ...
    def heartbeat(self, task_id: str, lease_token: str) -> None: ...
    def complete(self, task_id: str, lease_token: str, payload: dict) -> None: ...
    def fail(self, task_id: str, lease_token: str, payload: dict) -> None: ...
```

Create `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py`:
- poll VPS
- materialize inputs
- execute existing VASP logic locally
- request upload URLs
- upload artifacts
- report completion

Keep `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_worker.py` as the computation library, not the control-plane orchestrator.

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_hpc_pull_worker.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/hpc_pull_worker.py src/vasp_server/worker_client.py src/vasp_server/vasp_worker.py tests/test_hpc_pull_worker.py
git commit -m "feat: add hpc pull-worker runtime"
```

### Task 8: Rework Cancellation, Retry, And Lease Recovery

**Files:**
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/hpc_pull_worker.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_cancel_and_recovery.py`

**Step 1: Write the failing test**

```python
def test_cancel_requested_task_becomes_canceled_after_worker_ack():
    ...
    assert task.status == "canceled"
```

Add tests for:
- `cancel_requested` task is not marked `failed`
- worker issues `scancel` when a Slurm job id exists
- expired heartbeat returns task to `queued`
- retryable infrastructure failures increment `retry_count` and requeue

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_cancel_and_recovery.py`
Expected: FAIL because cancellation currently uses in-memory flags and falls into `failed`.

**Step 3: Write minimal implementation**

Refactor the state transitions in `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/task_manager/manager.py`:
- `request_cancel` sets `cancel_requested_at`
- worker polls cancel flag via heartbeat response or claim payload
- worker calls `scancel {scheduler_job_id}` when needed
- worker reports `cancel-ack`
- manager finalizes `status="canceled"`

Add a periodic VPS-side reconciler:
- `requeue_expired_leases()`
- `mark_orphaned_running_tasks()` only after timeout policy

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_cancel_and_recovery.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/task_manager/manager.py src/vasp_server/hpc_pull_worker.py tests/test_cancel_and_recovery.py
git commit -m "fix: make cancel and lease recovery work across distributed workers"
```

### Task 9: Migrate Public Submit And Status Endpoints To The New Internals

**Files:**
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/client.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/mcp_server.py`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_public_api_upgrade.py`

**Step 1: Write the failing test**

```python
def test_submit_task_only_creates_db_record_and_returns_task_id(client):
    resp = client.post("/vasp/scf-calculation", json=payload)
    assert resp.status_code == 200
    assert resp.json()["status"] == "queued"
```

Add tests for:
- submit endpoints no longer trigger local execution immediately
- status endpoint hydrates artifact URLs from DB-backed metadata
- MCP client still works without public API shape regressions

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_public_api_upgrade.py`
Expected: FAIL until API wiring uses the new orchestration flow.

**Step 3: Write minimal implementation**

Update `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/vasp_server/vasp_server_api.py` so public endpoints:
- validate request
- write `Task`
- return `queued`
- never call local worker execution

Update `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/client.py` and `/Users/ysl/Desktop/Code/MCP_IC_Tool/src/mcp_ic_tool/mcp_server.py` only as needed to preserve compatibility.

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_public_api_upgrade.py`
Expected: PASS

**Step 5: Commit**

```bash
git add src/vasp_server/vasp_server_api.py src/mcp_ic_tool/client.py src/mcp_ic_tool/mcp_server.py tests/test_public_api_upgrade.py
git commit -m "refactor: route public api through distributed control plane"
```

### Task 10: Add Deployment Gates, Backfill, And Cutover Controls

**Files:**
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/docs/runbooks/distributed-cutover.md`
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/scripts/backfill_artifacts.py`
- Create: `/Users/ysl/Desktop/Code/MCP_IC_Tool/scripts/requeue_stuck_tasks.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.env.example`
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_backfill_scripts.py`

**Step 1: Write the failing test**

```python
def test_backfill_artifacts_converts_local_artifact_rows_to_object_storage_rows():
    ...
    assert artifact.storage_backend == "oss"
```

Add tests for:
- backfill script can migrate legacy local paths
- requeue script only touches expired leases
- cutover flag disables old local artifact URLs

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_backfill_scripts.py`
Expected: FAIL because migration scripts do not exist.

**Step 3: Write minimal implementation**

Create:
- `/Users/ysl/Desktop/Code/MCP_IC_Tool/scripts/backfill_artifacts.py`
- `/Users/ysl/Desktop/Code/MCP_IC_Tool/scripts/requeue_stuck_tasks.py`
- `/Users/ysl/Desktop/Code/MCP_IC_Tool/docs/runbooks/distributed-cutover.md`

Document rollout phases:
- phase A: DB-backed queue, old local runner still enabled
- phase B: internal worker API enabled, single HPC worker canary
- phase C: object storage required for new tasks
- phase D: disable local executor, disable `/static` artifact serving

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/test_backfill_scripts.py`
Expected: PASS

**Step 5: Commit**

```bash
git add docs/runbooks/distributed-cutover.md scripts/backfill_artifacts.py scripts/requeue_stuck_tasks.py .env.example tests/test_backfill_scripts.py
git commit -m "docs: add rollout and backfill tooling for distributed execution"
```

### Task 11: End-To-End Verification Before Cutover

**Files:**
- Test: `/Users/ysl/Desktop/Code/MCP_IC_Tool/tests/e2e/test_distributed_flow.py`
- Modify: `/Users/ysl/Desktop/Code/MCP_IC_Tool/docs/runbooks/distributed-cutover.md`

**Step 1: Write the failing test**

```python
def test_end_to_end_structure_optimization_flow():
    """
    submit -> claim -> run -> upload -> complete -> query status
    """
```

Cover:
- structure optimization canary
- scf reusing upstream structure from object storage
- HTML report visible via signed URL
- cancel during running task

**Step 2: Run test to verify it fails**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/e2e/test_distributed_flow.py`
Expected: FAIL until all prior tasks are implemented.

**Step 3: Write minimal implementation**

Wire the test harness to a mocked object storage layer and a fake HPC worker loop. Update the cutover runbook with production acceptance criteria:
- no `queued` task older than 2 claim intervals
- no `leased` task older than lease timeout without requeue
- no public API response contains local filesystem path
- canary task type succeeds on one real HPC worker

**Step 4: Run test to verify it passes**

Run: `/Users/ysl/Desktop/Code/MCP_IC_Tool/.venv/bin/python -m pytest -q /Users/ysl/Desktop/Code/MCP_IC_Tool/tests/e2e/test_distributed_flow.py`
Expected: PASS

**Step 5: Commit**

```bash
git add tests/e2e/test_distributed_flow.py docs/runbooks/distributed-cutover.md
git commit -m "test: add end-to-end distributed execution coverage"
```

## Execution Notes

- Implement in this order. Do not start object-storage cutover before Task 2 and Task 3 are stable.
- Canary only one task type first: `structure_optimization`.
- Keep public API payloads backward-compatible while changing internals.
- Do not let HPC workers read PostgreSQL directly; they should only call internal worker APIs.
- Prefer `DATABASE_URL` over split `DB_*` variables during production rollout to avoid accidental SQLite fallback.

## Suggested Milestones

1. **Milestone A:** Database-backed queue and lease model are merged, but old local execution path can still be toggled on.
2. **Milestone B:** Internal worker API and one canary HPC pull worker are live in staging.
3. **Milestone C:** Artifacts flow through object storage and task status returns signed URLs only.
4. **Milestone D:** Legacy in-process execution and filesystem artifact serving are removed.
