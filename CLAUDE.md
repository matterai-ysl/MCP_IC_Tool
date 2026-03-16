# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture Overview

This is a VASP (Vienna Ab initio Simulation Package) computational chemistry service with two independently deployable layers:

**Layer 1 – VASP Server** (`src/vasp_server/`): A FastAPI backend that orchestrates VASP calculations. It accepts HTTP requests, manages a task queue backed by a PostgreSQL cloud database (Aliyun RDS), and spawns each calculation in its own daemon thread with a new asyncio event loop.

**Layer 2 – MCP Integration** (`src/mcp_ic_tool/`): A FastMCP server that exposes VASP workflows as LLM-callable tools. It is a thin HTTP client over Layer 1 — the MCP tools translate agent requests into REST calls to the VASP Server API.

### Data Flow

```
LLM Agent → MCP Server (port 8130, /mcp) → VASP Server API (port 8140) → VaspWorker thread → VASP binary on HPC
```

The MCP server and VASP server can run on different machines. `src/mcp_ic_tool/config.py` holds `VASP_SERVER_BASE_URL` pointing to the remote VASP API endpoint.

### Task Lifecycle (v2-lite)

1. API endpoint receives request → `TaskManager.submit_task()` writes a `Task` row (status=`queued`, analysis_status=`pending`)
2. A daemon thread starts immediately → creates `ExecutionAttempt` record → runs `VaspWorker` via fresh `asyncio.new_event_loop()`
3. `VaspWorker` fetches CIF, converts to POSCAR, generates VASP inputs, runs VASP binary via SLURM
4. Progress callbacks update `Task.progress` and `Task.progress_message` (NOT `error_message`)
5. On execution success → `ExecutionAttempt.status=succeeded` → auto-triggers `AnalysisRun`
6. Analysis extracts summary → registers `Artifact` records → writes `Task.result_summary` + `html_report_url`
7. Final state: `Task.status=completed`, `Task.analysis_status=completed`
8. If analysis fails: `Task.status=completed`, `Task.analysis_status=failed` — execution record preserved

### Internal Data Model (v2-lite)

- **Task** — primary object, backward-compatible. New fields: `analysis_status`, `result_summary`, `progress_message`, `input_snapshot`
- **ExecutionAttempt** — records each VASP run attempt (status: submitting→submitted→running→succeeded/runtime_failed/canceled)
- **AnalysisRun** — independent analysis phase (status: pending→running→completed/failed). Contains `summary` JSON
- **Artifact** — indexes all file outputs (OUTCAR, CONTCAR, HTML reports, etc.) with `storage_backend=local` and `storage_key`

### Calculation Types & Analyzers

Each task type has a dedicated analyzer module:
- `optimization_analyzer.py` – structure relaxation
- `scf_analyzer.py` – self-consistent field
- `dos_analyzer.py` – density of states
- `md_analyzer.py` – molecular dynamics (supports multi-temperature scan via `T_*K/` subdirectories)
- `bandgap_analyzer.py` – band gap extraction

`calc_type` parameter maps to INCAR templates defined in `Config.py`: `OXC` (oxide/sulfide), `ORC` (oxide reduction catalyst + vdW correction), `SSE`/`ECAT_OER`/`ECAT_HER` use OXC template. MD has its own separate template.

### Key Configuration

All server-side paths and INCAR templates are in `src/vasp_server/Config.py`:
- VASP binary: `/data/app/vasp/6.3.2-intel/bin/vasp_std`
- Pseudopotentials: `/data/home/ysl9527/software/psudopotential`
- Calculation output root: `/data/home/ysl9527/vasp_calculations/{user_id}/{task_id}/`

MCP client endpoint config is in `src/mcp_ic_tool/config.py` (`VASP_SERVER_BASE_URL`).

Database credentials are hardcoded in `src/vasp_server/task_manager/database.py` — override with env vars `DB_HOST`, `DB_PORT`, `DB_NAME`, `DB_USER`, `DB_PASSWORD`.

## Running the Services

### MCP Server (primary entry point for LLM agents)
```bash
# Single-port server: MCP at /mcp + file download service
python single_port_server.py                     # default port 8130
python single_port_server.py --port 9000 --host 0.0.0.0
```

### VASP Server API (backend computation server, typically on HPC)
```bash
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8140
```

### Install Dependencies
```bash
pip install -r requirements.txt
```

### Database Utilities
```bash
python inspect_tasks_db.py      # inspect task records
python debug_vasp_api.py        # debug API calls
```

## MCP Tools Exposed to Agents

- `submit_structure_optimization` – structure relaxation (requires `formula` XOR `cif_url`)
- `submit_scf_calculation` – SCF (requires one of `formula`/`cif_url`/`optimized_task_id`)
- `submit_dos_calculation` – DOS (requires one of `formula`/`cif_url`/`scf_task_id`)
- `submit_md_calculation` – MD with optional multi-temperature list (requires one of `formula`/`cif_url`/`scf_task_id`)
- `get_task_status` / `cancel_task` / `list_user_tasks` – task management
- `get_custom_incar_help` – returns INCAR parameter guide

`user_id` is read from the HTTP request header `user_id` by `get_user_id(ctx)` and defaults to `"123"` if absent.

## Adding a New Calculation Type

1. Create `src/vasp_server/{name}_analyzer.py` with analysis logic
2. Add INCAR template in `Config.py` → `base_incars` dict
3. Add Pydantic request/response schemas in `schemas.py`
4. Add API endpoint in `vasp_server_api.py`
5. Add `VaspWorker.run_{name}_calculation()` method in `vasp_worker.py`
6. Register in `manager.py` → `TASK_TYPE_TO_METHOD` and `TASK_TYPE_TO_ANALYSIS` dicts (no more if/elif chains)
7. Add MCP tool in `mcp_server.py` and input model in `models.py`

## Running Tests

```bash
uv venv .venv --python 3.12
uv pip install pytest pytest-asyncio sqlalchemy pydantic fastapi httpx
.venv/bin/python -m pytest tests/ -v
```

Tests use in-memory SQLite and mock VaspWorker — no PostgreSQL or VASP binary needed.
