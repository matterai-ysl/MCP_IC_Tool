# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Architecture Overview

This is a VASP (Vienna Ab initio Simulation Package) calculation service with a dual architecture:

1. **VASP Server API** (`src/vasp_server/`) - FastAPI backend that manages VASP computational tasks
2. **MCP Integration** (`src/mcp_ic_tool/`) - Model Context Protocol server for LLM agent integration

### Key Components

- **Task Management**: Centralized task queue system with SQLite database (`tasks.db`)
- **Analysis Modules**: Specialized analyzers for DOS, band gap, MD, SCF, and structure optimization
- **Configuration**: Centralized VASP parameter templates and paths (`src/vasp_server/Config.py`)

## Running the Services

### VASP Server API
```bash
# Run the main VASP computation server
python -m uvicorn src.vasp_server.vasp_server_api:app --host 0.0.0.0 --port 8000

# Single port development server
python single_port_server.py
```

### MCP Server
```bash
# The MCP server provides LLM agent integration
python -m src.mcp_ic_tool.mcp_server
```

## Development Commands

### Dependencies
```bash
# Install Python dependencies
pip install -r requirements.txt
```

### Database Management
```bash
# Inspect the task database
python inspect_tasks_db.py

# Debug VASP API calls
python debug_vasp_api.py
```

## Project Structure

- `src/vasp_server/` - Core VASP computational backend
  - `vasp_server_api.py` - Main FastAPI application
  - `task_manager/` - Task queue and database management
  - `*_analyzer.py` - Specialized computation modules (DOS, MD, SCF, etc.)
  - `Config.py` - VASP parameter templates and system paths
  - `schemas.py` - Pydantic models for API requests/responses

- `src/mcp_ic_tool/` - MCP integration layer
  - `mcp_server.py` - MCP server implementation
  - `client.py` - HTTP client for VASP API calls
  - `models.py` - Input validation models

- Test directories (`dos_test/`, `scf_test/`, `md_test/`) - Contain calculation examples and test cases

## Key Configuration

The system uses hardcoded paths for VASP binaries and pseudopotentials defined in `src/vasp_server/Config.py`:
- VASP executable: `/data/app/vasp/6.3.2-intel/bin/vasp_std`
- Pseudopotential path: `/data/home/ysl9527/software/psudopotential`

INCAR templates are provided for different calculation types (OXC, ORC, MD) with sensible defaults that can be customized via the `custom_incar` parameter in API requests.

## Task Types Supported

1. **Structure Optimization** - Geometric optimization with VASP
2. **SCF Calculations** - Self-consistent field calculations
3. **DOS Calculations** - Density of states analysis
4. **MD Simulations** - Molecular dynamics simulations

Each calculation type has its own analyzer module and predefined INCAR templates.