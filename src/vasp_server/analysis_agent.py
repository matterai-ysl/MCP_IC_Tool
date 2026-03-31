"""
VASP 数据分析 Agent
基于 LangGraph create_react_agent + 持久化 Python REPL + skills 库
"""
import importlib
import io
import logging
import os
import signal
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path
from typing import Any, Dict, List

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .analysis_skills import (
    extract_final_energy, extract_energy_per_atom, get_ionic_steps_energy,
    get_max_force, get_computation_time,
    get_final_structure, get_volume, get_lattice_params, get_formula,
    get_n_atoms, get_volume_change,
    is_converged, get_fermi_energy, get_band_gap,
    get_total_energy_from_vasprun, get_electronic_steps, get_magnetization,
    plot_energy_convergence, plot_dos, plot_band_structure, plot_forces_convergence,
)

logger = logging.getLogger(__name__)

_DEFAULT_MODEL_NAME = "qwen3.5-plus"

# ──────────────────────────────────────────────────────────────────────────────
#  System prompt
# ──────────────────────────────────────────────────────────────────────────────
_SYSTEM_PROMPT = """\
You are a VASP computational chemistry expert. Analyze VASP calculation outputs and answer questions precisely.

The work directory path is stored in the variable `work_dir`. Always start by calling list_files().

Pre-loaded variables (no import needed):
  work_dir  — path to VASP calculation directory (str)
  np        — numpy
  plt       — matplotlib.pyplot
  Path      — pathlib.Path

Pre-loaded skill functions:
  ── Energy ──────────────────────────────────────────────────────────────────
  extract_final_energy(work_dir)      → float | None   # final TOTEN from OUTCAR (eV)
  extract_energy_per_atom(work_dir)   → float | None   # eV/atom
  get_ionic_steps_energy(work_dir)    → List[float]    # energy per ionic step
  get_max_force(work_dir)             → float | None   # max |F| in final step (eV/Å)
  get_computation_time(work_dir)      → float | None   # total CPU time (s)

  ── Structure ────────────────────────────────────────────────────────────────
  get_final_structure(work_dir)       → pymatgen Structure | None
  get_volume(work_dir)                → float | None   # cell volume (Å³)
  get_lattice_params(work_dir)        → dict | None    # a,b,c,α,β,γ,volume
  get_formula(work_dir)               → str | None     # reduced chemical formula
  get_n_atoms(work_dir)               → int | None
  get_volume_change(work_dir)         → dict | None    # initial/final/change_pct

  ── Electronic structure ─────────────────────────────────────────────────────
  is_converged(work_dir)              → bool
  get_fermi_energy(work_dir)          → float | None   # eV
  get_band_gap(work_dir)              → dict | None    # energy_ev,direct,vbm_ev,cbm_ev
  get_total_energy_from_vasprun(work_dir) → float | None
  get_electronic_steps(work_dir)      → dict | None    # elec step count, convergence
  get_magnetization(work_dir)         → dict | None    # total magnetization

  ── Plots (return file path string) ──────────────────────────────────────────
  plot_energy_convergence(work_dir)   → str   # saves energy_convergence.png
  plot_forces_convergence(work_dir)   → str   # saves forces_convergence.png
  plot_dos(work_dir)                  → str   # saves dos_plot.png
  plot_band_structure(work_dir)       → str   # saves band_structure.png

Common VASP output files: OUTCAR, vasprun.xml, CONTCAR, POSCAR, DOSCAR, EIGENVAL, KPOINTS

Rules:
1. Call list_files() first to see what exists.
2. Use skill functions for standard analysis; write custom code for other needs.
3. Gracefully handle missing files — check existence before parsing.
4. For large files (CHGCAR, WAVECAR), never read them directly; use pymatgen APIs.
5. End your final answer with a JSON block summarising key numeric results, e.g.:
   ```json
   {"band_gap_eV": 1.23, "is_direct": false, "final_energy_eV": -123.45}
   ```
"""


# ──────────────────────────────────────────────────────────────────────────────
#  Persistent Python REPL with timeout
# ──────────────────────────────────────────────────────────────────────────────
class _SafeREPL:
    def __init__(self, initial_globals: Dict[str, Any], timeout: int = 60):
        self.timeout = timeout
        self.globals: Dict[str, Any] = {"__builtins__": __builtins__, **initial_globals}

    def run(self, code: str) -> str:
        out_buf = io.StringIO()
        err_buf = io.StringIO()

        def _alarm(signum, frame):
            raise TimeoutError(f"Code execution timed out ({self.timeout}s)")

        signal.signal(signal.SIGALRM, _alarm)
        signal.alarm(self.timeout)
        try:
            with redirect_stdout(out_buf), redirect_stderr(err_buf):
                exec(compile(code, "<agent>", "exec"), self.globals)  # noqa: S102
            output = out_buf.getvalue()
            err = err_buf.getvalue()
            if err:
                output += f"\n[stderr] {err[:400]}"
            return output or "(no output)"
        except TimeoutError as exc:
            return f"Error: {exc}"
        except Exception as exc:
            import traceback
            tb = traceback.format_exc()
            return f"Error: {type(exc).__name__}: {exc}\n{tb[-600:]}"
        finally:
            signal.alarm(0)


# ──────────────────────────────────────────────────────────────────────────────
#  Tool factory
# ──────────────────────────────────────────────────────────────────────────────
def _make_tools(work_dir: str) -> list:
    lc_tool = importlib.import_module("langchain_core.tools").tool
    repl = _SafeREPL(
        initial_globals={
            "work_dir": work_dir,
            "np": np,
            "plt": plt,
            "Path": Path,
            # ── energy ──
            "extract_final_energy": extract_final_energy,
            "extract_energy_per_atom": extract_energy_per_atom,
            "get_ionic_steps_energy": get_ionic_steps_energy,
            "get_max_force": get_max_force,
            "get_computation_time": get_computation_time,
            # ── structure ──
            "get_final_structure": get_final_structure,
            "get_volume": get_volume,
            "get_lattice_params": get_lattice_params,
            "get_formula": get_formula,
            "get_n_atoms": get_n_atoms,
            "get_volume_change": get_volume_change,
            # ── electronic ──
            "is_converged": is_converged,
            "get_fermi_energy": get_fermi_energy,
            "get_band_gap": get_band_gap,
            "get_total_energy_from_vasprun": get_total_energy_from_vasprun,
            "get_electronic_steps": get_electronic_steps,
            "get_magnetization": get_magnetization,
            # ── plot ──
            "plot_energy_convergence": plot_energy_convergence,
            "plot_forces_convergence": plot_forces_convergence,
            "plot_dos": plot_dos,
            "plot_band_structure": plot_band_structure,
        }
    )

    @lc_tool
    def execute_python(code: str) -> str:
        """Execute Python code in the VASP work directory.
        Variables persist between calls. work_dir, np, plt, Path and all
        skill functions are pre-loaded — use them directly without import."""
        return repl.run(code)

    @lc_tool
    def list_files() -> str:
        """List all files in the VASP work directory."""
        files = sorted(f.name for f in Path(work_dir).iterdir() if f.is_file())
        return "\n".join(files) if files else "(empty)"

    return [execute_python, list_files]


def _create_chat_model(model_name: str):
    normalized = model_name.strip()
    lowered = normalized.lower()

    if lowered.startswith("qwen"):
        api_key = os.environ.get("QWEN_API_KEY")
        base_url = os.environ.get("QWEN_BASE_URL")
        if not api_key:
            raise ValueError("QWEN_API_KEY 未设置，无法使用 Qwen 分析模型")
        if not base_url:
            raise ValueError("QWEN_BASE_URL 未设置，无法使用 Qwen 分析模型")

        chat_openai = importlib.import_module("langchain_openai").ChatOpenAI
        return chat_openai(
            model=normalized,
            temperature=0,
            api_key=api_key,
            base_url=base_url,
            max_tokens=4096,
        )

    api_key = os.environ.get("ANTHROPIC_API_KEY")
    if not api_key:
        raise ValueError("ANTHROPIC_API_KEY 未设置，无法使用 Claude 分析模型")

    chat_anthropic = importlib.import_module("langchain_anthropic").ChatAnthropic
    return chat_anthropic(
        model=normalized,
        temperature=0,
        api_key=api_key,
        max_tokens=4096,
    )


# ──────────────────────────────────────────────────────────────────────────────
#  Public API
# ──────────────────────────────────────────────────────────────────────────────
async def run_analysis(
    work_dir: str,
    question: str,
    model_name: str = _DEFAULT_MODEL_NAME,
) -> Dict[str, Any]:
    """
    在 work_dir 上运行 VASP 分析 Agent，返回分析结果。

    Returns:
        {
            "answer": str,             # Agent 的最终回答（含 JSON 摘要）
            "steps": int,              # 代码执行次数
            "generated_plots": list,   # 生成的 PNG 文件名列表
        }
    """
    tools = _make_tools(work_dir)
    model = _create_chat_model(model_name)
    create_react_agent = importlib.import_module("langgraph.prebuilt").create_react_agent
    agent = create_react_agent(model, tools, prompt=_SYSTEM_PROMPT)

    # Collect PNG files before the run to detect newly generated ones
    before = {f.name for f in Path(work_dir).glob("*.png")}

    result = await agent.ainvoke({
        "messages": [{"role": "user", "content": question}]
    })

    final_msg = result["messages"][-1].content

    # Count tool messages (= code executions)
    steps = sum(
        1 for m in result["messages"]
        if getattr(m, "type", None) == "tool"
    )

    after = {f.name for f in Path(work_dir).glob("*.png")}
    new_plots: List[str] = sorted(after - before)

    logger.info("Analysis agent finished: steps=%d new_plots=%s", steps, new_plots)
    return {
        "answer": final_msg,
        "steps": steps,
        "generated_plots": new_plots,
    }
