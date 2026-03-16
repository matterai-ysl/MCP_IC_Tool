from typing import Any, Dict, Optional

from fastmcp import FastMCP,Context

from .client import VaspAPIClient
from .config import mcp_config
from .models import StructOptInput, SCFInput, DOSInput, MDInput


mcp = FastMCP("VASP-MCP")
client = VaspAPIClient()


# Custom INCAR parameter usage guide
CUSTOM_INCAR_HELP = """
Custom INCAR Parameter Usage Guide:

The custom_incar parameter allows agents to directly specify VASP calculation INCAR parameters, providing maximum flexibility.

Common INCAR Parameter Examples:
- Electronic structure parameters:
  * EDIFF: Electronic convergence precision (e.g.: 1e-6, 1e-7)
  * NELM: Maximum electronic steps (e.g.: 100, 200)
  * ALGO: Algorithm selection (e.g.: "Fast", "Normal", "Very_Fast")
  * ISMEAR: Smearing method (e.g.: 0, -5)
  * SIGMA: Smearing parameter (e.g.: 0.05, 0.1)

- DOS calculation specific:
  * LORBIT: Orbital projection (e.g.: 11, 12)
  * NEDOS: Energy grid points (e.g.: 2000, 3000)
  * EMIN/EMAX: Energy range (e.g.: -20, 10)

- MD calculation specific:
  * SMASS: Thermostat mass (e.g.: 0, 1)
  * POTIM: Time step (e.g.: 0.5, 1.0)
  * ISYM: Symmetry (e.g.: 0 off, 1 on)
  * MDALGO: MD algorithm (e.g.: 1,2,3)

- Optimization calculation:
  * EDIFFG: Force convergence precision (e.g.: -0.01, -0.02)
  * IBRION: Ion movement method (e.g.: 1,2,3)
  * ISIF: Stress tensor (e.g.: 2,3,7)

Usage:
The custom_incar parameter accepts a dictionary format, with keys as INCAR parameter names and values as parameter values.
Parameter names are case-insensitive and will be automatically converted to uppercase.
Custom parameters will override default parameters with the same name.

Examples:
{"custom_incar": {"EDIFF": 1e-7, "NELM": 100, "ALGO": "Fast"}}
{"custom_incar": {"LORBIT": 11, "NEDOS": 3000, "ISMEAR": -5}}
{"custom_incar": {"SMASS": 1, "POTIM": 0.5, "ISYM": 0}}
"""
def get_user_id(ctx: Context) -> str:
    if ctx is not None:
        user_id = ctx.request_context.request.headers.get("user_id", None)  # type: ignore
    else:
        user_id = "123"
    return user_id #type: ignore

@mcp.tool()
async def submit_structure_optimization(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a structure optimization task and return task information.

    Parameters:
    - calc_type (required): Calculation type, options:
      * "OXC" - Oxide/sulfide solid electrolyte
      * "SSE" - Solid electrolyte (equivalent to OXC)
      * "ORC" - Oxide reduction catalyst
      * "ECAT_OER" - Oxygen evolution reaction catalyst
      * "ECAT_HER" - Hydrogen evolution reaction catalyst
    - formula (optional): Chemical formula, e.g. 'Li2O', 'LiFePO4', choose one between formula and cif_url
    - cif_url (optional): URL address of CIF file, choose one between formula and cif_url
    - spacegroup (optional): Space group symbol, e.g. "P1", "Fm-3m", only valid when using formula, used for structure filtering
    - max_energy_above_hull (optional): Maximum energy above hull (eV/atom), default 0.1, only valid when using formula, used for structure filtering
    - min_band_gap (optional): Minimum band gap (eV), only valid when using formula, used for structure filtering
    - max_band_gap (optional): Maximum band gap (eV), only valid when using formula, used for structure filtering
    - max_nsites (optional): Maximum number of atoms, only valid when using formula, used for structure filtering
    - min_nsites (optional): Minimum number of atoms, only valid when using formula, used for structure filtering
    - stable_only (optional): Only select stable materials, default true, only valid when using formula, used for structure filtering
    - selection_mode (optional): Selection mode, options "auto"/"stable"/"first", only valid when using formula, used for structure filtering
    - kpoint_density (optional): K-point density parameter, default 30.0
    - custom_incar (optional): Custom INCAR parameter dictionary for overriding or adding specific VASP calculation parameters

    Examples:
    submit_structure_optimization(calc_type="OXC", formula="Li2O", kpoint_density=25.0)
    submit_structure_optimization(calc_type="OXC", formula="Li2O", custom_incar={"EDIFF": 1e-7, "NELM": 100, "ALGO": "Fast"})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    print("submit_structure_optimization")
    print(params)
    payload = StructOptInput(**params).model_dump(mode="json", by_alias=True)
    if payload.get("user_id") is None:
        payload["user_id"] = "123"
    return await client.submit_structure_optimization(payload)


@mcp.tool()
async def submit_scf_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    optimized_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore

) -> Dict[str, Any]:
    """
    Submit a self-consistent field calculation task and return task information.

    Parameters:
    - calc_type (required): Calculation type, same as structure optimization
    - formula (optional): Chemical formula, e.g. 'Li2O', 'LiFePO4', choose one among formula, cif_url, and optimized_task_id
    - cif_url (optional): URL address of CIF file, choose one among formula, cif_url, and optimized_task_id
    - optimized_task_id (optional): Completed structure optimization task ID, choose one among formula, cif_url, and optimized_task_id
    - spacegroup (optional): Space group symbol, only valid when using formula
    - max_energy_above_hull (optional): Maximum energy above hull (eV/atom), default 0.1, only valid when using formula, used for structure filtering
    - min_band_gap (optional): Minimum band gap (eV), only valid when using formula, used for structure filtering
    - max_band_gap (optional): Maximum band gap (eV), only valid when using formula, used for structure filtering
    - max_nsites (optional): Maximum number of atoms, only valid when using formula, used for structure filtering
    - min_nsites (optional): Minimum number of atoms, only valid when using formula, used for structure filtering
    - stable_only (optional): Only select stable materials, default true, only valid when using formula, used for structure filtering
    - selection_mode (optional): Selection mode, options "auto"/"stable"/"first", only valid when using formula, used for structure filtering
    - kpoint_density (optional): K-point density parameter, default 30.0
    - precision (optional): Calculation precision, options "Normal"/"High"/"Accurate", default "Accurate"
    - custom_incar (optional): Custom INCAR parameter dictionary for overriding or adding specific VASP calculation parameters

    Examples:
    submit_scf_calculation(user_id="user123", calc_type="OXC", optimized_task_id="task_001", precision="High")
    submit_scf_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"ISMEAR": 0, "SIGMA": 0.05})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = SCFInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_scf(payload)


@mcp.tool()
async def submit_dos_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    kpoint_multiplier: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a density of states calculation task and return task information.

    Parameters:
    - calc_type (required): Calculation type, same as structure optimization
    - formula (optional): Chemical formula, e.g. 'Li2O', 'LiFePO4', choose one among formula, cif_url, and scf_task_id
    - cif_url (optional): URL address of CIF file, choose one among formula, cif_url, and scf_task_id
    - scf_task_id (optional): Completed self-consistent field calculation task ID, choose one among formula, cif_url, and scf_task_id
    - spacegroup (optional): Space group symbol, only valid when using formula
    - max_energy_above_hull (optional): Maximum energy above hull (eV/atom), default 0.1
    - min_band_gap (optional): Minimum band gap (eV)
    - max_band_gap (optional): Maximum band gap (eV)
    - max_nsites (optional): Maximum number of atoms
    - min_nsites (optional): Minimum number of atoms
    - stable_only (optional): Only select stable materials, default true
    - selection_mode (optional): Selection mode, options "auto"/"stable"/"first"
    - kpoint_density (optional): K-point density parameter, default 30.0
    - kpoint_multiplier (optional): K-point multiplier factor relative to optimization calculation, default 2.0
    - precision (optional): Calculation precision, options "Normal"/"High"/"Accurate", default "Accurate"
    - custom_incar (optional): Custom INCAR parameter dictionary for overriding or adding specific VASP calculation parameters

    Examples:
    submit_dos_calculation(user_id="user123", calc_type="OXC", scf_task_id="scf_001", kpoint_multiplier=3.0)
    submit_dos_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"LORBIT": 11, "NEDOS": 3000})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = DOSInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_dos(payload)


@mcp.tool()
async def submit_md_calculation(
    calc_type: str,
    formula: Optional[str] = None,
    cif_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    spacegroup: Optional[str] = None,
    max_energy_above_hull: Optional[float] = None,
    min_band_gap: Optional[float] = None,
    max_band_gap: Optional[float] = None,
    max_nsites: Optional[int] = None,
    min_nsites: Optional[int] = None,
    stable_only: Optional[bool] = None,
    selection_mode: Optional[str] = None,
    md_steps: Optional[int] = None,
    temperature: Optional[Any] = None,  # Can be float or List[float]
    time_step: Optional[float] = None,
    ensemble: Optional[str] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a molecular dynamics calculation task and return task information.

    Parameters:
    - calc_type (required): Calculation type, same as structure optimization
    - formula (optional): Chemical formula, e.g. 'Li2O', 'LiFePO4', choose one among formula, cif_url, and scf_task_id
    - cif_url (optional): URL address of CIF file, choose one among formula, cif_url, and scf_task_id
    - scf_task_id (optional): Completed self-consistent field calculation task ID, choose one among formula, cif_url, and scf_task_id
    - spacegroup (optional): Space group symbol, only valid when using formula
    - max_energy_above_hull (optional): Maximum energy above hull (eV/atom), default 0.1, only valid when using formula, used for structure filtering
    - min_band_gap (optional): Minimum band gap (eV), only valid when using formula, used for structure filtering
    - max_band_gap (optional): Maximum band gap (eV), only valid when using formula, used for structure filtering
    - max_nsites (optional): Maximum number of atoms, only valid when using formula, used for structure filtering
    - min_nsites (optional): Minimum number of atoms, only valid when using formula, used for structure filtering
    - stable_only (optional): Only select stable materials, default true, only valid when using formula, used for structure filtering
    - selection_mode (optional): Selection mode, options "auto"/"stable"/"first", only valid when using formula, used for structure filtering
    - md_steps (optional): Number of MD steps, default 1000
    - temperature (optional): Target temperature (K), supports single temperature or temperature list for multi-temperature scan
      * Single temperature: 300.0
      * Multi-temperature: [200.0, 300.0, 400.0, 500.0] - will create subtasks for each temperature
    - time_step (optional): Time step (fs), default 1.0
    - ensemble (optional): Ensemble type, options "NVT"/"NVE"/"NPT", default "NVT"
    - precision (optional): Calculation precision, options "Normal"/"High"/"Accurate", default "Normal"
    - custom_incar (optional): Custom INCAR parameter dictionary for overriding or adding specific VASP calculation parameters

    Multi-temperature MD calculation features:
    - Create independent subdirectories for each temperature point (e.g.: T_300K/, T_400K/)
    - All subtasks share the same task ID for easy management
    - Support independent success/failure status for subtasks
    - Automatically generate multi-temperature summary analysis report

    Examples:
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", md_steps=2000, temperature=350.0, ensemble="NPT")
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", temperature=[300.0, 400.0, 500.0], md_steps=1500)
    submit_md_calculation(user_id="user123", calc_type="OXC", formula="Li2O", custom_incar={"SMASS": 1, "POTIM": 0.5, "ISYM": 0})
    """
    params = {k: v for k, v in locals().items() if v is not None}
    params["user_id"] = get_user_id(ctx)
    payload = MDInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_md(payload)


@mcp.tool()
async def get_task_status(task_id: str,
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Query task status and result information.

    Parameters:
    - task_id (required): Task ID, returned when submitting the task

    Return information contains:
    - status: Task status ("queued"/"running"/"completed"/"failed"/"canceled")
    - progress: Progress percentage (0-100)
    - result_path: Result file path (when completed)
    - error_message: Error message (when failed)
    - result_data: Detailed result data (includes all calculation results)

    For multi-temperature MD tasks, result_data also includes:
    - multi_temperature_info: Multi-temperature subtask status summary
      * is_multi_temperature: Whether it is a multi-temperature calculation
      * total_subtasks: Total number of subtasks
      * completed_subtasks: Number of completed subtasks
      * failed_subtasks: Number of failed subtasks
      * subtask_status: Detailed status list for each temperature point

    Examples:
    get_task_status("task_001")
    """
    return await client.get_task_status(task_id, get_user_id(ctx))


@mcp.tool()
async def cancel_task(task_id: str, 
ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Cancel a running or queued task.

    Parameters:
    - task_id (required): Task ID to cancel

    Examples:
    cancel_task("task_001", "user123")
    """
    return await client.cancel_task(task_id, get_user_id(ctx))


@mcp.tool()
async def list_user_tasks(
ctx: Context = None #type: ignore
) -> Any:
    """
    List all tasks for the user.

    Parameters:

    Returns a task list, each task contains:
    - task_id: Task ID
    - task_type: Task type ("structure_optimization"/"scf_calculation"/"dos_calculation"/"md_calculation")
    - status: Task status
    - created_at: Creation time
    - updated_at: Update time

    Examples:
    list_user_tasks()
    """
    return await client.list_tasks(get_user_id(ctx))


@mcp.tool()
async def analyze_structure_optimization_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze structure optimization results independently — without submitting a new calculation.

    Provide one of:
    - task_id: ID of an existing completed structure optimization task
    - file_url: URL to a zip/tar.gz archive (.zip or .tar.gz) containing VASP output files

    Required files in the archive:
    - OUTCAR  (required — contains all convergence and energy data)
    Optional files (improve report quality if present):
    - POSCAR  (initial structure)
    - CONTCAR (optimized structure)
    Do NOT include CHGCAR/WAVECAR — they are large and not needed.

    Returns analysis summary (force_convergence, final_energy, etc.) and HTML report URL.

    Examples:
    analyze_structure_optimization_results(task_id="abc123")
    analyze_structure_optimization_results(file_url="https://example.com/opt_output.zip")
    """
    payload: Dict[str, Any] = {}
    if task_id:
        payload["task_id"] = task_id
        payload["user_id"] = get_user_id(ctx)
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_optimization(payload)


@mcp.tool()
async def analyze_scf_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze self-consistent field (SCF) calculation results independently.

    Provide one of:
    - task_id: ID of an existing completed SCF task
    - file_url: URL to a zip/tar.gz archive (.zip or .tar.gz) containing VASP output files

    Required files in the archive:
    - OUTCAR  (required — contains energy, Fermi level, band gap, convergence info)
    Optional files for Bader charge analysis (all three must be present together, or none):
    - POSCAR  (atom type and count info)
    - POTCAR  (valence electron counts per species; note: may have license restrictions)
    - ACF.dat (Bader charge output, generated by the external bader program — small file)
    If any of the three above is missing, Bader analysis is silently skipped.
    Do NOT include CHGCAR/WAVECAR/ELFCAR/DOSCAR — they can be hundreds of MB to GB and are not read.

    Returns analysis summary (convergence, total_energy, fermi_energy, band_gap) and HTML report URL.

    Examples:
    analyze_scf_results(task_id="abc123")
    analyze_scf_results(file_url="https://example.com/scf_output.zip")
    """
    payload: Dict[str, Any] = {}
    if task_id:
        payload["task_id"] = task_id
        payload["user_id"] = get_user_id(ctx)
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_scf(payload)


@mcp.tool()
async def analyze_dos_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze density of states (DOS) calculation results independently.

    Provide one of:
    - task_id: ID of an existing completed DOS task
    - file_url: URL to a zip/tar.gz archive (.zip or .tar.gz) containing VASP output files

    Required files in the archive:
    - vasprun.xml  (required — this analyzer is entirely based on pymatgen parsing vasprun.xml)
    Optional files:
    - OUTCAR  (optional — used for supplementary magnetic moment data)
    Do NOT include CHGCAR/WAVECAR/DOSCAR — not used by this analyzer.
    Note: vasprun.xml can itself be large (tens to hundreds of MB for dense k-point grids);
          only include it when file size is acceptable.

    Returns analysis summary (band_gap, fermi_energy, is_metal) and HTML report URL.

    Examples:
    analyze_dos_results(task_id="abc123")
    analyze_dos_results(file_url="https://example.com/dos_output.zip")
    """
    payload: Dict[str, Any] = {}
    if task_id:
        payload["task_id"] = task_id
        payload["user_id"] = get_user_id(ctx)
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_dos(payload)


@mcp.tool()
async def analyze_md_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze molecular dynamics (MD) calculation results independently.

    Provide one of:
    - task_id: ID of an existing completed MD task
    - file_url: URL to a zip/tar.gz archive (.zip or .tar.gz) containing VASP output files

    Required files in the archive:
    - XDATCAR  (required — atomic trajectory, used for MSD/diffusion/RDF analysis)
    Optional files (add more physical quantities if present):
    - OUTCAR       (temperature series, pressure, energy per step)
    - vasprun.xml  (supplementary: ionic steps, time step POTIM)
    - INCAR        (time step POTIM, used to convert steps to real time)
    Do NOT include CHGCAR/WAVECAR — not used.
    For multi-temperature MD, include the T_*K/ subdirectory structure inside the archive.

    Returns analysis summary and HTML report URL.
    Supports both single-temperature and multi-temperature MD results.

    Examples:
    analyze_md_results(task_id="abc123")
    analyze_md_results(file_url="https://example.com/md_output.zip")
    """
    payload: Dict[str, Any] = {}
    if task_id:
        payload["task_id"] = task_id
        payload["user_id"] = get_user_id(ctx)
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_md(payload)


@mcp.tool()
async def get_custom_incar_help() -> str:
    """
    Get detailed usage instructions and common parameter examples for custom INCAR parameters.

    This tool provides:
    - Usage methods for the custom_incar parameter
    - Detailed explanations of common INCAR parameters
    - Recommended parameters for different calculation types
    - Practical usage examples

    Calling this tool can help agents understand how to correctly use the custom_incar parameter to customize VASP calculations.
    """
    return CUSTOM_INCAR_HELP


# def run():
#     # 运行 FastMCP 服务器
#     mcp.run(host=mcp_config.host, port=mcp_config.port)


# if __name__ == "__main__":
#     run()


