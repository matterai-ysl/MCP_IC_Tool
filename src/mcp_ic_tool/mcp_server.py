import os
from typing import Any, Dict, List, Optional, Union

from fastmcp import FastMCP, Context

from .client import VaspAPIClient
from .models import StructOptInput, SCFInput, DOSInput, BandStructureInput, MDInput, NEBInput, PhononInput, CustomCalcInput


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
    user_id: str | None = None
    if ctx is not None:
        user_id = ctx.request_context.request.headers.get("user_id", None)  # type: ignore
    if not user_id:
        user_id = os.getenv("DEFAULT_USER_ID")
    if not user_id:
        raise ValueError("user_id is required: set the 'user_id' HTTP header or DEFAULT_USER_ID env var")
    return user_id

@mcp.tool()
async def submit_structure_optimization(
    cif_url: Optional[str] = None,
    queue_name: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a structure optimization task and return task information.

    Default INCAR parameters (override any of them via custom_incar):
      PREC=Normal, ENCUT=450, EDIFF=1E-4, EDIFFG=-0.02, NELM=100, NELMIN=2,
      ALGO=Normal, ISMEAR=0, SIGMA=0.05, IBRION=2, NSW=500, ISIF=2,
      POTIM=0.2, ISPIN=2, GGA=PE, LREAL=Auto, LWAVE=.FALSE., LCHARG=.FALSE.
    LDA+U 和 MAGMOM 会根据结构中的元素自动补齐。

    Parameters:
    - cif_url (required): URL address of the structure file
    - kpoint_density (optional): K-point density parameter, default 15.0
    - custom_incar (optional): Dict to override/add INCAR parameters.
      Examples: {"EDIFF": 1e-7, "ENCUT": 600}, {"IVDW": 12} to enable vdW correction

    Examples:
    submit_structure_optimization(cif_url="https://example.com/structure.cif")
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = StructOptInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_structure_optimization(payload)


@mcp.tool()
async def submit_scf_calculation(
    cif_url: Optional[str] = None,
    optimized_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a self-consistent field (SCF) calculation task and return task information.

    Default INCAR base is the same as optimization (PREC=Normal, ENCUT=450, etc.),
    then modified for SCF: NSW=0, IBRION=-1, LWAVE=.TRUE., LCHARG=.TRUE.
    Override any parameter via custom_incar.

    Parameters:
    - cif_url (optional): URL of structure file
    - optimized_task_id (optional): Completed structure optimization task ID
    - kpoint_density (optional): default 30.0
    - precision (optional): "Normal"/"High"/"Accurate", default "Accurate"
    - custom_incar (optional): Dict to override/add INCAR parameters.
      Examples: {"ISMEAR": -5, "SIGMA": 0.05}, {"IVDW": 12}

    Examples:
    submit_scf_calculation(optimized_task_id="task_001")
    submit_scf_calculation(cif_url="https://example.com/structure.cif", custom_incar={"ISMEAR": -5})
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = SCFInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_scf(payload)


@mcp.tool()
async def submit_dos_calculation(
    input_url: Optional[Union[str, List[str]]] = None,
    scf_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    kpoint_multiplier: Optional[float] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit a density of states (DOS) calculation task and return task information.

    Default INCAR for DOS: based on SCF settings, then modified with
    ISMEAR=-5 (tetrahedron), LORBIT=11, NEDOS=1500, ICHARG=11, NSW=0, IBRION=-1.
    K-point grid is multiplied by kpoint_multiplier (default 1.5x) relative to optimization.
    Override any parameter via custom_incar.

    Parameters:
    - input_url (optional): Unified external input source. Supported forms:
      1. Single structure file URL such as .cif, which runs single-point SCF + DOS.
      2. Single zip/tar.gz bundle URL.
      3. URL array such as [".../POSCAR", ".../POTCAR", ".../CHGCAR", ".../INCAR"].
         If CHGCAR is available, DOS directly reuses it; otherwise the system runs SCF + DOS.
    - scf_task_id (optional): Completed SCF task ID
    - kpoint_density (optional): default 20.0
    - kpoint_multiplier (optional): K-point multiplier vs optimization, default 1.5
    - precision (optional): "Normal"/"High"/"Accurate", default "High"
    - custom_incar (optional): Dict to override/add INCAR parameters.
      Examples: {"LORBIT": 12, "NEDOS": 3000}, {"EMIN": -20, "EMAX": 10}

    Examples:
    submit_dos_calculation(scf_task_id="scf_001", kpoint_multiplier=3.0)
    submit_dos_calculation(input_url="https://example.com/structure.cif", custom_incar={"LORBIT": 11, "NEDOS": 3000})
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = DOSInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_dos(payload)


@mcp.tool()
async def submit_band_structure_calculation(
    input_url: Optional[Union[str, List[str]]] = None,
    scf_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    line_density: Optional[int] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Submit a band structure calculation task and return task information.

    Band structure calculates the E(k) dispersion relation along high-symmetry
    k-paths in the Brillouin zone. This reveals whether the material is a metal,
    semiconductor, or insulator, and whether the band gap is direct or indirect.

    The calculation uses ICHARG=11 (read charge density from a prior SCF run),
    ISMEAR=0 (Gaussian smearing — NOT tetrahedron), LORBIT=11, NSW=0, IBRION=-1.
    K-points are generated as a line-mode path through high-symmetry points
    using pymatgen's HighSymmKpath.

    Parameters:
    - input_url (optional): Unified external input source. Supported forms:
      1. Single URL string, such as a .cif file or a zip/tar.gz bundle.
      2. URL array, such as [".../POSCAR", ".../POTCAR", ".../CHGCAR"].
      3. URL array with only structure inputs, such as [".../structure.cif", ".../POTCAR"],
         in which case the system automatically runs SCF seed -> Band.
    - scf_task_id (optional): Completed SCF task ID (recommended when available).
      This directly reuses the upstream CHGCAR and runs standard NSCF band structure.
    - kpoint_density (optional): K-point density for the internal SCF seed step, default 20.0
    - line_density (optional): Number of k-points per segment along high-symmetry path, default 20
    - precision (optional): "Normal"/"High"/"Accurate", default "High"
    - custom_incar (optional): Dict to override/add INCAR parameters.
      Note: when reusing CHGCAR, grid-related settings such as PREC are kept compatible with the source task/input.

    Examples:
    submit_band_structure_calculation(scf_task_id="scf_001")
    submit_band_structure_calculation(input_url="https://example.com/structure.cif", line_density=30)
    submit_band_structure_calculation(input_url=["https://example.com/POSCAR", "https://example.com/POTCAR", "https://example.com/CHGCAR"])
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = BandStructureInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_band_structure(payload)


@mcp.tool()
async def submit_md_calculation(
    input_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    md_steps: Optional[int] = None,
    temperature: Optional[Any] = None,
    time_step: Optional[float] = None,
    ensemble: Optional[str] = None,
    precision: Optional[str] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None #type: ignore
) -> Dict[str, Any]:
    """
    Submit molecular dynamics (MD) calculation task(s).

    Each task runs a single temperature. To scan multiple temperatures, pass a list
    (e.g. [300, 400, 500]) — one independent task is created per temperature.
    Use analyze_md_multi_results() afterwards to aggregate the results (Arrhenius, etc.).

    Default MD INCAR parameters (override via custom_incar):
      PREC=Normal (controlled by precision param), ISMEAR=0, SIGMA=0.1,
      IBRION=0, NSW=<md_steps>, POTIM=<time_step>, TEBEG/TEEND=<temperature>,
      SMASS=0 (Nose-Hoover), NBLOCK=1, ISYM=0, LCHARG=.FALSE., LWAVE=.FALSE.
    K-points: fixed 1x1x1 Gamma for MD.

    Parameters:
    - input_url (optional): URL of structure file or a single structure input bundle URL (zip/tar.gz)
    - scf_task_id (optional): Completed SCF task ID
    - md_steps (optional): Number of MD steps, default 1000
    - temperature (optional): Target temperature (K). Single float (300.0) or list of floats ([300, 400, 500])
    - time_step (optional): Time step (fs), default 1.0
    - ensemble (optional): "NVT"/"NVE"/"NPT", default "NVT"
    - precision (optional): "Normal"/"High"/"Accurate", default "Normal"
    - custom_incar (optional): Dict to override/add INCAR parameters.

    Returns:
    - Single temperature: {"task_id": "...", "status": "queued", "message": "..."}
    - Multiple temperatures: {"tasks": [{"temperature": 300, "task_id": "..."}, ...], "message": "..."}

    Examples:
    submit_md_calculation(input_url="https://example.com/structure.cif", md_steps=2000, temperature=350.0)
    submit_md_calculation(input_url="https://example.com/structure.cif", temperature=[300.0, 400.0, 500.0], md_steps=1500)
    """
    user_id = get_user_id(ctx)

    # 判断是否为多温度
    temperatures: list
    if temperature is not None and isinstance(temperature, list):
        temperatures = [float(t) for t in temperature]
    else:
        temperatures = [float(temperature) if temperature is not None else 300.0]

    # 构建共享参数（不含 temperature）
    base_params = {k: v for k, v in locals().items()
                   if v is not None and k not in ("ctx", "temperature", "temperatures", "user_id")}
    base_params["user_id"] = user_id

    if len(temperatures) == 1:
        # 单温度 — 直接提交，返回原有格式
        base_params["temperature"] = temperatures[0]
        payload = MDInput(**base_params).model_dump(mode="json", by_alias=True)
        return await client.submit_md(payload)
    else:
        # 多温度 — 逐个提交独立任务
        results = []
        for temp in temperatures:
            p = dict(base_params)
            p["temperature"] = temp
            payload = MDInput(**p).model_dump(mode="json", by_alias=True)
            resp = await client.submit_md(payload)
            results.append({"temperature": temp, "task_id": resp.get("task_id"), "status": resp.get("status")})
        return {
            "tasks": results,
            "message": f"已提交 {len(results)} 个MD任务（每个温度一个独立任务），"
                       f"全部完成后可使用 analyze_md_multi_results 进行聚合分析"
        }


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

    In distributed mode, `task_id` analysis may be dispatched back to the original HPC queue.
    In that case the response includes `analysis_task_id` and `status`, and you should poll
    `get_task_status(analysis_task_id)` until it completes.

    Returns analysis summary (band_gap, fermi_energy, is_metal) and HTML report URL
    when results are already available.

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
async def analyze_band_structure_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze band structure calculation results independently.

    Provide one of:
    - task_id: ID of an existing completed band structure task
    - file_url: URL to a zip/tar.gz archive (.zip or .tar.gz) containing VASP output files

    Required files in the archive:
    - vasprun.xml  (required — pymatgen parses this to extract band structure data)
    Optional files:
    - OUTCAR  (supplementary convergence and energy info)
    - EIGENVAL  (raw eigenvalues, used as fallback if vasprun.xml parsing fails)
    Do NOT include CHGCAR/WAVECAR — they are large and not needed for analysis.

    In distributed mode, `task_id` analysis may be dispatched back to the original HPC queue.
    In that case the response includes `analysis_task_id` and `status`, and you should poll
    `get_task_status(analysis_task_id)` until it completes.

    Returns analysis summary (band_gap, is_direct, vbm, cbm, fermi_energy) and HTML report URL
    with band structure plot when results are already available.

    Examples:
    analyze_band_structure_results(task_id="abc123")
    analyze_band_structure_results(file_url="https://example.com/band_output.zip")
    """
    payload: Dict[str, Any] = {}
    if task_id:
        payload["task_id"] = task_id
        payload["user_id"] = get_user_id(ctx)
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_band_structure(payload)


@mcp.tool()
async def analyze_md_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze a single MD task's results independently.

    Provide one of:
    - task_id: ID of an existing completed MD task
    - file_url: URL to a zip/tar.gz archive containing VASP output files

    Required files: XDATCAR (trajectory).
    Optional: OUTCAR, vasprun.xml, INCAR.

    For multi-temperature analysis across multiple tasks, use analyze_md_multi_results instead.
    In distributed mode, `task_id` analysis may be dispatched back to the original HPC queue.
    In that case the response includes `analysis_task_id` and `status`, and you should poll
    `get_task_status(analysis_task_id)` until it completes.

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
async def analyze_md_multi_results(
    task_ids: List[str] = [],
    ctx: Context = None  # type: ignore
) -> Dict[str, Any]:
    """
    Aggregate analysis across multiple completed MD tasks (e.g. multi-temperature scan).

    This tool collects results from several independent MD tasks (each at a different
    temperature) and performs cross-temperature analysis including:
    - Arrhenius plot (ln D vs 1/T) and activation energy fitting
    - Per-temperature diffusion coefficients and ionic conductivity
    - Comparative stability and RDF analysis

    Parameters:
    - task_ids (required): List of completed MD task IDs (at least 1, typically 3+)

    Workflow:
    1. Submit individual MD tasks: submit_md_calculation(input_url="https://example.com/structure.cif", temperature=[300, 400, 500])
    2. Wait for all tasks to complete: get_task_status(task_id) for each
    3. Run aggregation: analyze_md_multi_results(task_ids=["id1", "id2", "id3"])

    Examples:
    analyze_md_multi_results(task_ids=["task_001", "task_002", "task_003"])
    """
    payload: Dict[str, Any] = {
        "user_id": get_user_id(ctx),
        "task_ids": task_ids,
    }
    return await client.analyze_md_multi(payload)


@mcp.tool()
async def submit_neb_calculation(
    initial_cif_url: Optional[str] = None,
    initial_task_id: Optional[str] = None,
    final_cif_url: Optional[str] = None,
    final_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    n_images: Optional[int] = None,
    kpoint_density: Optional[float] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Submit a NEB (Nudged Elastic Band) transition state calculation.

    NEB finds the minimum energy path (MEP) and activation barrier between two known
    endpoint structures. Each intermediate image is relaxed simultaneously under NEB forces.

    Default NEB INCAR: ICHAIN=0, LCLIMB=.TRUE. (climbing image), IBRION=3, NSW=300,
    POTIM=0.5, EDIFF=1E-4, EDIFFG=-0.05, ISIF=2, ISYM=0, ENCUT=450.
    Intermediate images are generated by linear interpolation (pymatgen).

    Parameters:
    - initial_cif_url / initial_task_id: initial structure (choose one)
    - final_cif_url / final_task_id: final structure (choose one)
      Tip: use optimized task IDs for best results — ensures well-relaxed endpoints
    - n_images (optional): number of intermediate images excluding endpoints, default 5 (range 2-20)
      More images = smoother MEP but higher cost. 5-8 is typical.
    - kpoint_density (optional): K-point density, default 15.0
    - custom_incar (optional): override INCAR parameters, e.g. {"SPRING": -10, "EDIFFG": -0.03}

    Returns:
    - task_id: use get_task_status() to poll progress

    Examples:
    submit_neb_calculation(initial_task_id="opt_001", final_task_id="opt_002", n_images=7)
    submit_neb_calculation(initial_cif_url="...", final_cif_url="...", n_images=5)
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = NEBInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_neb(payload)


@mcp.tool()
async def submit_phonon_calculation(
    cif_url: Optional[str] = None,
    scf_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    kpoint_density: Optional[float] = None,
    displacement: Optional[float] = None,
    custom_incar: Optional[Dict[str, Any]] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Submit a phonon calculation at the Gamma point (finite displacement method, IBRION=6).

    Computes phonon modes at Gamma — all 3N frequencies for the unit cell. Used to:
    - Verify dynamic stability (no imaginary frequencies = stable)
    - Extract vibrational frequencies for IR/Raman comparison
    - Compute zero-point energy contributions

    Default INCAR: IBRION=6 (finite differences + symmetry), NFREE=2, POTIM=0.015 (displacement),
    EDIFF=1E-6, ENCUT=450, PREC=Normal, LREAL=Auto, ISMEAR=0, SIGMA=0.05.

    Parameters:
    - cif_url (optional): URL of structure file
    - scf_task_id (optional): Completed SCF/optimization task ID (reuses POSCAR+POTCAR)
    - kpoint_density (optional): K-point density, default 15.0
    - displacement (optional): finite displacement step in Å, default 0.015
      Smaller = more accurate but more noise; 0.01-0.02 is typical
    - custom_incar (optional): override INCAR, e.g. {"SIGMA": 0.05} for metals

    Returns:
    - task_id: use get_task_status() to poll. Result includes:
      n_modes, n_imaginary, dynamically_stable, max_imaginary_freq_cm1, max_real_freq_cm1

    Examples:
    submit_phonon_calculation(cif_url="https://example.com/structure.cif")
    submit_phonon_calculation(scf_task_id="scf_001", displacement=0.01)
    """
    params = {k: v for k, v in locals().items() if v is not None and k != "ctx"}
    params["user_id"] = get_user_id(ctx)
    payload = PhononInput(**params).model_dump(mode="json", by_alias=True)
    return await client.submit_phonon(payload)


@mcp.tool()
async def analyze_neb_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze NEB calculation results from a completed task or uploaded output files.

    Extracts energy profile across NEB images and computes:
    - forward_barrier_eV: activation energy (initial → TS)
    - backward_barrier_eV: reverse barrier (final → TS)
    - reaction_energy_eV: ΔE (initial → final)
    - ts_image_index: which image is the transition state
    - html_report_url: energy profile plot and data table

    Parameters:
    - task_id (optional): completed NEB task ID
    - file_url (optional): URL to zip/tar.gz containing numbered image directories (00/, 01/, ...)

    Examples:
    analyze_neb_results(task_id="neb_001")
    """
    payload: Dict[str, Any] = {"user_id": get_user_id(ctx)}
    if task_id:
        payload["task_id"] = task_id
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_neb(payload)


@mcp.tool()
async def analyze_phonon_results(
    task_id: Optional[str] = None,
    file_url: Optional[str] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze phonon calculation results from a completed task or uploaded output files.

    Parses Gamma-point phonon frequencies from OUTCAR (IBRION=5 or 6) and reports:
    - n_modes: total number of phonon modes (= 3 * N_atoms)
    - n_imaginary: number of imaginary frequency modes
    - dynamically_stable: True if no imaginary modes
    - max_imaginary_freq_cm1: largest imaginary frequency (negative sign convention)
    - max_real_freq_cm1: highest real frequency (cm⁻¹)
    - html_report_url: frequency distribution histogram

    Parameters:
    - task_id (optional): completed phonon task ID
    - file_url (optional): URL to zip/tar.gz containing OUTCAR

    Examples:
    analyze_phonon_results(task_id="phonon_001")
    """
    payload: Dict[str, Any] = {"user_id": get_user_id(ctx)}
    if task_id:
        payload["task_id"] = task_id
    if file_url:
        payload["file_url"] = file_url
    return await client.analyze_phonon(payload)


@mcp.tool()
async def submit_custom_calculation(
    cif_url: Optional[str] = None,
    from_task_id: Optional[str] = None,
    queue_name: Optional[str] = None,
    incar: Optional[Dict[str, Any]] = None,
    kpoint_density: Optional[float] = None,
    kpoint_mode: Optional[str] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Submit a generic custom VASP calculation with full INCAR control.

    Use this tool when the dedicated tools (SCF, DOS, band structure, NEB, phonon, MD) do not
    cover your calculation type. Typical use cases:

    Note: Bader charge analysis and ELF (electron localization function) are already supported
    by submit_scf_calculation + the SCF analyzer — do NOT use this tool for those.

    - **Hybrid functional (HSE06 single-point)**:
      incar={"ALGO": "All", "AEXX": 0.25, "HFSCREEN": 0.2, "LHFCALC": ".TRUE.",
             "NSW": 0, "IBRION": -1, "EDIFF": 1e-5, "ISMEAR": 0, "SIGMA": 0.05}

    - **DFPT dielectric constant / Born effective charges**:
      incar={"IBRION": 8, "LEPSILON": ".TRUE.", "EDIFF": 1e-8, "NSW": 1,
             "ISMEAR": 0, "SIGMA": 0.01, "LREAL": ".FALSE."}

    - **Wannier90 wave-function preparation**:
      incar={"NSW": 0, "IBRION": -1, "LWAVE": ".TRUE.", "LCHARG": ".TRUE.",
             "LWANNIER90": ".TRUE.", "EDIFF": 1e-8, "ISMEAR": 0}

    - **GW quasi-particle correction (G0W0)**:
      incar={"ALGO": "GW0", "NELM": 1, "NOMEGA": 50, "ENCUTGW": 200,
             "LWAVE": ".TRUE.", "NSW": 0, "IBRION": -1}

    - **Core-level spectroscopy (XPS/XAS core-hole)**:
      incar={"ISTART": 0, "ICHARG": 2, "NSW": 0, "IBRION": -1,
             "CLNT": 1, "CLN": 1, "CLL": 0, "CLZ": 1.0, "EDIFF": 1e-5}

    - **Large-cell single-point** (use kpoint_mode="gamma" for large supercells):
      incar={"NSW": 0, "IBRION": -1, "EDIFF": 1e-5, "ISMEAR": 0, "SIGMA": 0.1},
      kpoint_mode="gamma"

    - **Constrained magnetism / non-collinear SOC**:
      incar={"LSORBIT": ".TRUE.", "LNONCOLLINEAR": ".TRUE.", "NSW": 0,
             "IBRION": -1, "ISYM": -1, "EDIFF": 1e-6}

    Parameters:
    - cif_url / from_task_id: structure source (exactly one required).
      `from_task_id` reuses CONTCAR (or POSCAR) from any completed task.
    - incar: dict of VASP INCAR parameters. All keys are accepted — no blacklist.
      SYSTEM defaults to "CUSTOM" if not provided.
    - kpoint_mode: "mesh" (default, Monkhorst-Pack grid) or "gamma" (Gamma-only, 1×1×1).
    - kpoint_density: target k-point density for mesh mode (default 30.0).

    Returns task_id for status polling via get_task_status.
    """
    user_id = get_user_id(ctx)
    payload = CustomCalcInput(
        user_id=user_id,
        cif_url=cif_url,  # type: ignore
        from_task_id=from_task_id,
        queue_name=queue_name,
        incar=incar or {},
        kpoint_density=kpoint_density if kpoint_density is not None else 30.0,
        kpoint_mode=kpoint_mode or "mesh",
    ).model_dump(mode="json", exclude_none=True)
    return await client.submit_custom_calculation(payload)


@mcp.tool()
async def analyze_vasp_output(
    task_id: str,
    question: str,
    model: Optional[str] = None,
    ctx: Context = None,  # type: ignore
) -> Dict[str, Any]:
    """
    Analyze any completed VASP calculation using an AI agent with a persistent
    Python code interpreter and a pre-loaded skills library.

    The inner agent (LangGraph ReAct) can:
    - Call pre-built skill functions (get_band_gap, plot_dos, get_lattice_params, ...)
    - Write arbitrary Python using pymatgen / numpy / matplotlib
    - Persist variables across multiple code executions in one session
    - Generate plots saved as PNG files (URLs returned in response)

    Available skills injected into the Python environment:
      Energy:     extract_final_energy, extract_energy_per_atom,
                  get_ionic_steps_energy, get_max_force, get_computation_time
      Structure:  get_final_structure, get_volume, get_lattice_params,
                  get_formula, get_n_atoms, get_volume_change
      Electronic: is_converged, get_fermi_energy, get_band_gap,
                  get_total_energy_from_vasprun, get_electronic_steps, get_magnetization
      Plots:      plot_energy_convergence, plot_forces_convergence,
                  plot_dos, plot_band_structure

    Args:
        task_id: Any completed VASP task (optimization, SCF, DOS, band, MD, NEB, phonon, custom)
        question: Natural language question or analysis request, e.g.:
                  "What is the band gap? Is it direct or indirect?"
                  "Plot the energy and force convergence curves."
                  "Compare initial vs final volume. What changed?"
                  "Extract all ionic step energies and fit a convergence rate."
        model: Inner agent model (default: qwen3.5-plus).
               Claude models still work when ANTHROPIC_API_KEY is configured.

    In distributed mode this request is executed as an HPC-side analysis task. The response may
    therefore contain `analysis_task_id` and `status` first; poll `get_task_status()` on that
    task until it finishes, then read the final `answer`, `html_report_url`, and plot URLs.
    """
    user_id = get_user_id(ctx)
    payload: Dict[str, Any] = {
        "user_id": user_id,
        "task_id": task_id,
        "question": question,
    }
    if model:
        payload["model"] = model
    return await client.agent_analyze(payload)


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
