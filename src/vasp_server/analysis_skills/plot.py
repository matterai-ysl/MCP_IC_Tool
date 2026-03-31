"""可视化技能 — 保存 PNG 并返回文件路径"""
from pathlib import Path
from typing import Optional


def plot_energy_convergence(work_dir: str, filename: str = "energy_convergence.png") -> str:
    """绘制离子步能量收敛曲线，返回 PNG 路径"""
    from .energy import get_ionic_steps_energy
    import matplotlib.pyplot as plt

    energies = get_ionic_steps_energy(work_dir)
    if not energies:
        return "Error: no energy data in OUTCAR"

    fig, ax = plt.subplots(figsize=(7, 4))
    steps = list(range(1, len(energies) + 1))
    ax.plot(steps, energies, "b-o", markersize=4, linewidth=1.5)
    ax.set_xlabel("Ionic Step")
    ax.set_ylabel("Total Energy (eV)")
    ax.set_title("Energy Convergence")
    ax.grid(True, alpha=0.3)
    # Annotate final energy
    ax.annotate(f"{energies[-1]:.4f} eV", xy=(steps[-1], energies[-1]),
                xytext=(-40, 10), textcoords="offset points", fontsize=8)
    plt.tight_layout()
    out = str(Path(work_dir) / filename)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out


def plot_dos(work_dir: str, filename: str = "dos_plot.png",
             xlim: Optional[tuple] = (-6, 6)) -> str:
    """绘制总态密度图，返回 PNG 路径"""
    vr_path = Path(work_dir) / "vasprun.xml"
    if not vr_path.exists():
        return "Error: vasprun.xml not found"
    try:
        import matplotlib.pyplot as plt
        from pymatgen.io.vasp import Vasprun
        from pymatgen.electronic_structure.plotter import DosPlotter

        vr = Vasprun(str(vr_path), parse_dos=True)
        cdos = vr.complete_dos
        plotter = DosPlotter()
        plotter.add_dos("Total DOS", cdos)
        ax = plotter.get_plot(xlim=xlim)
        fig = ax.get_figure()
        out = str(Path(work_dir) / filename)
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        return out
    except Exception as e:
        return f"Error: {e}"


def plot_band_structure(work_dir: str, filename: str = "band_structure.png") -> str:
    """绘制能带结构图（需要 KPOINTS line-mode），返回 PNG 路径"""
    vr_path = Path(work_dir) / "vasprun.xml"
    kpts_path = Path(work_dir) / "KPOINTS"
    if not vr_path.exists():
        return "Error: vasprun.xml not found"
    try:
        import matplotlib.pyplot as plt
        from pymatgen.io.vasp import Vasprun
        from pymatgen.electronic_structure.plotter import BSPlotter

        vr = Vasprun(str(vr_path), parse_dos=False)
        bs = vr.get_band_structure(
            kpoints_filename=str(kpts_path) if kpts_path.exists() else None,
            line_mode=True,
        )
        plotter = BSPlotter(bs)
        ax = plotter.get_plot()
        fig = ax.get_figure()
        out = str(Path(work_dir) / filename)
        fig.savefig(out, dpi=120, bbox_inches="tight")
        plt.close(fig)
        return out
    except Exception as e:
        return f"Error: {e}"


def plot_forces_convergence(work_dir: str, filename: str = "forces_convergence.png") -> str:
    """绘制离子步最大力收敛曲线，返回 PNG 路径"""
    import re
    import matplotlib.pyplot as plt

    outcar = Path(work_dir) / "OUTCAR"
    if not outcar.exists():
        return "Error: OUTCAR not found"

    content = outcar.read_text(errors="ignore")
    blocks = list(re.finditer(r"TOTAL-FORCE.*?\n(.*?)(?=\n-{10})", content, re.DOTALL))
    if not blocks:
        return "Error: no force data found"

    max_forces = []
    for block in blocks:
        forces = []
        for line in block.group(1).strip().splitlines():
            parts = line.split()
            if len(parts) == 6:
                try:
                    fx, fy, fz = float(parts[3]), float(parts[4]), float(parts[5])
                    forces.append((fx**2 + fy**2 + fz**2) ** 0.5)
                except ValueError:
                    pass
        if forces:
            max_forces.append(max(forces))

    if not max_forces:
        return "Error: could not parse forces"

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(range(1, len(max_forces) + 1), max_forces, "r-o", markersize=4)
    ax.set_xlabel("Ionic Step")
    ax.set_ylabel("Max Force (eV/Å)")
    ax.set_title("Force Convergence")
    ax.axhline(0.05, color="gray", linestyle="--", alpha=0.6, label="0.05 eV/Å")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out = str(Path(work_dir) / filename)
    fig.savefig(out, dpi=120)
    plt.close(fig)
    return out
