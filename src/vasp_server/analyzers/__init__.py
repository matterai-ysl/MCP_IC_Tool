"""
VASP 分析器包 — 从各计算类型的 VASP 输出中提取数据并生成报告。

向后兼容导入：
    from src.vasp_server.analyzers import OUTCARAnalyzer, generate_optimization_report
    from src.vasp_server.analyzers import SCFAnalyzer, generate_scf_report
    from src.vasp_server.analyzers import PyMatGenDOSAnalyzer, generate_pymatgen_dos_report
    from src.vasp_server.analyzers import VASP_MDAnalyzer, generate_md_analysis_report
    from src.vasp_server.analyzers import BandGapAnalyzer, analyze_bandgap

DOS/MD/Bandgap 模块依赖 pymatgen 等重量级库，使用延迟导入避免启动时报错。
"""

import importlib as _importlib

# 轻量级模块 — 直接导入
from .base import BaseAnalyzer
from .html_base import BaseHTMLGenerator
from .optimization import OUTCARAnalyzer, OptimizationHTMLGenerator, generate_optimization_report
from .scf import SCFAnalyzer, SCFHTMLGenerator, generate_scf_report

__all__ = [
    # Base
    "BaseAnalyzer", "BaseHTMLGenerator",
    # Optimization
    "OUTCARAnalyzer", "OptimizationHTMLGenerator", "generate_optimization_report",
    # SCF
    "SCFAnalyzer", "SCFHTMLGenerator", "generate_scf_report",
    # DOS
    "PyMatGenDOSAnalyzer", "PyMatGenDOSHTMLGenerator", "generate_pymatgen_dos_report",
    # MD
    "VASP_MDAnalyzer", "MDHTMLReportGenerator", "generate_md_analysis_report",
    # Bandgap
    "BandGapAnalyzer", "analyze_bandgap",
]

# 延迟导入映射：名称 → (子模块, 属性名)
_LAZY_IMPORTS = {
    "PyMatGenDOSAnalyzer": (".dos", "PyMatGenDOSAnalyzer"),
    "PyMatGenDOSHTMLGenerator": (".dos", "PyMatGenDOSHTMLGenerator"),
    "generate_pymatgen_dos_report": (".dos", "generate_pymatgen_dos_report"),
    "VASP_MDAnalyzer": (".md", "VASP_MDAnalyzer"),
    "MDHTMLReportGenerator": (".md", "MDHTMLReportGenerator"),
    "generate_md_analysis_report": (".md", "generate_md_analysis_report"),
    "BandGapAnalyzer": (".bandgap", "BandGapAnalyzer"),
    "analyze_bandgap": (".bandgap", "analyze_bandgap"),
}


def __getattr__(name: str):
    if name in _LAZY_IMPORTS:
        submod_name, attr = _LAZY_IMPORTS[name]
        mod = _importlib.import_module(submod_name, __name__)
        val = getattr(mod, attr)
        # 缓存到模块命名空间，后续直接访问
        globals()[name] = val
        return val
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
