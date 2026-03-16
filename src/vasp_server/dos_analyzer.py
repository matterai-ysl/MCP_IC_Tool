"""向后兼容 — 实际实现已迁移至 analyzers.dos"""
from .analyzers.dos import (  # noqa: F401
    PyMatGenDOSAnalyzer,
    PyMatGenDOSHTMLGenerator,
    generate_pymatgen_dos_report,
)
