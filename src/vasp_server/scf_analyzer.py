"""向后兼容 — 实际实现已迁移至 analyzers.scf"""
from .analyzers.scf import (  # noqa: F401
    SCFAnalyzer,
    SCFHTMLGenerator,
    generate_scf_report,
)
