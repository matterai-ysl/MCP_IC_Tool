"""向后兼容 — 实际实现已迁移至 analyzers.optimization"""
from .analyzers.optimization import (  # noqa: F401
    OUTCARAnalyzer,
    OptimizationHTMLGenerator,
    generate_optimization_report,
)
