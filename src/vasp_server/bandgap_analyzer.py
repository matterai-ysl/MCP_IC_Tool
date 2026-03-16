"""向后兼容 — 实际实现已迁移至 analyzers.bandgap"""
from .analyzers.bandgap import (  # noqa: F401
    BandGapAnalyzer,
    analyze_bandgap,
)
