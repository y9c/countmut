"""
CountMut: Unified ultra-fast strand-aware mutation counter.

This package drives a C core (``backend/countmut_core``) through a thin Python
wrapper.  It supports strand-aware mutation counting (bisulfite NS/Zf/Yf tiers),
per-site base counting, and allele/VCF output, with mate-overlap deduplication.

Author: Ye Chang
"""

from .backend import run_backend
from .model import EngineConfig, FilterConfig, MutationConfig, StrandConfig
from .pipeline import PipelineResult, run_pipeline

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"

__all__ = [
    "run_backend",
    "run_pipeline",
    "PipelineResult",
    "EngineConfig",
    "FilterConfig",
    "MutationConfig",
    "StrandConfig",
]
