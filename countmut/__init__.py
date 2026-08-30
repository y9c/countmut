"""
CountMut: unified ultra-fast strand-aware mutation counter.

This package drives a C core (``backend/countmut_core``) through a thin Python
wrapper.  Both BAM-walk strategies -- read-walk and pileup -- are implemented in
C; the Python side only builds the command line and collects the output.  The
tool supports strand-aware mutation counting (bisulfite NS/Zf/Yf tiers),
per-site base counting, and allele/VCF output, with mate-overlap deduplication.

The original pure-Python implementation is preserved in :mod:`countmut.core`
for comparison and reference.

Author: Ye Chang
"""

from .backend import run_backend
from .model import EngineConfig, FilterConfig, MutationConfig, StrandConfig

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"

__all__ = [
    "run_backend",
    "EngineConfig",
    "FilterConfig",
    "MutationConfig",
    "StrandConfig",
]
