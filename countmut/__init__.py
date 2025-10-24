"""
CountMut: Ultra-fast strand-aware mutation counter.

This package provides a command-line tool and a library for counting mutations
from BAM pileup data, with a focus on bisulfite sequencing analysis.
"""

from .core import count_mutations
from .utils import format_duration, get_output_headers, write_output

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"

__all__ = [
    "count_mutations",
    "format_duration",
    "get_output_headers",
    "write_output",
]
