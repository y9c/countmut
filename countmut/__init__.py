"""
CountMut - Fast, parallel mutation counting from BAM pileup data

This package provides efficient mutation counting functionality with:
- Parallel processing using genomic windows
- Bisulfite conversion analysis
- Rich logging and progress tracking
- Modern CLI interface

Author: Ye Chang
Date: 2025-10-23
"""

from .cli import main
from .core import count_mutations
from .utils import format_duration, get_output_headers, write_output

__author__ = "Ye Chang"
__email__ = "yech1990@gmail.com"

__all__ = [
    "count_mutations",
    "format_duration",
    "get_output_headers",
    "main",
    "write_output",
]
