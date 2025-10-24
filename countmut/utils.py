#!/usr/bin/env python3
"""
Utility functions for CountMut.

This module provides common utility functions used across the package.

Author: Ye Chang
Date: 2025-10-23
"""

import os


def format_duration(sec: float) -> str:
    """
    Format duration in seconds to human-readable string.

    Args:
        sec: Duration in seconds

    Returns:
        Formatted duration string (e.g., "12.34s", "5m 23s", "2h 15m 30s")
    """
    if sec < 60:
        return f"{sec:.2f}s"
    if sec < 3600:
        minutes = int(sec // 60)
        seconds = int(round(sec % 60))
        return f"{minutes}m {seconds}s"
    hours = int(sec // 3600)
    rem = sec % 3600
    minutes = int(rem // 60)
    seconds = int(round(rem % 60))
    return f"{hours}h {minutes}m {seconds}s"


def get_output_headers(save_rest: bool = False) -> list[str]:
    """Return the output column headers based on whether to include 'other' counts.

    The column headers are defined as follows:
    - u: unconverted (reference base)
    - m: mutation (mutation base only)
    - o: other bases (only with save_rest)

    Count categories (x0, x1, x2):
    - x0 (low quality): Bases failing quality filters (trim region, max-sub, min-mapq, min-baseq)
    - x1 (insufficient conversion): Bases from reads with insufficient conversion efficiency (high Zf or low Yf)
    - x2 (high conversion): Bases from reads with high conversion efficiency (low Zf and high Yf)
    """
    headers = ["chrom", "pos", "strand", "motif"]
    base_cols = ["u0", "u1", "u2", "m0", "m1", "m2"]
    if save_rest:
        base_cols += ["o0", "o1", "o2"]
    headers.extend(base_cols)
    return headers


def write_output(
    results: list[list],
    output_file: str | None = None,
    save_rest: bool = False,
) -> None:
    """
    Write results to file or stdout.

    Args:
        results: List of result rows to write
        output_file: Path to output file (if None, prints to stdout)
        save_rest: Whether to include additional statistics columns
    """
    if not results:
        # If there are no results, write an empty file with headers
        if output_file:
            os.makedirs(os.path.dirname(output_file), exist_ok=True)
            with open(output_file, "w") as out_f:
                out_f.write("\t".join(get_output_headers(save_rest)) + "\n")
        else:
            print("\t".join(get_output_headers(save_rest)))
        return

    # headers are always generated without alt_ref/alt_mut for TSV output
    headers = get_output_headers(save_rest)

    if output_file:
        # Ensure output directory exists
        os.makedirs(os.path.dirname(output_file), exist_ok=True)

        with open(output_file, "w") as out_f:
            # Write header
            out_f.write("\t".join(headers) + "\n")

            # Write data
            for row in results:
                out_f.write("\t".join(map(str, row)) + "\n")
    else:
        # Print to stdout
        print("\t".join(headers))

        for row in results:
            print("\t".join(map(str, row)))
