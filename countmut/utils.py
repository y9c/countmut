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
    """
    Get the appropriate output headers based on save_rest parameter.

    Args:
        save_rest: Whether to include additional statistics columns

    Returns:
        List of column headers
        - u0, u1, u2: unconverted (reference base) counts
        - m0, m1, m2: mutation (mutation base only) counts
        - o0, o1, o2: other bases counts (only with save_rest)
    """
    if save_rest:
        return [
            "chrom",
            "pos",
            "strand",
            "motif",
            "u0",
            "u1",
            "u2",
            "m0",
            "m1",
            "m2",
            "o0",
            "o1",
            "o2",
        ]
    return [
        "chrom",
        "pos",
        "strand",
        "motif",
        "u0",
        "u1",
        "u2",
        "m0",
        "m1",
        "m2",
    ]


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
