"""Shared helper: drive the C countmut core from tests.

Both the read-walk and pileup engines live in the C binary
(``backend/countmut_core``); tests therefore run the binary directly instead of
importing Python engines (there are none -- Python only wraps C).
"""

import subprocess

from countmut.backend import ensure_backend

_BINARY = None


def binary() -> str:
    global _BINARY
    if _BINARY is None:
        p = ensure_backend()
        if p is None:
            raise RuntimeError("C backend not built; run `make` in backend/")
        _BINARY = str(p)
    return _BINARY


def _num(v):
    try:
        return int(v)
    except (ValueError, TypeError):
        return v


def run_c(
    bam,
    fa,
    *,
    mode="mutation",
    region=None,
    engine="read-walk",
    extra=None,
    threads=2,
    binary_path=None,
):
    """Run the C core and return (header_cols, rows) with numeric cells parsed."""
    cmd = [
        binary_path or binary(),
        "--bam",
        bam,
        "--fa",
        fa,
        "--out",
        "-",
        "--mode",
        mode,
        "--engine",
        engine,
        "--min-mapq",
        "0",
        "--min-baseq",
        "0",
        "--trim-fragment-start",
        "0",
        "--trim-fragment-end",
        "0",
        "--threads",
        str(threads),
    ]
    if region:
        cmd += ["--region", region]
    cmd += list(extra or [])
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        raise AssertionError(f"C core failed ({r.returncode}): {r.stderr}")
    lines = [ln for ln in r.stdout.splitlines() if ln and not ln.startswith("#")]
    if not lines:
        return [], []
    header = lines[0].split("\t")
    rows = [[_num(x) for x in ln.split("\t")] for ln in lines[1:]]
    return header, rows


def norm_rows(rows):
    """Rows -> sorted tuples of (str, int, ...) for comparison."""
    return sorted(tuple(r) for r in rows)
