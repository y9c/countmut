#!/usr/bin/env python3
"""
Python wrapper around the C countmut core.

All computation happens in ``backend/countmut_core`` (both the read-walk and
the pileup engine are implemented in C).  This module only: builds the command
line from the config dataclasses, makes sure indices exist, shells out to the
binary and collects a small summary.  There is deliberately no Python counting
engine: results come from C alone.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

from .model import EngineConfig, FilterConfig, StrandConfig

_PACKAGE_DIR = Path(__file__).resolve().parent
# Binary shipped inside the installed package (bundled in the wheel).
SHIPPED_BINARY = _PACKAGE_DIR / "_core" / "countmut_core"
# Binary built from the source checkout (dev / `pip install -e .`).
BACKEND_DIR = _PACKAGE_DIR.parent / "backend"
SRC_BINARY = BACKEND_DIR / "countmut_core"


# ---------------------------------------------------------------------------
# locating / building the binary
# ---------------------------------------------------------------------------
def find_binary() -> Path | None:
    for candidate in (SHIPPED_BINARY, SRC_BINARY):
        if candidate.exists():
            return candidate
    return None


def build_binary() -> Path | None:
    """Compile the backend with `make` (source checkout); return the binary or None."""
    if not (BACKEND_DIR / "Makefile").exists():
        return None
    try:
        subprocess.run(
            ["make"],
            cwd=BACKEND_DIR,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"[countmut] backend build failed: {exc.stderr.decode()}\n")
        return None
    return SRC_BINARY if SRC_BINARY.exists() else None


def ensure_backend() -> Path | None:
    binary = find_binary()
    if binary is None:
        binary = build_binary()
    return binary


# ---------------------------------------------------------------------------
# config -> CLI args
# ---------------------------------------------------------------------------
def _build_cmd(
    binary: Path,
    *,
    samfile: str,
    reference: str,
    output: str | None,
    fcfg: FilterConfig,
    scfg: StrandConfig,
    ecfg: EngineConfig,
) -> list[str]:
    cmd = [
        str(binary),
        "--bam",
        samfile,
        "--fa",
        reference,
        "--out",
        output or "-",
        "--engine",
        {"read-walk": "read-walk", "pileup": "pileup", "auto": "auto"}[ecfg.engine],
        "--threads",
        str(ecfg.threads or min(os.cpu_count() or 1, 8)),
        "--max-depth",
        str(fcfg.max_depth),  # 0 = unlimited
    ]
    if ecfg.count_indels:
        cmd.append("--count-indels")
    if ecfg.strandless:
        cmd.append("--strandless")
    if scfg.process in ("forward", "reverse"):
        cmd += ["--strand", scfg.process]
    if ecfg.read_expr:
        cmd += ["--read-expr", ecfg.read_expr]
    if ecfg.pile_expr:
        cmd += ["--pile-expr", ecfg.pile_expr]
    if ecfg.output_expr:
        cmd += ["--output-expr", ecfg.output_expr]
    if ecfg.fmt_header:
        cmd += ["--fmt-header", ecfg.fmt_header]
    if ecfg.motif_pad:
        cmd += ["--motif-pad", str(ecfg.motif_pad)]
    if ecfg.region:
        cmd += ["--region", ecfg.region]
    if ecfg.vcf:
        cmd.append("--vcf")
    # --min-depth / --min-allele-support are expressed via -p (depth >= N, g >= N)
    if ecfg.verbose:
        cmd.append("--verbose")
    return cmd


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------
def _forward_stderr(stream) -> None:
    """Forward a child process' stderr lines to our stderr (--verbose)."""
    try:
        for line in stream:
            sys.stderr.write(line)
            sys.stderr.flush()
    except (ValueError, OSError):
        pass


def run_backend(
    samfile: str,
    reference: str,
    output: str | None = None,
    *,
    fcfg: FilterConfig | None = None,
    scfg: StrandConfig | None = None,
    ecfg: EngineConfig | None = None,
) -> dict:
    """Run the C core; returns a summary dict.

    The C binary must exist (it is auto-built by :func:`ensure_backend`); there
    is no pure-Python counting fallback.
    """
    fcfg = fcfg or FilterConfig()
    scfg = scfg or StrandConfig()
    ecfg = ecfg or EngineConfig()

    # ensure indices exist
    if not os.path.exists(samfile + ".bai"):
        import pysam

        pysam.index(samfile)
    if not os.path.exists(reference + ".fai"):
        import pysam

        pysam.faidx(reference)

    binary = ensure_backend()
    if binary is None:
        raise RuntimeError(
            "countmut core binary not found and could not be built "
            f"(looked in {BACKEND_DIR}). Run `make` there."
        )

    cmd = _build_cmd(
        binary,
        samfile=samfile,
        reference=reference,
        output=output,
        fcfg=fcfg,
        scfg=scfg,
        ecfg=ecfg,
    )
    start = time.time()
    stdout = ""
    if ecfg.verbose:
        # stream the C core's stderr (progress) to our stderr in real time
        proc = subprocess.Popen(
            cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
        )
        import threading

        fwd = threading.Thread(target=_forward_stderr, args=(proc.stderr,), daemon=True)
        fwd.start()
        out, _ = proc.communicate()
        fwd.join(timeout=2)
        stdout = out or ""
    else:
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stdout = proc.stdout
    elapsed = time.time() - start

    if output is None and proc.returncode == 0:
        # stream the TSV from stdout
        sys.stdout.write(stdout)
        sys.stdout.flush()

    rows = 0
    if output:
        with open(output) as fh:
            rows = sum(1 for _ in fh) - 1
    else:
        rows = len(stdout.splitlines()) - 1

    if proc.returncode != 0:
        sys.stderr.write(f"[countmut] backend error: {proc.stderr}\n")
        return {
            "backend": "c",
            "success": False,
            "error": proc.stderr or "unknown",
            "total_sites": rows,
            "total_depth": rows,
            "elapsed": elapsed,
            "output": output,
        }

    return {
        "backend": "c",
        "success": True,
        "total_sites": rows,
        "total_depth": rows,
        "elapsed": elapsed,
        "output": output,
        "error": None,
    }
