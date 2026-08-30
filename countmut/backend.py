#!/usr/bin/env python3
"""
Bridge between the Python wrapper and the C countmut backend.

The C binary (``backend/countmut_core``) does ALL the computation.  This module
only: builds the config, makes sure indices exist, shells out to the binary, and
collects the output for a rich summary.  If the binary is not built, it falls
back to the pure-Python pipeline (``countmut.pipeline``) so the tool always runs.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

from .model import EngineConfig, FilterConfig, MutationConfig, StrandConfig

_PACKAGE_DIR = Path(__file__).resolve().parent
BACKEND_DIR = _PACKAGE_DIR.parent / "backend"
BINARY = BACKEND_DIR / "countmut_core"


# ---------------------------------------------------------------------------
# locating / building the binary
# ---------------------------------------------------------------------------
def find_binary() -> Path | None:
    if BINARY.exists():
        return BINARY
    return None


def build_binary() -> bool:
    """Compile the backend with `make`; return True on success."""
    if not (BACKEND_DIR / "Makefile").exists():
        return False
    try:
        subprocess.run(
            ["make"], cwd=BACKEND_DIR, check=True,
            stdout=subprocess.DEVNULL, stderr=subprocess.PIPE,
        )
    except subprocess.CalledProcessError as exc:
        sys.stderr.write(f"[countmut] backend build failed: {exc.stderr.decode()}\n")
        return False
    return BINARY.exists()


def ensure_backend() -> Path | None:
    binary = find_binary()
    if binary is None:
        binary = build_binary()
    return binary


# ---------------------------------------------------------------------------
# config -> CLI args
# ---------------------------------------------------------------------------
def _build_cmd(binary: Path, *, samfile: str, reference: str, output: str | None,
               fcfg: FilterConfig, mcfg: MutationConfig | None, scfg: StrandConfig,
               ecfg: EngineConfig) -> list[str]:
    cmd = [
        str(binary), "--bam", samfile, "--fa", reference,
        "--out", output or "-",
        "--mode", ecfg.mode,
        "--engine", {"read-walk": "read-walk", "pileup": "pileup", "auto": "auto"}[ecfg.engine],
        "--min-mapq", str(fcfg.min_mapq),
        "--min-baseq", str(fcfg.min_baseq),
        "--max-sub", str(fcfg.max_sub),
        "--trim-start", str(fcfg.trim_start),
        "--trim-end", str(fcfg.trim_end),
        "--threads", str(ecfg.threads or min(os.cpu_count() or 1, 8)),
        "--max-depth", str(fcfg.max_depth or 8000),
    ]
    if fcfg.max_unc is not None:
        cmd += ["--max-unc", str(fcfg.max_unc)]
    if fcfg.min_con is not None:
        cmd += ["--min-con", str(fcfg.min_con)]
    if ecfg.count_indels:
        cmd.append("--count-indels")
    if ecfg.split_strand or ecfg.mode != "mutation":
        cmd.append("--split-strand")
    cmd += ["--strand", {"both": "both", "forward": "forward", "reverse": "reverse"}[scfg.process]]
    if ecfg.region:
        cmd += ["--region", ecfg.region]
    if ecfg.mode == "mutation" and mcfg is not None:
        cmd += ["--ref-base", mcfg.ref_base, "--mut-base", mcfg.mut_base,
                "--pad", str(mcfg.pad)]
        if mcfg.save_rest:
            cmd.append("--save-rest")
        if mcfg.ref_base2:
            cmd += ["--ref-base2", mcfg.ref_base2, "--mut-base2", mcfg.mut_base2 or "T"]
    if ecfg.mode == "allele":
        cmd += ["--min-allele-support", str(ecfg.min_allele_support)]
        if ecfg.vcf:
            cmd.append("--vcf")
        if ecfg.min_depth:
            cmd += ["--min-depth", str(ecfg.min_depth)]
    return cmd


# ---------------------------------------------------------------------------
# public entry point
# ---------------------------------------------------------------------------
def run_backend(
    samfile: str,
    reference: str,
    output: str | None = None,
    *,
    fcfg: FilterConfig | None = None,
    mcfg: MutationConfig | None = None,
    scfg: StrandConfig | None = None,
    ecfg: EngineConfig | None = None,
) -> dict:
    """Run the C backend; returns a summary dict.

    Returns ``{"backend": "c", ...}`` on success, or ``{"backend": "python"}``
    if it had to fall back (with the unified pipeline result attached).
    """
    fcfg = fcfg or FilterConfig()
    scfg = scfg or StrandConfig()
    ecfg = ecfg or EngineConfig()

    # ---- string expression filters need the Python engine (C can't eval) ----
    if ecfg.read_expr or ecfg.pile_expr:
        from .pipeline import run_pipeline
        start = time.time()
        res = run_pipeline(samfile, reference, output, fcfg=fcfg, mcfg=mcfg,
                           scfg=scfg, ecfg=ecfg)
        return {"backend": "python", "success": res.success, "total_sites": res.total_sites,
                "total_depth": res.total_depth, "elapsed": time.time() - start,
                "error": res.error, "output": output, "note": "expression filter"}

    # ensure indices exist
    if not os.path.exists(samfile + ".bai"):
        import pysam
        pysam.index(samfile)
    if not os.path.exists(reference + ".fai"):
        import pysam
        pysam.faidx(reference)

    binary = ensure_backend()
    if binary is None:
        # graceful fall back to pure Python
        from .pipeline import run_pipeline
        start = time.time()
        res = run_pipeline(samfile, reference, output, fcfg=fcfg, mcfg=mcfg,
                           scfg=scfg, ecfg=ecfg)
        return {"backend": "python", "success": res.success, "total_sites": res.total_sites,
                "total_depth": res.total_depth, "elapsed": time.time() - start,
                "error": res.error, "output": output}

    cmd = _build_cmd(binary, samfile=samfile, reference=reference, output=output,
                     fcfg=fcfg, mcfg=mcfg, scfg=scfg, ecfg=ecfg)
    start = time.time()
    proc = subprocess.run(cmd, capture_output=True, text=True)
    elapsed = time.time() - start

    if output is None and proc.returncode == 0:
        # stream the TSV from stdout (already printed by the driver)
        sys.stdout.write(proc.stdout)
        sys.stdout.flush()

    # count output rows for the summary (read from output or captured stdout)
    rows = 0
    if output:
        with open(output) as fh:
            rows = sum(1 for _ in fh) - 1
    else:
        rows = len(proc.stdout.splitlines()) - 1

    if proc.returncode != 0:
        sys.stderr.write(f"[countmut] backend error: {proc.stderr}\n")
        return {"backend": "c", "success": False, "error": proc.stderr or "unknown",
                "total_sites": rows, "total_depth": rows, "elapsed": elapsed, "output": output}

    return {"backend": "c", "success": True, "total_sites": rows, "total_depth": rows,
            "elapsed": elapsed, "output": output, "error": None}
