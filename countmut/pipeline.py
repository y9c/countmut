#!/usr/bin/env python3
"""
Unified countmut pipeline.

Ties configuration, engine selection, region binning, parallel execution and
output together.  Both engines (read-walk and pileup) are driven here, so the
rest of the code-base never has to know which "way" produced a result.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from typing import Any

import pysam

from .core import read_fasta_index
from .engine_pileup import pileup_region
from .engine_readwalk import readwalk_region
from .formatter import (
    BASE_HEADER,
    BASE_HEADER_INDEL,
    BASE_HEADER_SPLIT,
    BASE_HEADER_SPLIT_INDEL,
    allele_rows,
    base_rows,
    mutation_header,
    mutation_rows,
)
from .model import EngineConfig, FilterConfig, MutationConfig, StrandConfig

# ---------------------------------------------------------------------------
# globals per worker process (opened once by the pool initializer)
# ---------------------------------------------------------------------------
_WORKER_SAM = None
_WORKER_REF = None
_WORKER_FCFG = None
_WORKER_MCFG = None
_WORKER_SCFG = None
_WORKER_ECFG = None
_WORKER_STRANDS = ("+", "-")
_WORKER_BISULFITE = False


def _init_worker(
    samfile: str,
    reffile: str,
    fcfg: FilterConfig,
    mcfg: MutationConfig | None,
    scfg: StrandConfig,
    ecfg: EngineConfig,
    strand_process: str,
    has_bisulfite: bool,
) -> None:
    global \
        _WORKER_SAM, \
        _WORKER_REF, \
        _WORKER_FCFG, \
        _WORKER_MCFG, \
        _WORKER_SCFG, \
        _WORKER_ECFG
    global _WORKER_STRANDS, _WORKER_BISULFITE
    if _WORKER_SAM is None:
        _WORKER_SAM = pysam.AlignmentFile(samfile, "rb")
    _WORKER_REF = pysam.FastaFile(reffile)
    _WORKER_FCFG = fcfg
    _WORKER_MCFG = mcfg
    _WORKER_SCFG = scfg
    _WORKER_ECFG = ecfg
    _WORKER_STRANDS = ("+", "-") if strand_process == "both" else (strand_process,)
    _WORKER_BISULFITE = has_bisulfite


def _worker(args: tuple) -> dict[str, Any]:
    """Process a single bin and return its output rows + counters."""
    chrom, start, end = args
    try:
        mode = _WORKER_ECFG.mode
        from .expression import compile_pile_pred, compile_read_pred

        read_pred = (
            compile_read_pred(_WORKER_ECFG.read_expr)
            if _WORKER_ECFG.read_expr
            else None
        )
        pile_pred = (
            compile_pile_pred(_WORKER_ECFG.pile_expr)
            if _WORKER_ECFG.pile_expr
            else None
        )
        if _WORKER_ECFG.engine == "read-walk" or (
            _WORKER_ECFG.engine == "auto" and mode == "mutation"
        ):
            cols = readwalk_region(
                _WORKER_SAM,
                _WORKER_REF,
                chrom,
                start,
                end,
                _WORKER_FCFG,
                _WORKER_MCFG,
                mode=mode,
                strand_process=_WORKER_SCFG.process,
                has_bisulfite_tags=_WORKER_BISULFITE,
                read_pred=read_pred,
                pile_pred=pile_pred,
            )
        else:
            cols = pileup_region(
                _WORKER_SAM,
                _WORKER_REF,
                chrom,
                start,
                end,
                _WORKER_FCFG,
                _WORKER_MCFG,
                mode=mode,
                strand_process=_WORKER_SCFG.process,
                has_bisulfite_tags=_WORKER_BISULFITE,
                read_pred=read_pred,
                pile_pred=pile_pred,
            )

        rows = _render(cols)
        depth = sum(c.total_depth() for c in cols)
        return {
            "region": f"{chrom}:{start}-{end}",
            "rows": rows,
            "depth": depth,
            "sites": len(cols),
            "success": True,
            "error": None,
        }
    except Exception as exc:  # pragma: no cover - surfaced to caller
        return {
            "region": f"{chrom}:{start}-{end}",
            "rows": [],
            "depth": 0,
            "sites": 0,
            "success": False,
            "error": str(exc),
        }


def _render(cols) -> list[list]:
    mode = _WORKER_ECFG.mode
    strands = _WORKER_STRANDS
    if mode == "mutation":
        return list(mutation_rows(cols, _WORKER_MCFG, strands))
    if mode == "base":
        return list(
            base_rows(
                cols,
                split_strand=_WORKER_SCFG.split,
                count_indels=_WORKER_ECFG.count_indels,
                strands=strands,
                min_depth=_WORKER_ECFG.min_depth,
            )
        )
    # allele table / VCF
    return list(
        allele_rows(
            cols,
            min_support=_WORKER_ECFG.min_allele_support,
            min_depth=_WORKER_ECFG.min_depth,
            vcf=_WORKER_ECFG.vcf,
        )
    )


# ---------------------------------------------------------------------------
# public result object
# ---------------------------------------------------------------------------
@dataclass
class PipelineResult:
    header: list[str]
    rows: list[list] = field(default_factory=list)
    total_sites: int = 0
    total_depth: int = 0
    elapsed: float = 0.0
    error: str | None = None
    success: bool = True


def pick_engine(mode: str, requested: str) -> str:
    """Resolve the 'auto' engine selection based on mode."""
    if requested in ("read-walk", "pileup"):
        return requested
    # auto
    return "read-walk" if mode == "mutation" else "pileup"


def _detect_bisulfite(sam: pysam.AlignmentFile) -> bool:
    """Peek at the first reads to see if Yf/Zf bisulfite tags exist."""
    count = 0
    for read in sam.fetch():
        count += 1
        if read.has_tag("Yf") or read.has_tag("Zf"):
            return True
        if count >= 1000:
            break
    return False


def _build_bins(
    bam_chroms: list[str],
    ref_lengths: dict[str, int],
    bin_size: int,
    region: str | None,
) -> list[tuple[str, int, int]]:
    bins: list[tuple[str, int, int]] = []
    if region and ":" in region and "-" in region:
        chrom, pos_range = region.split(":")
        s, e = map(int, pos_range.split("-"))
        bins.append((chrom, s - 1, e))
        return bins
    if region:
        raise ValueError(f"invalid region '{region}'; use 'chr:start-end'")
    for chrom in bam_chroms:
        length = ref_lengths.get(chrom, 0)
        if length == 0:
            continue
        for bstart in range(0, length, bin_size):
            bins.append((chrom, bstart, min(bstart + bin_size, length)))
    return bins


def run_pipeline(
    samfile: str,
    reference: str,
    output: str | None = None,
    *,
    fcfg: FilterConfig | None = None,
    mcfg: MutationConfig | None = None,
    scfg: StrandConfig | None = None,
    ecfg: EngineConfig | None = None,
    verbose: bool = False,
) -> PipelineResult:
    """Run the unified pipeline and return a :class:`PipelineResult`.

    This is the single entry point used by the CLI.  It decides the engine,
    bins the genome, processes bins in parallel, aggregates the rows, and writes
    them (or returns them for stdout streaming).
    """
    start = time.time()
    fcfg = fcfg or FilterConfig()
    scfg = scfg or StrandConfig()
    ecfg = ecfg or EngineConfig()

    # ---- indices -------------------------------------------------------------
    if not os.path.exists(samfile + ".bai"):
        pysam.index(samfile)
    if not os.path.exists(reference + ".fai"):
        pysam.faidx(reference)

    threads = ecfg.threads or min(os.cpu_count() or 1, 8)

    # ---- build bins ----------------------------------------------------------
    ref_lengths = read_fasta_index(reference)
    with pysam.AlignmentFile(samfile, "rb") as sam:
        bam_chroms = list(sam.references)
        has_bisulfite = _detect_bisulfite(sam) if ecfg.mode == "mutation" else False
    bins = _build_bins(bam_chroms, ref_lengths, ecfg.bin_size or 10_000, ecfg.region)

    # ---- header --------------------------------------------------------------
    header = _header_for(mode=ecfg.mode, mcfg=mcfg, scfg=scfg, ecfg=ecfg)

    # ---- execute -------------------------------------------------------------
    all_rows: list[list] = []
    total_sites = 0
    total_depth = 0
    # Chromosome ordering from the BAM header for deterministic output.
    chrom_order = {c: i for i, c in enumerate(bam_chroms)}

    with ProcessPoolExecutor(
        max_workers=threads,
        initializer=_init_worker,
        initargs=(
            samfile,
            reference,
            fcfg,
            mcfg,
            scfg,
            ecfg,
            scfg.process,
            has_bisulfite,
        ),
    ) as pool:
        futures = [pool.submit(_worker, b) for b in bins]
        for fut in as_completed(futures):
            res = fut.result()
            if not res["success"]:
                if verbose:
                    sys.stderr.write(
                        f"[countmut] region {res['region']}: {res['error']}\n"
                    )
                continue
            all_rows.extend(res["rows"])
            total_sites += res["sites"]
            total_depth += res["depth"]

    # ---- sort + write --------------------------------------------------------
    all_rows.sort(key=lambda r: (chrom_order.get(r[0], 1 << 30), r[1]))
    _write(
        header,
        all_rows,
        output,
        vcf=ecfg.mode == "allele" and ecfg.vcf,
    )
    return PipelineResult(
        header=header,
        rows=all_rows,
        total_sites=total_sites,
        total_depth=total_depth,
        elapsed=time.time() - start,
        success=True,
    )


def _header_for(mode: str, mcfg, scfg, ecfg) -> list[str]:
    if mode == "mutation":
        return mutation_header(mcfg.save_rest if mcfg else False)
    if mode == "base":
        if scfg.split and ecfg.count_indels:
            return BASE_HEADER_SPLIT_INDEL
        if scfg.split:
            return BASE_HEADER_SPLIT
        if ecfg.count_indels:
            return BASE_HEADER_INDEL
        return BASE_HEADER
    if ecfg.vcf:
        return [
            "##fileformat=VCFv4.2",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
        ]
    return ["chrom", "pos", "ref", "depth", "ref_count", "alt", "alt_count"]


def _write(
    header: list[str], rows: list, output: str | None, *, vcf: bool = False
) -> None:
    if vcf:
        # header holds whole VCF lines; rows are pre-rendered VCF records.
        lines: list[str] = list(header)
        lines += [r if isinstance(r, str) else "\t".join(map(str, r)) for r in rows]
    else:
        lines = ["\t".join(header), *("\t".join(map(str, r)) for r in rows)]
    text = "\n".join(lines) + "\n"
    if output:
        os.makedirs(os.path.dirname(os.path.abspath(output)), exist_ok=True)
        with open(output, "w") as f:
            f.write(text)
    else:
        sys.stdout.write(text)
