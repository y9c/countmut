#!/usr/bin/env python3
"""
Configuration for the countmut C core.

The heavy lifting (both the read-walk and pileup engines) lives in the C
backend (``backend/countmut_core``); this module only holds the dataclasses the
Python wrapper uses to build the C command line.  Semantics are distilled from
the original countmut (bisulfite NS/Zf/Yf tiers, biological strand, trim) and
perbase / pbr / minipileup (mate-aware overlap dedup, base counts, indels).

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass
class FilterConfig:
    """Only the computational depth cap is a real setting here.

    All read-level filtering and trimming belongs in the ``-e`` read expression
    (``mapq``, ``bq``, ``tag('NS')/tag('Zf')/tag('Yf')/tag('NM')``,
    ``dist5``/``dist3``/``qpos``/``flags``) and all site-level filtering in the
    ``-p`` pile expression (``depth``, ``ref``, ``a/c/g/t/n``, ``pos``, ...) --
    not as dedicated flags.
    """

    max_depth: int = 0  # per-position depth cap (pileup engine; 0 = unlimited)


@dataclass
class MutationConfig:
    """Which substitution to count and how many reference bases to report."""

    ref_base: str = "A"
    mut_base: str = "G"
    ref_base2: str | None = None  # alternative ref base for Yc/Zc tagging
    mut_base2: str | None = None
    pad: int = 15  # motif half-window
    save_rest: bool = False  # also emit o0/o1/o2 (other bases)

    def __post_init__(self) -> None:
        # Case-normalize like the C core (driver uppercases --ref-base/--mut-base).
        self.ref_base = self.ref_base.upper()
        self.mut_base = self.mut_base.upper()
        if self.ref_base2 is not None:
            self.ref_base2 = self.ref_base2.upper()
        if self.mut_base2 is not None:
            self.mut_base2 = self.mut_base2.upper()


@dataclass
class StrandConfig:
    """Which biological strand(s) to report."""

    process: str = "both"  # 'both' | 'forward' | 'reverse'
    strandless: bool = False  # True = sum both strands into one row (base/allele)


@dataclass
class EngineConfig:
    """Selects the BAM walk strategy (both are implemented in the C core)."""

    engine: str = "auto"  # 'auto' | 'read-walk' | 'pileup'
    mode: str = "mutation"  # 'auto' | 'mutation' | 'base' | 'allele'
    bin_size: int = 10_000
    threads: int | None = None
    region: str | None = None
    save_rest: bool = False
    # Lua filter expressions (evaluated inside the C core, pbr-style)
    read_expr: str | None = None  # -e read-level filter (per aligned base)
    pile_expr: str | None = None  # -p site-level filter (per site)
    # allele / output
    vcf: bool = False  # allele mode: emit VCF
    count_indels: bool = False
    strandless: bool = False  # base/allele: collapse '+'/'-' into one row
    verbose: bool = False  # real-time per-region progress on stderr
