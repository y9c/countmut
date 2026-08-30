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
    """Read- and base-level filters (mirrors countmut/perbase/minipileup)."""

    min_mapq: int = 0
    min_baseq: int = 20
    max_sub: int = 1  # max substitutions (NS tag); read-level
    max_unc: int | None = 3  # max unconverted (Zf tag); read-level, None=ignore
    min_con: int | None = 1  # min converted (Yf tag); read-level, None=ignore
    trim_start: int = 2  # bases trimmed from fragment 5' end
    trim_end: int = 2  # bases trimmed from fragment 3' end
    include_flags: int = 0  # SAM flags that must be set
    exclude_flags: int = 1796  # UNMAP|SECONDARY|QCFAIL|DUP (samtools default)
    max_depth: int = 0  # per-position depth cap (0 = unlimited)


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
    split: bool = True  # emit separate '+'/'-' rows (False = sum both strands)


@dataclass
class EngineConfig:
    """Selects the BAM walk strategy (both are implemented in the C core)."""

    engine: str = "auto"  # 'auto' | 'read-walk' | 'pileup'
    mode: str = "mutation"  # 'mutation' | 'base' | 'allele'
    bin_size: int = 10_000
    threads: int | None = None
    region: str | None = None
    save_rest: bool = False
    # Lua filter expressions (evaluated inside the C core, pbr-style)
    read_expr: str | None = None  # -e read-level filter (per aligned base)
    pile_expr: str | None = None  # -p site-level filter (per site)
    # allele / vcf options
    min_allele_support: int = 1
    min_allele_fraction: float = 0.0
    min_strand_support: int = 0
    min_allele_depth: int = 0
    min_mean_depth: int = 0
    min_depth: int = 0  # min site depth to report (base/allele)
    vcf: bool = False  # allele mode: emit VCF
    count_indels: bool = False
    split_strand: bool = True
    verbose: bool = False  # real-time per-region progress on stderr
    report_reference_bases: bool = False  # -k flanking window (pbr)
    flanking: int = 0
