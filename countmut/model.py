#!/usr/bin/env python3
"""
Shared counting model + configuration for the unified countmut core.

This module is the single source of truth for *what* is counted. Both engines
(``engine_readwalk`` and ``engine_pileup``) produce :class:`SiteColumn` objects
with exactly this shape, so the two "ways" of walking a BAM (read-by-read vs
pileup-by-position) yield *identical* output.  The semantics here are distilled
line-by-line from:

* ``countmut`` (read-walk, bisulfite NS/Zf/Yf tiers, biological strand, motifs)
* ``perbase`` / ``pbr`` (mate-aware overlap dedup, PileupPosition: a/c/g/t/n,
  ins/del/ref_skip/fail, depth, per-strand)
* ``minipileup`` / ``mpileup`` / ``cpup`` (allele ordering, read- and base-level
  filters, strand switching)

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable

# ---------------------------------------------------------------------------
# DNA helpers
# ---------------------------------------------------------------------------
DNA_COMPLEMENT = str.maketrans("ATGCNatgcn", "TACGNtacgn")
DNA_ALPHABET = ("A", "C", "G", "T", "N")


def reverse_complement(seq: str) -> str:
    """Return the reverse complement of a DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def complement(base: str) -> str:
    """Return the complement of a single base."""
    return base.translate(DNA_COMPLEMENT)


# ---------------------------------------------------------------------------
# Quality / conversion categories
# ---------------------------------------------------------------------------
# Mutation-mode quality tiers, in output order (x0, x1, x2).
LOW_QUALITY = "low_quality"        # x0 -- fails base-level filter (trim/qual/sub)
INSUFFICIENT = "insufficient"      # x1 -- passes quality, fails read conversion filter
HIGH_CONVERSION = "high"           # x2 -- passes quality + conversion filter
MUTATION_CATEGORIES = (LOW_QUALITY, INSUFFICIENT, HIGH_CONVERSION)

# Single category used in base-count / allele modes (no conversion tiers).
BASE_CATEGORY = "base"

CATEGORY_COLUMNS = {
    LOW_QUALITY: "0",
    INSUFFICIENT: "1",
    HIGH_CONVERSION: "2",
}


# ---------------------------------------------------------------------------
# Configuration dataclasses
# ---------------------------------------------------------------------------
@dataclass
class FilterConfig:
    """Read- and base-level filters (mirrors countmut/perbase/minipileup)."""

    min_mapq: int = 0
    min_baseq: int = 20
    max_sub: int = 1            # max substitutions (NS tag); read-level
    max_unc: int | None = 3     # max unconverted (Zf tag); read-level, None=ignore
    min_con: int | None = 1     # min converted (Yf tag); read-level, None=ignore
    trim_start: int = 2         # bases trimmed from fragment 5' end
    trim_end: int = 2           # bases trimmed from fragment 3' end
    include_flags: int = 0      # SAM flags that must be set
    exclude_flags: int = 4 | 256 | 512  # unmapped, secondary, duplicate
    max_depth: int = 0          # per-position depth cap (0 = unlimited)


@dataclass
class MutationConfig:
    """Which substitution to count and how many reference bases to report."""

    ref_base: str = "A"
    mut_base: str = "G"
    ref_base2: str | None = None  # alternative ref base for Yc/Zc tagging
    mut_base2: str | None = None
    pad: int = 15                 # motif half-window
    save_rest: bool = False        # also emit o0/o1/o2 (other bases)


@dataclass
class StrandConfig:
    """Which biological strand(s) to report."""

    process: str = "both"  # 'both' | 'forward' | 'reverse'
    split: bool = True     # emit separate '+'/'-' rows (False = sum both strands)


@dataclass
class EngineConfig:
    """Selects the BAM walk strategy."""

    engine: str = "auto"          # 'auto' | 'read-walk' | 'pileup'
    mode: str = "mutation"        # 'mutation' | 'base' | 'allele'
    bin_size: int = 10_000
    threads: int | None = None
    region: str | None = None
    save_rest: bool = False
    # allele / vcf options
    min_allele_support: int = 1
    min_allele_fraction: float = 0.0
    min_strand_support: int = 0
    min_allele_depth: int = 0
    min_mean_depth: int = 0
    min_depth: int = 0              # min site depth to report (base/allele)
    vcf: bool = False               # allele mode: emit VCF
    report_reference_bases: bool = False  # -k flanking window (pbr) -> also motif
    flanking: int = 0
    count_indels: bool = False
    split_strand: bool = True
    read_expr: str | None = None    # pbr -e read-string filter
    pile_expr: str | None = None    # pbr -p pileup-string filter


# ---------------------------------------------------------------------------
# SiteColumn: the neutral per-site accumulator emitted by both engines
# ---------------------------------------------------------------------------
@dataclass
class SiteColumn:
    """Per-site, per-strand base counts (+ special events) for one genomic position.

    Both engines fill exactly this object, so downstream formatting is shared.
    ``counts[strand][category][base]`` holds observed bases binned by strand and
    quality/conversion category.  ``ins/del/ref_skip/fail`` are strand-keyed
    counters for the events that do not contribute a base.
    """

    chrom: str
    pos: int                                   # 0-based genomic coordinate
    ref_base: str                              # reference allele (forward strand)
    motif: str = ""                            # pad-flanked reference window
    counts: dict[str, dict[str, Counter]] = field(default_factory=dict)
    ins: dict[str, int] = field(default_factory=dict)
    deletes: dict[str, int] = field(default_factory=dict)
    ref_skip: dict[str, int] = field(default_factory=dict)
    fail: dict[str, int] = field(default_factory=dict)

    # lifecycle ------------------------------------------------------------
    @classmethod
    def make(
        cls,
        chrom: str,
        pos: int,
        ref_base: str,
        *,
        motif: str = "",
        categories: Iterable[str] = MUTATION_CATEGORIES,
        strands: Iterable[str] = ("+", "-"),
    ) -> "SiteColumn":
        """Create a SiteColumn pre-populated with zeroed counters."""
        col = cls(chrom=chrom, pos=pos, ref_base=ref_base, motif=motif)
        for s in strands:
            col.counts[s] = {cat: Counter() for cat in categories}
            col.ins[s] = 0
            col.deletes[s] = 0
            col.ref_skip[s] = 0
            col.fail[s] = 0
        return col

    # mutation --------------------------------------------------------------
    def add_observation(self, strand: str, base: str, category: str) -> None:
        """Record one observed base under the given strand + category."""
        base = base.upper()
        bucket = self.counts.setdefault(strand, {}).get(category)
        if bucket is None:
            bucket = self.counts.setdefault(strand, {})[category] = Counter()
        bucket[base] += 1

    # base/allele mode special events ---------------------------------------
    def add_indel(self, strand: str, kind: str, n: int = 1) -> None:
        target = self.ins if kind == "ins" else self.deletes
        target[strand] = target.get(strand, 0) + n

    def add_ref_skip(self, strand: str, n: int = 1) -> None:
        self.ref_skip[strand] = self.ref_skip.get(strand, 0) + n

    def add_fail(self, strand: str, n: int = 1) -> None:
        self.fail[strand] = self.fail.get(strand, 0) + n

    # accessors -------------------------------------------------------------
    def strand_depth(self, strand: str) -> int:
        """Number of observations that produced a base on a strand."""
        total = 0
        for cat in self.counts.get(strand, {}).values():
            total += sum(cat.values())
        return total

    def total_depth(self) -> int:
        return sum(self.strand_depth(s) for s in self.counts)

    def is_empty(self) -> bool:
        return self.total_depth() <= 0 and all(
            not self.ref_skip.get(s, 0) and not self.deletes.get(s, 0)
            for s in self.counts
        )

    def has_data(self, strand: str | None = None) -> bool:
        """True if this column (or one strand of it) actually observed anything."""
        strands = [strand] if strand is not None else list(self.counts)
        for s in strands:
            for cat in self.counts.get(s, {}).values():
                if sum(cat.values()) > 0:
                    return True
            if self.ins.get(s, 0) or self.deletes.get(s, 0) \
               or self.ref_skip.get(s, 0) or self.fail.get(s, 0):
                return True
        return False
