#!/usr/bin/env python3
"""
Rendering of SiteColumn objects into the various countmut outputs.

Several renderers are provided, one per ``--mode``:

* ``mutation`` -- countmut's strand-aware bisulfite substitution table
                   (``chrom pos strand motif u0 u1 u2 m0 m1 m2 [o0 o1 o2]``)
* ``base``     -- per-site depth / base counts (perbase / mpileup style)
* ``allele``   -- per-sample allele counts + optional VCF (minipileup style)

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

from collections.abc import Iterable

from .model import (
    HIGH_CONVERSION,
    INSUFFICIENT,
    LOW_QUALITY,
    MutationConfig,
    SiteColumn,
)

MUTATION_HEADER = [
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
MUTATION_HEADER_REST = MUTATION_HEADER + ["o0", "o1", "o2"]

BASE_HEADER = ["chrom", "pos", "ref", "depth", "a", "c", "g", "t", "n"]
BASE_HEADER_SPLIT = ["chrom", "pos", "strand", "ref", "depth", "a", "c", "g", "t", "n"]
BASE_HEADER_INDEL = BASE_HEADER + ["ins", "del", "ref_skip", "fail"]
BASE_HEADER_SPLIT_INDEL = BASE_HEADER_SPLIT + ["ins", "del", "ref_skip", "fail"]


# ---------------------------------------------------------------------------
# mutation mode
# ---------------------------------------------------------------------------
def mutation_header(save_rest: bool) -> list[str]:
    return MUTATION_HEADER_REST if save_rest else MUTATION_HEADER


def mutation_row(col: SiteColumn, mcfg: MutationConfig, strand: str) -> list:
    ref_base = mcfg.ref_base
    mut_base = mcfg.mut_base
    counts = col.counts.get(strand, {})
    lq = counts.get(LOW_QUALITY, {})
    insuf = counts.get(INSUFFICIENT, {})
    high = counts.get(HIGH_CONVERSION, {})

    def pick(cat: dict, base: str) -> int:
        return int(cat.get(base, 0))

    def others(cat: dict, ref: int, mut: int) -> int:
        return max(0, sum(cat.values()) - ref - mut)

    u0, m0 = pick(lq, ref_base), pick(lq, mut_base)
    u1, m1 = pick(insuf, ref_base), pick(insuf, mut_base)
    u2, m2 = pick(high, ref_base), pick(high, mut_base)
    o0, o1, o2 = others(lq, u0, m0), others(insuf, u1, m1), others(high, u2, m2)

    # All engines store bases in reference-forward orientation (BAM SEQ is
    # reference-forward regardless of read strand), so u/m counts are strand-
    # agnostic and the motif is the reference-forward window for BOTH strands.
    # (Reverse-complementing it for '-' would contradict the reference-forward
    # bases we are counting.)
    motif = col.motif

    row = [col.chrom, col.pos + 1, strand, motif, u0, u1, u2, m0, m1, m2]
    if mcfg.save_rest:
        row += [o0, o1, o2]
    return row


def mutation_rows(
    columns: Iterable[SiteColumn], mcfg: MutationConfig, strands=("+", "-")
):
    """Yield TSV rows for mutation mode following countmut's inclusion rule."""
    for col in columns:
        for strand in strands:
            counts = col.counts.get(strand, {})
            high = counts.get(HIGH_CONVERSION, {})
            insuf = counts.get(INSUFFICIENT, {})
            # countmut only emits a site if any (insufficient or high) ref+mut base was seen
            if (
                insuf.get(mcfg.ref_base, 0)
                + insuf.get(mcfg.mut_base, 0)
                + high.get(mcfg.ref_base, 0)
                + high.get(mcfg.mut_base, 0)
                > 0
            ):
                yield mutation_row(col, mcfg, strand)


# ---------------------------------------------------------------------------
# base mode
# ---------------------------------------------------------------------------
def _base_counts(col: SiteColumn, strand: str | None) -> dict:
    """Aggregate base counts, optionally restricted to one strand."""
    agg = {"a": 0, "c": 0, "g": 0, "t": 0, "n": 0}
    strands = [strand] if strand is not None else list(col.counts)
    for s in strands:
        for cat in col.counts.get(s, {}).values():
            for base, count in cat.items():
                key = base.lower()
                if key in agg:
                    agg[key] += count
                else:
                    agg["n"] += count
    return agg


def _special(col: SiteColumn, strand: str | None, kind: str) -> int:
    strands = [strand] if strand is not None else list(col.counts)
    table = {
        "ins": col.ins,
        "del": col.deletes,
        "ref_skip": col.ref_skip,
        "fail": col.fail,
    }
    return sum(table[kind].get(s, 0) for s in strands)


def base_rows(
    columns: Iterable[SiteColumn],
    *,
    split_strand: bool = False,
    count_indels: bool = False,
    strands=("+", "-"),
    min_depth: int = 0,
):
    """Yield TSV rows for base-count mode."""
    for col in columns:
        ref = col.ref_base.upper()
        if split_strand:
            for strand in strands:
                if not col.has_data(strand):
                    continue
                agg = _base_counts(col, strand)
                depth = sum(agg.values())
                if min_depth > 0 and depth < min_depth:
                    continue
                row = [
                    col.chrom,
                    col.pos + 1,
                    strand,
                    ref,
                    depth,
                    agg["a"],
                    agg["c"],
                    agg["g"],
                    agg["t"],
                    agg["n"],
                ]
                if count_indels:
                    row += [
                        _special(col, strand, "ins"),
                        _special(col, strand, "del"),
                        _special(col, strand, "ref_skip"),
                        _special(col, strand, "fail"),
                    ]
                yield row
        else:
            if not col.has_data(None):
                continue
            agg = _base_counts(col, None)
            depth = sum(agg.values())
            if min_depth > 0 and depth < min_depth:
                continue
            row = [
                col.chrom,
                col.pos + 1,
                ref,
                depth,
                agg["a"],
                agg["c"],
                agg["g"],
                agg["t"],
                agg["n"],
            ]
            if count_indels:
                row += [
                    _special(col, None, "ins"),
                    _special(col, None, "del"),
                    _special(col, None, "ref_skip"),
                    _special(col, None, "fail"),
                ]
            yield row


# ---------------------------------------------------------------------------
# allele mode
# ---------------------------------------------------------------------------
def allele_counts(col: SiteColumn) -> list[tuple[str, int]]:
    """Return [(allele, support)] sorted by descending support for one site."""
    agg = _base_counts(col, None)
    pairs = [(b, agg[b.lower()]) for b in ("A", "C", "G", "T", "N")]
    pairs = [(b, n) for b, n in pairs if n > 0]
    pairs.sort(key=lambda x: (-x[1], x[0]))
    return pairs


def allele_rows(
    columns: Iterable[SiteColumn],
    ref_allele: str | None = None,
    min_support: int = 1,
    min_depth: int = 0,
    vcf: bool = False,
):
    """Yield rows matching the C backend.

    Table mode: ``chrom pos ref depth ref_count alt alt_count`` (alt='.' if no
    alt found).  VCF mode (``vcf=True``): a VCFv4.2 record per line, identical
    to the C core's ``chrom pos . ref alt . PASS . GT:AD 0/1:ref,alt``.
    """
    for col in columns:
        ref = (ref_allele or col.ref_base.upper()).upper()
        agg = _base_counts(col, None)
        depth = sum(agg.values())
        if min_depth > 0 and depth < min_depth:
            continue
        if depth <= 0:
            continue
        ref_n = agg.get(ref.lower(), 0)
        best, bn = ".", 0
        for b in ("A", "C", "G", "T"):
            if b != ref and b != "N":
                n = agg.get(b.lower(), 0)
                if n > bn:
                    bn, best = n, b
        if bn < min_support:
            best, bn = ".", 0
        if vcf:
            if best == ".":
                continue
            yield (
                f"{col.chrom}\t{col.pos + 1}\t.\t{ref}\t{best}\t.\tPASS"
                f"\t.\tGT:AD\t0/1:{ref_n},{bn}"
            )
        else:
            yield [col.chrom, col.pos + 1, ref, depth, ref_n, best, bn]
