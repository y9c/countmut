#!/usr/bin/env python3
"""
Read-walk countmut engine (the "non-pileup" way).

Iterates a BAM *read by read* (like countmut's core), walks each read's aligned
query/reference pairs, and only touches the target sites of interest.  This is
the most efficient strategy when you care about a limited set of positions
(e.g. every reference base equal to ``ref_base`` in bisulfite analysis).

It fills identical :class:`SiteColumn` objects to :mod:`countmut.engine_pileup`,
so the two "ways" produce interchangeable output.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import pysam

from . import reads
from .model import (
    BASE_CATEGORY,
    DNA_COMPLEMENT,
    FilterConfig,
    HIGH_CONVERSION,
    LOW_QUALITY,
    MUTATION_CATEGORIES,
    MutationConfig,
    SiteColumn,
)


def _target_sites(reference, chrom, start, end, ref_base: str | None) -> set[int] | None:
    if ref_base is None:
        return None
    seq = reference.fetch(chrom, max(start, 0), end)
    return {start + i for i, b in enumerate(seq) if b.upper() == ref_base}


def readwalk_region(
    sam: pysam.AlignmentFile,
    reference: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    fcfg: FilterConfig,
    mcfg: MutationConfig | None = None,
    mode: str = "mutation",
    strand_process: str = "both",
    has_bisulfite_tags: bool = False,
    read_pred=None,  # optional per-base read predicate (-e)
    pile_pred=None,  # optional per-site pileup predicate (-p)
) -> list[SiteColumn]:
    """Count one region (half-open [start,end)) by walking reads directly."""
    is_mutation = mode == "mutation"
    targets = _target_sites(reference, chrom, start, end, mcfg.ref_base if is_mutation else None)

    pad = mcfg.pad if mcfg else 0
    left = reference.fetch(chrom, max(start - pad, 0), start)
    left = "N" * (pad - len(left)) + left if len(left) < pad else left
    right = reference.fetch(chrom, end, end + pad).ljust(pad, "N")
    ext_seq = left + reference.fetch(chrom, start, end) + right

    # best[(pos, qname)] = (key, strand, base, qual, category)
    best: dict[tuple[int, str], tuple] = {}
    for read in sam.fetch(chrom, start, end):
        strand = reads.actual_strand(read)
        if strand_process == "forward" and strand != "+":
            continue
        if strand_process == "reverse" and strand != "-":
            continue
        if reads.read_fail_reason(read, fcfg, has_bisulfite_tags) is not None:
            continue

        qs = read.query_sequence
        qq = read.query_qualities
        if not qs:
            continue
        qlen = len(qs)
        qname = read.query_name
        mapq = read.mapping_quality
        is_read1 = read.is_read1

        for qpos, ref_pos in read.get_aligned_pairs(matches_only=True):
            if qpos is None or ref_pos is None:
                continue
            if targets is not None and ref_pos not in targets:
                continue
            if not (start <= ref_pos < end):
                continue
            if qpos >= qlen:
                continue
            if not reads.is_internal(qpos, qlen, strand, fcfg.trim_start, fcfg.trim_end):
                continue
            if read_pred is not None and not read_pred(read, qpos):
                continue

            # BAM stores SEQ in reference-forward orientation (verified on real
            # aligner BAMs), so the base is already reference-forward here.
            base = qs[qpos].upper()
            qual = int(qq[qpos]) if qq is not None else 0

            if is_mutation:
                category = HIGH_CONVERSION if qual >= fcfg.min_baseq else LOW_QUALITY
            else:
                category = BASE_CATEGORY

            key = reads.obs_key(mapq, is_read1, qual)
            cur = best.get((ref_pos, qname))
            if cur is None or key > cur[0]:
                best[(ref_pos, qname)] = (key, strand, base, qual, category)

    # Flush winners into per-position buckets.
    pos_buckets: dict[int, dict[str, dict[str, dict[str, int]]]] = {}
    for (ref_pos, qname), (_key, strand, base, qual, category) in best.items():
        cats = pos_buckets.setdefault(ref_pos, {})
        strands = cats.setdefault(category, {})
        bases = strands.setdefault(strand, {})
        bases[base] = bases.get(base, 0) + 1

    columns: list[SiteColumn] = []
    for ref_pos in sorted(pos_buckets):
        if targets is not None and ref_pos not in targets:
            continue
        ref_base = reference.fetch(chrom, ref_pos, ref_pos + 1).upper()
        col = SiteColumn.make(
            chrom, ref_pos, ref_base,
            categories=MUTATION_CATEGORIES if is_mutation else (BASE_CATEGORY,),
        )
        if is_mutation:
            idx = ref_pos - start + pad
            col.motif = ext_seq[idx - pad: idx + pad + 1]
        for category, strands in pos_buckets[ref_pos].items():
            for strand, bases in strands.items():
                for base, count in bases.items():
                    for _ in range(count):
                        col.add_observation(strand, base, category)
        if pile_pred is not None and not pile_pred(col):
            continue
        columns.append(col)
    return columns
