#!/usr/bin/env python3
"""
Pileup-based countmut engine.

Walks a BAM *position by position* with pysam's bundled htslib pileup engine
(the same walk used by minipileup / perbase / pbr / mpileup).  At each position
it groups the covering alignments by query name, collapses overlapping mates via
:func:`countmut.reads.obs_key`, and fills a :class:`countmut.model.SiteColumn`.

This is "way one" (pileup-based).  The read-walk engine fills identical
SiteColumn objects, making the two BAM-walk strategies interchangeable.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import pysam

from . import reads
from .model import (
    BASE_CATEGORY,
    HIGH_CONVERSION,
    LOW_QUALITY,
    MUTATION_CATEGORIES,
    FilterConfig,
    MutationConfig,
    SiteColumn,
)


def _target_sites(
    reference: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    ref_base: str | None,
) -> set[int] | None:
    """0-based positions in [start,end) whose ref base matches ``ref_base``.

    Returns None in base-count/allele mode (every covered position is a target).
    """
    if ref_base is None:
        return None
    seq = reference.fetch(chrom, max(start, 0), end)
    return {start + i for i, b in enumerate(seq) if b.upper() == ref_base}


def pileup_region(
    sam: pysam.AlignmentFile,
    reference: pysam.FastaFile,
    chrom: str,
    start: int,
    end: int,
    fcfg: FilterConfig,
    mcfg: MutationConfig | None = None,
    mode: str = "mutation",
    has_bisulfite_tags: bool = False,
    read_pred=None,  # optional per-base read predicate (-e)
    pile_pred=None,  # optional per-site pileup predicate (-p)
) -> list[SiteColumn]:
    """Count one region (half-open [start,end)) via a position-by-position pileup.

    Returns :class:`SiteColumn` objects for every position we care about.
    """
    is_mutation = mode == "mutation"
    targets = _target_sites(
        reference, chrom, start, end, mcfg.ref_base if is_mutation else None
    )
    categories = MUTATION_CATEGORIES if is_mutation else (BASE_CATEGORY,)

    # Reference window padded by `pad` on each side for motif extraction.
    pad = mcfg.pad if mcfg else 0
    left = reference.fetch(chrom, max(start - pad, 0), start)
    left = "N" * (pad - len(left)) + left if len(left) < pad else left
    right = reference.fetch(chrom, end, end + pad).ljust(pad, "N")
    ext_seq = left + reference.fetch(chrom, start, end) + right

    columns: list[SiteColumn] = []
    try:
        pileups = sam.pileup(
            chrom,
            max(start, 0),
            end,
            truncate=True,
            max_depth=fcfg.max_depth or 8000,
            min_base_quality=0,
        )
    except TypeError:  # older pysam without all kwargs
        pileups = sam.pileup(chrom, max(start, 0), end, truncate=True)

    for p in pileups:
        pos = p.reference_pos
        if not (start <= pos < end):
            continue
        if is_mutation and pos not in targets:
            continue
        ref_base = reference.fetch(chrom, pos, pos + 1).upper()

        col = SiteColumn.make(chrom, pos, ref_base, categories=categories)
        if is_mutation:
            idx = pos - start + pad
            col.motif = ext_seq[idx - pad : idx + pad + 1]

        # Collapse overlapping mates.  CRITICAL canonical rule (consistent with
        # engine_readwalk): a read is only a dedup candidate if it passes the
        # read-level filters AND its base here is internal (post-trim).  Deletions,
        # ref-skips and filter-failures are tallied per-read and are NOT bases.
        best: dict[str, tuple] = {}
        for pr in p.pileups:
            rec = pr.alignment
            strand = reads.actual_strand(rec)
            if reads.read_fail_reason(rec, fcfg, has_bisulfite_tags) is not None:
                col.add_fail(strand)
                continue
            qpos = pr.query_position
            if qpos is None:
                if pr.is_refskip:
                    col.add_ref_skip(strand)
                elif pr.is_del:
                    col.add_indel(strand, "del")
                continue
            qlen = len(rec.query_sequence) if rec.query_sequence is not None else 0
            if not reads.is_internal(
                qpos, qlen, strand, fcfg.trim_start, fcfg.trim_end
            ):
                continue  # trimmed: not counted, not a candidate
            if read_pred is not None and not read_pred(rec, qpos):
                continue  # -e expression filter rejected this base
            qual = (
                int(rec.query_qualities[qpos]) if rec.query_qualities is not None else 0
            )
            key = reads.obs_key(rec.mapping_quality, rec.is_read1, qual)
            _remember(best, key, pr, rec, qpos)

        for _key, (pr, rec, qpos) in best.values():
            strand = reads.actual_strand(rec)
            # BAM stores SEQ in reference-forward orientation, so the base is
            # already reference-forward here (verified on real aligner BAMs).
            base = rec.query_sequence[qpos].upper()
            qual = (
                int(rec.query_qualities[qpos]) if rec.query_qualities is not None else 0
            )
            if is_mutation:
                category = HIGH_CONVERSION if qual >= fcfg.min_baseq else LOW_QUALITY
            else:
                category = BASE_CATEGORY
            col.add_observation(strand, base, category)

        if pile_pred is not None and not pile_pred(col):
            continue
        columns.append(col)

    return columns


def _remember(best: dict, key: tuple, pr, rec, qpos: int | None) -> None:
    """Keep the observation with the largest preference key per qname."""
    cur = best.get(rec.query_name)
    if cur is None or key > cur[0]:
        best[rec.query_name] = (key, (pr, rec, qpos))
