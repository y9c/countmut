#!/usr/bin/env python3
"""
Shared read-handling helpers for the countmut engines.

Both the read-walk engine and the pileup engine must classify reads *identically*:

* biological strand (``actual_strand``) -- from countmut
* which query positions are "internal" after 5'/3' trimming (fragment orient) --
  from countmut
* which observation a fragment contributes at a site when mates overlap --
  collapse via ``obs_key`` (from perbase/pbr mate-aware dedup, refined to tie-break
  by base-quality like countmut)

Keeping these in one place is what makes the two BAM-walk strategies emit
identical output.

Author: Ye Chang
Date: 2026-08-30
"""

from __future__ import annotations

import pysam

from .model import (
    FilterConfig,
)


# Mapping quality / first-in-pair / base-quality preference tuple.
# We keep the observation with the *largest* key.  Higher MAPQ wins, then a
# read1 (first-in-pair) wins over a read2, then higher base-quality wins.
def obs_key(mapq: int, is_read1: bool, base_qual: int) -> tuple:
    return (int(mapq), 1 if is_read1 else 0, int(base_qual))


def actual_strand(read: pysam.AlignedSegment) -> str:
    """Biological strand ('+' / '-') of a read (countmut semantics).

    Paired-end: read1 forward = '+', read2 reverse-complemented = '+'.
    """
    if read.is_paired:
        if read.is_read1:
            return "+" if not read.is_reverse else "-"
        return "+" if read.is_reverse else "-"
    return "+" if not read.is_reverse else "-"


def is_internal(
    query_pos: int, query_len: int, strand: str, trim_start: int, trim_end: int
) -> bool:
    """Is this query position inside the trimmed fragment (countmut semantics)?"""
    if strand == "+":
        return query_pos >= trim_start and query_len - query_pos > trim_end
    return query_pos >= trim_end and query_len - query_pos > trim_start


def read_fail_reason(
    read: pysam.AlignedSegment,
    fcfg: FilterConfig,
    has_bisulfite_tags: bool,
) -> str | None:
    """Return a skip reason if the whole read should be dropped, else None.

    Filter order mirrors countmut: unmapped/duplicate/secondary flags, then MAPQ,
    then NS (substitutions), then conversion tags (Zf/Yf).  Conversion checking is
    only applied when the BAM carries the bisulfite tags.
    """
    if read.is_unmapped:
        return "unmapped"
    if read.is_duplicate:
        return "duplicate"
    # NOTE: supplementary reads are kept (same as the C core / samtools default
    # exclflags: UNMAP|SECONDARY|QCFAIL|DUP).  Older countmut only dropped
    # SECONDARY too, so this keeps all backends byte-consistent.
    if read.is_secondary:
        return "secondary"
    if fcfg.exclude_flags and (read.flag & fcfg.exclude_flags):
        return "excluded_flags"
    if fcfg.include_flags and (read.flag & fcfg.include_flags) != fcfg.include_flags:
        return "excluded_flags"
    if read.mapping_quality < fcfg.min_mapq:
        return "mapq"
    # NS tag (number of substitutions) -- only if present
    if fcfg.max_sub is not None and read.has_tag("NS"):
        if read.get_tag("NS") > fcfg.max_sub:
            return "mismatch"
    if has_bisulfite_tags and (
        fcfg.max_unc is not None
        and read.has_tag("Zf")
        and fcfg.min_con is not None
        and read.has_tag("Yf")
    ):
        zf = read.get_tag("Zf")
        yf = read.get_tag("Yf")
        if not (zf <= fcfg.max_unc and yf >= fcfg.min_con):
            return "conversion"
    return None
