"""Regression tests for correctness bugs found during the unified-core audit.

Each test encodes a bug that was confirmed on real output and then fixed:

* ``test_minus_row_uses_forward_motif`` -- the '-' mutation row reverse-
  complemented the motif while counting *reference-forward* bases, i.e. the
  motif base and the counted u/m bases disagreed.  Both strands now show the
  reference-forward window.
* ``test_readwalk_indel_parity`` -- the read-walk engine ignored deletions /
  ref-skips / failures (and omitted those positions entirely), diverging from
  pileup / C.
* ``test_pileup_strand_process`` -- the pileup engine ignored the forward /
  reverse strand gate, so (unsplit) base output counted both strands.
* ``test_min_depth_filters`` -- ``--min-depth`` was a no-op in pure Python.
* ``test_allele_mode`` -- allele header/row shape + ``min_allele_support``.
* ``test_mutation_config_case`` -- lowercase ``--ref-base a`` etc. must work
  like the C core (uppercased).

The reference is intentionally *non-palindromic* (``ACGT``*15) so a reverse-
complemented motif is distinguishable from the forward one.
"""

import os

import pysam
import pytest

from countmut.model import EngineConfig, FilterConfig, MutationConfig, StrandConfig
from countmut.pipeline import run_pipeline

SEQ60 = "ACGT" * 15  # 0-based: 0=A 1=C 2=G 3=T 4=A 5=C ...


def _mk(header, qname, flag, start, seq, mapq=60, cigar=None):
    a = pysam.AlignedSegment(header)
    a.query_name = qname
    a.flag = flag
    a.reference_id = 0
    a.reference_start = start
    a.mapping_quality = mapq
    a.cigarstring = cigar or f"{len(seq)}M"
    a.query_sequence = seq
    a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
    return a


def _write_bam(root, reads):
    fa = os.path.join(root, "ref.fa")
    bam = os.path.join(root, "test.bam")
    with open(fa, "w") as f:
        f.write(">chr1\n" + SEQ60 + "\n")
    pysam.faidx(fa)
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 60}]}
    )
    reads = [r(hdr) for r in reads]
    reads.sort(key=lambda r: (r.reference_start, r.flag))
    with pysam.AlignmentFile(bam, "wb", header=hdr) as out:
        for r in reads:
            out.write(r)
    pysam.index(bam)
    return bam, fa


@pytest.fixture(scope="module")
def motif_data(tmp_path_factory):
    """plus (fwd) / minus (rev) / minus_mut (rev, C->T at 0-based 1) reads."""
    root = tmp_path_factory.mktemp("mfix")
    return _write_bam(
        str(root),
        [
            lambda h, q="plus": _mk(h, q, 0, 0, "ACGTACGTAC"),
            lambda h, q="minus": _mk(h, q, 16, 0, "ACGTACGTAC"),
            lambda h, q="minus_mut": _mk(h, q, 16, 0, "ATGTACGTAC"),
        ],
    )


@pytest.fixture(scope="module")
def indel_data(tmp_path_factory):
    """D = 5M2D5M (2bp deletion @[35,37)), S = 5M3N4M (3bp ref-skip)."""
    root = tmp_path_factory.mktemp("ind")
    return _write_bam(
        str(root),
        [
            lambda h, q="D": _mk(h, q, 0, 30, "A" * 10, cigar="5M2D5M"),
            lambda h, q="S": _mk(h, q, 0, 30, "A" * 9, cigar="5M3N4M"),
        ],
    )


@pytest.fixture(scope="module")
def fcfg():
    return FilterConfig(min_mapq=0, min_baseq=0, trim_start=0, trim_end=0)


def _run(bam, fa, fcfg, mode, engine, region, split=True, **eckw):
    mcfg = (
        MutationConfig(ref_base="C", mut_base="T", pad=1, save_rest=True)
        if mode == "mutation"
        else None
    )
    ecfg = EngineConfig(mode=mode, engine=engine, region=region, threads=1, **eckw)
    return run_pipeline(
        bam, fa, output=None, fcfg=fcfg, mcfg=mcfg,
        scfg=StrandConfig(process="both", split=split), ecfg=ecfg,
    )


# ---------------------------------------------------------------------------
# BUG: '-' mutation row motif reverse-complemented + reference-forward bases
# ---------------------------------------------------------------------------
def test_minus_row_uses_forward_motif(motif_data, fcfg):
    bam, fa = motif_data
    res = _run(bam, fa, fcfg, "mutation", "pileup", "chr1:1-8")
    by = {(r[1], r[2]): r for r in res.rows}
    # position 2 (1-based; 0-based 1 = 'C', motif ACG with pad 1)
    plus = by[(2, "+")]
    minus = by[(2, "-")]
    assert plus[3] == "ACG"
    # the '-' row must show the SAME forward motif (reference-forward bases)
    assert minus[3] == "ACG", f"- row motif {minus[3]!r} must be the forward ACG"
    # minus read: C (unconverted -> u2), minus_mut read: T (-> m2)
    assert (minus[6], minus[9]) == (1, 1), minus  # u2, m2
    assert (plus[6], plus[9]) == (1, 0), plus


def test_minus_motif_identical_both_engines(motif_data, fcfg):
    bam, fa = motif_data
    rows = {}
    for engine in ("pileup", "read-walk"):
        res = _run(bam, fa, fcfg, "mutation", engine, "chr1:1-8")
        rows[engine] = sorted(map(tuple, res.rows))
    assert rows["pileup"] == rows["read-walk"]
    for r in rows["pileup"]:
        assert r[2] in ("+", "-")
        assert r[3] == "ACG", f"Unexpected motif {r[3]!r} on strand {r[2]}"


# ---------------------------------------------------------------------------
# BUG: read-walk ignored deletions / ref-skips (and omitted those positions)
# ---------------------------------------------------------------------------
def test_readwalk_indel_parity(indel_data, fcfg):
    bam, fa = indel_data
    out = {}
    for engine in ("pileup", "read-walk"):
        res = _run(bam, fa, fcfg, "base", engine, "chr1:31-45", split=False,
                   count_indels=True)
        ix = {h: i for i, h in enumerate(res.header)}
        out[engine] = {
            r[ix["pos"]]: (r[ix["del"]], r[ix["ref_skip"]], r[ix["fail"]])
            for r in res.rows
        }
    assert out["pileup"] == out["read-walk"]
    # 1-based 36/37 = deletion [35,37) + ref-skip [35,38); 38 = ref-skip only
    assert out["pileup"][36] == (1, 1, 0)
    assert out["pileup"][37] == (1, 1, 0)
    assert out["pileup"][38] == (0, 1, 0)


# ---------------------------------------------------------------------------
# BUG: pileup engine ignored strand_process
# ---------------------------------------------------------------------------
def test_pileup_strand_process(motif_data, fcfg):
    bam, fa = motif_data
    depths = {}
    for engine in ("read-walk", "pileup"):
        ecfg = EngineConfig(mode="base", engine=engine, region="chr1:1-10", threads=1)
        res = run_pipeline(
            bam, fa, output=None, fcfg=fcfg,
            scfg=StrandConfig(process="forward", split=False), ecfg=ecfg,
        )
        depths[engine] = next(r[3] for r in res.rows if r[1] == 2)  # depth at pos2
    # only the forward 'plus' read counts -> depth exactly 1 (not 3)
    assert depths == {"read-walk": 1, "pileup": 1}


# ---------------------------------------------------------------------------
# BUG: --min-depth was a no-op
# ---------------------------------------------------------------------------
def test_min_depth_filters(motif_data, fcfg):
    bam, fa = motif_data
    for split in (False, True):
        res = _run(bam, fa, fcfg, "base", "pileup", "chr1:1-25", split=split,
                   min_depth=1000)
        assert res.rows == [], f"min_depth=1000 should drop all base rows (split={split})"
    res = _run(bam, fa, fcfg, "allele", "pileup", "chr1:1-25", min_depth=1000)
    assert res.rows == [], "min_depth=1000 should drop all allele rows"


# ---------------------------------------------------------------------------
# BUG: allele header/row shape + min_allele_support no-op
# ---------------------------------------------------------------------------
def test_allele_mode_shape_and_support(motif_data, fcfg):
    bam, fa = motif_data
    res = _run(bam, fa, fcfg, "allele", "pileup", "chr1:1-25", min_allele_support=100)
    assert res.header == ["chrom", "pos", "ref", "depth", "ref_count", "alt", "alt_count"]
    assert all(len(r) == 7 for r in res.rows)
    # with support=100 no alt allele should satisfy the bar
    assert all(r[5] == "." and r[6] == 0 for r in res.rows)


def test_mutation_config_case(motif_data, fcfg):
    bam, fa = motif_data
    mcfg = MutationConfig(ref_base="c", mut_base="t", pad=1)
    assert (mcfg.ref_base, mcfg.mut_base) == ("C", "T")
    res = run_pipeline(
        bam, fa, output=None, fcfg=fcfg, mcfg=mcfg,
        scfg=StrandConfig(process="both"), ecfg=EngineConfig(mode="mutation"),
    )
    assert len(res.rows) > 0
