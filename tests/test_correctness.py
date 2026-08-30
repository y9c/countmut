"""Regression tests for correctness bugs found during the unified-core audit.

The only implementation is the C core (both the read-walk and pileup engines
live in ``backend/countmut_core``); these tests drive the binary directly.

The reference is intentionally non-palindromic (``ACGT``*15) so a
reverse-complemented motif is distinguishable from the forward one.

Covered regressions:
* ``test_minus_row_uses_forward_motif`` -- the '-' mutation-row motif must be
  the reference-forward window (it used to be reverse-complemented while the
  bases were reference-forward).
* ``test_indel_parity`` -- both engines tally deletions / ref-skips /
  failures identically (base mode + ``--count-indels``).
* ``test_strand_gate`` -- ``--strand forward/reverse`` filters reads, not just
  output rows.
* ``test_min_depth`` -- ``--min-depth`` actually filters base/allele rows.
* ``test_allele_mode`` -- header/row shape + ``min_allele_support``.
* ``test_mutation_config_case`` -- lowercase ``--ref-base a`` works.
* ``test_expr_e`` / ``test_expr_p`` / ``test_expr_read_equals_pileup`` --
  `-e` / `-p` Lua filters work in C, in both engines.
"""

import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c_runner import run_c

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


# ---------------------------------------------------------------------------
# BUG: '-' mutation-row motif reverse-complemented with reference-forward bases
# ---------------------------------------------------------------------------
def test_minus_row_uses_forward_motif(motif_data):
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
        mode="mutation",
        engine="pileup",
        region="chr1:1-8",
        extra=["--ref-base", "C", "--mut-base", "T", "--pad", "1", "--save-rest"],
    )
    by = {(r[1], r[2]): r for r in rows}
    plus = by[(2, "+")]
    minus = by[(2, "-")]
    assert plus[3] == "ACG"
    assert minus[3] == "ACG", f"- row motif {minus[3]!r} must be the forward ACG"
    # minus read: C -> u2; minus_mut read: T -> m2 (reference-forward)
    assert (minus[6], minus[9]) == (1, 1), minus  # u2, m2
    assert (plus[6], plus[9]) == (1, 0), plus


def test_minus_motif_identical_both_engines(motif_data):
    bam, fa = motif_data
    rw = run_c(
        bam,
        fa,
        mode="mutation",
        engine="read-walk",
        region="chr1:1-8",
        extra=["--ref-base", "C", "--mut-base", "T", "--pad", "1", "--save-rest"],
    )
    pl = run_c(
        bam,
        fa,
        mode="mutation",
        engine="pileup",
        region="chr1:1-8",
        extra=["--ref-base", "C", "--mut-base", "T", "--pad", "1", "--save-rest"],
    )
    assert rw == pl
    for r in pl[1]:
        assert r[3] == "ACG", f"Unexpected motif {r[3]!r} on strand {r[2]}"


# ---------------------------------------------------------------------------
# BUG: read-walk ignored deletions / ref-skips (diverged from pileup)
# ---------------------------------------------------------------------------
def test_indel_parity(indel_data):
    bam, fa = indel_data
    results = {}
    for engine in ("pileup", "read-walk"):
        _h, rows = run_c(
            bam,
            fa,
            mode="base",
            engine=engine,
            region="chr1:31-45",
            extra=["--split-strand", "--count-indels"],
        )
        # header: chrom pos strand ref depth a c g t n ins del ref_skip fail
        per = {}
        for r in rows:
            per.setdefault(r[1], [0, 0, 0, 0])
            for i in range(4):  # ins, del, ref_skip, fail at idx 10..13
                per[r[1]][i] += int(r[10 + i])
        results[engine] = {p: tuple(v) for p, v in per.items()}
    assert results["pileup"] == results["read-walk"]
    # 36/37 = deletion + ref-skip; 38 = ref-skip only
    assert results["pileup"][36] == (0, 1, 1, 0)
    assert results["pileup"][37] == (0, 1, 1, 0)
    assert results["pileup"][38] == (0, 0, 1, 0)


# ---------------------------------------------------------------------------
# BUG: strand gate must filter reads (not just rows)
# ---------------------------------------------------------------------------
def test_strand_gate_depths(motif_data):
    bam, fa = motif_data
    fwd = run_c(
        bam,
        fa,
        mode="base",
        engine="pileup",
        region="chr1:1-10",
        extra=["--strand", "forward"],
    )
    # only the forward 'plus' read -> depth 1 at pos 2 (1-based)
    pos2 = [r for r in fwd[1] if r[1] == 2]
    assert pos2 and pos2[0][4] == 1, pos2


# ---------------------------------------------------------------------------
# BUG: --min-depth was a no-op
# ---------------------------------------------------------------------------
def test_min_depth(motif_data):
    bam, fa = motif_data
    for mode, xtra in (("base", ["--split-strand"]), ("allele", [])):
        _h, rows = run_c(
            bam,
            fa,
            mode=mode,
            engine="pileup",
            region="chr1:1-25",
            extra=xtra + ["--min-depth", "1000"],
        )
        assert rows == [], f"min_depth=1000 should drop all {mode} rows"


# ---------------------------------------------------------------------------
# BUG: allele header/row shape + min_allele_support no-op
# ---------------------------------------------------------------------------
def test_allele_mode_shape_and_support(motif_data):
    bam, fa = motif_data
    header, rows = run_c(
        bam,
        fa,
        mode="allele",
        engine="pileup",
        region="chr1:1-25",
        extra=["--min-allele-support", "100"],
    )
    assert header == ["chrom", "pos", "ref", "depth", "ref_count", "alt", "alt_count"]
    assert all(len(r) == 7 for r in rows)
    assert all(r[5] == "." and r[6] == 0 for r in rows)


# ---------------------------------------------------------------------------
# case-insensitive --ref-base / --mut-base
# ---------------------------------------------------------------------------
def test_mutation_config_case(motif_data):
    bam, fa = motif_data
    lo = run_c(
        bam,
        fa,
        mode="mutation",
        engine="pileup",
        region="chr1:1-8",
        extra=["--ref-base", "c", "--mut-base", "t", "--pad", "1"],
    )
    up = run_c(
        bam,
        fa,
        mode="mutation",
        engine="pileup",
        region="chr1:1-8",
        extra=["--ref-base", "C", "--mut-base", "T", "--pad", "1"],
    )
    assert lo == up
    assert len(lo[1]) > 0


# ---------------------------------------------------------------------------
# -e / -p Lua filters in C (both engines)
# ---------------------------------------------------------------------------
def test_expr_e_read_filter(motif_data):
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
        mode="base",
        engine="pileup",
        region="chr1:1-10",
        extra=["--read-expr", "flags & 16 == 0"],
    )  # keep forward only
    for r in rows:
        assert r[2] == "+", f"reverse read leaked with -e: {r}"


def test_expr_p_pile_filter(motif_data):
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
        mode="base",
        engine="pileup",
        region="chr1:1-10",
        extra=["--pile-expr", "ref == 'A'"],
    )  # only A-reference sites
    assert rows, "expected some A-reference sites"
    assert all(r[3] == "A" for r in rows), "a non-A site leaked through -p"


def test_expr_identical_both_engines(motif_data):
    bam, fa = motif_data
    for expr in ("mapq >= 30", "qname ~= 'minus' and flags & 16 != 0", "bq >= 20"):
        rw = run_c(
            bam,
            fa,
            mode="mutation",
            engine="read-walk",
            region="chr1:1-8",
            extra=[
                "--ref-base",
                "C",
                "--mut-base",
                "T",
                "--pad",
                "1",
                "--save-rest",
                "--read-expr",
                expr,
            ],
        )
        pl = run_c(
            bam,
            fa,
            mode="mutation",
            engine="pileup",
            region="chr1:1-8",
            extra=[
                "--ref-base",
                "C",
                "--mut-base",
                "T",
                "--pad",
                "1",
                "--save-rest",
                "--read-expr",
                expr,
            ],
        )
        assert rw == pl, f"engines diverged with -e {expr!r}"


def test_expr_invalid_syntax(motif_data):
    bam, fa = motif_data
    # C validates the Lua up-front and exits 2 on a syntax error.
    with pytest.raises(AssertionError):
        run_c(
            bam,
            fa,
            mode="mutation",
            engine="pileup",
            extra=["--read-expr", "(( not lua !!"],
        )


def test_expr_file_like_usage(motif_data, tmp_path):
    """-p with a ref base + depth filter (countmut doc example style)."""
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
        mode="mutation",
        engine="pileup",
        region="chr1:1-8",
        extra=[
            "--ref-base",
            "C",
            "--mut-base",
            "T",
            "--pad",
            "1",
            "--pile-expr",
            "ref == 'C' and g >= 0",
        ],
    )
    assert rows
    for r in rows:
        assert r[1] in (2, 6)  # only C-reference sites
