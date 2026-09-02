"""Regression tests for the -e group router (countmut >= 0.2.2).

The read expression is now a group router: ``nil``/``false`` drops the base,
``true`` routes it to group 0, and an integer ``0..3`` routes it to that
group.  Per-group counts surface in ``-o`` row templates as ``{a.0}`` …
``{n.3}`` (plain ``{a}`` stays the per-strand total over all groups).

The fixture mirrors the bisulfite A->G pipeline router (2 groups):

    ([NS] <= 1) and (([Yf] >= 1 and [Zf] <= 3 and bq >= 20
                      and qpos >= 2 and qlen - qpos > 2) and 1 or 0)

group 1 = high-conversion bases, group 0 = everything else that passes the
hard NS gate (low quality / read-end positions); NS > 1 drops the read.

Minus-strand convention (0.2.1): a minus read's stored SEQ is reference-
forward, so the reference-frame base is complement(stored[qpos]) with qpos in
CIGAR order (pos - POS0, no flip).  A minus read storing T reports A.
"""

import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c_runner import run_c

ROUTER = (
    "([NS] <= 1) and (([Yf] >= 1 and [Zf] <= 3 and bq >= 20 "
    "and qpos >= 2 and qlen - qpos > 2) and 1 or 0)"
)
TPL = "{chrom}\t{pos+1}\t{strand}\t{a.0}\t{a.1}\t{g.0}\t{g.1}"
TPL_MOTIF = "{chrom}\t{pos+1}\t{strand}\t{motif}\t{a.0}\t{a.1}\t{g.0}\t{g.1}"


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("cmr")
    fa = str(root / "ref.fa")
    bam = str(root / "test.bam")
    ref = "GGGGAAAACCCC"  # A at 0-based 4..7 (1-based 5..8)
    with open(fa, "w") as f:
        f.write(">chrT\n" + ref + "\n")
    pysam.faidx(fa)

    hdr = {
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": "chrT", "LN": len(ref)}],
    }
    header = pysam.AlignmentHeader.from_dict(hdr)

    def mk(qname, flag, start, seq, qual, tags):
        a = pysam.AlignedSegment(header)
        a.query_name = qname
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = 60
        a.cigarstring = f"{len(seq)}M"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array(qual)
        for tag in tags:
            a.set_tag(*tag)
        a.next_reference_id = 0
        a.next_reference_start = start
        return a

    Q40 = "I" * 9  # Phred 40
    Q10 = "+" * 9  # Phred 10
    ok = [("NS", 0), ("Yf", 1), ("Zf", 0)]
    reads = [
        # plus reads, 9M at 0-based 0..8 (1-based 1..9)
        # p_hi: A at 0-based 4,5,6 (qpos 4,5,6 pass the gate) -> group 1
        mk("p_hi", 0, 0, "GGGGAAATT", Q40, ok),
        # p_lo: same bases but bq < 20 -> group 0
        mk("p_lo", 0, 0, "GGGGAAATT", Q10, ok),
        # p_g: G at 0-based 4, A at 5,6 -> group 1
        mk("p_g", 0, 0, "GGGGGAATT", Q40, ok),
        # p_x: NS > 1 -> dropped entirely
        mk("p_x", 0, 0, "GGGGAAATT", Q40, [("NS", 2), ("Yf", 1), ("Zf", 0)]),
        # minus reads, 9M at 0-based 4..12; stored T -> reference-frame A.
        # qpos = pos - 4, so group 1 (qpos 2..6) = 0-based 6..10.
        mk("m_hi", 16, 4, "TTTTTTTTT", Q40, ok),
        # m_lo: no Yf -> group 0 at every base
        mk("m_lo", 16, 4, "TTTTTTTTT", Q40, [("NS", 0), ("Yf", 0), ("Zf", 0)]),
        # m_x: NS > 1 -> dropped entirely
        mk("m_x", 16, 4, "TTTTTTTTT", Q40, [("NS", 2), ("Yf", 1), ("Zf", 0)]),
    ]
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(bam)
    return bam, fa, ref


HDR = "chrom\tpos\tstrand\ta0\ta1\tg0\tg1"


def test_router_per_group_counts(data):
    bam, fa, ref = data
    for engine in ("pileup", "read-walk"):
        h, rows = run_c(
            bam,
            fa,
            engine=engine,
            extra=[
                "--read-expr",
                ROUTER,
                "--pile-expr",
                "ref == 'A' and (a + g) > 0",
                "--output-expr",
                TPL,
                "--fmt-header",
                HDR,
            ],
        )
        got = {(str(r[1]), str(r[2])): (int(r[3]), int(r[4]), int(r[5]), int(r[6]))
               for r in rows}
        exp = {
            ("5", "+"): (1, 1, 0, 1),   # p_lo A g0; p_hi A g1; p_g G g1
            ("5", "-"): (2, 0, 0, 0),   # m_hi qpos0->g0, m_lo g0 (drops excluded)
            ("6", "+"): (1, 2, 0, 0),   # p_lo A g0; p_hi+p_g A g1
            ("6", "-"): (2, 0, 0, 0),   # m_hi qpos1->g0, m_lo g0
            ("7", "+"): (1, 2, 0, 0),
            ("7", "-"): (1, 1, 0, 0),   # m_hi qpos2->g1, m_lo g0
            ("8", "+"): (0, 0, 0, 0),   # plus reads store T here (site kept by minus)
            ("8", "-"): (1, 1, 0, 0),   # m_hi qpos3->g1, m_lo g0
        }
        assert got == exp, f"[{engine}] group counts mismatch:\ngot {got}\nexp {exp}"


def test_engines_identical(data):
    bam, fa, ref = data
    extra = ["--read-expr", ROUTER, "--pile-expr", "ref == 'A'",
             "--output-expr", TPL, "--fmt-header", HDR]
    a = run_c(bam, fa, engine="pileup", extra=extra)
    b = run_c(bam, fa, engine="read-walk", extra=extra)
    assert a == b, "engines diverged under the group router"


def test_plain_totals_sum_groups(data):
    """{a}/{g} stay the per-strand total over all groups."""
    bam, fa, ref = data
    h, rows = run_c(
        bam,
        fa,
        engine="pileup",
        extra=[
            "--read-expr",
            ROUTER,
            "--pile-expr",
            "ref == 'A'",
            "--output-expr",
            "{pos+1}\t{strand}\t{a}\t{a.0}\t{a.1}\t{g}\t{g.0}\t{g.1}",
            "--fmt-header",
            "pos\tstrand\ta\ta0\ta1\tg\tg0\tg1",
        ],
    )
    for r in rows:
        pos, s, a, a0, a1, g, g0, g1 = r
        assert int(a) == int(a0) + int(a1), (pos, s, r)
        assert int(g) == int(g0) + int(g1), (pos, s, r)


def test_drop_semantics(data):
    """NS > 1 drops the whole read: p_x/m_x must contribute nothing."""
    bam, fa, ref = data
    # Without the NS gate the same reads would add to the counts.
    no_ns = ROUTER.replace("([NS] <= 1) and ", "")
    h1, rows1 = run_c(bam, fa, engine="pileup",
                      extra=["--read-expr", no_ns, "--pile-expr", "ref == 'A'",
                             "--output-expr", "{pos+1}\t{strand}\t{a.0 + a.1}",
                             "--fmt-header", "pos\tstrand\ta"])
    h2, rows2 = run_c(bam, fa, engine="pileup",
                      extra=["--read-expr", ROUTER, "--pile-expr", "ref == 'A'",
                             "--output-expr", "{pos+1}\t{strand}\t{a.0 + a.1}",
                             "--fmt-header", "pos\tstrand\ta"])
    d1 = {(str(r[0]), str(r[1])): int(r[2]) for r in rows1}
    d2 = {(str(r[0]), str(r[1])): int(r[2]) for r in rows2}
    # p_x adds 1 A at 0-based 4..6 (pos 5..7, +); m_x adds 1 A at 0-based
    # 4..12 (pos 5..8, -) via the no-NS gate.
    for pos, s, add in (("5", "+", 1), ("6", "+", 1), ("7", "+", 1),
                        ("8", "+", 0), ("5", "-", 1), ("6", "-", 1),
                        ("7", "-", 1), ("8", "-", 1)):
        assert d1[(pos, s)] == d2[(pos, s)] + add, (pos, s, d1, d2)


def test_boolean_backward_compat(data):
    """A boolean -e still behaves as a keep/drop filter (group 0 only)."""
    bam, fa, ref = data
    h, rows = run_c(
        bam, fa, engine="pileup",
        extra=["--read-expr", "bq >= 30",
               "--output-expr", "{pos+1}\t{strand}\t{a}\t{a.0}\t{a.1}",
               "--fmt-header", "pos\tstrand\ta\ta0\ta1"],
    )
    by = {(str(r[0]), str(r[1])): (int(r[2]), int(r[3]), int(r[4])) for r in rows}
    # bq >= 30 keeps the Q40 reads (p_hi, p_g, p_x, m_hi, m_lo, m_x) and
    # drops p_lo (Q10).  A bare boolean result routes every kept base to
    # group 0 -> a.1 == 0 everywhere and a == a.0.
    for (pos, s), (a, a0, a1) in by.items():
        assert a1 == 0, (pos, s)
        assert a == a0, (pos, s)
    # pos5+: p_hi A + p_x A (p_g is G) -> a=2 ; pos6+: p_hi, p_g, p_x all A.
    assert by[("5", "+")] == (2, 2, 0)
    assert by[("6", "+")] == (3, 3, 0)


def test_motif_window_and_revcomp(data):
    """--motif-pad 2: 5-mer window; minus rows get the reverse complement."""
    bam, fa, ref = data
    h, rows = run_c(
        bam, fa, engine="pileup",
        extra=["--read-expr", ROUTER, "--pile-expr", "ref == 'A'",
               "--motif-pad", "2", "--output-expr", TPL_MOTIF,
                       "--fmt-header", "chrom\tpos\tstrand\tmotif\ta0\ta1\tg0\tg1"],
    )
    exp_motif = {
        ("5", "+"): "GGAAA", ("5", "-"): "TTTCC",
        ("6", "+"): "GAAAA", ("6", "-"): "TTTTC",
        ("7", "+"): "AAAAC", ("7", "-"): "GTTTT",
        ("8", "+"): "AAACC", ("8", "-"): "GGTTT",
    }
    for r in rows:
        pos, s, motif = str(r[1]), str(r[2]), r[3]
        assert motif == exp_motif[(pos, s)], (pos, s, motif)
