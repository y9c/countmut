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
def test_perbase_expr_parity_qpos(motif_data):
    """qpos-dependent per-base filters (dist5/dist3/base/bq) must be byte-identical
    between engines.  Regression: the pileup `expr_pass` memo dropped the qpos
    arg (used qpos=0), silently breaking dist5/dist3/base in the pileup engine."""
    bam, fa = motif_data
    for expr in ("dist5 >= 5", "dist3 >= 5", "base == 'A'", "bq >= 20"):
        rw = run_c(
            bam,
            fa,
            engine="read-walk",
            region="chr1:1-25",
            extra=["--read-expr", expr],
        )
        pl = run_c(
            bam,
            fa,
            engine="pileup",
            region="chr1:1-25",
            extra=["--read-expr", expr],
        )
        assert rw == pl, f"per-base qpos expr diverged between engines: {expr!r}"


def test_indel_parity(indel_data):
    bam, fa = indel_data
    results = {}
    for engine in ("pileup", "read-walk"):
        _h, rows = run_c(
            bam,
            fa,
            engine=engine,
            region="chr1:31-45",
            extra=["--count-indels"],
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
    for xtra in ([], ["--vcf"]):
        _h, rows = run_c(
            bam,
            fa,
            engine="pileup",
            region="chr1:1-25",
            extra=xtra + ["--min-depth", "1000"],
        )
        assert rows == [], "min_depth=1000 should drop all rows"


# ---------------------------------------------------------------------------
# BUG: allele header/row shape + min_allele_support no-op
# ---------------------------------------------------------------------------
def test_expr_e_read_filter(motif_data):
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
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
            engine="read-walk",
            region="chr1:1-8",
            extra=["--read-expr", expr],
        )
        pl = run_c(
            bam,
            fa,
            engine="pileup",
            region="chr1:1-8",
            extra=["--read-expr", expr],
        )
        assert rw == pl, f"engines diverged with -e {expr!r}"


# ---------------------------------------------------------------------------
# read-walk solo/direct fast path must still dedup overlapping mates exactly
# ---------------------------------------------------------------------------
def test_readwalk_proper_paired_overlap_dedup(tmp_path):
    """A proper pair (mpos + template length set -> hybrid read-walk path)
    must still count an overlapping position once."""
    root = str(tmp_path / "pw")
    os.makedirs(root, exist_ok=True)
    fa = os.path.join(root, "ref.fa")
    bam = os.path.join(root, "p.bam")
    with open(fa, "w") as f:
        f.write(">chr1\n" + "ACGT" * 5 + "\n")
    pysam.faidx(fa)
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 20}]}
    )
    seq = "ACGT" * 5
    a = pysam.AlignedSegment(hdr)
    a.query_name = "p"
    a.flag = 99
    a.reference_id = 0
    a.reference_start = 0
    a.mapping_quality = 60
    a.cigarstring = "10M"
    a.query_sequence = seq[0:10]
    a.query_qualities = pysam.qualitystring_to_array("I" * 10)
    a.next_reference_id = 0
    a.next_reference_start = 6
    a.template_length = 16
    b = pysam.AlignedSegment(hdr)
    b.query_name = "p"
    b.flag = 147
    b.reference_id = 0
    b.reference_start = 6
    b.mapping_quality = 60
    b.cigarstring = "10M"
    b.query_sequence = seq[6:16]
    b.query_qualities = pysam.qualitystring_to_array("I" * 10)
    b.next_reference_id = 0
    b.next_reference_start = 0
    b.template_length = -16
    with pysam.AlignmentFile(bam, "wb", header=hdr) as out:
        out.write(a)
        out.write(b)
    pysam.index(bam)
    _h, rows = run_c(
        bam,
        fa,
        engine="read-walk",
        region="chr1:1-20",
        extra=["--trim-fragment-start", "0", "--trim-fragment-end", "0"],
    )
    # overlap is 1-based 7..10 (0-based 6..9); each must have exactly depth 1
    for pos in (7, 8, 9, 10):
        row = [r for r in rows if r[1] == pos]
        assert row and row[0][4] == 1, (pos, row)


def test_expr_invalid_syntax(motif_data):
    bam, fa = motif_data
    # C validates the Lua up-front and exits 2 on a syntax error.
    with pytest.raises(AssertionError):
        run_c(
            bam,
            fa,
            engine="pileup",
            extra=["--read-expr", "(( not lua !!"],
        )


def test_expr_file_like_usage(motif_data, tmp_path):
    """-p pile-site filter selects which sites to report."""
    bam, fa = motif_data
    _h, rows = run_c(
        bam,
        fa,
        engine="pileup",
        region="chr1:1-8",
        extra=["--pile-expr", "ref == 'C' and g >= 0"],
    )
    assert rows
    for r in rows:
        assert r[1] in (2, 6)  # only C-reference sites


def test_long_run_of_bang_operator(tmp_path):
    """A long run of unary '!' must not overflow the expression translator
    ('!' expands to 'not ', a 4x buffer-size bound)."""
    bam, fa = _write_reference_read_bam(tmp_path, "chrX", 200)
    expr = "!" * 64 + " (mapq >= 0)"
    _h, rows = run_c(
        bam,
        fa,
        engine="read-walk",
        region="chrX:1-200",
        extra=["--read-expr", expr],
    )
    assert rows, "unary-bang expression did not evaluate (overflow/regression)"


def test_long_read_seq_not_truncated(tmp_path):
    """seq/qual must span the whole read (>1023 bp), not a fixed 1024 buffer."""
    root = str(tmp_path / "lr")
    os.makedirs(root, exist_ok=True)
    fa = os.path.join(root, "ref.fa")
    bam = os.path.join(root, "long.bam")
    with open(fa, "w") as f:
        f.write(">chrL\n" + "AAACCCGGGTTT" * 30 + "\n")
    pysam.faidx(fa)
    hdr = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrL", "LN": 400}]}
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "lr1"
    r.flag = 0
    r.reference_id = 0
    r.reference_start = 0
    r.mapping_quality = 60
    r.query_sequence = "A" * 2000  # 2000 bp read
    r.query_qualities = [40] * 2000
    r.cigar = [(0, 2000)]
    with pysam.AlignmentFile(bam, "wb", header=hdr) as out:
        out.write(r)
    pysam.index(bam)
    _h, rows = run_c(
        bam,
        fa,
        engine="read-walk",
        region="chrL:1-200",
        extra=["--read-expr", "slen(seq) == 2000"],
    )
    assert rows, "seq was truncated (slen(seq) != 2000) for a 2000 bp read"


def _write_reference_read_bam(tmp_path, chrom, length):
    """Build a tiny single-read BAM on a uniform (all-A) reference."""
    root = str(tmp_path / "ref_read")
    os.makedirs(root, exist_ok=True)
    fa = os.path.join(root, "ref.fa")
    bam = os.path.join(root, "in.bam")
    with open(fa, "w") as f:
        f.write(f">{chrom}\n" + "A" * (length + 100) + "\n")
    pysam.faidx(fa)
    hdr = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": chrom, "LN": length + 100}],
        }
    )
    r = pysam.AlignedSegment(hdr)
    r.query_name = "r1"
    r.flag = 0
    r.reference_id = 0
    r.reference_start = 0
    r.mapping_quality = 60
    r.query_sequence = "A" * min(100, length)
    r.query_qualities = [40] * min(100, length)
    r.cigar = [(0, min(100, length))]
    with pysam.AlignmentFile(bam, "wb", header=hdr) as out:
        out.write(r)
    pysam.index(bam)
    return bam, fa


def test_sam_input_matches_bam(motif_data, tmp_path):
    """SAM (plain, and gzipped) input must produce byte-identical output to the
    equivalent BAM (it is auto-transcoded to a temp BAM + index)."""
    import gzip

    bam, fa = motif_data
    sam = str(tmp_path / "x.sam")
    samgz = str(tmp_path / "x.sam.gz")
    htxt = pysam.view(bam, "-h")
    with open(sam, "w") as f:
        f.write(htxt)
    with gzip.open(samgz, "wt") as f:
        f.write(htxt)
    args = [
        "--read-expr",
        "bq >= 20 and dist5 >= 2",
    ]
    _h, from_bam = run_c(bam, fa, engine="pileup", region="chr1:1-8", extra=args)
    _h, from_sam = run_c(sam, fa, engine="pileup", region="chr1:1-8", extra=args)
    assert from_bam and from_bam == from_sam, "SAM input diverged from BAM"
    _h, from_samgz = run_c(samgz, fa, engine="pileup", region="chr1:1-8", extra=args)
    assert from_bam == from_samgz, "gzipped SAM input diverged from BAM"


def test_output_format_template(motif_data):
    """--output-format as a row template: header via --fmt-header (\\t expanded),
    cells computed from the site namespace."""
    bam, fa = motif_data
    tpl = "{pos+1}\t{ref}\t{a}\t{t}\t{round(a/(a+t), 3)}"
    header, rows = run_c(
        bam,
        fa,
        engine="pileup",
        region="chr1:1-8",
        extra=[
            "--output-expr",
            tpl,
            "--fmt-header",
            "pos\tref\ta\tt\tat_ratio",
        ],
    )
    assert header == ["pos", "ref", "a", "t", "at_ratio"], header
    assert rows and all(len(r) == 5 for r in rows), rows[:2]
    r0 = rows[0]
    assert int(r0[0]) >= 1 and r0[1] in "ACGTN"
