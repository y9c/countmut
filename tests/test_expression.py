"""Tests for the samtools-style filter-expression engine (-e / -p)."""

import pysam
import pytest

from countmut.expression import compile_pile_pred, compile_read_pred
from countmut.engine_pileup import pileup_region
from countmut.model import FilterConfig


@pytest.fixture
def read():
    H = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 100}]})
    r = pysam.AlignedSegment(H)
    r.query_name = "q"
    r.flag = 0x400  # DUP
    r.reference_id = 0
    r.reference_start = 10
    r.mapping_quality = 30
    r.cigarstring = "50M"
    r.query_sequence = "A" * 50
    r.query_qualities = pysam.qualitystring_to_array("I" * 50)
    r.set_tag("RG", "sampleA")
    r.set_tag("NM", 3)
    return r


READ_EXPRS = [
    ("mapq >= 20", True),
    ("flag & 16 == 0", True),          # not reverse
    ("flag & DUP != 0", True),         # is duplicate
    ("[NM] >= 3", True),
    ("[NM] < 2", False),
    ("tag('RG') == 'sampleA'", True),
    ("exists([RG])", True),
    ("qname == 'q'", True),
    ("pos >= 10", True),
    ("qlen == 50", True),
    ("mapq >= 30 && qlen == 50", True),
    ("mapq >= 40 || qlen < 5", False),
    ("not (mapq >= 40)", True),
    ("avg(qual) > 20", True),
    ("length(seq) == 50", True),
    ("mapq >= 20 and [NM] >= 3", True),  # 'and' alias
]


@pytest.mark.parametrize("expr,expected", READ_EXPRS)
def test_read_expressions(read, expr, expected):
    assert compile_read_pred(expr)(read, 0) is expected


def test_backslash_and_or_not_aliases(read):
    assert compile_read_pred("mapq >= 20 or mapq < 5")(read, 0) is True
    assert compile_read_pred("not (flag & 16)")(read, 0) is True  # ! binds tighter, so use parens


def test_bad_expression_raises():
    with pytest.raises(ValueError):
        compile_read_pred("mapq >>= 5")


def test_pile_expressions(data_bam):
    bam, fa = data_bam
    cols = pileup_region(pysam.AlignmentFile(bam, "rb"), pysam.FastaFile(fa),
                         "chr1", 0, 80, FilterConfig(min_mapq=0, min_baseq=0,
                                                     trim_start=0, trim_end=0),
                         mode="base")
    assert any(compile_pile_pred("depth >= 1")(c) for c in cols)
    assert compile_pile_pred("pos >= 5 and ref == 'A'")(cols[0])


@pytest.fixture
def data_bam(tmp_path_factory):
    root = tmp_path_factory.mktemp("expr")
    fa = str(root / "ref.fa")
    bam = str(root / "t.bam")
    open(fa, "w").write(">chr1\n" + "A" * 80 + "\n")
    pysam.faidx(fa)
    H = pysam.AlignmentHeader.from_dict(
        {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 80}]})
    a = pysam.AlignedSegment(H)
    a.query_name = "f"; a.flag = 0; a.reference_id = 0; a.reference_start = 5
    a.mapping_quality = 60; a.cigarstring = "25M"; a.query_sequence = "A" * 25
    a.query_qualities = pysam.qualitystring_to_array("I" * 25)
    with pysam.AlignmentFile(bam, "wb", header=H) as out:
        out.write(a)
    pysam.index(bam)
    return bam, fa
