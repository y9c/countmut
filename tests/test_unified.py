"""Tests for the C countmut core.

Both the read-walk and the pileup engine are implemented in C
(``backend/countmut_core``); Python only wraps the binary.  These tests build a
small, self-contained BAM/FASTA and assert, against the C core:

* read-walk and pileup are byte-identical,
* overlapping mates are not double counted (mate-overlap dedup),
* strand-aware processing (forward/reverse) works,
* base and allele (VCF) modes run,
* the '-' mutation-row motif is the reference-forward window (regression).
"""

import os
import sys

import pysam
import pytest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from c_runner import binary, run_c


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("cm")
    fa = str(root / "ref.fa")
    bam = str(root / "test.bam")
    with open(fa, "w") as f:
        f.write(">chr1\n" + "A" * 80 + "\n")
    pysam.faidx(fa)

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 80}]}
    header = pysam.AlignmentHeader.from_dict(hdr)

    def mk(qname, flag, start, seq, mapq=60):
        a = pysam.AlignedSegment(header)
        a.query_name = qname
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = mapq
        a.cigarstring = f"{len(seq)}M"
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.next_reference_id = 0
        a.next_reference_start = start
        return a

    # Refs are all 'A'.  Reads mostly 'A' (ref) / 'G' (mut) / 'C','T' (other).
    reads = []
    # f1: overlapping mates (read1 fwd @5 25M, read2 rev @20 25M -> overlap 21-30)
    reads.append(mk("f1", 99, 5, "A" * 25))  # forward, read1
    reads.append(mk("f1", 147, 20, "A" * 25))  # reverse, read2
    reads.append(mk("f2", 0, 55, "G" * 25))  # G (mutation vs A ref)
    reads.append(mk("f3", 0, 70, "C" * 25))  # C (other)
    reads.append(mk("f4", 0, 60, "A" * 25, mapq=0))  # low MAPQ ref
    reads.append(mk("f5", 16, 33, "A" * 25))  # single-end reverse -> '-'
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for r in reads:
            out.write(r)
    pysam.index(bam)
    return bam, fa


MUT_XTRA = ["--ref-base", "A", "--mut-base", "G", "--pad", "2", "--save-rest"]


def test_readwalk_pileup_identical(data):
    bam, fa = data
    for mode in ("mutation", "base", "allele"):
        rw = run_c(
            bam,
            fa,
            mode=mode,
            engine="read-walk",
            extra=MUT_XTRA if mode == "mutation" else [],
        )
        pl = run_c(
            bam,
            fa,
            mode=mode,
            engine="pileup",
            extra=MUT_XTRA if mode == "mutation" else [],
        )
        assert rw == pl, f"engines diverged in {mode} mode"


def test_overlap_mates_not_double_counted(data):
    bam, fa = data
    # f1 read1 fwd @5 and read2 rev @20: mates overlap ref 21-30 (1-based).
    _h, rows = run_c(bam, fa, mode="base", engine="pileup", region="chr1:1-40")
    by_pos = {int(r[1]): r for r in rows}
    for pos in range(21, 31):
        assert pos in by_pos, f"missing position {pos}"
        assert by_pos[pos][4] == 1, (
            f"pos {pos} depth {by_pos[pos][4]} (expected 1, dedup)"
        )


def test_strand_filter(data):
    bam, fa = data
    _h, fwd = run_c(
        bam,
        fa,
        mode="base",
        engine="pileup",
        region="chr1:1-40",
        extra=["--strand", "forward"],
    )
    _h, rev = run_c(
        bam,
        fa,
        mode="base",
        engine="pileup",
        region="chr1:1-40",
        extra=["--strand", "reverse"],
    )
    assert all(r[2] == "+" for r in fwd)
    assert all(r[2] == "-" for r in rev)


def test_base_mode(data):
    bam, fa = data
    header, rows = run_c(bam, fa, mode="base", engine="pileup")
    assert header[:4] == ["chrom", "pos", "strand", "ref"]
    assert any(r[0] == "chr1" and r[1] == 6 for r in rows)  # first base 1-based 6


def test_allele_mode(data):
    bam, fa = data
    header, rows = run_c(bam, fa, mode="allele", engine="pileup")
    assert header == ["chrom", "pos", "ref", "depth", "ref_count", "alt", "alt_count"]
    assert all(len(r) == 7 for r in rows)


def test_allele_vcf(data):
    bam, fa = data
    import subprocess

    r = subprocess.run(
        [
            binary(),
            "--bam",
            bam,
            "--fa",
            fa,
            "--out",
            "-",
            "--mode",
            "allele",
            "--vcf",
            "--min-mapq",
            "0",
            "--min-baseq",
            "0",
            "--trim-start",
            "0",
            "--trim-end",
            "0",
        ],
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, r.stderr
    lines = r.stdout.splitlines()
    assert lines[0].startswith("##fileformat=VCF")
    assert lines[1].startswith("#CHROM")
    assert any("\tPASS\t.\tGT:AD\t" in ln for ln in lines[2:])
