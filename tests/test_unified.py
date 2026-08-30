"""Unified countmut engine + backend tests.

These build a small, self-contained BAM/FASTA with overlapping mates and reverse
reads, then assert:

* read-walk and pileup produce identical mutation counts,
* the C backend (if built) matches the Python engines,
* overlapping mates are not double counted (dedup),
* strand-aware filtering works,
* base and allele (VCF) modes run.
"""

import os
import subprocess
import sys

import pysam
import pytest

from countmut.backend import ensure_backend
from countmut.model import EngineConfig, FilterConfig, MutationConfig, StrandConfig
from countmut.pipeline import run_pipeline


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    root = tmp_path_factory.mktemp("cm")
    fa = str(root / "ref.fa")
    bam = str(root / "test.bam")
    with open(fa, "w") as f:
        f.write(">chr1\n" + "A" * 80 + "\n")
    pysam.faidx(fa)

    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 80}]}
    H = pysam.AlignmentHeader.from_dict(hdr)

    def mk(qname, flag, start, seq, mapq=60):
        a = pysam.AlignedSegment(H)
        a.query_name = qname
        a.flag = flag
        a.reference_id = 0
        a.reference_start = start
        a.mapping_quality = mapq
        a.cigarstring = "%dM" % len(seq)
        a.query_sequence = seq
        a.query_qualities = pysam.qualitystring_to_array("I" * len(seq))
        a.next_reference_id = 0
        a.next_reference_start = start
        return a

    # Refs are all 'A'.  Reads mostly 'A' (ref) / 'G' (mut) / 'C','T' (other).
    reads = []
    # f1: overlapping mates (read1 fwd @5 25M, read2 rev @20 25M -> overlap 21-30)
    reads.append(mk("f1", 99, 5, "A" * 25))          # forward, read1
    reads.append(mk("f1", 147, 20, "A" * 25))        # reverse, read2
    # separate non-overlapping fragments (kept clear of f1's overlap region)
    reads.append(mk("f2", 0, 55, "G" * 25))          # G (mutation vs A ref)
    reads.append(mk("f3", 0, 70, "C" * 25))          # C (other)
    reads.append(mk("f4", 0, 60, "A" * 25, mapq=0))  # low MAPQ ref
    # a genuinely reverse (minus-strand) fragment for strand filtering
    reads.append(mk("f5", 16, 33, "A" * 25))         # single-end reverse -> '-'
    reads.sort(key=lambda r: r.reference_start)
    with pysam.AlignmentFile(bam, "wb", header=H) as out:
        for r in reads:
            out.write(r)
    pysam.index(bam)
    return bam, fa


@pytest.fixture(scope="module")
def fcfg():
    return FilterConfig(min_mapq=0, min_baseq=0, trim_start=0, trim_end=0)


TEST_MCFG = dict(ref_base="A", mut_base="G", pad=2, save_rest=True)


def _mutation_rows(bam, fa, engine, fcfg):
    res = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, mcfg=MutationConfig(**TEST_MCFG),
        scfg=StrandConfig(process="both", split=True),
        ecfg=EngineConfig(mode="mutation", engine=engine, bin_size=1000, threads=2),
    )
    return sorted(map(tuple, res.rows))


def test_readwalk_pileup_identical(data, fcfg):
    bam, fa = data
    rw = _mutation_rows(bam, fa, "read-walk", fcfg)
    pl = _mutation_rows(bam, fa, "pileup", fcfg)
    assert rw == pl
    assert len(rw) > 0


def test_overlap_mates_not_double_counted(data, fcfg):
    bam, fa = data
    # f1 read1 fwd @5 and read2 rev @20: mates overlap ref 21-30 (0-based).
    # Combined-strand depth in that overlap must be exactly 1 (not 2).
    res = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, scfg=StrandConfig(process="both", split=False),
        ecfg=EngineConfig(mode="base", engine="pileup", bin_size=1000, threads=2,
                          region="chr1:1-40"),
    )
    by_pos = {r[1]: r for r in res.rows if r[0] == "chr1"}
    # positions 21-30 (1-based) are inside the fragment-1 mate overlap
    for pos in range(21, 31):
        assert pos in by_pos, f"missing position {pos}"
        assert by_pos[pos][3] == 1, f"pos {pos} depth {by_pos[pos][3]} (expected 1, dedup)"


def test_strand_filter(data, fcfg):
    bam, fa = data
    fwd = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, scfg=StrandConfig(process="forward", split=True),
        ecfg=EngineConfig(mode="base", engine="pileup", bin_size=1000,
                          threads=2, region="chr1:1-40"),
    )
    assert all(r[2] == "+" for r in fwd.rows)
    rev = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, scfg=StrandConfig(process="reverse", split=True),
        ecfg=EngineConfig(mode="base", engine="pileup", bin_size=1000,
                          threads=2, region="chr1:1-40"),
    )
    assert all(r[2] == "-" for r in rev.rows)


def test_base_mode(data, fcfg):
    bam, fa = data
    res = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, scfg=StrandConfig(process="both", split=False),
        ecfg=EngineConfig(mode="base", engine="pileup", bin_size=1000, threads=2),
    )
    assert res.header[:3] == ["chrom", "pos", "ref"]
    # reads start at ref 5 (f1 read1) -> 1-based position 6 is the first base
    assert any(r[0] == "chr1" and r[1] == 6 for r in res.rows)


def test_allele_vcf_runs(data, fcfg):
    bam, fa = data
    res = run_pipeline(
        samfile=bam, reference=fa, output=None,
        fcfg=fcfg, ecfg=EngineConfig(mode="allele", engine="pileup",
                                     bin_size=1000, threads=2),
    )
    # fixed-width rows: chrom,pos,ref,depth,ref_count,alt,alt_count
    assert all(len(r) == 7 for r in res.rows)
    assert any(r[0] == "chr1" for r in res.rows)


@pytest.mark.skipif(ensure_backend() is None, reason="C backend not built")
def test_c_backend_matches_python(data, fcfg):
    bam, fa = data
    py = _mutation_rows(bam, fa, "read-walk", fcfg)
    binary = ensure_backend()
    out = os.path.join(os.path.dirname(bam), "c_mut.tsv")
    cmd = [str(binary), "--bam", bam, "--fa", fa, "--out", out,
           "--mode", "mutation", "--ref-base", "A", "--mut-base", "G",
           "--pad", "2", "--save-rest", "--min-mapq", "0", "--min-baseq", "0",
           "--trim-start", "0", "--trim-end", "0", "--threads", "2"]
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    with open(out) as fh:
        header = fh.readline().split("\t")
        assert header[:4] == ["chrom", "pos", "strand", "motif"]
        c = [l.rstrip("\n").split("\t") for l in fh]
    # Convert to ints FIRST, then sort numerically (string sort would put "10" < "6")
    c_norm = sorted(tuple([r[0], int(r[1]), r[2], r[3], *map(int, r[4:])]) for r in c)
    py_norm = sorted(tuple([r[0], int(r[1]), r[2], r[3], *map(int, r[4:])]) for r in py)
    assert c_norm == py_norm
