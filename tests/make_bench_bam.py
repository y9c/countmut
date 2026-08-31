#!/usr/bin/env python3
"""Build a reproducible bimodal benchmark BAM for countmut.

Reproduces the shape that matters for performance: a long, shallow body of
transcripts plus one or two ultra-deep "rRNA-hotspot" contigs, ~100 bp
paired-end reads (proper FR pairs, reverse reads stored reference-forward).

Usage:
    python tests/make_bench_bam.py --out /path/to/bench [options]
"""

import argparse
import os
import random

import pysam


def main():
    ap = argparse.ArgumentParser(description="build a bimodal benchmark BAM")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--shallow", type=int, default=3000, help="shallow transcripts")
    ap.add_argument("--deep", type=int, default=2, help="deep hotspot contigs")
    ap.add_argument(
        "--deep-depth", type=int, default=3000, help="per-position depth at hotspots"
    )
    ap.add_argument("--read-len", type=int, default=100)
    ap.add_argument("--insert", type=int, default=300)
    ap.add_argument(
        "--tags", action="store_true", help="add NS:i0 tags (for -e tag() examples)"
    )
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()
    rng = random.Random(args.seed)

    os.makedirs(args.out, exist_ok=True)
    fa = os.path.join(args.out, "bench.fa")
    bam = os.path.join(args.out, "bench.bam")

    # (name, sequence, target read1 depth)
    targets = []
    for i in range(args.shallow):
        name = f"PCG{i:05d}"
        seq = "".join(rng.choice("ACGT") for _ in range(rng.randint(500, 2000)))
        targets.append((name, seq, 1))
    for i in range(args.deep):
        name = "rRNA-18S" if i == 0 else f"rRNA-{18 + (i * 10) % 11}S"
        seq = "".join(rng.choice("ACGT") for _ in range(1500))
        targets.append((name, seq, args.deep_depth))

    with open(fa, "w") as f:
        for name, seq, _d in targets:
            f.write(f">{name}\n{seq}\n")
    pysam.faidx(fa)

    hdr = pysam.AlignmentHeader.from_dict(
        {
            "HD": {"VN": "1.6", "SO": "coordinate"},
            "SQ": [{"SN": t[0], "LN": len(t[1])} for t in targets],
        }
    )

    recs = []
    n_reads = 0
    for tid, (name, seq, depth) in enumerate(targets):
        # number of fragments so read1 depth ~= target:
        #   fragments/len[reads-per-bp] * read_len reads-per-position = depth
        nfrag = max(1, int(depth * len(seq) / args.read_len))
        for f in range(nfrag):
            s1 = rng.randrange(0, len(seq) - args.read_len + 1)
            s2 = s1 + args.insert - args.read_len
            for start, flag, is_read1 in ((s1, 99, True), (s2, 147, False)):
                if start < 0 or start + args.read_len > len(seq):
                    continue
                qseq = seq[start : start + args.read_len]
                # BAM stores SEQ reference-forward: for a (reverse-strand) read
                # the stored sequence is the reference span itself, with SEQ[0]
                # at the leftmost reference position (NOT the observed molecule;
                # flag bits 16 alone carry the strand).  Reverse-complementing
                # here would store the molecule form and make every reverse-read
                # base disagree with the reference.
                r = pysam.AlignedSegment(hdr)
                r.query_name = f"{name}_{f}_{1 if is_read1 else 2}"
                r.flag = flag
                r.reference_id = tid
                r.reference_start = start
                r.mapping_quality = 60
                r.query_sequence = qseq
                r.query_qualities = [40] * args.read_len
                r.cigar = [(0, args.read_len)]
                r.next_reference_id = tid
                r.next_reference_start = s2 if is_read1 else s1
                r.template_length = args.insert if is_read1 else -args.insert
                if args.tags:
                    r.set_tag("NS", 0, value_type="i")
                recs.append(r)
    n_reads = len(recs)

    recs.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(bam, "wb", header=hdr) as f:
        for r in recs:
            f.write(r)
    pysam.index(bam)

    nbases = sum(len(r.query_sequence) for r in recs)
    print(f"wrote {n_reads:,} reads / {nbases:,} read-bases  ->  {bam}")
    print(
        f"reference: {len(targets):,} contigs, {sum(len(t[1]) for t in targets):,} bp  ->  {fa}"
    )


if __name__ == "__main__":
    main()
