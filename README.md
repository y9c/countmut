# CountMut

[![Pypi Releases](https://img.shields.io/pypi/v/countmut.svg)](https://pypi.python.org/pypi/countmut)
[![Downloads](https://img.shields.io/pepy/dt/countmut)](https://pepy.tech/project/countmut)
[![Development Status](https://img.shields.io/badge/status-alpha-orange.svg)](https://github.com/y9c/countmut)

> **Ultra-fast strand-aware mutation counter for bisulfite sequencing analysis**

CountMut is a high-performance tool for counting mutations from bisulfite sequencing BAM files (BS-seq, CAM-seq, GLORI-seq, eTAM-seq). It features parallel processing, quality-based mate overlap deduplication, and optimized file I/O for maximum speed.

## Features

- 🚀 **Ultra-Fast**: Call mutation without pileup reads
- 🧬 **Bisulfite Support**: NS, Zf, Yf tag filtering for conversion analysis
- 🎯 **Accurate**: Quality-based mate overlap deduplication prevents double-counting
- ⚡ **Parallel**: Multi-threaded genomic window processing
- 🔧 **Flexible**: Configurable filtering, strand-specific processing, auto-indexing

## Installation

```bash
pip install countmut
```

## Quick Start

```bash
# Basic usage - auto-creates indices if needed
countmut -i input.bam -r reference.fa -o mutations.tsv

# Count T→C mutations (common in bisulfite sequencing)
countmut -i input.bam -r reference.fa -o mutations.tsv --ref-base T --mut-base C

# With custom threads and filtering
countmut -i input.bam -r reference.fa -o mutations.tsv -t 8 --max-unc 5 --min-con 2
```

## Options

**Input/Output**
```bash
-i, --input PATH       Input BAM file (coordinate-sorted) [required]
-r, --reference PATH   Reference FASTA file [required]
-o, --output PATH      Output TSV file (default: stdout)
-f, --force            Overwrite output without prompting
```

**Mutation Analysis**
```bash
--ref-base TEXT        Reference base to count from [default: A]
--mut-base TEXT        Mutation base to count [default: G]
--strand TEXT          Strand: both/forward/reverse [default: both]
--region TEXT          Genomic region (e.g., 'chr1:1000000-2000000')
```

**Performance**
```bash
-t, --threads INTEGER  Number of parallel threads [default: auto]
-b, --bin-size INTEGER Genomic bin size in bp [default: 10000]
```

**Alternative Mutation Tagging**
```bash
--ref-base2 TEXT       Alternative reference base for tagging (e.g., 'C')
--mut-base2 TEXT       Alternative mutation base for tagging (e.g., 'T')
--output-bam PATH      Output BAM with alternative tags (Yc, Zc)
```

**Quality Filters**
```bash
--min-baseq INTEGER    Min base quality (Phred score) [default: 20]
--min-mapq INTEGER     Min mapping quality (MAPQ) [default: 0]
--max-sub INTEGER      Max substitutions (NS tag) [default: 1]
--trim-start INTEGER   Trim N bases from read 5' end (fragment orientation) [default: 2]
--trim-end INTEGER     Trim N bases from read 3' end (fragment orientation) [default: 2]
--max-unc INTEGER      Max unconverted (Zf tag) [default: 3]
--min-con INTEGER      Min converted (Yf tag) [default: 1]
```

**Output Records**
```bash
-p, --pad INTEGER      Motif window half-size [default: 15]
-s, --save-rest        Include other bases (o0, o1, o2 columns)
```

> **Note**: BAM files must have **NS**, **Zf**, and **Yf** tags (essential for bisulfite analysis).
> Indices (.bai, .fai) are created automatically if missing.

## Output Format

TSV file with the following columns:

| Column | Description |
|--------|-------------|
| `chrom` | Chromosome name |
| `pos` | Genomic position (1-based) |
| `strand` | Strand (+ or -) |
| `motif` | Sequence context (2×pad+1 bp window) |
| `u0`, `u1`, `u2` | **Unconverted** (reference base) counts |
| `m0`, `m1`, `m2` | **Mutation** (mutation base only) counts |
| `o0`, `o1`, `o2` | **Other bases** counts (with `--save-rest`) |

**Count categories** (x0, x1, x2):
- **x0 (low quality)**: Bases failing quality filters (trim region, max-sub, min-mapq, min-baseq)
- **x1 (insufficient conversion)**: Bases from reads with insufficient conversion efficiency (high Zf or low Yf)
- **x2 (high conversion)**: Bases from reads with high conversion efficiency (low Zf and high Yf)

### Example Output

Without `--save-rest`:

| `chrom` | `pos` | `strand` | `motif` | `u0` | `u1` | `u2` | `m0` | `m1` | `m2` |
|---------|--------|----------|---------|------|------|------|------|------|------|
| `chr1`  | `10000` | `+`      | `AAG`   | `10` | `5`  | `2`  | `0`  | `0`  | `0`  |
| `chr1`  | `10001` | `-`      | `TTC`   | `10` | `5`  | `2`  | `0`  | `0`  | `0`  |

With `--save-rest`:

| `chrom` | `pos` | `strand` | `motif` | `u0` | `u1` | `u2` | `m0` | `m1` | `m2` | `o0` | `o1` | `o2` |
|---------|--------|----------|---------|------|------|------|------|------|------|------|------|------|
| `chr1`  | `10000` | `+`      | `AAG`   | `10` | `5`  | `2`  | `0`  | `0`  | `0`  | `1`  | `2`  | `3`  |
| `chr1`  | `10001` | `-`      | `TTC`   | `10` | `5`  | `2`  | `0`  | `0`  | `0`  | `1`  | `2`  | `3`  |


&nbsp;

<p align="center">
  <img
    src="https://raw.githubusercontent.com/y9c/y9c/master/resource/footer_line.svg?sanitize=true"
  />
</p>
<p align="center">
  Copyright &copy; 2025-present
  <a href="https://github.com/y9c" target="_blank">Chang Y</a>
</p>
<p align="center">
  <a href="https://github.com/y9c/countmut/blob/main/LICENSE">
    <img src="https://img.shields.io/static/v1.svg?style=for-the-badge&label=License&message=MIT&logoColor=d9e0ee&colorA=282a36&colorB=c678dd" />
  </a>
</p>
