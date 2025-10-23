# CountMut 🧬

[![CI](https://github.com/y9c/countmut/actions/workflows/ci.yml/badge.svg)](https://github.com/y9c/countmut/actions/workflows/ci.yml)
[![PyPI version](https://badge.fury.io/py/countmut.svg)](https://badge.fury.io/py/countmut)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Ultra-fast strand-aware mutation counter**

CountMut counts mutations from bisulfite sequencing /  BAM files with parallel processing, quality-based mate overlap deduplication, and optimized file I/O.

## Features

- 🚀 **Ultra-Fast**: Direct FASTA index reading, shared file handles, BGZF multi-threading
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

## Key Options

```bash
Required:
  -i, --input PATH           Input BAM file
  -r, --reference PATH       Reference FASTA file

Output:
  -o, --output PATH          Output TSV file (default: stdout)

Mutation:
  --ref-base TEXT            Reference base [default: A]
  --mut-base TEXT            Mutation base [default: G]
  --strand [both|forward|reverse]  Strand processing [default: both]
  --region TEXT              Specific region (e.g., 'chr1:1000000-2000000')

Performance:
  -t, --threads INTEGER      Number of threads [default: auto]
  -b, --bin-size INTEGER     Genomic bin size [default: 10000]

Filtering (Bisulfite):
  --pad INTEGER              Motif window padding [default: 15]
  --trim-start INTEGER       Trim 5' bases [default: 2]
  --trim-end INTEGER         Trim 3' bases [default: 2]
  --max-unc INTEGER          Max unconverted (Zf) [default: 3]
  --min-con INTEGER          Min converted (Yf) [default: 1]
  --max-sub INTEGER          Max substitutions (NS) [default: 1]
```

**Note**: BAM must have NS, Zf, and Yf tags for bisulfite analysis.

## Output Format

TSV file with columns:
- `chrom`, `pos`, `strand`, `motif` - Position and context
- `u0`, `d0` - Drop counts (trimmed/unmapped bases)
- `u1`, `d1` - Clean counts (converted reads)
- `u2`, `d2` - Unconverted counts

## Performance

Real-world test (3 rRNA genes, 17k reads):
- **Before optimization**: 16.7s
- **After optimization**: 4.5s (3.7x faster)
- **Startup**: 0.2s vs 3.2s (16x faster FASTA index loading)

## Python API

```python
from countmut import count_mutations

success = count_mutations(
    samfile="input.bam",
    reffile="reference.fa",
    output_file="mutations.tsv",
    ref_base="T",
    mut_base="C",
    threads=8,
    max_unc=3,
    min_con=1,
    max_sub=1
)
```

## Requirements

- Python 3.10+
- pysam, rich, click, rich-click, numpy
- BAM files: coordinate-sorted with `.bai` index (auto-created if missing)
- FASTA files: with `.fai` index (auto-created if missing)

## Development

```bash
# Clone and install
git clone https://github.com/y9c/countmut.git
cd countmut
uv sync --dev

# Run tests
uv run pytest tests/ -v

# Lint
uv run ruff check .
```

## License

MIT License - see [LICENSE](LICENSE) file.

## Citation

If you use CountMut in your research, please cite:

```
Ye Chang. (2025). CountMut: Ultra-fast parallel mutation counter for bisulfite sequencing.
GitHub: https://github.com/y9c/countmut
```

## Contact

- **Author**: Ye Chang
- **Email**: yech1990@gmail.com
- **Issues**: https://github.com/y9c/countmut/issues

---

**Version 0.0.1** | **MIT License**
