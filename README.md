# CountMut

> **Unified ultra-fast strand-aware mutation counter** — C backend, Python wrapper.

CountMut counts base/substitution ratios from BAM files with **fast C core** and
a thin Python wrapper. It fuses the two classic ways of walking a BAM:

- **pileup-based** (`bam_mplp_auto` / pysam pileup) — fast, sees indels/ref-skips, general.
- **read-walk** (countmut's "no pileup") — walk reads directly, only touch the target sites.

Both produce **identical** output, and the tool can process whole genomes in
parallel (threads).

## Why it's fast & correct

- The hot loop (BAM read, pileup, per-(site,strand) base counting, mate-overlap
  dedup, quality/conversion classification) is in **C** (`backend/countmut_core`,
  built on the self-contained htslib subset from lh3/minipileup).
- **Strand-aware** (countmut biological-strand rule for paired-end reads).
- **Paired-end overlap dedup**: at an overlapping position a fragment is counted
  once, choosing the best mate by `(mapq, read1, base-qual)` — the thing
  minipileup gets wrong.
- **Parallel**: divides the genome into bins and processes them across threads
  (`--threads`).
- **Memory-clean** (verified under AddressSanitizer).

## Install

```bash
pip install -e .
# or, to prebuild the C core:
make backend
```

## Quick start

```bash
# strand-aware A->G mutation count (bisulfite / m6A style)
countmut -i in.bam -r ref.fa -o mut.tsv --ref-base A --mut-base G

# per-site base counts (perbase/mpileup style)
countmut -i in.bam -r ref.fa --mode base -o depth.tsv

# alleles -> VCF (minipileup style)
countmut -i in.bam -r ref.fa --mode allele --vcf -o allele.vcf
```

## Modes

| Mode | Output |
|------|--------|
| `mutation` | `chrom pos strand motif u0 u1 u2 m0 m1 m2 [o0 o1 o2]` (strand-aware substitution table) |
| `base` | `chrom pos [strand] ref depth a c g t n [ins del ref_skip fail]` |
| `allele` | `chrom pos ref depth ref_count alt alt_count`, or VCF with `--vcf` |

## Filtering with expressions (`-e` / `-p`)

Filtering is done with **samtools-style filter expressions** — there are no
separate `--min-mapq`/`--min-baseq`/`--trim-*` flags; write them as expressions
instead.

* `-e, --expression <STR>` — per-base **read** filter (samtools SAM fields).
* `-p, --pile-expression <STR>` — per-**site** filter (pileup fields).

Grammar is the samtools `filter=STRING` expression language (C-style precedence,
`&&`/`||`/`!`, bit fields, tags, regex). See `docs/filter_grammar.md`.

```bash
# keep high-quality, non-5prime, properly paired reads
countmut -i x.bam -r ref.fa -e "mapq >= 20 && bq >= 20 && dist5 >= 2 && flag & PROPER_PAIR"

# restrict to one RG group
countmut -i x.bam -r ref.fa -e "tag('RG') == 'sampleA'"

# report only A-reference sites with depth >= 5 and > 2 G alleles
countmut -i x.bam -r ref.fa -p "ref == 'A' && depth >= 5 && g > 2"
```

Read variables: `mapq`, `flag` (+ `flag.dup`, `flag.unmap`, ...), `qname`, `pos`,
`endpos`, `pnext`, `rname`, `mrname`, `tlen`, `qlen`, `rlen`, `ncigar`, `seq`,
`qual`, `sclen`, `hclen`, `bq`, `dist5`/`dist3`, `strand`, `[NM]`/`[RG]` tags,
`avg(qual)`, `exists([NM])`, `sqrt(mapq)`, ...

Site variables: `depth`, `pos`, `ref`, `a c g t n`, `ins`, `del`, `ref_skip`, `fail`.

> When `-e`/`-p` is given, counting runs on the Python engine (the C core cannot
> evaluate strings). Without expressions, the fast C backend is used.

## Engine selection

`--engine {auto|read-walk|pileup}` (default `auto`):

- `auto` → read-walk for `mutation` (targeted sites), pileup for `base`/`allele`.
- `read-walk` / `pileup` → force a strategy.

## Options

```
-i/--input, -r/--reference, -o/--output
--mode {mutation,base,allele}   --engine {auto,read-walk,pileup}
--region, --threads/-t
--ref-base, --mut-base, --pad, --save-rest
--split-strand, --count-indels, --min-depth, --min-allele-support, --vcf
-e/--expression, -p/--pile-expression
```

## Design

```
countmut/
  cli.py               rich CLI (routes to the backend)
  backend.py           builds/loads the C binary; calls it; Python fallback
  pipeline.py          region binning + parallel dispatch
  model.py             FilterConfig / MutationConfig / StrandConfig / EngineConfig
  engine_readwalk.py   read-walk engine
  engine_pileup.py     pileup engine
  formatter.py         TSV/VCF renderers
  expression.py        samtools-style filter-expression engine
backend/
  countmut_core.c      the computation core (htslib subset)
  countmut_core_main.c CLI wrapper
  Makefile             builds the `countmut_core` binary
```

Both engines fill the same `SiteColumn` (per-site, per-strand base counts), so
the two "ways" are interchangeable; the C backend implements the pileup engine.

## References this tool learns from

- [minipileup](https://github.com/lh3/minipileup) — pileup walk, filters, allele counting
- [perbase](https://github.com/sstadick/perbase) / [pbr](https://github.com/brentp/pbr) — mate-aware overlap dedup, base counts
- [countmut](https://github.com/y9c/countmut) — biological strand, bisulfite NS/Zf/Yf tiers
- [mpileup](https://github.com/y9c/mpileup) / [cpup](https://github.com/y9c/cpup) — base-count output
- samtools `--input-fmt-option filter=STRING` — the expression grammar

## License

MIT
