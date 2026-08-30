# CountMut

Fast, strand-aware **base / mutation / allele counting** from BAM files.
C core (htslib-subset + embedded Lua), thin Python CLI, both walk strategies
(read-walk and pileup) produce byte-identical output.

```bash
pip install -e .
```

## Quick start

```bash
# C→T conversion rate at target sites          -> mutation table
countmut -i in.bam -r ref.fa -o mut.tsv --ref-base C --mut-base T

# every base, per strand                       -> base table
countmut -i in.bam -r ref.fa -o depth.tsv

# alleles as VCF                               -> VCF
countmut -i in.bam -r ref.fa --vcf -o allele.vcf
```

There is **no `--mode` flag** — the output follows the inputs:

| you give | output columns |
|---|---|
| nothing | `chrom pos strand ref depth a c g t n` |
| `--ref-base C --mut-base T` | `chrom pos strand motif u0 u1 u2 m0 m1 m2 mutation_rate` |
| `--vcf` | VCF `GT:AD` |

`mutation_rate` = converted / (converted + unconverted). Output is **per
strand** by default (like the mutation view); `--strandless` merges the two.

## Filtering with `-e` / `-p`

Read-level QC and trimming are **expressions**, not flags. `-e` filters reads
per base, `-p` filters sites.

```bash
# quality + keep away from read ends
countmut -i x -r ref -o out -e "mapq >= 20 and bq >= 20 and dist5 >= 2"

# one sample (RG tag)
countmut -i x -r ref -o out -e "tag('RG') == 'sampleA'"

# samtools-style: low mismatch, not a PCR dup, read1 only
countmut -i x -r ref -o out -e "[NM] <= 3 and not (flag.dup ~= 0) and flag.read1 != 0"

# site-level: only well-covered sites, ≥2 G reads
countmut -i x -r ref -o out -p "depth >= 5 and g >= 2"
```

**You'll use ~10 variables 90% of the time:** `mapq`, `bq` (base quality),
`flags`, `qpos` (position in read), `dist5`/`dist3` (distance to read ends),
`base`/`ref`, `tag('XX')`, `rname`.

**Gotchas (read these once):**
- **Missing tags** — if a read has no `NM`, then `tag('NM')` is nothing, and
  `tag('NM') <= 3` *errors* (that read is dropped). Guard with `exists('NM')`.
- **Per-base fields** — only `qpos, bq, base, ref, dist5, dist3` change per
  base. Everything else (`mapq`, `flags`, tags, `rname`, …) is evaluated once
  per read and costs essentially nothing.
- A **syntax error exits with code 2** — a typo is never silently ignored.

Full grammar: [`docs/filter_grammar.md`](docs/filter_grammar.md)
· Exhaustive reference: [`docs/expression_reference.md`](docs/expression_reference.md)

## Engine

`--engine {auto|read-walk|pileup}` (default `auto`) — `auto` uses read-walk for
the targeted mutation view and pileup otherwise. Either way the output is
identical.

## Options

```
-i/--input  -r/--reference  -o/--output
--region  --threads/-t  --engine
--ref-base  --mut-base  --pad  --save-rest
--strandless  --count-indels  --vcf   # + --min-depth, --min-allele-support
-e/--expression  -p/--pile-expression
```

## Why it's fast

- read-constant filters (`mapq >= 20`, `[NM] <= 3`) run **once per read** in
  both engines — ≈ free, even on whole genomes
- per-base filters (~55–140 ns/base) materialize **only** the fields you use
- both engines are C; `--threads` splits the genome across workers

## License

MIT
