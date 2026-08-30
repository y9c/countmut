# CountMut

At each site in a modification assay you want the same small number: how many
reads still show the reference base, how many show the conversion, and what
that is as a rate. Existing tools made this hard — a flag for every QC idea,
two BAM-walking strategies that disagreed on deep overlapping sites, and
read-level filters priced once per aligned position instead of once per read.

CountMut collapses this into one design: QC and trimming are expressions (the
samtools grammar, in C), the two walks provably agree, and read-level rules
run once per read. On a deep rRNA transcriptome (784 k reads / 90 Mb) a
genome-wide `mapq >= 20` adds ~0.5 s, the per-base filter was cut ~3×, and
the mutation view ends with a `mutation_rate`.

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

There is no `--mode` flag: the output follows what you asked for. Bare runs
give the per-strand base composition; adding a reference/mutation pair gives
the conversion view with a `mutation_rate` column; `--vcf` gives an allele
VCF. Output is per strand by default, and `--strandless` merges the two.

## Filtering with one expression instead of ten flags

In RNA there is no genomic mutation to count: the "converted" base is a
modification read out through reverse transcription, so the conversion rate
reports modification level rather than a variant. Whatever your sample, the
QC and trimming live in one expression language — read-level rules are `-e`
expressions, site-level rules `-p`, in the samtools `filter=` grammar,
evaluated inside the C core. The old `--min-mapq` / `--min-baseq` /
`--trim-*` flags are gone.

```bash
# quality, and not on the error-prone read ends
countmut -i x -r ref -o out -e "mapq >= 20 and bq >= 20 and dist5 >= 2"

# one sample
countmut -i x -r ref -o out -e "tag('RG') == 'sampleA'"

# samtools-style: low mismatch, not a PCR duplicate, read 1 only
countmut -i x -r ref -o out -e "[NM] <= 3 and not (flag.dup ~= 0) and flag.read1 != 0"

# site-level: only well-covered sites with ≥2 G reads
countmut -i x -r ref -o out -p "depth >= 5 and g >= 2"
```

Most filters use roughly ten variables — `mapq`, `bq` (base quality), `flags`,
`qpos` (position in the read), `dist5`/`dist3` (distance to the read ends),
`base`/`ref`, `tag('XX')`, and `rname`. A couple of things are worth knowing.
A missing tag is nothing: comparing `tag('NM')` *errors* and drops that read,
so guard with `exists('NM')` when tags are optional. Only the six per-base
values (`qpos`, `bq`, `base`, `ref`, `dist5`, `dist3`) cost anything to
evaluate; everything else runs once per read and is essentially free. A syntax
error stops the run (exit code 2), so a typo can never silently change your
numbers.

The full grammar is in
[`docs/filter_grammar.md`](docs/filter_grammar.md), with an exhaustive
reference in [`docs/expression_reference.md`](docs/expression_reference.md).

## Engines and options

Two BAM-walking strategies live in the C core. `--engine auto` (default) uses
the read walk for the targeted mutation view and the pileup walk otherwise;
both produce identical output, so the choice only affects speed. The remaining
options are few: input/reference/output, `--region`, `--threads/-t`,
`--engine`, `--ref-base`, `--mut-base`, `--pad`, `--save-rest`,
`--strandless`, `--count-indels`, `--vcf` (+ `--min-depth`,
`--min-allele-support`), and `-e`/`-p`.

## Input formats

**BAM** (indexed, fast, threaded) and **SAM** (plain or gzipped — auto-transcoded
to a temp BAM + index, same output as the equivalent BAM) are both supported,
detected automatically. **CRAM** is not read by this self-contained core; convert
first: `samtools view -b in.cram -o out.bam` — for a CRAM with an embedded
reference, that conversion also works without a separate FASTA.

## Performance

![scaling + filter overhead](docs/perf-scaling.png)

Measured on a bimodal benchmark (232 k reads, 23.2 M read-bases, deep rRNA-style
hotspots): read-walk mutation 2.14 s @1 thread → 0.92 s @16; pileup base
1.96 s → 0.69 s (dynamic work queue keeps deep hotspots from serializing).
Read-constant filters run once per read (≈ free); only true per-base filters
`qpos/bq/base/ref/dist5/dist3` add cost (~0.3 s @8 threads for `bq and dist5`).
`tests/make_bench_bam.py` regenerates the fixture; `scripts/plot_perf.py` re-renders
this figure.

## License

MIT