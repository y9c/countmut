# CountMut

Counting a modification assay — bisulfite DNA, or an RNA modification library —
always comes down to the same small question: at each site, how many reads
still show the reference base and how many show the conversion, and what is
that as a rate? The tools that answered it made the small thing hard. Quality
control meant another flag for every idea. Two BAM-walking strategies existed
and quietly disagreed exactly where the data gets hard — the deep, heavily
overlapped sites. And a filter as cheap as `mapq >= 20` was priced once per
aligned position instead of once per read, so something that should be free
could silently add minutes to every whole-genome run.

CountMut removes all of that with one design. QC and trimming stopped being
flags and became expressions — the samtools filter grammar, evaluated in C.
The two walks fill the same per-site count table, so they provably agree and
the engine choice is only about speed. The parser reads a filter once, sees
which variables it actually touches, and runs every read-level rule exactly
once per read; only the six values that genuinely change per base are
evaluated per base.

It shows in the numbers. On a deep rRNA transcriptome (784 k reads, 90 Mb), a
genome-wide `mapq >= 20` costs about half a second on top of an ~11 s run,
the one filter that must run per base was cut roughly threefold (per-base
evaluation ~180 → ~55 ns), and both engines output byte-identical rows. The
mutation view ends with a `mutation_rate` column — the number that actually
belongs in the paper.

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

## License

MIT
