# CountMut

CountMut is a C-speed, strand-aware counter for BAM files: point it at a BAM
and a reference and it tells you, at every site, what the reads actually
show. Give it a conversion pair and it turns that into a mutation or
conversion rate.

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

## DNA or RNA?

The counting itself is agnostic — it compares every read base against a
reference and a target base — but what the ratio means depends on your sample.
In **DNA**, `--ref-base C --mut-base T` is the classic bisulfite or
cytosine-damage readout: a genuine C→T change on one strand. In **RNA** there
is no genomic mutation to find; the transcript itself is the reference, and a
modified base (for example m5C, in the chemistries that detect it) makes
reverse transcription misread it — a modified C reads out as a U that the
sequencer records as T. So the column labelled "converted" is really the
**modification signal**, and the `mutation_rate` is a readout of how much of
the RNA at that site carries the modification — not a variant frequency.

Because RNA is single-stranded, strand matters differently than in DNA. The
transcript is copied from only one strand of the genome, and a modification
rides on that one RNA strand, so the two biological strands are *not*
symmetrical — which is precisely why countmut keeps the strand per row, and
why a per-strand filter such as `-e "strand == 1"` is a meaningful question in
an RNA library.

## Filtering with one expression instead of ten flags

Read-level QC and trimming are expressions (`-e`), site-level rules are
expressions too (`-p`), in the samtools `filter=` grammar, evaluated inside the
C core. The old `--min-mapq` / `--min-baseq` / `--trim-*` flags are gone.

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

## Why it's designed this way

CountMut answers a simple question from modification and damage assays: at a
given site, how many reads show the reference base versus the converted one,
and what is that as a rate? Earlier tools made it harder than it needed to be —
they buried QC under dozens of filter flags, forced a choice between two
BAM-walking strategies that disagreed on paired-end overlaps, and re-priced a
cheap read-level filter at every base of every read, which quietly dominates
deep ribosomal hotspots.

The design fixes all three at once. QC and trimming are one expression
language evaluated in the C core, so there are no filter flags to grow. Both
walking strategies share a single count structure and are byte-identical, so
the engine is purely a speed choice. And the parser decides at compile time
which variables an expression touches: read-constant filters run once per read
(a whole-genome `mapq >= 20` costs about half a second), while the six
per-base values stay accurate but cheap. The mutation view ends with a
`mutation_rate`, and read QC stays honest because BAMs already store reverse
reads reference-forward.

## License

MIT
