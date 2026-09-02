# Expression reference (the details)

Reference notes for `-e`/`-p` filter expressions.  For a quick guide see
[`filter_grammar.md`](filter_grammar.md); this is the complete spec.

- `-e, --expression <STR>` — **read** filter **and group router**, per aligned
  base (or once per read, see *Read-constant* below).
- `-p, --pile-expression <STR>` — **site** filter, once per reported site.

Return-value semantics:

| `-e` returns          | effect                                             |
|-----------------------|----------------------------------------------------|
| `nil` / `false`       | base dropped (the old filter behavior)             |
| `true`                | base counted into **group 0** (old filter behavior)|
| integer `0` … `3`     | base counted into that group                       |
| anything else         | base dropped + warning on stderr                   |

Groups surface in `--output-format` templates as per-group cells `{a.0}` …
`{n.3}` (e.g. a 2-group A→G view: `{a.0}\t{a.1}\t{g.0}\t{g.1}`); the plain
`{a}` … `{n}` cells stay the per-strand total over all groups.  The `-p` site
filter stays a plain boolean predicate.  A **syntax error is fatal (exit code
2)**; a runtime error just rejects that item.

## Language

The grammar is **samtools `filter=STRING`**: C-style precedence, bit fields,
tags, regex, string-stat helpers, evaluated with Lua 5.4 semantics.  Both
source styles are fine:

```
countmut -i x.bam -r ref.fa -e "mapq >= 20 and bq >= 20"
countmut -i x.bam -r ref.fa -e "return read.mapq >= 20"      # pbr style
```

Bare predicates are auto-wrapped in `return (...)`; a source already starting
with `return` is used verbatim.  `-p` is the same but with the pile namespace.

### Operator sugar

| as written        | Lua meaning |
|-------------------|-------------|
| `!=`              | `~=` (not equal) |
| `&&`, `and`       | `and`       |
| `||`, `or`        | `or`        |
| unary `!`         | `not`       |
| `A =~ "re"`       | `re_match(A, "re")` |
| `A !~ "re"`       | `not re_match(A, "re")` |

### Tag and flag shorthand

| as written          | rewrites to                  | notes |
|---------------------|------------------------------|-------|
| `[XX]`              | `tag('XX')`                  | value of aux tag `XX` (`Z`/`i`/`f`; `nil` if absent) |
| `exists([XX])`      | `exists('XX')`               | `true` iff the tag is present |
| `flag.read1`        | `flags & 64`                 | any `flag.<bit>` keyword |

Bit keywords: `flag.paired proper_pair unmap munmap reverse mreverse read1
read2 secondary qcfail dup supplementary`.

Symbolic constants: `PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE MREVERSE READ1
READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY` — e.g. `flags & (SECONDARY|DUP) == 0`.

## Read namespace (`-e`)

Every name is a flat global **and** a `read.<name>` field.  Two cost classes:

### Per-base (evaluated for every aligned base)

| name | meaning |
|------|---------|
| `qpos` / `QPOS` | 0-based position inside the read |
| `bq` / `BQ` / `baseq` | base (Phred) quality at the current position, `-1` if none |
| `base` | aligned base at the current position (`A C G T N` / IUPAC) |
| `ref` | reference base at the current position |
| `dist5` / `DIST5` / `distance_from_5prime` | bases to the fragment 5′ end (reverse reads measure from their 3′ end) |
| `dist3` / `DIST3` / `distance_from_3prime` | bases to the fragment 3′ end |

### Read-constant (evaluated once per read)

| name | meaning |
|------|---------|
| `mapq` / `MAPQ` / `read.mapping_quality` | mapping quality |
| `flag` / `flags` / `FLAGS` | SAM flag bitfield |
| `strand` / `STRAND` | `+1` forward / `-1` reverse (biological, paired-aware) |
| `pos` / `POS` / `read.start` | 1-based reference start |
| `endpos` / `stop` | 1-based reference end (CIGAR span) |
| `qlen` / `length` / `LEN` | read length (query bases) |
| `rlen` | reference length consumed (CIGAR) |
| `ncigar` | number of CIGAR operations |
| `tid` / `refid` | reference index |
| `mtid` | mate reference index |
| `mpos` / `pnext` | 1-based mate position (0 if unmapped) |
| `tlen` / `insert_size` | TLEN (template length, signed) |
| `sclen` | soft-clipped bases |
| `hclen` | hard-clipped bases |
| `n_indel` / `indel_count` | number of insertions+deletions in the CIGAR |
| `soft_clips_5_prime` / `soft_clips_3_prime` | soft-clip length at each end |
| `is_reverse` / `is_paired` / `r1` / `r2` | boolean flag tests (1/0) |
| `qname` | read name |
| `rname` | reference name |
| `mrname` / `rnext` | mate reference name (`""` if none) |
| `seq` / `sequence` | read sequence, uppercase, already reference-forward |
| `qual` | raw quality bytes as a string (Phred = byte − 33 on standard BAMs) |
| `library` | the `LB` aux tag (`""` if absent) |

### Functions

| function | meaning |
|----------|---------|
| `tag('XX')` | aux tag value (`nil` if absent) |
| `exists('XX')` | `true` iff tag `XX` is present (also `exists([XX])`) |
| `re_match(s, "re")` / `re_find(s, "re")` | POSIX ERE match anywhere in `s` |
| `slen(s)` / `length(s)` | string length |
| `smin(s)` / `min(s)`, `smax(s)` / `max(s)` | min / max byte value in a string |
| `savg(s)` / `avg(s)` | mean byte value (`NaN` on empty) — so `avg(qual) - 33` is mean Phred |
| `n5(n)` / `n3(n)` | fraction of `N`s in the first / last `n` read bases |
| `read.n_proportion_5_prime(n)` / `read.n_proportion_3_prime(n)` | same, pbr-style |
| `sqrt(x)` `log(x)` `exp(x)` `pow(x,y)` | math helpers (Lua `math.*` also available) |

### Missing-tag semantics

`tag('XX')` is `nil` for an absent tag.  `nil` is false in a boolean context,
but **comparing or doing arithmetic with `nil` raises an error** (that read is
rejected).  Guard presence:

```
exists('NM') and tag('NM') <= 3          # skip reads without NM
exists('NM') ~= true or tag('NM') <= 3   # ...or tolerate absence
```

## Pile namespace (`-p`, per site)

`depth` (total base depth, both strands), `pos` (0-based), `ref`/`ref_base`,
`a c g t n` (= `A C G T N`), `ins`, `del`, `ref_skip`, `fail`; in the mutation
view the reference window is also exposed as `motif`.

When `-e` routes bases into groups, the pile-level `a c g t n` are the totals
over **all groups** (both strands); per-group cells (`a.0` … `n.3`) are only
available in `-o` row templates, where they are per-strand.

## Performance

Two cost classes decide how fast a filter is:

**Read-constant** (anything not in the per-base list): detected at compile
time, evaluated **once per read** in both engines — read-walk in the read loop
(and a failed read is skipped before any base is walked), pileup via a
per-slot memo.  Measured full genome (784 k reads, 90.5 Mb): `mapq >= 20`
≈ +0.5–0.7 s on an ~11 s run (≈ free); a rejecting filter can be *faster*
than baseline because whole reads are skipped.  A regex filter
(`rname =~ 'rRNA.*'`), previously ~93 s on the old per-position path, runs in
~12 s.

**Per-base** (only `qpos, bq, base, ref, dist5, dist3`): evaluated per aligned
base at **~55–140 ns/base** (measured: full-genome `bq >= 20` adds +1.6 s
read-walk / +4.3 s pileup; `bq and dist5 and base` +3.9 s / +10.9 s).  Kept
low because the engine materializes only the fields the expression references,
computes `seq`/CIGAR lazily, and skips the `read.`/`pile.` tables entirely for
bare predicates (`read.bq >= 20` costs ~1.7× more than `bq >= 20`).

**Engine choice** also matters: at ultra-deep hotspots pileup is the faster
walk (`-e "bq >= 20"`: 4.3 s vs read-walk 18 s); read-walk wins when a filter
rejects many reads early (whole reads are skipped before their bases are
walked).

## Examples

```bash
# unique, high-quality, not near the 5' end
countmut -i x.bam -r ref.fa -e "mapq >= 20 and bq >= 20 and dist5 >= 2"

# forward-strand, properly-paired
countmut -i x.bam -r ref.fa -e "strand == 1 and flags & 2 != 0"

# RG group
countmut -i x.bam -r ref.fa -e "tag('RG') == 'sampleA'"

# samtools-style: few mismatches, not a PCR duplicate, read 1 only
countmut -i x.bam -r ref.fa -e "[NM] <= 3 and not (flag.dup ~= 0) and flag.read1 != 0"

# rRNA / snRNA drop-out (regex)
countmut -i x.bam -r ref.fa -e "rname =~ 'rRNA.*'"

# A-reference sites with depth >= 5 and > 2 G alleles
countmut -i x.bam -r ref.fa -p "ref == 'A' and depth >= 5 and g > 2"

# group router: bisulfite A->G, 2 groups (g1 = high-conversion bases,
# g0 = low-quality / read-end bases that pass the hard NS gate); per-strand
# 2-group view via the template
countmut -i x.bam -r ref.fa \
  -e "([NS] <= 1) and (([Yf] >= 1 and [Zf] <= 3 and bq >= 20 and qpos >= 2 and qlen - qpos > 2) and 1 or 0)" \
  --motif-pad 15 \
  --fmt-header "chrom\tpos\tstrand\tmotif\tu0\tu1\tm0\tm1" \
  --output-format "{chrom}\t{pos+1}\t{strand}\t{motif}\t{a.0}\t{a.1}\t{g.0}\t{g.1}"
```

Both engines apply the filters identically and emit byte-identical output with
`-e` / `-p` set.
