# CountMut filter expressions (`-e` / `-p`)

`-e` / `-p` are evaluated **inside the C core** by an embedded Lua 5.4 state
(the approach used by pbr).  There is no Python expression engine and no
per-base Python fallback: the fast C backend evaluates the expressions itself.

**CountMut has no dedicated filter/trim CLI flags.**  All read-level QC and
trimming is a `-e` expression, all site-level QC is a `-p` expression:

| old flag / concept            | expression                                     |
|-------------------------------|------------------------------------------------|
| `min_mapq`                    | `-e "mapq >= N"`                               |
| `min_baseq` (quality QC)      | `-e "bq >= N"`                                 |
| `max_sub` (NS tag)            | `-e "tag('NS') <= N"`                          |
| `max_unc` / `min_con` (Yf/Zf) | `-e "tag('Zf') <= N and tag('Yf') >= N"`       |
| mismatches (NM tag)           | `-e "tag('NM') <= N"`                          |
| include/exclude SAM flags     | `-e "flags & N == N"`, `-e "flags & N == 0"`   |
| fragment 5'/3' trim           | `-e "dist5 >= N and dist3 >= N"`               |
| read R1 3'-end trim           | `-e "not (flags & 64 ~= 0 and qpos >= length - N)"` |
| read R2 5'-start trim         | `-e "not (flags & 128 ~= 0 and qpos < N)"`     |
| `min_depth`                   | `-p "depth >= N"`                              |
| `min_allele_support`          | `-p "g >= N"` (or a/c/t), `-p "depth >= N and g > 0"` |
| basic strand selection        | `-e "strand == 1"` (fwd) / `-e "strand == -1"` (rev) |

Only genuinely structural options stay as flags: input/reference/output,
`--engine`, `--region`, `--threads`, `--ref-base`/`--mut-base`/`--pad`/`--save-rest`,
`--strandless`, `--count-indels`, `--max-depth`, `--vcf`, `--verbose`.
(The output view — mutation / base / allele — is inferred: mutation when
`--ref-base` & `--mut-base` are both given, allele when `--vcf`, else base;
there is no `--mode` flag.)

| Flag | Scope     | Evaluated when                |
|------|-----------|-------------------------------|
| `-e, --expression` | read-level | for every aligned read base (or once per read, see below) |
| `-p, --pile-expression` | site-level | once per reported site (before output) |

An expression that evaluates to `false`/`nil` excludes that base (read) or
omits that site (pile).  Runtime errors reject the item; a **syntax error is
fatal** (exit code 2) so a typo is never silently ignored.

## Language

The grammar is the **samtools `filter=STRING` expression language** (C-style
precedence, bit fields, tags, regex, string-stat helpers), evaluated with Lua
5.4 semantics.  Two source styles are accepted:

```
countmut -i x.bam -r ref.fa -e "mapq >= 20 and bq >= 20"
countmut -i x.bam -r ref.fa -e "return read.mapq >= 20"      # pbr style
```

Bare predicates (the default countmut style) are auto-wrapped in
`return (...)`; a source that already begins with `return` is used verbatim.

### Operator sugar (commonly used Python / samtools spellings)

| as written        | Lua meaning |
|-------------------|-------------|
| `!=`              | `~=` (not equal) |
| `&&`, `and`       | `and`       |
| `||`, `or`        | `or`        |
| unary `!`         | `not`       |
| `A =~ "regexp"`   | `re_match(A, "regexp")` |
| `A !~ "regexp"`   | `not re_match(A, "regexp")` |

### Tag and flag shorthand (samtools idioms)

| as written          | rewrites to                  | notes |
|---------------------|------------------------------|-------|
| `[XX]`              | `tag('XX')`                  | value of aux tag `XX` (`Z`/`i`/`f`; `nil` if absent) |
| `exists([XX])`      | `exists('XX')`               | `true` iff the tag is present |
| `flag.read1`        | `flags & 64`                 | any `flag.<bit>` keyword |
| `flags & READ1 != 0`| direct Lua                   | symbolic constants below |

Bit keywords: `flag.paired proper_pair unmap munmap reverse mreverse read1
read2 secondary qcfail dup supplementary`.

Symbolic constants: `PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE MREVERSE READ1
READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY` (the numeric SAM flag bits), e.g.
`flags & (SECONDARY|DUP) == 0`.

## Read namespace (`-e`)

Flat globals and the corresponding `read.<name>` fields.  Values that are
constant over a read are computed **once per read**; values that change per
aligned base are computed per base.

### Per-base (evaluated for every aligned base)

| name | meaning |
|------|---------|
| `qpos` / `QPOS` | 0-based position inside the read |
| `bq` / `BQ` / `baseq` | base (Phred) quality at the current position, `-1` if none |
| `base` | aligned base at the current position (`A C G T N` / IUPAC) |
| `ref` | reference base at the current position |
| `dist5` / `DIST5` / `distance_from_5prime` | bases to the fragment 5′ end (reverse-strand reads are measured from the read 3′ end) |
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
| `mrname` / `rnext` | mate reference name ("" if none) |
| `seq` / `sequence` | read sequence (uppercase bases) |
| `qual` | raw quality bytes as a string (ASCII; Phred = byte − 33) |
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

> **Gotcha — missing tags.** `tag('XX')` is `nil` when the tag is absent.  In a
> boolean context `nil` is false, but comparing / doing arithmetic with `nil`
> raises an error (which rejects the read).  Guard presence first:
> `exists('NM') and tag('NM') <= 3`, or write the whole filter to tolerate it
> (`exists('NM') ~= true or tag('NM') <= 3`).

### Performance

Filters that reference **only read-constant values** (`mapq`, `flags`, `tlen`,
`tag('XX')`, `rname`, …) are detected at compile time and evaluated **once per
read** — in **both** engines:
- read-walk evaluates them once in the read loop (and skips the whole read when
  it fails, before walking any base);
- pileup memoises them by pileup slot, so a read spanning ~100 positions still
  only pays one evaluation.

Measured on the deep rRNA transcriptome fixture (784 k reads, 90.5 Mb):
- read-constant `-e` adds ≈0–0.1× full-genome (e.g. `mapq >= 20` ≈ +0.5–0.7 s
  on an ~11 s read-walk mutation run; rejecting filters can be *faster* than
  baseline because the read is skipped before its bases are walked);
- an expensive read-constant filter (regex `rname =~ 'rRNA.*'`) went from ~93 s
  (the old per-position path) to ~12 s via the once-per-read memo;
- only expressions that touch `qpos`/`bq`/`base`/`ref`/`dist5`/`dist3` are
  evaluated per aligned base.  That cost is **~55–140 ns/base** (measured, full
  genome: read-walk mutation +1.6 s, pileup base +4.3 s for `bq >= 20`; +3.9 s /
  +10.9 s for `bq and dist5 and base`).  It is kept low by
  - materialising **only the fields the expression references**,
  - computing `seq`/CIGAR stats **lazily** (only when the expr needs them — no
    full read-sequence decode per eval for `bq`),
  - populating the `read.`/`pile.` table **only when the expression uses dotted
    access** — a bare predicate like `bq >= 20` skips it entirely (`read.bq`
    costs ~1.7× more).
  Read-walk mutation also evaluates ~3× fewer bases than pileup (it only visits
  positions whose reference base matches `--ref-base`), so per-base filters are
  cheapest in that configuration.

Pileup is the faster engine at ultra-deep hotspots (mutation + `-e "bq >= 20"`:
4.3 s vs read-walk 18 s); read-walk wins where the filter rejects many reads
early.

## Pile namespace (`-p`, per site)

`depth` (total base depth, both strands), `pos` (0-based), `ref`/`ref_base`,
`a c g t n` (= `A C G T N`), `ins`, `del`, `ref_skip`, `fail`; and in mutation
mode the reference window as `motif`.

## Examples

```
# unique, high-quality, not near the 5' end
countmut -i x.bam -r ref.fa -e "mapq >= 20 and bq >= 20 and dist5 >= 2"

# keep forward-strand, properly-paired reads
countmut -i x.bam -r ref.fa -e "strand == 1 and flags & 2 != 0"

# restrict to an RG group
countmut -i x.bam -r ref.fa -e "tag('RG') == 'sampleA'"

# samtools-style: few mismatches, not a PCR duplicate, read 1 only
countmut -i x.bam -r ref.fa -e "[NM] <= 3 and not (flag.dup ~= 0) and flag.read1 != 0"

# regex on the reference name (rRNA / snRNA drop-out)
countmut -i x.bam -r ref.fa -e "rname =~ 'rRNA.*'"

# report only A-reference sites with depth >= 5 and > 2 G alleles
countmut -i x.bam -r ref.fa -p "ref == 'A' and depth >= 5 and g > 2"
```

Both engines (read-walk and pileup) apply the filters identically and are
byte-identical with `-e` / `-p` set.
