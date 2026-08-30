# CountMut filter expressions (`-e` / `-p`)

`-e` / `-p` are evaluated **inside the C core** by an embedded Lua 5.4 state
(the approach used by pbr).  There is no Python expression engine.

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

Only genuinely structural options stay as flags: input/reference/output, `--mode`,
`--engine`, `--region`, `--threads`, `--ref-base`/`--mut-base`/`--pad`/`--save-rest`,
`--split-strand`, `--count-indels`, `--max-depth`, `--vcf`, `--verbose`.

| Flag | Scope     | Evaluated when                |
|------|-----------|-------------------------------|
| `-e, --expression` | read-level | once per aligned base (after trim, before overlap dedup) |
| `-p, --pile-expression` | site-level | once per reported site (before output) |

An expression that evaluates to `false`/`nil` excludes that base (read) or
omits that site (pile).  Runtime errors reject the item; a **syntax error is
fatal** (exit code 2) so a typo is never silently ignored.

Both **styles** are accepted:

```
countmut -i x.bam -r ref.fa -e "MAPQ >= 20 and bq >= 20"
countmut -i x.bam -r ref.fa -e "return read.mapq >= 20"
```

Bare predicates (the original countmut style) are auto-wrapped in
`return (...)`; pbr-style chunks with an explicit `return` work as-is.
For convenience the Python/samtools operator spellings are translated to Lua:
`!=` → `~=`, `&&` → `and`, `||` → `or`, unary `!` → `not`.

## Read namespace (`-e`, per aligned base)

Flat globals (and the `read` table): `mapq`/`MAPQ`, `bq`/`baseq`/`BQ`
(base quality at the current position, `-1` if none), `flags`/`flag`/`FLAGS`,
`strand`/`STRAND` (`+1`/`-1`, biological, paired-aware), `qname`, `rname`,
`pos`/`POS` (1-based start), `endpos`, `qlen`/`length`/`LEN`, `rlen`, `ncigar`,
`qpos`/`QPOS` (0-based position in the read), `dist5`/`DIST5` and
`dist3`/`DIST3` (bases to the fragment 5′/3′ ends).  `tag('XX')` reads an aux
tag (`exists('XX')` tests presence), `'Z'`/`'i'`/`'f'` types supported.

SAM flag constants are pre-defined: `PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE
MREVERSE READ1 READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY` -- e.g.
`flags & (SECONDARY|DUP) == 0`.

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

# report only A-reference sites with depth >= 5 and > 2 G alleles
countmut -i x.bam -r ref.fa -p "ref == 'A' and depth >= 5 and g > 2"
```

Both engines (read-walk and pileup) apply the filters identically and are
byte-identical with `-e` / `-p` set.
