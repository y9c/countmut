# Filter expressions — the 30-second guide

`-e` filters **reads** (per aligned base), `-p` filters **sites**. There are no
old-school filter/trim flags — everything QC-like is an expression, evaluated
in the fast C backend (embedded Lua). Both engines apply them identically.

## Quick reference: what you used before

| old flag / idea | write this |
|---|---|
| min mapq | `-e "mapq >= 20"` |
| min base quality | `-e "bq >= 20"` |
| NS / Yf / Zf tags | `-e "tag('NS') <= 3"`, `-e "tag('Zf') <= 3 and tag('Yf') >= 3"` |
| mismatches (NM) | `-e "[NM] <= 3"` |
| include/exclude SAM flags | `-e "flags & 2 != 0"`, `-e "flags & (SECONDARY\|DUP) == 0"` |
| trim fragment 5′/3′ | `-e "dist5 >= 2 and dist3 >= 2"` |
| trim read ends | `-e "not (flag.read1 ~= 0 and qpos >= length - 10)"` |
| strand | `-e "strand == 1"` (fwd) / `-e "strand == -1"` (rev) |
| min depth | `-p "depth >= 20"` |
| min allele support | `-p "g >= 2"` (or `a`/`c`/`t`) |

## The 10 variables you need most

| var | means |
|---|---|
| `mapq` | mapping quality |
| `bq` | base (Phred) quality at this position |
| `flags` | SAM flag bitfield (`flag.dup`, `flag.reverse`, … also work) |
| `qpos` | 0-based position inside the read |
| `dist5` / `dist3` | bases to the read 5′ / 3′ end |
| `base` / `ref` | aligned base / reference base here |
| `tag('XX')` | any aux tag value (`[XX]` is the short form: `[NM] <= 3`) |
| `rname` | reference name |

Bare names work (`mapq >= 20`); pbr-style `read.mapq` works too. `!=`, `&&`,
`||`, `!` are accepted alongside Lua's `~=`, `and`, `or`, `not`; `=~`/`!~` do
regex (`rname =~ 'rRNA.*'`).

## Gotchas

- **Missing tags.** If a read has no `NM`, then `tag('NM')` is nothing and
  `tag('NM') <= 3` errors → that read is dropped. Test first:
  `exists('NM') and tag('NM') <= 3`.
- **Per-base fields are the only ones that cost.** `qpos, bq, base, ref,
  dist5, dist3` are evaluated per aligned base (~tens of ns); everything else
  runs once per read and is essentially free.
- **Syntax errors are fatal** (exit code 2) so a typo can't silently pass.
- `avg(qual)` returns the *raw stored byte* mean (Phred = that − 33); on some
  BAMs qual is stored as Phred+33 and on others as plain Phred.

## Examples

```bash
# high-quality, away from read ends, forward-strand & properly paired
countmut -i x -r ref -o out \
  -e "mapq >= 20 and bq >= 20 and dist5 >= 2 and strand == 1 and flags & 2 != 0"

# drop rRNA/snRNA contigs (regex)
countmut -i x -r ref -o out -e "rname =~ 'rRNA.*'"

# report only A-reference sites, depth ≥ 5, > 2 G alleles
countmut -i x -r ref -o out -p "ref == 'A' and depth >= 5 and g > 2"
```

Need every variable, function, or the performance model?
See [`expression_reference.md`](expression_reference.md).
