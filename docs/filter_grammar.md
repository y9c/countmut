# CountMut filter-grammar proposal (`-e` / `-p`)

> **Status: proposal for review.** Two string filters only. The grammar is a
> **strict samtools-style predicate** over a flat namespace of scalar values
> (no object/method fields, no `exec`/`import`/assignment), with samtools-style
> variable names. Implemented as a *safe Python-expression subset*.

---

## 1. The two filters

| Flag | Scope | Evaluated when |
|------|-------|----------------|
| `-e, --expression <STR>` | read-level | once per aligned base |
| `-p, --pile-expression <STR>` | site-level | once per reported site |

The expression must evaluate to a boolean. A `false` read is excluded; a `false`
site is omitted. With neither flag, the fast C path runs with no extra filtering.

---

## 2. Values & operators

* **number** (int/float), **string** (`'A'` or `"sample"`), **boolean** (`true`/`false`).
* comparison: `==` `!=` `<` `<=` `>` `>=` (also `in` / `not in`)
* logical: `&&` `||` `!` (samtools style); `and` `or` `not` are accepted as aliases
* arithmetic: `+` `-` `*` `/` `//` `%`
* bitwise (samtools flag tests): `&` `|` `^`
* regex: `=~` / `!~`
* grouping: `( ... )`

There are **no** separate `--min-mapq` / `--min-baseq` / `--trim-*` flags —
express them with `-e` / `-p` instead (e.g. `-e "bq >= 20 and dist5 >= 2"`).

---

## 3. Namespace — read variables (`-e`)

Flat, scalar, samtools-style names (lower-case and UPPER-case aliases accepted).

| Name | Type | Description |
|------|------|-------------|
| `mapq` / `MAPQ` | int | read mapping quality |
| `baseq` / `BQ` / `MIN_BQ` | int | base quality (Phred) at the current position (`-1` if none) |
| `flags` / `FLAGS` | int | raw SAM flag bitmask |
| `strand` / `STRAND` | int | `+1` forward, `-1` reverse (biological/paired-aware) |
| `qname` / `QNAME` | str | query name |
| `length` / `LEN` | int | read length |
| `dist5` / `DIST5` / `distance_from_5prime` | int | bases from fragment 5′ end |
| `dist3` / `DIST3` / `distance_from_3prime` | int | bases from fragment 3′ end |
| `tag('XX')` | value | BAM aux tag (string/int); missing → `None` |

Flag masks may be written symbolically: `flags & UNMAP == 0`, `flags & (SECONDARY|DUP) == 0`,
`flags & PAIRED != 0`, or numerically (`flags & 4 == 0`).

---

## 4. Namespace — site variables (`-p`)

| Name | Type | Description |
|------|------|-------------|
| `depth` / `DEPTH` | int | observed base depth |
| `pos` / `POS` | int | 0-based genomic position |
| `ref` / `REF` / `ref_base` | str | reference allele |
| `a c g t n` / `A C G T N` | int | per-residue base counts |
| `ins` / `INS`, `del` / `DEL`, `ref_skip`, `fail` | int | indels / ref-skips / filter-failures |

---

## 5. Examples (samtools style)

    # unique, high-quality, not at the 5' end
    countmut -i x.bam -r ref.fa -e "MAPQ >= 20 and BQ >= 20 and dist5 >= 2"

    # keep only properly-paired reads
    countmut -i x.bam -r ref.fa -e "flags & 2 != 0 and flags & UNMAP == 0"

    # restrict to an RG group
    countmut -i x.bam -r ref.fa -e "tag('RG') == 'sampleA'"

    # report only A-reference sites with depth >= 5 and > 2 G alleles
    countmut -i x.bam -r ref.fa -p "ref == 'A' and depth >= 5 and g > 2"

---

## 6. Safety & notes

* Evaluated via `ast` under a restricted node allow-list; errors/`NameError` →
  `false` (rejected), never raised; values are read-only.
* Symbolic flag constants available: `PAIRED PROPER_PAIR UNMAP MUNMAP REVERSE
  MREVERSE READ1 READ2 SECONDARY QCFAIL DUP SUPPLEMENTARY`.
* When `-e` or `-p` is set, processing uses the Python engine (expressions are
  not evaluated in the C core, which stays the zero-overhead path otherwise).

