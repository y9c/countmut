# Legacy vs unified countmut comparison

Dated: 2026-08-30.  Compares the original pure-Python `count_mutations`
(`countmut/core.py`) with the new unified tool (C core + engines) on a synthetic
bisulfite BAM (C→T, `NS`/`Zf`/`Yf` tags, 300 bp, 80 forward-only reads).

## Headline result

On forward-only reads the **total mutation counts agree 100%** (1306 = 1306 over
every site), but the legacy tool **puts them in the wrong columns**.

```
site chr1:4        legacy: u1=0 u2=0 m1=1 m2=0     new: u1=0 u2=0 m1=0 m2=1
```

`count_mutations` emits high-conversion counts in the **`x1`** columns and
leaves `x2` empty, whereas its own README documents `x1`=insufficient,
`x2`=high. The unified tool emits high-conversion in `x2` (consistent with the
documentation). This is a **column-mapping inversion** in the legacy code, not a
counting difference.

## Where the legacy tool diverges (and the unified tool improves)

1. **Column semantics inverted.** In `core.py` the code writes
   `u1/m1 = high_conversion_count` and `u2/m2 = insufficient_conversion_count`,
   but `insufficient_conversion_count` is never populated (conversion-failed
   reads are skipped before categorisation), so `x2` is always 0 and all passing
   bases land in `x1`. README says the opposite. The unified tool follows the
   documented intent.

2. **Reverse-strand reads mishandled.** Countmut reverse-complements
   `query_sequence` for `-` strand. Because BAM stores `SEQ` already in
   reference-forward orientation, this flips the base to the wrong residue, so
   reverse-strand contributions are dropped/mis-counted. The unified tool keeps
   `bam_seqi`/`query_sequence` as reference-forward (no complement) and reports
   both strands. To stay consistent with those reference-forward base counts the
   motif is also rendered reference-forward for **both** strands (no reverse
   complement on the `-` row). On the mixed-strand dataset legacy reported
   **111** sites vs **206** for the unified tool (the missing 95 are
   reverse-strand sites).

3. **Mandatory bisulfite tags.** Legacy `get_tag("NS"/"Zf"/"Yf")` raises and
   skips the read if the tag is absent, so generic alignments produce nothing.
   The unified tool auto-detects whether the BAM carries the tags and only then
   enforces the conversion filter.

4. **Dedup tie-break.** Legacy keeps the highest base-quality observation per
   `(ref_pos, qname)`; the unified tool keeps `(mapq, read1, base-qual)`. The
   unified rule is the mate-aware rule from perbase/pbr plus a base-quality
   tie-break, and it is applied *before* trimming (a usable mate's base is kept).

## Test data

`/tmp/opencode/bis/` — generated with pysam. Note pysam does not reverse-complement
on write, so these tests write the aligned (reference-forward) sequence explicitly,
mimicking bwa output.

## Re-run

```
count_mutations(...)             # legacy: writes legacy.tsv
countmut_core --bam test.bam --fa ref.fa --mode conversion \
  --engine read-walk --ref-base C --mut-base T --pad 15 ...   # new (C core)
```
