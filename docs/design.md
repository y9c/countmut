# Design: pain points, improvements, and the smart solutions

## The pain points

1. **Filter-flag explosion.** Old-style tools had a knob for everything —
   `--min-mapq`, `--min-baseq`, `--min-depth`, `--min-allele-support`,
   `--max-sub`, `--trim-*`… Dozens of orthogonal flags to wire, test, and
   remember; every new QC idea meant a new flag.
2. **Two engines that didn't agree.** Pileup sees indels/ref-skips; read-walk
   sees reads directly. Classic tools pick one and hope — and overlap-dedup
   (counting a fragment twice at the overlap, or keeping the wrong mate) is
   exactly where minipileup-style tools go wrong.  Users couldn't trust either
   walk.
3. **Slow per-base filtering.** Base-level QC (quality, trim positions,
   distance-to-read-ends) naively re-evaluates per base *and* re-evaluates
   read-level filters hundreds of times per read (once per pileup position).
   At 50–100k× rRNA hotspots that is brutal — and whole-genome read filters
   that cost "once per position" instead of "once per read" were the silent
   killer.
4. **Unreadable conversion-assay output.** m5C / TAPS / bisulfite-like views
   dumped raw tiered counts (`u0 u1 u2 m0 m1 m2`) with no rate, and strand
   reality was confusing (split? merged?).
5. **Redundant "modes".** `mutation` vs `base` looked like two tools, but base
   output already contained everything mutation reports — choosing a mode
   forced the user to decide before seeing the data.

## The improvements

- **One expression language replaces the flag zoo.** Every read-level QC and
  trim is a `-e` expression, every site-level QC a `-p` expression — the full
  samtools `filter=` grammar, evaluated inside the C core (embedded Lua 5.4).
  No new QC idea needs a new flag; it needs a one-liner.
- **Both engines are interchangeable.** Two walks, one shared `[strand][tier]
  [base]` count structure, byte-identical output — verified continuously by
  the test suite (including qpos-dependent per-base filters).
- **Read-level QC is ~free.** Read-constant filters run **once per read** in
  both engines.  A whole-genome `mapq >= 20` adds ≈ +0.5–0.7 s.
- **Per-base filters made ~3× faster** (~55 ns/base vs ~180 ns/base), then a
  real bug found and fixed along the way: the pileup engine was silently using
  the wrong read position, corrupting `dist5`/`dist3`/`base` on reverse reads.
- **Views unified.** No `--mode` flag; the output follows the inputs, and the
  mutation view now carries `mutation_rate` with per-strand-by-default output
  (`--strandless` to collapse).

## The smart solutions

1. **One counting core, two walks, choose by shape.** Both engines fill the
   *same* `[strand][tier][base]` array, so `auto` picks read-walk for sparse
   targeted sites and pileup for dense modes — without ever changing the
   answer.
2. **Compile-time read-constant vs per-base classification.** The engine
   parses the expression once, knows exactly which variables it touches, and
   routes it to the once-per-read path or the per-base path.  That single idea
   is why "QC is free, trimming is fast".
3. **A pileup-slot memo that survives 100k× depth.** Keyed by the pileup slot
   pointer with a pos/qlen/qname verify — cheap integer lookups, correct
   against buffer recycling, and it doesn't collapse at ultra-deep sites the
   way a fixed-size (LRU) cache does.
4. **Lazy everything.** Only the fields an expression references are computed
   and materialized; `bq >= 20` never decodes the read sequence; the
   `read.`/`pile.` tables are skipped entirely for bare predicates.
5. **Stored SEQ is reference-forward.** BAM stores reverse-mapped reads
   reverse-complemented, so base counting never needs per-read complementing —
   one less thing to get wrong (verified against real data).
6. **Embedded Lua in the C hot loop.** Full expression power without a
   Python-process penalty, and pbr-style `return …` works next to bare
   predicates.
