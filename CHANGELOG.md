# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.6] - 2026-08-30

### Changed
- **All filtering and trimming are expressions (`-e`/`-p`), not flags.**
  The countmut CLI removed its filter flags (`--min-depth`, `--min-allele-support`)
  and never exposes read filters/trims; every read-level QC lives in `-e`
  (`mapq`, `bq`, `tag('NS'/'Zf'/'Yf'/'NM')`, `flags & ...`, `dist5`/`dist3`/
  `qpos` for trimming) and every site-level QC in `-p` (`depth`, `ref`,
  `a/c/g/t/n`, ...).  See `docs/filter_grammar.md` for the complete mapping.
- Mutation quality-QC moved to `-e "bq >= N"`: all bases that pass `-e` land in
  tier `x2` (the internal quality boundary is 0; `x0`/`x1` stay 0 unless you
  opt into tiered reports).  The C core keeps the low-level flags for internal
  use/tests but defaults to "no filtering" (`max_sub/unc/con = -1`, no trim).
- `FilterConfig`/`EngineConfig` slimmed to structural settings only.

## [0.1.5] - 2026-08-30

### Performance
- **read-walk engine: solo-vs-overlap counting.**  Only reads whose mate can
  also cover a base go through the `(pos,qname)` overlap-dedup hash; bases that
  only one mate can cover (solo bases, single-end reads, mates on another
  contig) are counted **directly into the site**.  This removes the huge
  accumulating dedup hash at deep sites, where previously every fragment ×
  target base inserted an entry (~10M entries at a 50k-deep rRNA site).
  Overlapping mates are still deduped exactly; unknown mate geometry (mpos<0 /
  TLEN unusable) keeps the exact hash path for those reads.
- **Effect on the deep PE transcriptome BAM** (784k reads, rRNA hotspots
  ~50–100k depth):
  - full-genome mutation read-walk: **74 s -> 5.9 s** (pileup: 3.4 s)
  - deep rRNA region read-walk: **60.7 s -> 4.3 s** (old binary: 60.7 s)
- Read-walk now also matches pileup on sparse BAMs (82k reads / 41k
  transcripts: read-walk 0.84 s vs pileup 0.90 s).
- Verified: read-walk == pileup byte-identical on all fixtures/modes; the
  overlap-dedup depth stays exact (10 fragments -> depth 10); full test suite
  (46) passes; output is unchanged for the intended overlapping-mate case.

## [0.1.4] - 2026-08-30

### Fixed / changed
- **`--max-depth` is now unlimited (0) by default** — the pileup engine counts
  *all* reads.  Previously the wrapper silently used 8000 and, worse, htslib's
  own 8000 pileup ceiling truncated very deep sites regardless.  Now wired
  through `bam_mplp_set_maxcnt` (0 = unlimited), `--max-depth N` caps the
  per-position depth for users who want it, and the CLI exposes `--max-depth`.
- **Many-contig runs**: output no longer uses one temp file *per region* (which
  exhausted the file-descriptor limit on transcriptomes with 40k+ contigs);
  each worker now writes its own contiguous region slice to one temp file and
  slices are concatenated in order.
- **`--verbose` progress**: throttled to ~100 updates, shows completion % and
  the process's current RSS (MB) in real time; the CLI now streams the C
  progress to stderr live (instead of swallowing it).
- New `scripts/countmut_monitor.py` — a runtime memory monitor that samples the
  process-tree RSS of any countmut run and reports peak + a timeline.

### Verified
- On a 40k-transcript RNA-Seq PE BAM with ~63k-read (100k-fragment) rRNA
  hotspots: default depth is now the true fragment count (51,084 at a hotspot),
  `--max-depth 100` caps at 100; full-genome mutation runs (2.1M rows) complete
  and are byte-identical across runs.

## [0.1.3] - 2026-08-30

### Performance (results unchanged -- byte-identical output, no logic change)
- **read-walk engine**: indel-free reads jump straight to the (sorted) target
  positions instead of walking every aligned base; the `(pos,qname)` overlap
  dedup now keys on integer qname ids (one copy per unique qname instead of a
  `strdup` per entry) and uses a single hash probe per base.
- **read-level filter memoization**: the NS/Zf/Yf checks are computed once per
  read in the pileup engine (was: once per read per pileup position).
- Reference sequence pre-uppercased once per contig (removes per-base
  `toupper`); the mutation motif window is built once per site instead of once
  per strand.
- Benchmark (500 kb chr, 40x, 200k reads), median of 3:
  - mutation: **0.64s -> 0.36s (1.8x)** single-thread, 0.108s -> 0.072s (8 thr)
  - base: **0.82s -> 0.69s (1.2x)** single-thread.
- Verified byte-identical to the previous binary on mutation/base/allele,
  1 & 8 threads, with and without `-e`/`-p`, plus the full regression suite and
  the cpup / legacy countmut external comparisons.

## [0.1.2] - 2026-08-30

### Changed
- **Architecture: C is the only counting implementation.**  The Python
  read-walk and pileup engines (`engine_readwalk.py`, `engine_pileup.py`) and
  their orchestration (`pipeline.py`) were removed.  `--engine read-walk` and
  `--engine pileup` are both fully implemented in `backend/countmut_core`
  (the read-walk engine is new), and Python is now a thin wrapper that only
  builds the command line (`backend.py`) and drives the CLI.
- **`-e` / `-p` expressions moved into C via embedded Lua 5.4** (pbr-style),
  replacing the Python expression engine.  The original countmut bare-predicate
  style (`mapq >= 20`) still works (auto-wrapped in `return (...)`) and the
  `!=`/`&&`/`||`/`!` spellings are translated to Lua (`~=`/`and`/`or`/`not`).
  A syntax error is fatal (exit 2).  See `docs/filter_grammar.md`.
- Python fallbacks removed: if the C binary is missing, `run_backend` raises
  instead of silently using a Python engine.
- `-p` `depth` is the total base depth (both strands), matching perbase/pbr.

### Fixed
- `min_allele_support` was only honoured in the VCF allele output; it now also
  filters the plain allele table.
- Validation vs external references: `samtools mpileup | cpup` 0/1480
  base-count mismatches on forward-only data; legacy `count_mutations`
  agrees on all shared sites (modulo the documented x1/x2 column inversion).

## [0.1.1] - 2026-08-30

### Fixed
- **`-` strand mutation rows were internally inconsistent**: the motif was
  reverse-complemented while the u/m counts were reference-forward, so the
  motif's reference base disagreed with the counted bases.  Both strands now
  show the reference-forward window, matching perbase/pbr/mpileup and the
  reference-forward bases all three backends already counted.
- **read-walk engine ignored deletions / ref-skips / failures**, so base mode
  with `--count-indels` (and the emitted position set) diverged from the pileup
  engine and C core.  It now tallies `del`/`ref_skip`/`fail` per read and emits
  those positions, matching pileup/C exactly.
- **pileup engine ignored the strand gate** (`--strand forward/reverse`), so
  unsplit base output counted both strands; it now skips wrong-strand reads
  like the read-walk engine and C core.
- **`--min-depth` was a no-op**: applied in the Python base/allele formatters and
  in the C base-mode (split) branch.
- **Python allele mode**: fixed the header (was 5 columns "alleles_and_counts"
  while rows have 7) and wired `min_allele_support`; added VCF output parity
  with the C core (`--vcf` in the Python fallback path).
- **C region parsing was off by one**: `--region chr1:1-60` was treated as
  `[1,60)` (0-based, dropping the first position) instead of samtools' 1-based
  `[0,60)`; matches the Python wrapper now.
- **supplementary reads**: Python now keeps them (same as the C core / samtools
  default `UNMAP|SECONDARY|QCFAIL|DUP`), fixing a C/Python divergence.
- `MutationConfig` bases are case-normalized (lowercase `--ref-base a` works),
  matching the C driver.
- `FilterConfig.exclude_flags` default fixed to 1796 (`UNMAP|SECONDARY|QCFAIL|DUP`).

### Tests
- Added `tests/test_correctness.py` (non-palindromic reference) covering the
  `-`-motif consistency, read-walk indel parity, pileup strand gate, `min-depth`,
  allele mode shape/support and base-case handling.

## [0.1.0] - 2026-08-30

### Added
- Unified tool: `--mode {mutation, base, allele}` and `--engine {auto, read-walk, pileup}`.
- **C backend** (`backend/countmut_core`): pileup walk, per-(site, strand) base counts,
  mate-overlap dedup, quality/conversion classification; built on a self-contained htslib
  subset from [lh3/minipileup](https://github.com/lh3/minipileup); auto-built on first use;
  multicore (`--threads`).
- **samtools-style filter expressions** `-e/--expression` (read) and `-p/--pile-expression`
  (site) — the samtools `filter=STRING` grammar (`&&`/`||`/`!`, `flag.dup`, `[NM]`, `=~`,
  `avg(qual)`, `exists([NM])`); syntax documented in `docs/filter_grammar.md`.
- `allele` mode with `--vcf` (minipileup style); BED include/exclude (`-b`/`-x`);
  samtools `--incl-flags`/`--excl-flags`/`--input-fmt-option`.
- Tests: `tests/test_unified.py`, `tests/test_expression.py`.

### Changed
- Removed granular numeric filter flags (`--min-mapq`, `--min-baseq`, `--max-sub`,
  `--max-unc`, `--min-con`, `--trim-*`) in favour of `-e`/`-p` expressions.
- Reverse-strand base handling fixed (no more double-complement); trim applied before
  mate dedup (a usable mate's base is kept). See `docs/legacy_comparison.md`.

### Correctness
- Byte-identical output between the C backend and the Python read-walk/pileup engines;
  memory-clean under AddressSanitizer.

## [0.0.8] - 2025-01-27

### Added
- **Region Chunking for High-Density Reads**: Automatic splitting of regions with too many reads
  - New `--max-reads-per-chunk` option (default: 100,000) to control chunking threshold
  - Regions exceeding the threshold are automatically split into smaller chunks
  - Chunks are processed in parallel for improved performance on high-density regions (e.g., rRNA genes)
  - Optimized BAM file handle reuse during chunking phase

### Performance
- Improved handling of regions with extremely high read density (e.g., >10M reads)
- Parallel processing of chunked regions reduces lagging on dense genomic regions

## [0.0.7] - 2025-01-27

### Fixed
- Fixed `tag_read_with_alternative_mutations` to properly return the modified read

## [0.0.6] - 2025-10-24
### Fixed
- Resolved `UndefinedName` errors by correctly initializing counters and result lists.
- Fixed `ProcessPoolExecutor` shutdown issue (`cannot schedule new futures after shutdown`).
- Ensured `Progress` bar displays after worker warmup, improving perceived startup time.
### Changed
- Implemented per-process shard files for parallel BAM output, centralizing BAM writing and reducing temporary file creation.
- Refined worker initialization (`_init_worker`) to use `os.getpid()` for worker IDs and removed redundant `worker_id` passing.
- Added `countmut_` prefix to temporary directories created by `tempfile.TemporaryDirectory`.
### Added
- Introduced worker warmup phase to reduce initial processing delay.
- `_warmup_worker` and `_worker_shutdown_task` functions for explicit worker lifecycle management.

## [0.0.5] - 2025-10-24

### Added
- **Alternative Mutation Tagging**: New feature to count a secondary set of mutations
  - `--ref-base2` and `--mut-base2` flags to specify alternative mutations (e.g., C->T)
  - Adds `Yc` (alt-mut) and `Zc` (alt-ref) tags to each read
  - Corrects `NS` tag by subtracting alternative mutation counts
- **BAM Output**: `--output-bam` flag to save the newly tagged BAM file
  - If not specified, a temporary BAM is used for counting and then deleted
  - Enables seamless integration with other bioinformatic tools

### Performance
- **Thread-Safe Tagging**: Fully parallelized read processing and BAM writing
  - Each worker writes to a temporary BAM file
  - Merged and sorted in parallel for maximum speed

### Changed
- The main `count_mutations` function now orchestrates the tagging and counting process
- CLI updated with new "Alternative Mutation Tagging" option group

## [0.0.4] - 2025-10-24

### Fixed
- **Critical bug fix: Strand-aware trimming** - Trimming now correctly respects fragment orientation
  - Forward strand (+): `trim_start` from beginning (5'), `trim_end` from end (3')
  - Reverse strand (-): `trim_start` from END (5'), `trim_end` from BEGINNING (3')
  - Previous behavior incorrectly applied uniform trimming regardless of strand
- Fixed off-by-one error in trim logic: `>` changed to `>=` for correct base counting
  - `trim_start=2` now correctly skips 2 bases (not 3)

### Changed
- Enhanced code comments explaining fragment orientation vs genomic mapping
- Clarified documentation that trimming is based on **fragment orientation** (5'→3' of molecule)

### Why This Matters
For bisulfite sequencing, enzymatic/adapter artifacts occur at the molecular ends of DNA fragments:
- Reverse-strand fragments have 5' molecular end at the END of the read sequence
- Trimming must respect fragment biology, not genome coordinates
- Critical for accurate analysis of eTAM-seq, GLORI-seq, CAM-seq, and BS-seq data

## [0.0.3] - 2025-10-24

### Changed
- **Major internal refactoring for clarity**: Renamed internal variables for better code readability
  - `is_converted` → `passes_conversion_filter` (clearer boolean logic)
  - `drop_count` → `low_quality_count` (failing quality filters)
  - `clean_count` → `high_conversion_count` (high conversion efficiency)
  - `unc_count` → `insufficient_conversion_count` (poor conversion efficiency)
- **Reorganized CLI option groups** for better UX:
  - Combined Required/Output into "Input/Output Options"
  - Merged Quality and Bisulfite filters into single "Quality Filters" group
  - Added "Output Records" group for pad and save-rest
- **Enhanced documentation**:
  - Updated README feature description: "Call mutation without pileup reads"
  - Improved count category descriptions with accurate terminology
  - Added `-p` short flag for `--pad` option
  - Fixed `--save-rest` description to show correct column names (o0,o1,o2)
  - Added markdown table for output format
  - Better blockquote styling for header

### Fixed
- Corrected misleading variable names that confused conversion status
- Updated all documentation to reflect accurate filtering logic

### Documentation
- x0 (low quality): Bases failing quality filters
- x1 (high conversion): Reads with high conversion efficiency (low Zf and high Yf)
- x2 (insufficient conversion): Reads with poor conversion efficiency (high Zf or low Yf)

## [0.0.2] - 2025-10-24

### Added
- **Organized CLI option groups** for better readability and user experience
  - Required Options (input, reference)
  - Output Options (output, save-rest, force)
  - Mutation Analysis (ref-base, mut-base, strand, region)
  - Performance Options (threads, bin-size)
  - Quality Filters (min-baseq, min-mapq, max-sub, trim-start, trim-end)
  - Bisulfite Filters (max-unc, min-con)
  - Sequence Context (pad)
  - Help & Version
- Enhanced CLI help text with clearer descriptions
- Improved output format documentation in README

### Changed
- Renamed `--min-base-qual` to `--min-baseq` for consistency
- Added `--min-mapq` option for filtering by mapping quality (default: 0)
- Reorganized filtering parameters into logical groups
- Updated CLI docstring to reflect latest features and output format
- Enhanced README with detailed explanation of output columns (drop/clean/unconverted)

### Fixed
- Query base extraction now correctly uses `read.query_sequence[query_pos]` instead of reference base
- Simplified `get_aligned_pairs` to use `with_seq=False` for better performance
- Fixed tuple unpacking for modern pysam (always 2 elements when `with_seq=False`)

### Performance
- Using `get_aligned_pairs(matches_only=True, with_seq=False)` for optimal performance
- No reference sequence construction overhead
- Direct query position to base mapping

## [0.0.1] - 2025-10-23

### Added
- Initial release of CountMut
- Fast, parallel mutation counting from BAM pileup data
- Bisulfite conversion analysis support (NS, Zf, Yf tags)
- Multi-threaded genomic window processing
- Rich CLI interface with progress tracking
- Support for custom reference/mutation bases
- Region-specific analysis capabilities
- Memory-efficient streaming processing
- Strand-specific processing (both/forward/reverse)
- Configurable filtering parameters (pad, trim, thresholds)
- Quality-based mate overlap deduplication
- Comprehensive test suite
- Complete documentation and examples

### Performance Optimizations
- ⚡ 7x faster FASTA index loading via direct .fai reading
- 🚀 3.7x overall speedup through multiple optimizations
- 💾 Shared file handles per worker process (eliminates redundant open/close)
- 🔧 BGZF multi-threading support for parallel BAM decompression
- 🎯 Smart read iteration with get_aligned_pairs(matches_only=True, with_seq=True)
- 🔍 Early region skipping for empty regions
- 📊 Removed granular timing overhead from hot paths

### Features
- **Parallel Processing**: Multi-threaded genomic window processing for maximum speed
- **Thread Safety**: Completely thread-safe implementation using ProcessPoolExecutor
- **Memory Efficient**: Streaming processing prevents memory overflow on large files
- **Bisulfite Analysis**: Built-in conversion detection and analysis
- **Rich Output**: Beautiful CLI with progress bars and detailed statistics
- **Flexible Configuration**: Support for custom bases, regions, and processing parameters
- **Error Handling**: Robust error handling and resource cleanup
- **Performance Optimized**: Optimized for high-performance processing of large BAM files

### Technical Details
- Built with modern Python tooling (uv, ruff, rich, click)
- Thread-safe parallel processing using ProcessPoolExecutor
- Efficient memory usage with streaming BAM file processing
- Comprehensive error handling and resource cleanup
- Extensive test coverage including thread safety tests
- Performance benchmarking and optimization
- GitHub Actions for continuous integration and deployment

### Dependencies
- Python 3.10+
- pysam >= 0.21.0
- rich >= 13.0.0
- click >= 8.0.0
- rich-click >= 1.6.0
- numpy >= 1.21.0

### Installation
```bash
pip install countmut
```

### Basic Usage
```bash
# Count A→G mutations (default)
countmut input.bam reference.fa

# Count T→C mutations (common in bisulfite sequencing)
countmut input.bam reference.fa --ref-base T --mut-base C

# Save to file with custom parameters
countmut input.bam reference.fa -o mutations.tsv --bin-size 5000 --threads 8
```

### Performance
- Typical performance on modern workstation:
  - 1 GB BAM: ~2 minutes (8 threads)
  - 5 GB BAM: ~8 minutes (8 threads)
  - 10 GB BAM: ~15 minutes (8 threads)
- Parallel processing provides 2-4x speedup over sequential processing
- Memory usage scales linearly with bin size, not file size

### Thread Safety
- Completely thread-safe implementation
- Each worker process opens its own file handles
- No shared state between workers
- Proper error handling and resource cleanup
- Extensive testing for concurrent access patterns

### Repository
- GitHub: https://github.com/y9c/countmut
- PyPI: https://pypi.org/project/countmut/
- Documentation: https://github.com/y9c/countmut#readme
- Issues: https://github.com/y9c/countmut/issues

### License
MIT License - see LICENSE file for details

### Author
Ye Chang (yech1990@gmail.com)

---

**Note**: This is the initial release of CountMut. Future releases will include additional features, performance improvements, and bug fixes based on user feedback and requirements.