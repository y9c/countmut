# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

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