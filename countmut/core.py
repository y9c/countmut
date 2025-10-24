#!/usr/bin/env python3
"""
Core mutation counting logic for CountMut.

This module provides the main functionality for counting mutations from BAM pileup data,
including bisulfite conversion analysis and parallel processing capabilities.

Author: Ye Chang
Date: 2025-10-23
"""

import atexit
import logging
import os
import shutil
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

import pysam
from rich.progress import Progress, TimeElapsedColumn, TimeRemainingColumn

from .utils import get_output_headers, write_output

# Set up logger
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.setLevel(logging.WARNING)

# DNA complement mapping for reverse complement
DNA_COMPLEMENT = str.maketrans("ATGCNatgcn", "TACGNtacgn")


# Determine the biological strand for a read
def determine_actual_strand(read: pysam.AlignedSegment) -> str:
    """Return '+' or '-' as the biological strand for the read.

    For paired-end reads: read1 forward = '+', read2 reverse = '+' (mirror logic),
    otherwise '-' respectively. For single-end, use read.is_reverse.
    """
    try:
        if read.is_paired:
            if read.is_read1:
                return "+" if not read.is_reverse else "-"
            # read2: reverse complemented indicates '+'
            return "+" if read.is_reverse else "-"
        # single-end
        return "+" if not read.is_reverse else "-"
    except Exception:
        return "+"


# Global per-process file handles initialized once per worker process
_GLOBAL_SAM = None
_GLOBAL_REF = None


def _init_worker(samfile_path: str, reffile_path: str) -> None:
    """Open BAM/FASTA once per process to avoid repeated open/close overhead and enable BGZF threads."""
    global _GLOBAL_SAM, _GLOBAL_REF
    _GLOBAL_SAM = pysam.AlignmentFile(samfile_path, "rb")
    try:
        # Enable multi-threaded BGZF decompression for faster reading
        _GLOBAL_SAM.set_threads(2)
    except Exception:
        pass
    _GLOBAL_REF = pysam.FastaFile(reffile_path)

    import atexit

    def _cleanup() -> None:
        try:
            if _GLOBAL_SAM is not None:
                _GLOBAL_SAM.close()
        except Exception:
            pass
        try:
            if _GLOBAL_REF is not None:
                _GLOBAL_REF.close()
        except Exception:
            pass

    atexit.register(_cleanup)


def reverse_complement(seq: str) -> str:
    """Get reverse complement of DNA sequence."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def read_fasta_index(fasta_path: str) -> dict[str, int]:
    """
    Quickly read chromosome lengths from FASTA index file (.fai).

    This is much faster than opening the full FASTA file with pysam.FastaFile(),
    especially for large reference files with many contigs.

    Args:
        fasta_path: Path to FASTA file (will look for .fai)

    Returns:
        Dictionary mapping chromosome names to their lengths
    """
    fai_path = fasta_path + ".fai"
    chrom_lengths = {}
    try:
        with open(fai_path) as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split("\t")
                    if len(parts) >= 2:
                        chrom_name = parts[0]
                        chrom_length = int(parts[1])
                        chrom_lengths[chrom_name] = chrom_length
    except FileNotFoundError as e:
        raise FileNotFoundError(
            f"FASTA index file not found: {fai_path}. "
            f"Please index your reference with: samtools faidx {fasta_path}"
        ) from e
    return chrom_lengths


def get_motif(seq: str, start: int, end: int) -> str:
    """Extract motif from sequence with padding."""
    seq_len = len(seq)
    if start < 0:
        left = "N" * abs(start)
    else:
        left = ""
    if end > len(seq):
        right = "N" * (end - seq_len)
    else:
        right = ""
    return left + seq[max(0, start) : min(seq_len, end)] + right


def parse_region_worker(args: tuple) -> dict[str, Any]:
    """
    Worker function for parsing a single genomic region.

    This function is designed to be completely thread-safe by:
    1. Opening its own file handles
    2. Not sharing any global state
    3. Returning results instead of writing to shared files
    4. Proper error handling and cleanup
    5. Processing both strands in one pass for efficiency
    """
    try:
        overall_start = time.time()
        # Unpack arguments
        (
            samfile_path,
            region_chrom,
            region_start,
            region_end,
            strand_option,  # both/forward/reverse
            ref_base,
            mut_base,
            ref_base2,
            mut_base2,
            temp_bam_path,
            save_rest,
            pad,
            trim_start,
            trim_end,
            max_unc,
            min_con,
            max_sub,
            min_baseq,
            min_mapq,
            worker_id,
        ) = args

        # Setup for optional BAM writing
        outfile = None
        if temp_bam_path:
            with pysam.AlignmentFile(samfile_path, "rb") as infile:
                outfile = pysam.AlignmentFile(temp_bam_path, "wb", header=infile.header)

        # Require per-process global handles (must be initialized by pool initializer)
        if _GLOBAL_SAM is None or _GLOBAL_REF is None:
            raise RuntimeError("Worker not initialized: missing global readers")
        samfile = _GLOBAL_SAM
        reffile = _GLOBAL_REF
        counts = []

        # Determine strand processing options
        process_forward_only = strand_option.lower() == "forward"
        process_reverse_only = strand_option.lower() == "reverse"

        # Parameters are passed from caller (CLI-configurable)

        # Check if chromosome exists in reference file
        contig_exists = region_chrom in reffile.references
        if not contig_exists:
            return {
                "worker_id": worker_id,
                "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
                "counts": [],
                "success": False,
                "error": f"invalid contig `{region_chrom}`",
                "reads": 0,
                "timings": {"total": time.time() - overall_start},
            }

        # Get target sequence and sites
        target_seq = reffile.fetch(region_chrom, region_start, region_end)
        target_sites_set = {
            i
            for i, b in zip(range(region_start, region_end), target_seq, strict=False)
            if b.upper() == ref_base
        }
        target_sites_list = sorted(target_sites_set)

        if not target_sites_set:
            return {
                "worker_id": worker_id,
                "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
                "counts": [],
                "success": True,
                "error": None,
                "reads": 0,
                "timings": {"total": time.time() - overall_start},
            }

        # Quick check: skip regions without any reads
        # This is much faster than doing full pileup processing
        try:
            # Count reads in the region (very fast operation)
            read_count = samfile.count(region_chrom, region_start, region_end)
            if read_count == 0:
                return {
                    "worker_id": worker_id,
                    "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
                    "counts": [],
                    "success": True,
                    "error": None,
                    "skipped": True,
                    "reason": "no_reads",
                    "reads": 0,
                    "timings": {"total": time.time() - overall_start},
                }
        except Exception:
            # If count fails, continue with normal processing
            pass

        # Add padding to target sequence
        extend_pad = 20
        target_seq_left = reffile.fetch(
            region_chrom, max(region_start - extend_pad, 0), region_start
        ).rjust(extend_pad, "N")
        target_seq_right = reffile.fetch(
            region_chrom, region_end, region_end + extend_pad
        ).ljust(extend_pad, "N")
        extended_target_seq = target_seq_left + target_seq + target_seq_right

        # Process reads directly (much faster than pileup)
        # Initialize counters for each target position and strand
        position_data = {}
        for pos in target_sites_list:
            position_data[pos] = {
                "+": {
                    "high_conversion_count": Counter(),
                    "insufficient_conversion_count": Counter(),
                    "low_quality_count": Counter(),
                },
                "-": {
                    "high_conversion_count": Counter(),
                    "insufficient_conversion_count": Counter(),
                    "low_quality_count": Counter(),
                },
            }

        # Process all reads in the region
        total_reads = 0
        skipped_wrong_strand = 0
        skipped_unmapped_dup_secondary = 0
        skipped_missing_tags = 0
        skipped_no_sequence = 0
        processed_reads = 0

        # Track best observation per (ref_pos, query_name) to avoid double counting overlapping mates
        # Value: (strand, base, qual, is_internal, is_mapped, is_converted)
        best_obs: dict[tuple[int, str], tuple[str, str, int, bool, bool, bool]] = {}

        reads_to_process = samfile.fetch(region_chrom, region_start, region_end)
        if outfile:
            # If we are writing a BAM, we must iterate over a copy,
            # as modifying tags can affect iteration over the original.
            reads_to_process = list(reads_to_process)

        for read in reads_to_process:
            total_reads += 1
            try:
                # If tagging is enabled, count alternative mutations and update tags
                if ref_base2 and mut_base2:
                    alt_ref_count = (
                        read.query_sequence.upper().count(ref_base2)
                        if read.query_sequence
                        else 0
                    )
                    alt_mut_count = (
                        read.query_sequence.upper().count(mut_base2)
                        if read.query_sequence
                        else 0
                    )

                    read.set_tag("Zc", alt_ref_count, "i")
                    read.set_tag("Yc", alt_mut_count, "i")

                    if read.has_tag("NS"):
                        ns_val = read.get_tag("NS")
                        read.set_tag("NS", max(0, ns_val - alt_mut_count), "i")

                # Determine which strand to process based on read orientation
                actual_strand = determine_actual_strand(read)

                # Skip if we don't want to process this strand
                if process_forward_only and actual_strand != "+":
                    skipped_wrong_strand += 1
                    continue
                if process_reverse_only and actual_strand != "-":
                    skipped_wrong_strand += 1
                    continue

                # Skip reads with deletions or reference skips
                if read.is_unmapped or read.is_duplicate or read.is_secondary:
                    skipped_unmapped_dup_secondary += 1
                    continue

                # Get read properties (avoid exceptions on missing tags)
                if not (
                    read.has_tag("NS") and read.has_tag("Zf") and read.has_tag("Yf")
                ):
                    skipped_missing_tags += 1
                    continue

                # Check mapping quality
                passes_mapq_filter = read.mapping_quality >= min_mapq

                # Check mismatch filter
                ns = read.get_tag("NS")
                passes_mismatch_filter = ns <= max_sub

                # Check conversion filter (reads passing are sufficiently converted)
                zf = read.get_tag("Zf")
                yf = read.get_tag("Yf")
                passes_conversion_filter = (zf <= max_unc) and (yf >= min_con)

                # Process each position in the read
                query_sequence = read.query_sequence
                if not query_sequence:
                    skipped_no_sequence += 1
                    continue
                query_qualities = read.query_qualities or []

                # Mark this read as successfully processed
                processed_reads += 1

                # Iterate via reference positions (fast path) and filter
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    # matches_only=True: only aligned positions (no insertions/deletions with None)
                    if query_pos is None or ref_pos is None:
                        continue
                    if ref_pos not in target_sites_set:
                        continue
                    if query_pos >= len(query_sequence):
                        continue

                    # Check if position is internal (not in trimmed regions)
                    # Trim based on actual strand orientation (fragment 5' to 3')
                    # For forward strand (+): 5' is at start (low query_pos), 3' is at end (high query_pos)
                    # For reverse strand (-): 5' is at end (high query_pos), 3' is at start (low query_pos)
                    if actual_strand == "+":
                        # Forward: trim_start from 5' (beginning), trim_end from 3' (end)
                        is_internal = (
                            query_pos >= trim_start
                            and len(query_sequence) - query_pos > trim_end
                        )
                    else:
                        # Reverse: trim_start from 5' (end), trim_end from 3' (beginning)
                        is_internal = (
                            query_pos >= trim_end
                            and len(query_sequence) - query_pos > trim_start
                        )
                    # Get query base from read sequence (base_char is reference base, not query)
                    query_base = query_sequence[query_pos].upper()
                    base_qual = (
                        int(query_qualities[query_pos])
                        if query_qualities and query_pos < len(query_qualities)
                        else 0
                    )

                    # Check base quality
                    passes_baseq_filter = base_qual >= min_baseq

                    if actual_strand == "-":
                        query_base = query_base.translate(DNA_COMPLEMENT)
                    key = (ref_pos, read.query_name)
                    prev = best_obs.get(key)
                    if (prev is None) or (base_qual > prev[2]):
                        best_obs[key] = (
                            actual_strand,
                            query_base,
                            base_qual,
                            bool(is_internal),
                            bool(passes_mismatch_filter),
                            bool(passes_mapq_filter),
                            bool(passes_baseq_filter),
                            bool(passes_conversion_filter),
                        )

            except (KeyError, AttributeError) as e:
                # Skip reads with missing tags or invalid data
                logger.debug(f"Skipping read due to missing tag: {e}")
                if outfile:
                    outfile.write(read)  # Write unmodified read on error
                continue

            if outfile:
                outfile.write(read)

        # Apply best observations to position counts (deduplicated across overlapping mates)
        for (ref_pos, _qname), (
            strand_symbol,
            query_base,
            _q,
            is_internal,
            passes_mismatch_filter,
            passes_mapq_filter,
            passes_baseq_filter,
            passes_conversion_filter,
        ) in best_obs.items():
            # Skip positions that are not in our target sites (safety)
            if ref_pos not in position_data:
                continue

            # Count in low_quality if ANY quality filter fails
            if not (
                is_internal
                and passes_mismatch_filter
                and passes_mapq_filter
                and passes_baseq_filter
            ):
                position_data[ref_pos][strand_symbol]["low_quality_count"][
                    query_base
                ] += 1
            # Count in high_conversion/insufficient_conversion only if ALL quality filters pass
            else:
                if passes_conversion_filter:
                    # High conversion efficiency (pass quality + conversion filters)
                    position_data[ref_pos][strand_symbol]["high_conversion_count"][
                        query_base
                    ] += 1
                else:
                    # Insufficient conversion efficiency (pass quality but fail conversion)
                    position_data[ref_pos][strand_symbol][
                        "insufficient_conversion_count"
                    ][query_base] += 1

        # Calculate total skipped
        total_skipped = (
            skipped_wrong_strand
            + skipped_unmapped_dup_secondary
            + skipped_missing_tags
            + skipped_no_sequence
        )

        # Log read processing statistics
        if total_reads > 0:
            logger.info(
                f"📊 Region {region_chrom}:{region_start}-{region_end}:{strand_option} - "
                f"Total: {total_reads}, Processed: {processed_reads}, Skipped: {total_skipped}"
            )
            if total_skipped > 0:
                logger.debug(
                    f"   Skipped details: wrong_strand={skipped_wrong_strand}, "
                    f"unmapped/dup/secondary={skipped_unmapped_dup_secondary}, "
                    f"missing_tags={skipped_missing_tags}, no_sequence={skipped_no_sequence}"
                )

        # Process each target position for each strand
        for pos in target_sites_list:
            for strand_symbol in ["+", "-"]:
                # Skip if we don't want to process this strand
                if process_forward_only and strand_symbol != "+":
                    continue
                if process_reverse_only and strand_symbol != "-":
                    continue

                high_conversion_count = position_data[pos][strand_symbol][
                    "high_conversion_count"
                ]
                insufficient_conversion_count = position_data[pos][strand_symbol][
                    "insufficient_conversion_count"
                ]
                low_quality_count = position_data[pos][strand_symbol][
                    "low_quality_count"
                ]

                # Get motif
                motif = extended_target_seq[
                    (pos - region_start - pad + extend_pad) : (
                        pos - region_start + pad + extend_pad + 1
                    )
                ]
                if strand_symbol == "-":
                    motif = reverse_complement(motif)

                # Calculate counts
                # u = unconverted (reference base), m = mutation (mutation base only), o = others
                u0 = low_quality_count[ref_base]
                m0 = low_quality_count[mut_base]
                o0 = low_quality_count.total() - u0 - m0
                u1 = high_conversion_count[ref_base]
                m1 = high_conversion_count[mut_base]
                o1 = high_conversion_count.total() - u1 - m1
                u2 = insufficient_conversion_count[ref_base]
                m2 = insufficient_conversion_count[mut_base]
                o2 = insufficient_conversion_count.total() - u2 - m2

                if u1 + m1 + u2 + m2 > 0:
                    site_info = [region_chrom, pos + 1, strand_symbol, motif]
                    if save_rest:
                        counts.append(site_info + [u0, u1, u2, m0, m1, m2, o0, o1, o2])
                    else:
                        counts.append(site_info + [u0, u1, u2, m0, m1, m2])

        if outfile:
            outfile.close()

        return {
            "worker_id": worker_id,
            "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
            "counts": counts,
            "success": True,
            "error": None,
            "reads": processed_reads,
            "total_reads": total_reads,
            "skipped_reads": total_skipped,
            "temp_bam_path": temp_bam_path if outfile else None,
            "timings": {"total": time.time() - overall_start},
        }
    except Exception as e:
        logger.error(f"Error in worker {worker_id}: {e}")
        return {
            "worker_id": worker_id,
            "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
            "counts": [],
            "success": False,
            "error": str(e),
            "reads": 0,
            "timings": {"total": time.time() - overall_start},
        }
    finally:
        # Readers are shared per process and closed by initializer cleanup
        pass


_temp_files_to_clean = set()


def cleanup_temp_files():
    for f in _temp_files_to_clean:
        try:
            os.remove(f)
        except OSError:
            pass
    _temp_files_to_clean.clear()


atexit.register(cleanup_temp_files)


def count_mutations(
    samfile: str,
    reffile: str,
    output_file: str | None = None,
    output_bam: str | None = None,
    ref_base: str = "A",
    mut_base: str = "G",
    ref_base2: str | None = None,
    mut_base2: str | None = None,
    bin_size: int = 10_000,
    threads: int | None = None,
    save_rest: bool = False,
    region: str | None = None,
    strand: str = "both",
    pad: int = 15,
    trim_start: int = 2,
    trim_end: int = 2,
    max_unc: int = 3,
    min_con: int = 1,
    max_sub: int = 1,
    min_baseq: int = 20,
    min_mapq: int = 0,
) -> bool:
    """
    Count mutations from BAM pileup data with parallel processing.

    This function is completely thread-safe and optimized for performance:
    1. Uses ProcessPoolExecutor for true parallelism
    2. Each worker opens its own file handles
    3. No shared state between workers
    4. Efficient memory usage with streaming
    5. Proper error handling and cleanup

    Args:
        samfile: Path to input BAM file
        reffile: Path to reference FASTA file
        output_file: Path to output TSV file (if None, prints to stdout)
        ref_base: Reference base to count (default: 'A')
        mut_base: Mutation base to count (default: 'G')
        bin_size: Size of genomic bins for processing (default: 10,000)
        threads: Number of parallel threads (default: auto-detect)
        save_rest: Whether to save additional statistics (default: False)
        region: Genomic region to process (e.g., 'chr1:1000000-2000000')
        pad: Motif half-window padding around site (default: 15)
        trim_start: Number of bases to trim from read 5' end (fragment orientation) (default: 2)
        trim_end: Number of bases to trim from read 3' end (fragment orientation) (default: 2)
        max_unc: Max unconverted threshold (Zf) to consider converted (default: 3)
        min_con: Min converted threshold (Yf) to consider converted (default: 1)
        max_sub: Max substitutions (NS) to consider mapped (default: 1)

    Returns:
        True if successful, False otherwise
    """
    start_time = time.time()

    # Convert bases to uppercase once at the beginning
    ref_base = ref_base.upper()
    mut_base = mut_base.upper()

    # Handle alternative mutation tagging
    tagging_enabled = ref_base2 and mut_base2
    temp_dir = None
    if tagging_enabled:
        ref_base2 = ref_base2.upper()
        mut_base2 = mut_base2.upper()
        temp_dir = tempfile.mkdtemp()
        print(
            f"🏷️ Alternative mutation tagging enabled. Temporary directory: {temp_dir}"
        )

    # Validate input files exist
    if not os.path.exists(samfile):
        print(f"❌ Input BAM file '{samfile}' does not exist!")
        return False
    if not os.path.exists(reffile):
        print(f"❌ Reference file '{reffile}' does not exist!")
        return False

    # Check and create BAM index if needed
    bam_index = samfile + ".bai"
    if not os.path.exists(bam_index):
        print(f"📇 BAM index not found. Creating index: {bam_index}")
        try:
            pysam.index(samfile)
            print("✅ BAM index created successfully")
        except Exception as e:
            print(f"❌ Failed to create BAM index: {e}")
            return False

    # Check and create FASTA index if needed
    fasta_index = reffile + ".fai"
    if not os.path.exists(fasta_index):
        print(f"📇 FASTA index not found. Creating index: {fasta_index}")
        try:
            pysam.faidx(reffile)
            print("✅ FASTA index created successfully")
        except Exception as e:
            print(f"❌ Failed to create FASTA index: {e}")
            return False

    # Set default threads
    if threads is None:
        threads = min(os.cpu_count() or 1, 8)

    # Validate base parameters
    valid_bases = {"A", "T", "G", "C"}
    if ref_base not in valid_bases:
        print(
            f"❌ Invalid reference base '{ref_base}'. Must be one of: {', '.join(valid_bases)}"
        )
        return False
    if mut_base not in valid_bases:
        print(
            f"❌ Invalid mutation base '{mut_base}'. Must be one of: {', '.join(valid_bases)}"
        )
        return False

    try:
        print("📖 Creating genomic bins...")

        # Read FASTA index file directly (much faster than opening full FASTA)
        ref_chrom_lengths = read_fasta_index(reffile)

        # Get BAM header
        samfile_open = pysam.AlignmentFile(samfile, "rb")
        bam_chroms = list(samfile_open.references)
        samfile_open.close()

        print(
            f"🔍 Filtering chromosomes: {len(bam_chroms)} in BAM, {len(ref_chrom_lengths)} in reference"
        )

        bin_list = []

        if region:
            # Parse region specification
            if ":" in region and "-" in region:
                chrom, pos_range = region.split(":")
                start, end = map(int, pos_range.split("-"))
                # Check if chromosome exists in BAM
                if chrom not in bam_chroms:
                    print(f"❌ Chromosome '{chrom}' not found in BAM file!")
                    print(f"Available chromosomes: {', '.join(sorted(bam_chroms))}")
                    return False
                # Convert from 1-based to 0-based coordinates for pysam
                bin_list = [(chrom, start - 1, end)]
            else:
                print(f"❌ Invalid region format: {region}. Use 'chr1:1000000-2000000'")
                return False
        else:
            # Process only chromosomes present in BAM; query lengths from reference
            valid_chroms = []
            missing_in_ref = []
            for chrom in bam_chroms:
                if chrom not in ref_chrom_lengths:
                    missing_in_ref.append(chrom)
                    continue
                chrom_length = ref_chrom_lengths[chrom]
                valid_chroms.append(chrom)
                bin_start = 0
                while bin_start < chrom_length:
                    bin_end = min(bin_start + bin_size, chrom_length)
                    bin_list.append((chrom, bin_start, bin_end))
                    bin_start += bin_size

            if missing_in_ref:
                print(
                    f"⚠️  {len(missing_in_ref)} BAM chromosomes not found in reference"
                )
                if len(missing_in_ref) <= 10:
                    print(f"   Missing: {', '.join(missing_in_ref)}")
                else:
                    print(
                        f"   Missing: {', '.join(missing_in_ref[:10])} ... and {len(missing_in_ref) - 10} more"
                    )

            print(f"✅ Processing {len(valid_chroms)} valid chromosomes")

        print(
            f"✅ Created {len(bin_list)} bins across {len({b[0] for b in bin_list})} chromosomes"
        )

        # Determine which strands to process
        process_both_strands = strand.lower() == "both"
        process_forward_only = strand.lower() == "forward"
        process_reverse_only = strand.lower() == "reverse"

        if not any([process_both_strands, process_forward_only, process_reverse_only]):
            print(
                f"❌ Invalid strand option '{strand}'. Must be 'both', 'forward', or 'reverse'"
            )
            return False

        # Use all regions for now (pre-filtering can be added later)
        print("🔍 Using all regions for processing...")
        total_skipped = 0
        filtered_bin_list = bin_list
        print(f"✅ Processing {len(filtered_bin_list)} regions")

        # Prepare worker arguments - now process both strands in one worker
        worker_args = []
        for i, (chrom, bin_start, bin_end) in enumerate(filtered_bin_list):
            temp_bam_path = (
                os.path.join(temp_dir, f"temp_{i}.bam") if tagging_enabled else None
            )
            worker_args.append(
                (
                    samfile,
                    chrom,
                    bin_start,
                    bin_end,
                    strand,  # Pass the strand option to worker
                    ref_base,
                    mut_base,
                    ref_base2,
                    mut_base2,
                    temp_bam_path,
                    save_rest,
                    pad,
                    trim_start,
                    trim_end,
                    max_unc,
                    min_con,
                    max_sub,
                    min_baseq,
                    min_mapq,
                    i,
                )
            )

        # Process all regions with optimal parallelism
        total_processed = 0
        total_counts = 0
        total_skipped = 0
        total_reads = 0
        all_results = []
        all_timings: list[dict[str, float]] = []
        worker_temp_bams = []

        # Write header immediately if outputting to stdout
        if output_file is None:
            headers = get_output_headers(save_rest)
            print("\t".join(headers))
            # Flush to ensure immediate output
            import sys

            sys.stdout.flush()

        print(f"🚀 Processing {len(worker_args)} regions with {threads} threads...")
        print(
            f"📊 Strand processing: {strand} ({'2 strands' if strand.lower() == 'both' else '1 strand'})"
        )

        with Progress(
            "[progress.description]{task.description}",
            "[progress.percentage]{task.percentage:>3.0f}%",
            "[cyan]{task.completed}/{task.total} regions",
            "[green]{task.fields[counts]:,} mutations",
            "[magenta]{task.fields[reads]:,} reads",
            TimeElapsedColumn(),
            TimeRemainingColumn(),
            expand=False,
        ) as progress:
            task = progress.add_task(
                "🔄 Processing regions...", total=len(worker_args), counts=0, reads=0
            )

            # Use ProcessPoolExecutor for optimal parallelism
            # This automatically handles load balancing and memory management
            with ProcessPoolExecutor(
                max_workers=threads,
                initializer=_init_worker,
                initargs=(samfile, reffile),
            ) as executor:
                # Submit all tasks at once for maximum parallelism
                future_to_args = {
                    executor.submit(parse_region_worker, args): args
                    for args in worker_args
                }

                # Process results as they complete and stream output
                for future in as_completed(future_to_args):
                    result = future.result()
                    total_processed += 1

                    if result["success"]:
                        # Check if region was skipped
                        if result.get("skipped", False):
                            total_skipped += 1
                        else:
                            total_counts += len(result["counts"])
                        # Accumulate reads and timing
                        total_reads += result.get("reads", 0)
                        if result.get("timings"):
                            all_timings.append(result["timings"])

                        # Stream results immediately if outputting to stdout
                        if output_file is None and result["counts"]:
                            for row in result["counts"]:
                                print("\t".join(map(str, row)))
                                sys.stdout.flush()
                        else:
                            # Collect for file output
                            all_results.extend(result["counts"])

                        if result.get("temp_bam_path"):
                            worker_temp_bams.append(result["temp_bam_path"])

                    else:
                        logger.warning(
                            f"Failed to process region {result['region']}: {result['error']}"
                        )

                    # Update progress
                    progress.update(
                        task, advance=1, counts=total_counts, reads=total_reads
                    )

        # Write results to file if specified
        if output_file:
            print("📝 Writing results to file...")
            # Sort results by chromosome and position
            all_results.sort(key=lambda x: (x[0], x[1]))
            write_output(all_results, output_file, save_rest)

        # Merge, sort, and index temporary BAM files if tagging was enabled
        if tagging_enabled and worker_temp_bams:
            print(f"Merging {len(worker_temp_bams)} temporary BAM files...")

            final_tagged_bam = (
                output_bam
                or tempfile.NamedTemporaryFile(delete=False, suffix=".bam").name
            )

            try:
                # Use pysam cat for robust merging
                pysam.cat("-o", final_tagged_bam, *worker_temp_bams)

                print("Sorting and indexing final BAM...")
                sorted_bam_path = final_tagged_bam + ".sorted"
                pysam.sort("-@", str(threads), "-o", sorted_bam_path, final_tagged_bam)
                shutil.move(sorted_bam_path, final_tagged_bam)
                pysam.index(final_tagged_bam)

                if not output_bam:
                    _temp_files_to_clean.add(final_tagged_bam)
                    _temp_files_to_clean.add(final_tagged_bam + ".bai")

                print(f"✅ Final tagged BAM created: {final_tagged_bam}")

            except Exception as e:
                print(f"❌ Failed during BAM processing: {e}")
            finally:
                if temp_dir:
                    shutil.rmtree(temp_dir)

        # Print summary
        elapsed_time = time.time() - start_time
        print("✅ Processing completed!")
        print(f"   Regions processed: {total_processed}")
        print(f"   Regions skipped (no reads): {total_skipped}")
        print(f"   Total mutations found: {total_counts}")
        print(f"   Time elapsed: {elapsed_time:.2f}s")
        print(f"   Processing rate: {total_processed / elapsed_time:.1f} regions/sec")
        if total_skipped > 0:
            print(f"   ⚡ Performance boost: Skipped {total_skipped} empty regions!")
        if all_timings:
            # Calculate average total time
            total_time = sum(t.get("total", 0) for t in all_timings)
            n = len(all_timings)
            avg_time_ms = (total_time * 1000 / n) if n > 0 else 0
            print(f"   ⏱️ Average per-window time: {avg_time_ms:.1f} ms")

        return True

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        print(f"❌ Processing failed: {e}")
        return False
    finally:
        cleanup_temp_files()
