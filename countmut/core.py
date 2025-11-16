#!/usr/bin/env python3
"""
Core mutation counting logic for CountMut.

This module provides the main functionality for counting mutations from BAM pileup data,
including bisulfite conversion analysis and parallel processing capabilities.

Author: Ye Chang
Date: 2025-10-23
"""

import atexit
import glob
import logging
import os
import sys
import tempfile
import time
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import Any

import pysam
from rich.console import Console
from rich.progress import (
    Progress,
    SpinnerColumn,
    TimeElapsedColumn,
    TimeRemainingColumn,
)

from .bam_tags import (
    calculate_alternative_mutations_in_region,
    tag_read_with_alternative_mutations,
)
from .utils import get_output_headers, write_output

# Set up logger
logger = logging.getLogger(__name__)
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter("%(asctime)s - %(message)s", datefmt="%H:%M:%S")
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
    if read.is_paired:
        if read.is_read1:
            return "+" if not read.is_reverse else "-"
        # read2: reverse complemented indicates '+'
        return "+" if read.is_reverse else "-"
    # single-end
    return "+" if not read.is_reverse else "-"


# Global per-process file handles initialized once per worker process
_GLOBAL_SAM = None
_GLOBAL_REF = None


# New globals for worker-specific BAM writer
_WORKER_SHARD_PATH = None
_WORKER_WRITER = None
_WORKER_ID = None

# Special value to identify shutdown tasks
# _SHUTDOWN_SENTINEL = "__SHUTDOWN__"


def _worker_shutdown_task():
    """Special task to close worker files before process termination."""
    global _WORKER_WRITER
    try:
        if _WORKER_WRITER is not None:
            _WORKER_WRITER.close()
            _WORKER_WRITER = None
    except Exception:
        pass
    return {"shutdown": True}


def _warmup_worker():
    """A dummy worker task to force ProcessPoolExecutor workers to initialize."""
    return True


def _init_worker(samfile_path: str, reffile_path: str, shard_dir: str) -> None:
    """Open BAM/FASTA once per process to avoid repeated open/close overhead and enable BGZF threads."""
    global _GLOBAL_SAM, _GLOBAL_REF, _WORKER_SHARD_PATH, _WORKER_WRITER, _WORKER_ID
    # Only create reader once per worker process
    if _GLOBAL_SAM is None:
        _GLOBAL_SAM = pysam.AlignmentFile(samfile_path, "rb")
    try:
        # Enable multi-threaded BGZF decompression for faster reading
        _GLOBAL_SAM.set_threads(2)
    except Exception:
        pass
    _GLOBAL_REF = pysam.FastaFile(reffile_path)

    # Create a unique shard path for this worker
    pid = os.getpid()
    _WORKER_SHARD_PATH = os.path.join(shard_dir, f"shard_{pid}.bam")
    _WORKER_WRITER = None  # Lazily initialized
    _WORKER_ID = pid  # Store the worker ID (using PID)

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
        try:
            if _WORKER_WRITER is not None:
                _WORKER_WRITER.close()
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

    This function is designed to be completely thread-safe and optimized by:
    1. Using per-process global file handles initialized by the pool's initializer
    2. Not sharing any mutable global state between workers
    3. Returning results instead of writing to a shared file/queue
    4. Processing both strands in one pass for efficiency
    """
    try:
        overall_start = time.time()
        # Unpack arguments
        (
            region_chrom,
            region_start,
            region_end,
            strand_option,  # both/forward/reverse
            ref_base,
            mut_base,
            ref_base2,
            mut_base2,
            save_rest,
            pad,
            trim_start,
            trim_end,
            max_unc,
            min_con,
            max_sub,
            min_baseq,
            min_mapq,
        ) = args

        # Get worker-specific BAM writer
        global _WORKER_SHARD_PATH, _WORKER_WRITER, _GLOBAL_SAM, _WORKER_ID
        worker_id = _WORKER_ID  # Retrieve from global

        # Setup for optional BAM writing to the worker's shard file
        if _WORKER_SHARD_PATH and _WORKER_WRITER is None:
            _WORKER_WRITER = pysam.AlignmentFile(
                _WORKER_SHARD_PATH, "wb", header=_GLOBAL_SAM.header
            )
        outfile = _WORKER_WRITER

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
        skipped_unmapped = 0  # New: skipped due to unmapped
        skipped_duplicate = 0  # New: skipped due to duplicate
        skipped_secondary = 0  # New: skipped due to secondary
        skipped_mismatch_filter = 0  # New: skipped due to mismatch filter
        skipped_mapq_filter = 0  # New: skipped due to mapq filter
        skipped_conversion_filter = 0  # New: skipped due to conversion filter
        skipped_missing_tags = 0
        skipped_no_sequence = 0
        processed_reads = 0
        # Removed: reads_with_no_valid_bases = 0 # Counter no longer needed

        # Initialize is_skipped and total_skipped_reads for the worker result
        is_skipped = False
        total_skipped_reads = 0

        # Track best observation per (ref_pos, query_name) to avoid double counting overlapping mates
        # Value: (strand, base, qual, is_internal, passes_baseq_filter, passes_conversion_filter)
        best_obs: dict[
            tuple[int, str],
            tuple[str, str, int, bool, bool, bool],
        ] = {}

        reads_to_process = samfile.fetch(region_chrom, region_start, region_end)
        if outfile:
            # If we are writing a BAM, we must iterate over a copy,
            # as modifying tags can affect iteration over the original.
            reads_to_process = list(reads_to_process)

        for read in reads_to_process:
            total_reads += 1
            try:
                actual_strand = determine_actual_strand(read)
                # If tagging is enabled, count alternative mutations and update tags
                if ref_base2 and mut_base2:
                    alt_ref_count, alt_mut_count = (
                        calculate_alternative_mutations_in_region(
                            read,
                            ref_base2,
                            mut_base2,
                            region_start,
                            region_end,
                            target_seq,
                            actual_strand,
                            reverse_complement,
                        )
                    )

                    read = tag_read_with_alternative_mutations(
                        read, alt_ref_count, alt_mut_count
                    )

                # Skip if we don't want to process this strand
                if process_forward_only and actual_strand != "+":
                    skipped_wrong_strand += 1
                    if outfile:
                        outfile.write(read)  # Write unmodified read on error
                    continue
                if process_reverse_only and actual_strand != "-":
                    skipped_wrong_strand += 1
                    if outfile:
                        outfile.write(read)  # Write unmodified read on error
                    continue

                # Skip reads that are unmapped, duplicate, or secondary
                if read.is_unmapped:
                    skipped_unmapped += 1
                    if outfile:
                        outfile.write(read)
                    continue
                if read.is_duplicate:
                    skipped_duplicate += 1
                    if outfile:
                        outfile.write(read)
                    continue
                if read.is_secondary:
                    skipped_secondary += 1
                    if outfile:
                        outfile.write(read)
                    continue

                # Get read properties (assume NS, Zf, and Yf tags exist and are correct)
                # Check mapping quality (read-level filter)
                if read.mapping_quality < min_mapq:
                    skipped_mapq_filter += 1
                    if outfile:
                        outfile.write(read)  # Write unmodified read
                    continue

                # Check mismatch filter (read-level filter)
                # Assume NS tag exists and is correct
                ns = read.get_tag("NS")
                if ns > max_sub:
                    skipped_mismatch_filter += 1
                    if outfile:
                        outfile.write(read)
                    continue

                # Check conversion filter (read-level filter)
                # Assume Zf and Yf tags exist and are correct
                zf = read.get_tag("Zf")
                yf = read.get_tag("Yf")
                if not ((zf <= max_unc) and (yf >= min_con)):
                    skipped_conversion_filter += 1
                    if outfile:
                        outfile.write(read)
                    continue

                # Process each position in the read
                query_sequence = read.query_sequence
                if not query_sequence:
                    skipped_no_sequence += 1
                    if outfile:
                        outfile.write(read)  # Write unmodified read on error
                    continue
                query_qualities = read.query_qualities or []

                # Mark this read as successfully processed (considered for base-level analysis)
                processed_reads += 1

                # Iterate via reference positions (fast path) and filter
                for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
                    if query_pos is None or ref_pos is None:
                        continue
                    if ref_pos not in target_sites_set:
                        continue
                    if query_pos >= len(query_sequence):
                        continue

                    # Check if position is internal (not in trimmed regions)
                    # Trim based on actual strand orientation (fragment 5' to 3')
                    if actual_strand == "+":
                        is_internal = (
                            query_pos >= trim_start
                            and len(query_sequence) - query_pos > trim_end
                        )
                    else:
                        is_internal = (
                            query_pos >= trim_end
                            and len(query_sequence) - query_pos > trim_start
                        )

                    if not is_internal:
                        continue

                    query_base = query_sequence[query_pos].upper()
                    base_qual = (
                        int(query_qualities[query_pos])
                        if query_qualities and query_pos < len(query_qualities)
                        else 0
                    )

                    # Check base quality
                    passes_baseq_filter = base_qual >= min_baseq
                    if not passes_baseq_filter:
                        continue

                    if actual_strand == "-":
                        query_base = query_base.translate(DNA_COMPLEMENT)
                    key = (ref_pos, read.query_name)
                    prev = best_obs.get(key)

                    # If any of the remaining filters fail, this base is considered dropped (low_quality)
                    if not (
                        is_internal
                        and passes_baseq_filter  # Only these two are now base-level filters
                    ):
                        # Even if dropped, we still keep the best observation to correctly populate low_quality_count
                        if (prev is None) or (base_qual > prev[2]):
                            best_obs[key] = (
                                actual_strand,
                                query_base,
                                base_qual,
                                False,  # is_internal (false if dropped)
                                False,  # passes_baseq_filter (false if dropped)
                                False,  # passes_conversion_filter (false if dropped at base level)
                            )
                        continue  # Skip further processing for this base as it's low quality

                    # If all quality filters pass, proceed to store the best observation
                    if (prev is None) or (base_qual > prev[2]):
                        best_obs[key] = (
                            actual_strand,
                            query_base,
                            base_qual,
                            bool(is_internal),
                            bool(passes_baseq_filter),
                            True,  # passes_conversion_filter (always true at this point, since read-level filter passed)
                        )
                    # Removed: read_contributes_to_counts = True # Flag no longer needed

            except (KeyError, AttributeError) as e:
                # Skip reads with missing tags or invalid data
                logger.debug(f"Skipping read due to missing tag: {e}")
                if outfile:
                    outfile.write(read)  # Write unmodified read on error
                continue

            # Removed: After processing all positions for a read, update processed_reads and total_skipped_reads
            # Removed: if read_contributes_to_counts:
            # Removed: processed_reads += 1
            # Removed: elif not read_contributes_to_counts: # If read didn't contribute valid bases
            # Removed: # This means it passed initial filters but all its bases were dropped.
            # Removed: reads_with_no_valid_bases += 1

            if outfile:
                outfile.write(
                    read
                )  # Write the (potentially modified) read to the shard

        # Apply best observations to position counts (deduplicated across overlapping mates)
        for (ref_pos, _qname), (
            strand_symbol,
            query_base,
            _q,
            is_internal,  # Keep this, used below
            passes_baseq_filter,  # Keep this, used below
            passes_conversion_filter,  # Keep this, used below for counting
        ) in best_obs.items():
            # Skip positions that are not in our target sites (safety)
            if ref_pos not in position_data:
                continue

            # Count in low_quality if ANY quality filter fails
            if not (is_internal and passes_baseq_filter):
                position_data[ref_pos][strand_symbol]["low_quality_count"][
                    query_base
                ] += 1
            # Count in high_conversion/insufficient_conversion only if ALL base-level quality filters pass
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

        # Calculate total skipped reads (including those with no valid bases)
        total_skipped_reads = (
            skipped_wrong_strand
            + skipped_unmapped
            + skipped_duplicate
            + skipped_secondary
            + skipped_mismatch_filter  # New: Include reads skipped by mismatch filter
            + skipped_mapq_filter  # New: Include reads skipped by mapq filter
            + skipped_conversion_filter  # New: Include reads skipped by conversion filter
            + skipped_missing_tags
            + skipped_no_sequence
            # Removed: + reads_with_no_valid_bases # Include reads that passed initial filters but had no valid bases
        )
        # A region is considered skipped if no reads were successfully processed through all filters
        if processed_reads == 0:
            is_skipped = True

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

        # Return results to the main process
        return {
            "worker_id": worker_id,
            "success": True,
            "error": None,
            "region": f"{region_chrom}:{region_start}-{region_end}:{strand_option}",
            "counts": counts,  # Use the 'counts' list generated in the worker
            "reads": processed_reads,
            "total_reads": total_reads,  # Total reads for this region
            "skipped": is_skipped,  # Correctly reflect if the region was skipped
            "skipped_reads": total_skipped_reads,  # Total reads skipped in this region
            "skipped_wrong_strand": skipped_wrong_strand,
            "skipped_unmapped": skipped_unmapped,  # New: individual skipped count
            "skipped_duplicate": skipped_duplicate,  # New: individual skipped count
            "skipped_secondary": skipped_secondary,  # New: individual skipped count
            "skipped_mismatch_filter": skipped_mismatch_filter,  # New: individual skipped count
            "skipped_mapq_filter": skipped_mapq_filter,  # New: individual skipped count
            "skipped_conversion_filter": skipped_conversion_filter,  # New: individual skipped count
            "skipped_missing_tags": skipped_missing_tags,
            "skipped_no_sequence": skipped_no_sequence,
            # Removed: "reads_with_no_valid_bases": reads_with_no_valid_bases, # New: Detailed skipped count
            "temp_bam_path": _WORKER_SHARD_PATH,  # Path to the temporary BAM file for this worker
            "timings": {
                "total": time.time() - overall_start
            },  # Use the already constructed timings dictionary
        }

    except Exception as e:
        # Ensure worker_id is always available for logging
        worker_id_for_log = _WORKER_ID if _WORKER_ID is not None else "unknown"
        logger.error(f"Error in worker {worker_id_for_log}: {e}")
        logger.error(f"❌ Processing failed: {e}")
        return {
            "worker_id": worker_id_for_log,
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


def _split_region_if_needed(
    samfile: pysam.AlignmentFile | str,
    chrom: str,
    region_start: int,
    region_end: int,
    max_reads_per_chunk: int = 100_000,
    min_chunk_size: int = 1_000,
) -> list[tuple[str, int, int]]:
    """
    Split a region into smaller chunks if it has too many reads.
    
    Args:
        samfile: Either a pysam.AlignmentFile object (for reuse) or path to BAM file
        chrom: Chromosome name
        region_start: Start position (0-based)
        region_end: End position (0-based)
        max_reads_per_chunk: Maximum reads per chunk before splitting
        min_chunk_size: Minimum chunk size in base pairs
        
    Returns:
        List of (chrom, start, end) tuples for chunks
    """
    try:
        # Handle both file path and open file object
        if isinstance(samfile, str):
            samfile_open = pysam.AlignmentFile(samfile, "rb")
            read_count = samfile_open.count(chrom, region_start, region_end)
            samfile_open.close()
        else:
            # Reuse the provided file handle
            read_count = samfile.count(chrom, region_start, region_end)
        
        # If read count is below threshold, return original region
        if read_count <= max_reads_per_chunk:
            return [(chrom, region_start, region_end)]
        
        # Calculate number of chunks needed
        num_chunks = (read_count + max_reads_per_chunk - 1) // max_reads_per_chunk
        region_size = region_end - region_start
        chunk_size = max(region_size // num_chunks, min_chunk_size)
        
        # Split region into chunks
        chunks = []
        chunk_start = region_start
        while chunk_start < region_end:
            chunk_end = min(chunk_start + chunk_size, region_end)
            chunks.append((chrom, chunk_start, chunk_end))
            chunk_start = chunk_end
        
        logger.info(
            f"📦 Split region {chrom}:{region_start}-{region_end} "
            f"({read_count:,} reads) into {len(chunks)} chunks"
        )
        return chunks
    except Exception as e:
        logger.warning(f"Failed to count reads for {chrom}:{region_start}-{region_end}: {e}")
        # Return original region if counting fails
        return [(chrom, region_start, region_end)]


def count_mutations(
    samfile: str,
    reference: str,
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
    force: bool = False,
    strand: str = "both",
    pad: int = 15,
    trim_start: int = 2,
    trim_end: int = 2,
    max_unc: int = 3,
    min_con: int = 1,
    max_sub: int = 1,
    min_baseq: int = 20,
    min_mapq: int = 0,
    verbose: bool = False,
    max_reads_per_chunk: int = 100_000,
) -> dict:
    """
    Count mutations from BAM pileup data with parallel processing.

    This function is completely thread-safe and optimized for performance:
    1. Uses ProcessPoolExecutor for true parallelism
    2. Each worker opens its own file handles
    3. No shared state between workers
    4. Efficient memory usage with streaming
    5. Proper error handling and cleanup
    6. Automatically chunks regions with too many reads for parallel processing

    Args:
        samfile: Path to input BAM file
        reference: Path to reference FASTA file
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
        max_reads_per_chunk: Maximum reads per chunk before splitting region (default: 100,000)

    Returns:
        True if successful, False otherwise
    """
    start_time = time.time()

    # Configure logging level based on verbose flag
    if verbose:
        logger.setLevel(logging.DEBUG)  # Changed to DEBUG
    else:
        logger.setLevel(logging.WARNING)

    # Convert bases to uppercase once at the beginning
    ref_base = ref_base.upper()
    mut_base = mut_base.upper()

    # Handle alternative mutation tagging
    tagging_enabled = ref_base2 and mut_base2
    if tagging_enabled:
        ref_base2 = ref_base2.upper()
        mut_base2 = mut_base2.upper()

    # All validation and setup steps prior to main processing loop
    try:
        # Check and create BAM index if needed
        bam_index = samfile + ".bai"
        if not os.path.exists(bam_index):
            logger.info(f"📇 BAM index not found. Creating index: {bam_index}")
            try:
                pysam.index(samfile)
                logger.info("✅ BAM index created successfully")
            except Exception as e:
                logger.error(f"❌ Failed to create BAM index: {e}")
                return False

        # Check and create FASTA index if needed
        fasta_index = reference + ".fai"
        if not os.path.exists(fasta_index):
            logger.info(f"📇 FASTA index not found. Creating index: {fasta_index}")
            try:
                pysam.faidx(reference)
                logger.info("✅ FASTA index created successfully")
            except Exception as e:
                logger.error(f"❌ Failed to create FASTA index: {e}")
                return False

        # Set default threads
        if threads is None:
            threads = min(os.cpu_count() or 1, 8)

        # Validate base parameters
        valid_bases = {"A", "T", "G", "C"}
        if ref_base not in valid_bases:
            logger.error(
                f"❌ Invalid reference base '{ref_base}'. Must be one of: {', '.join(valid_bases)}"
            )
            return False
        if mut_base not in valid_bases:
            logger.error(
                f"❌ Invalid mutation base '{mut_base}'. Must be one of: {', '.join(valid_bases)}"
            )
            return False

        logger.info("🔍 Validating BAM file...")
        # Validate input files exist
        if not os.path.exists(samfile):
            logger.error(f"❌ Input BAM file '{samfile}' does not exist!")
            return False
        if not os.path.exists(reference):
            logger.error(f"❌ Reference file '{reference}' does not exist!")
            return False
        logger.info("✅ BAM file validation passed")

        logger.info("🚀 Starting mutation counting...")

        with tempfile.TemporaryDirectory(prefix="countmut_") as temp_dir:
            if tagging_enabled:
                logger.info(
                    f"🏷️ Alternative mutation tagging enabled. Temporary directory: {temp_dir}"
                )

            logger.info("📖 Creating genomic bins...")

            # Read FASTA index file directly (much faster than opening full FASTA)
            ref_chrom_lengths = read_fasta_index(reference)

            # Get BAM header
            samfile_open = pysam.AlignmentFile(samfile, "rb")
            bam_chroms = list(samfile_open.references)
            samfile_open.close()

            logger.info(
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
                        logger.error(f"❌ Chromosome '{chrom}' not found in BAM file!")
                        logger.error(
                            f"Available chromosomes: {', '.join(sorted(bam_chroms))}"
                        )
                        return False
                    # Convert from 1-based to 0-based coordinates for pysam
                    bin_list = [(chrom, start - 1, end)]
                else:
                    logger.error(
                        f"❌ Invalid region format: {region}. Use 'chr1:1000000-2000000'"
                    )
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
                    logger.warning(
                        f"⚠️  {len(missing_in_ref)} BAM chromosomes not found in reference"
                    )
                    if len(missing_in_ref) <= 10:
                        logger.warning(f"   Missing: {', '.join(missing_in_ref)}")
                    else:
                        logger.warning(
                            f"   Missing: {', '.join(missing_in_ref[:10])} ... and {len(missing_in_ref) - 10} more"
                        )

                logger.info(f"✅ Processing {len(valid_chroms)} valid chromosomes")

            logger.info(
                f"✅ Created {len(bin_list)} bins across {len({b[0] for b in bin_list})} chromosomes"
            )

            # Determine which strands to process
            process_both_strands = strand.lower() == "both"
            process_forward_only = strand.lower() == "forward"
            process_reverse_only = strand.lower() == "reverse"

            if not any(
                [process_both_strands, process_forward_only, process_reverse_only]
            ):
                logger.error(
                    f"❌ Invalid strand option '{strand}'. Must be 'both', 'forward', or 'reverse'"
                )
                return False

            # Use all regions for now (pre-filtering can be added later)
            logger.info("🔍 Using all regions for processing...")
            total_skipped = 0
            filtered_bin_list = bin_list
            logger.info(f"✅ Processing {len(filtered_bin_list)} regions")
            
            # Check for regions with too many reads and split them into chunks
            logger.info("🔍 Checking for regions with high read density...")
            chunked_bin_list = []
            regions_split = 0
            # Open BAM file once for all read count checks
            samfile_for_counting = pysam.AlignmentFile(samfile, "rb")
            try:
                for chrom, bin_start, bin_end in filtered_bin_list:
                    chunks = _split_region_if_needed(
                        samfile_for_counting, chrom, bin_start, bin_end, max_reads_per_chunk
                    )
                    if len(chunks) > 1:
                        regions_split += 1
                    chunked_bin_list.extend(chunks)
            finally:
                samfile_for_counting.close()
            
            if regions_split > 0:
                logger.info(
                    f"📦 Split {regions_split} regions with high read density "
                    f"into {len(chunked_bin_list)} total chunks"
                )
            else:
                logger.info("✅ No regions needed splitting")
            
            # Use chunked bin list for processing
            filtered_bin_list = chunked_bin_list
            logger.info(f"✅ Processing {len(filtered_bin_list)} regions/chunks")

            # Prepare worker arguments - now process both strands in one worker
            worker_args = []
            for i, (chrom, bin_start, bin_end) in enumerate(filtered_bin_list):
                worker_args.append(
                    (
                        chrom,
                        bin_start,
                        bin_end,
                        strand,  # Pass the strand option to worker
                        ref_base,
                        mut_base,
                        ref_base2,
                        mut_base2,
                        save_rest,
                        pad,
                        trim_start,
                        trim_end,
                        max_unc,
                        min_con,
                        max_sub,
                        min_baseq,
                        min_mapq,
                    )
                )

            # Initialize counters and result lists
            total_processed = 0
            total_counts = 0
            total_skipped = 0  # Regions skipped (no reads)
            total_raw_reads_all_workers = (
                0  # New: Accumulate raw reads fetched from all workers
            )
            total_reads_processed_all_workers = (
                0  # Accumulate processed reads from all workers
            )
            total_reads_skipped_all_workers = (
                0  # Accumulate skipped reads from all workers
            )
            # Detailed skipped read counts across all workers
            total_skipped_wrong_strand_agg = 0
            total_skipped_unmapped_agg = 0  # New: Aggregated unmapped reads
            total_skipped_duplicate_agg = 0  # New: Aggregated duplicate reads
            total_skipped_secondary_agg = 0  # New: Aggregated secondary reads
            total_skipped_mismatch_filter_agg = (
                0  # New: Aggregated mismatch filter skipped reads
            )
            total_skipped_mapq_filter_agg = (
                0  # New: Aggregated mapq filter skipped reads
            )
            total_skipped_conversion_filter_agg = (
                0  # New: Aggregated conversion filter skipped reads
            )
            total_skipped_missing_tags_agg = 0
            total_skipped_no_sequence_agg = 0
            # Removed: total_reads_with_no_valid_bases_agg = 0 # Accumulator no longer needed
            all_results = []
            all_timings: list[dict[str, float]] = []

            # Write header immediately if outputting to stdout
            if output_file is None:
                headers = get_output_headers(save_rest)
                sys.stdout.write(
                    "\t".join(headers) + "\nNo valid bases"
                )  # Use sys.stdout.write for raw output
                sys.stdout.flush()

            logger.info(
                f"🚀 Processing {len(worker_args)} regions with {threads} threads..."
            )
            logger.info(
                f"📊 Strand processing: {strand} ({'2 strands' if strand.lower() == 'both' else '1 strand'})"
            )

            # Use ProcessPoolExecutor for optimal parallelism
            # This automatically handles load balancing and memory management
            with ProcessPoolExecutor(
                max_workers=threads,
                initializer=_init_worker,
                initargs=(samfile, reference, temp_dir),
            ) as executor:
                # Initialize rich console once for consistent output
                console = Console()

                # Use a simple Progress bar for warmup phase
                with (
                    Progress(
                        SpinnerColumn(),
                        "[progress.description]{task.description}",
                        "[cyan]{task.completed}/{task.total} workers warmed up",
                        TimeElapsedColumn(),
                        expand=False,
                        console=console,  # Added console for transient behavior
                        transient=True,  # Added transient=True to hide the bar after completion
                    ) as warmup_progress
                ):
                    warmup_task = warmup_progress.add_task(
                        "🚀 Warming up workers...", total=threads
                    )

                    warmup_futures = [
                        executor.submit(_warmup_worker) for _ in range(threads)
                    ]
                    for future in as_completed(warmup_futures):
                        future.result()  # Wait for each worker to initialize
                        warmup_progress.update(warmup_task, advance=1)
                # The warmup progress bar will now hide automatically.
                logger.info(
                    "✅ Workers warmed up. Submitting main tasks..."
                )  # This message will appear after the bar hides

                # Use a Live context to keep the progress bar at the bottom
                with Progress(
                    SpinnerColumn(),
                    "[progress.description]{task.description}",
                    "[progress.percentage]{task.percentage:>3.0f}%",
                    "[cyan]{task.completed}/{task.total} regions",
                    "[green]{task.fields[counts]:,} mutations",
                    "[magenta]{task.fields[reads]:,} reads",
                    TimeElapsedColumn(),
                    TimeRemainingColumn(),
                    expand=False,
                    console=console,  # Ensure console is passed for transient behavior
                    transient=True,  # Make the main progress bar transient
                ) as progress:
                    task = progress.add_task(
                        "🔄 Processing regions...",
                        total=len(worker_args),
                        counts=0,
                        reads=0,
                    )

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
                            if result.get("skipped", False):
                                total_skipped += 1  # Count of regions skipped
                            else:
                                total_counts += len(result["counts"])

                            # Accumulate raw reads, processed reads and skipped reads from worker results
                            total_raw_reads_all_workers += result.get("total_reads", 0)
                            total_reads_processed_all_workers += result.get("reads", 0)
                            total_reads_skipped_all_workers += result.get(
                                "skipped_reads", 0
                            )

                            # Accumulate detailed skipped read counts
                            total_skipped_wrong_strand_agg += result.get(
                                "skipped_wrong_strand", 0
                            )
                            total_skipped_unmapped_agg += result.get(
                                "skipped_unmapped", 0
                            )  # New: Accumulate unmapped
                            total_skipped_duplicate_agg += result.get(
                                "skipped_duplicate", 0
                            )  # New: Accumulate duplicate
                            total_skipped_secondary_agg += result.get(
                                "skipped_secondary", 0
                            )  # New: Accumulate secondary
                            total_skipped_mismatch_filter_agg += result.get(
                                "skipped_mismatch_filter", 0
                            )  # New: Accumulate mismatch filter skipped
                            total_skipped_mapq_filter_agg += result.get(
                                "skipped_mapq_filter", 0
                            )  # New: Accumulate mapq filter skipped
                            total_skipped_conversion_filter_agg += result.get(
                                "skipped_conversion_filter", 0
                            )  # New: Accumulate conversion filter skipped
                            total_skipped_missing_tags_agg += result.get(
                                "skipped_missing_tags", 0
                            )
                            total_skipped_no_sequence_agg += result.get(
                                "skipped_no_sequence", 0
                            )
                            # Removed: total_reads_with_no_valid_bases_agg += result.get("reads_with_no_valid_bases", 0) # Accumulate new counter

                            if result.get("timings"):
                                all_timings.append(result["timings"])

                            # Stream results immediately if outputting to stdout
                            if output_file is None and result["counts"]:
                                for row in result["counts"]:
                                    sys.stdout.write("\t".join(map(str, row)) + "\n")
                                    sys.stdout.flush()
                            else:
                                # Collect for file output
                                all_results.extend(result["counts"])

                        else:
                            logger.warning(
                                f"Failed to process region {result['region']}: {result['error']}"
                            )

                        # Update progress with total mutations and reads processed
                        progress.update(
                            task,
                            advance=1,
                            counts=total_counts,
                            reads=total_reads_processed_all_workers,
                        )

                    # Explicitly shut down the executor after all results are collected
                    executor.shutdown(wait=True)

                logger.info(
                    f"📊 Processed {total_processed} regions, {total_counts:,} sites, "
                    f"{total_reads_processed_all_workers:,} reads."
                )

            # Write results to file if specified
            if output_file:
                logger.info("📝 Writing results to file...")
                # Sort results by chromosome and position
                all_results.sort(key=lambda x: (x[0], x[1]))
                write_output(all_results, output_file, save_rest)

            # Merge, sort, and index temporary BAM files if tagging was enabled
            _final_bam_processing(
                tagging_enabled, output_bam, threads, temp_dir, samfile
            )

            # Print summary
            elapsed_time = time.time() - start_time
            logger.info(f"✅ Processing completed! (⏱️ {elapsed_time:.2f}s)")

            # Return key statistics for CLI to display in a final panel
            return {
                "success": True,
                "total_processed_regions": total_processed,
                "total_skipped_regions": total_skipped,
                "total_raw_reads": total_raw_reads_all_workers,
                "total_reads_processed": total_reads_processed_all_workers,
                "total_reads_skipped": total_reads_skipped_all_workers,
                "total_mutations_found": total_counts,
                "elapsed_time": elapsed_time,
                "total_skipped_wrong_strand_agg": total_skipped_wrong_strand_agg,
                "total_skipped_unmapped_agg": total_skipped_unmapped_agg,  # New: return aggregated unmapped
                "total_skipped_duplicate_agg": total_skipped_duplicate_agg,  # New: return aggregated duplicate
                "total_skipped_secondary_agg": total_skipped_secondary_agg,  # New: return aggregated secondary
                "total_skipped_mismatch_filter_agg": total_skipped_mismatch_filter_agg,  # New: return aggregated mismatch filter skipped
                "total_skipped_mapq_filter_agg": total_skipped_mapq_filter_agg,  # New: return aggregated mapq filter skipped
                "total_skipped_conversion_filter_agg": total_skipped_conversion_filter_agg,  # New: return aggregated conversion filter skipped
                "total_skipped_missing_tags_agg": total_skipped_missing_tags_agg,
                "total_skipped_no_sequence_agg": total_skipped_no_sequence_agg,
                # Removed: "total_reads_with_no_valid_bases_agg": total_reads_with_no_valid_bases_agg, # Add to return dict
            }

    except Exception as e:
        logger.error(f"Error during processing: {e}")
        logger.error(f"❌ Processing failed: {e}")
        return {"success": False, "error": str(e)}
    finally:
        # Readers are shared per process and closed by initializer cleanup
        pass


def _final_bam_processing(
    tagging_enabled: bool,
    output_bam: str | None,
    threads: int,
    temp_dir: str,
    samfile: str,
) -> None:
    """
    Handles merging, sorting, and indexing of temporary BAM files.
    """
    if tagging_enabled and output_bam:
        # Glob all shard files from the temporary directory
        shard_files = glob.glob(os.path.join(temp_dir, "shard_*.bam"))
        if not shard_files:
            logger.warning("No temporary BAM files found to merge.")
            return

        logger.info(f"Merging {len(shard_files)} temporary BAM files...")
        final_tagged_bam = output_bam

        try:
            # Use pysam cat for robust merging
            pysam.cat("-o", final_tagged_bam, *shard_files)

            logger.info("Sorting and indexing final BAM...")
            sorted_bam_path = tempfile.NamedTemporaryFile(
                delete=False, suffix=".bam", dir=os.path.dirname(output_bam)
            ).name
            pysam.sort("-@", str(threads), "-o", sorted_bam_path, final_tagged_bam)
            os.replace(sorted_bam_path, final_tagged_bam)
            pysam.index(final_tagged_bam)

            logger.info(f"✅ Final tagged BAM created: {final_tagged_bam}")

        except Exception as e:
            logger.error(f"❌ Failed during final BAM processing: {e}")
        finally:
            # The temporary directory is automatically cleaned up by the context manager
            pass
