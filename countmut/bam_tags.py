#!/usr/bin/env python3
"""
BAM tag operations for CountMut.

This module provides functionality for reading and writing BAM tags,
specifically for tracking alternative mutations with Yc and Zc tags.

Author: Ye Chang
Date: 2025-01-27
"""

import os
import tempfile

import pysam


def calculate_alternative_mutations_in_region(
    read: pysam.AlignedSegment,
    ref_base2: str,
    mut_base2: str,
    region_start: int,
    region_end: int,
    target_seq: str,
    actual_strand: str,
    reverse_complement_func,
) -> tuple[int, int]:
    """
    Calculate alternative mutations (ref_base2 -> mut_base2) for a read within a region.

    This function counts alternative mutations by comparing the reference sequence
    at each position with the query sequence, handling both forward and reverse strands.

    Args:
        read: BAM read segment
        ref_base2: Alternative reference base to count
        mut_base2: Alternative mutation base to count
        region_start: Start position of the region (0-based)
        region_end: End position of the region (0-based)
        target_seq: Reference sequence for the region
        actual_strand: Biological strand ('+' or '-')
        reverse_complement_func: Function to reverse complement a sequence

    Returns:
        Tuple of (alt_ref_count, alt_mut_count) for alternative mutations
    """
    alt_ref_count = 0
    alt_mut_count = 0

    if not read.query_sequence:
        return alt_ref_count, alt_mut_count

    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if not (region_start <= ref_pos < region_end):
            continue

        ref_base_at_pos = target_seq[ref_pos - region_start].upper()
        query_base_at_pos = read.query_sequence[query_pos].upper()

        if actual_strand == "+":
            if ref_base_at_pos == ref_base2:
                if query_base_at_pos == ref_base2:
                    alt_ref_count += 1
                elif query_base_at_pos == mut_base2:
                    alt_mut_count += 1
        else:  # Reverse strand
            if ref_base_at_pos == reverse_complement_func(ref_base2):
                query_base_revcomp = reverse_complement_func(query_base_at_pos)
                if query_base_revcomp == ref_base2:
                    alt_ref_count += 1
                elif query_base_revcomp == mut_base2:
                    alt_mut_count += 1

    return alt_ref_count, alt_mut_count


def tag_read_with_alternative_mutations(
    read: pysam.AlignedSegment,
    alt_ref_count: int,
    alt_mut_count: int,
) -> pysam.AlignedSegment:
    """
    Tag a read with alternative mutation counts (Yc, Zc) and adjust NS tag.

    Assumes NS tag exists and is correct.

    Args:
        read: BAM read segment to tag
        alt_ref_count: Count of alternative reference bases (Zc tag)
        alt_mut_count: Count of alternative mutation bases (Yc tag)

    Returns:
        The tagged BAM read segment
    """
    read.set_tag("Yc", alt_mut_count, "i")
    read.set_tag("Zc", alt_ref_count, "i")

    # Assume NS tag exists and is correct
    ns_val = read.get_tag("NS")
    read.set_tag("NS", max(0, ns_val - alt_mut_count), "i")
    
    return read


def count_alternative_mutations(
    read: pysam.AlignedSegment,
    ref_base2: str,
    mut_base2: str,
    min_baseq: int = 20,
    trim_start: int = 2,
    trim_end: int = 2,
) -> tuple[int, int]:
    """
    Count alternative mutations (ref_base2 -> mut_base2) in a read.

    Args:
        read: BAM read segment
        ref_base2: Alternative reference base to count
        mut_base2: Alternative mutation base to count
        min_baseq: Minimum base quality threshold
        trim_start: Bases to trim from read start
        trim_end: Bases to trim from read end

    Returns:
        Tuple of (ref_count, mut_count) for alternative mutations
    """
    if not read.query_sequence:
        return 0, 0

    query_sequence = read.query_sequence.upper()
    query_qualities = read.query_qualities or []

    # Determine actual strand
    if read.is_paired:
        if read.is_read1:
            actual_strand = "+" if not read.is_reverse else "-"
        else:
            actual_strand = "+" if read.is_reverse else "-"
    else:
        actual_strand = "+" if not read.is_reverse else "-"

    # DNA complement mapping
    dna_complement = str.maketrans("ATGCNatgcn", "TACGNtacgn")

    ref_count = 0
    mut_count = 0

    # Process each aligned position
    for query_pos, ref_pos in read.get_aligned_pairs(matches_only=True):
        if query_pos is None or ref_pos is None:
            continue
        if query_pos >= len(query_sequence):
            continue

        # Check if position is internal (not in trimmed regions)
        if actual_strand == "+":
            is_internal = (
                query_pos >= trim_start and len(query_sequence) - query_pos > trim_end
            )
        else:
            is_internal = (
                query_pos >= trim_end and len(query_sequence) - query_pos > trim_start
            )

        if not is_internal:
            continue

        # Get query base and quality
        query_base = query_sequence[query_pos].upper()
        base_qual = (
            int(query_qualities[query_pos])
            if query_qualities and query_pos < len(query_qualities)
            else 0
        )

        if base_qual < min_baseq:
            continue

        # Apply reverse complement for reverse strand
        if actual_strand == "-":
            query_base = query_base.translate(dna_complement)

        # Count alternative mutations
        if query_base == ref_base2:
            ref_count += 1
        elif query_base == mut_base2:
            mut_count += 1

    return ref_count, mut_count


def update_read_tags(
    read: pysam.AlignedSegment,
    ref_base2: str,
    mut_base2: str,
    min_baseq: int = 20,
    trim_start: int = 2,
    trim_end: int = 2,
) -> pysam.AlignedSegment:
    """
    Update BAM read with alternative mutation tags (Yc, Zc) and adjust NS tag.

    Args:
        read: BAM read segment
        ref_base2: Alternative reference base
        mut_base2: Alternative mutation base
        min_baseq: Minimum base quality threshold
        trim_start: Bases to trim from read start
        trim_end: Bases to trim from read end

    Returns:
        Updated BAM read segment
    """
    # Count alternative mutations
    ref_count, mut_count = count_alternative_mutations(
        read, ref_base2, mut_base2, min_baseq, trim_start, trim_end
    )

    # Create a copy of the read
    updated_read = read.copy()

    # Add Yc and Zc tags
    updated_read.set_tag("Yc", ref_count, "i")  # Alternative reference count
    updated_read.set_tag("Zc", mut_count, "i")  # Alternative mutation count

    # Adjust NS tag by subtracting alternative mutations
    if read.has_tag("NS"):
        current_ns = read.get_tag("NS")
        # Subtract alternative mutations from NS (they're now tracked separately)
        new_ns = max(0, current_ns - mut_count)
        updated_read.set_tag("NS", new_ns, "i")

    return updated_read


def process_bam_with_alternative_mutations(
    input_bam: str,
    output_bam: str,
    ref_base2: str,
    mut_base2: str,
    threads: int = 4,
    min_baseq: int = 20,
    trim_start: int = 2,
    trim_end: int = 2,
    region: str | None = None,
) -> bool:
    """
    Process BAM file to add alternative mutation tags (Yc, Zc) and update NS tags.

    Args:
        input_bam: Path to input BAM file
        output_bam: Path to output BAM file
        ref_base2: Alternative reference base
        mut_base2: Alternative mutation base
        threads: Number of threads for processing
        min_baseq: Minimum base quality threshold
        trim_start: Bases to trim from read start
        trim_end: Bases to trim from read end
        region: Genomic region to process (optional)

    Returns:
        True if successful, False otherwise
    """
    try:
        # Open input BAM file
        with pysam.AlignmentFile(input_bam, "rb") as infile:
            # Create temporary BAM file for writing
            temp_bam = tempfile.NamedTemporaryFile(
                suffix=".bam", delete=False, mode="wb"
            )
            temp_bam_path = temp_bam.name
            temp_bam.close()

            # Create output BAM file with same header
            with pysam.AlignmentFile(
                temp_bam_path, "wb", header=infile.header
            ) as outfile:
                # Process reads
                reads_processed = 0
                reads_updated = 0

                # Determine which reads to process
                if region:
                    reads = infile.fetch(region=region)
                else:
                    reads = infile.fetch()

                for read in reads:
                    reads_processed += 1

                    # Skip unmapped, duplicate, or secondary reads
                    if read.is_unmapped or read.is_duplicate or read.is_secondary:
                        outfile.write(read)
                        continue

                    # Update read with alternative mutation tags
                    updated_read = update_read_tags(
                        read, ref_base2, mut_base2, min_baseq, trim_start, trim_end
                    )

                    outfile.write(updated_read)
                    reads_updated += 1

                    if reads_processed % 10000 == 0:
                        print(
                            f"Processed {reads_processed} reads, updated {reads_updated}"
                        )

        # Sort the temporary BAM file
        print(f"Sorting BAM file with {threads} threads...")
        pysam.sort("-@", str(threads), "-o", output_bam, temp_bam_path)

        # Index the output BAM file
        print("Indexing output BAM file...")
        pysam.index(output_bam)

        # Clean up temporary file
        os.unlink(temp_bam_path)

        print(
            f"Successfully processed {reads_processed} reads, updated {reads_updated}"
        )
        return True

    except Exception as e:
        print(f"Error processing BAM file: {e}")
        # Clean up temporary file if it exists
        if "temp_bam_path" in locals() and os.path.exists(temp_bam_path):
            os.unlink(temp_bam_path)
        return False


def get_alternative_mutation_stats(
    bam_file: str,
    ref_base2: str,
    mut_base2: str,
    region: str | None = None,
) -> dict[str, int]:
    """
    Get statistics about alternative mutations in a BAM file.

    Args:
        bam_file: Path to BAM file
        ref_base2: Alternative reference base
        mut_base2: Alternative mutation base
        region: Genomic region to analyze (optional)

    Returns:
        Dictionary with mutation statistics
    """
    stats = {
        "total_reads": 0,
        "reads_with_yc": 0,
        "reads_with_zc": 0,
        "total_yc": 0,
        "total_zc": 0,
        "reads_with_alternative_mutations": 0,
    }

    try:
        with pysam.AlignmentFile(bam_file, "rb") as infile:
            reads = infile.fetch(region=region) if region else infile.fetch()

            for read in reads:
                stats["total_reads"] += 1

                if read.has_tag("Yc"):
                    stats["reads_with_yc"] += 1
                    stats["total_yc"] += read.get_tag("Yc")

                if read.has_tag("Zc"):
                    stats["reads_with_zc"] += 1
                    stats["total_zc"] += read.get_tag("Zc")

                if read.has_tag("Yc") and read.has_tag("Zc"):
                    yc = read.get_tag("Yc")
                    zc = read.get_tag("Zc")
                    if yc > 0 or zc > 0:
                        stats["reads_with_alternative_mutations"] += 1

    except Exception as e:
        print(f"Error reading BAM file: {e}")

    return stats
