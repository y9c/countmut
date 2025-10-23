#!/usr/bin/env python3
"""
Basic usage examples for CountMut.

This script demonstrates how to use CountMut for mutation counting
with different parameters and configurations.

Author: Ye Chang
Date: 2025-10-23
"""

import os
import tempfile

from countmut import count_mutations
from tests.test_utils import create_test_bam, create_test_fasta


def create_example_bam(output_path: str, num_reads: int = 1000):
    """Create an example BAM file for testing."""
    create_test_bam(output_path, num_reads=num_reads, chrom_length=100000, read_length=100)


def create_example_fasta(output_path: str, length: int = 100000):
    """Create an example FASTA file for testing."""
    create_test_fasta(output_path, length=length)


def example_basic_usage():
    """Example 1: Basic mutation counting."""
    print("🧬 Example 1: Basic mutation counting")
    print("=" * 50)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create example files
        bam_path = os.path.join(tmp_dir, "example.bam")
        fasta_path = os.path.join(tmp_dir, "example.fa")
        output_path = os.path.join(tmp_dir, "mutations.tsv")

        create_example_bam(bam_path, num_reads=1000)
        create_example_fasta(fasta_path, length=100000)

        # Count mutations
        success = count_mutations(
            samfile=bam_path,
            reffile=fasta_path,
            output_file=output_path,
            ref_base="A",
            mut_base="G",
            bin_size=5000,
            threads=4
        )

        if success:
            print("✅ Successfully counted mutations")
            print(f"📄 Results saved to: {output_path}")

            # Show first few lines of output
            with open(output_path) as f:
                lines = f.readlines()
                print(f"📊 Found {len(lines)-1} mutation sites")
                print("First few lines:")
                for line in lines[:5]:
                    print(f"  {line.strip()}")
        else:
            print("❌ Failed to count mutations")


def example_bisulfite_analysis():
    """Example 2: Bisulfite sequencing analysis."""
    print("\n🧬 Example 2: Bisulfite sequencing analysis")
    print("=" * 50)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create example files
        bam_path = os.path.join(tmp_dir, "bisulfite.bam")
        fasta_path = os.path.join(tmp_dir, "reference.fa")
        output_path = os.path.join(tmp_dir, "bisulfite_mutations.tsv")

        create_example_bam(bam_path, num_reads=2000)
        create_example_fasta(fasta_path, length=100000)

        # Count C→T conversions (unmethylated cytosines)
        success = count_mutations(
            samfile=bam_path,
            reffile=fasta_path,
            output_file=output_path,
            ref_base="C",
            mut_base="T",
            bin_size=10000,
            threads=4,
            save_rest=True  # Save additional statistics
        )

        if success:
            print("✅ Successfully analyzed bisulfite data")
            print(f"📄 Results saved to: {output_path}")
        else:
            print("❌ Failed to analyze bisulfite data")


def example_region_specific():
    """Example 3: Region-specific analysis."""
    print("\n🧬 Example 3: Region-specific analysis")
    print("=" * 50)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create example files
        bam_path = os.path.join(tmp_dir, "region.bam")
        fasta_path = os.path.join(tmp_dir, "reference.fa")
        output_path = os.path.join(tmp_dir, "region_mutations.tsv")

        create_example_bam(bam_path, num_reads=1500)
        create_example_fasta(fasta_path, length=100000)

        # Analyze specific region
        success = count_mutations(
            samfile=bam_path,
            reffile=fasta_path,
            output_file=output_path,
            ref_base="A",
            mut_base="G",
            region="chr1:20000-80000",  # Specific region
            bin_size=5000,
            threads=4
        )

        if success:
            print("✅ Successfully analyzed region chr1:20000-80000")
            print(f"📄 Results saved to: {output_path}")
        else:
            print("❌ Failed to analyze region")


def example_performance_comparison():
    """Example 4: Performance comparison."""
    print("\n🧬 Example 4: Performance comparison")
    print("=" * 50)

    with tempfile.TemporaryDirectory() as tmp_dir:
        # Create example files
        bam_path = os.path.join(tmp_dir, "perf.bam")
        fasta_path = os.path.join(tmp_dir, "reference.fa")

        create_example_bam(bam_path, num_reads=5000)
        create_example_fasta(fasta_path, length=200000)

        # Test different thread counts
        thread_counts = [1, 2, 4, 8]

        for threads in thread_counts:
            import time
            start_time = time.time()

            success = count_mutations(
                samfile=bam_path,
                reffile=fasta_path,
                output_file=os.path.join(tmp_dir, f"perf_{threads}.tsv"),
                ref_base="A",
                mut_base="G",
                bin_size=10000,
                threads=threads
            )

            elapsed_time = time.time() - start_time

            if success:
                print(f"✅ {threads} threads: {elapsed_time:.2f}s")
            else:
                print(f"❌ {threads} threads: Failed")


def main():
    """Run all examples."""
    print("🧬 CountMut Examples")
    print("=" * 50)

    try:
        example_basic_usage()
        example_bisulfite_analysis()
        example_region_specific()
        example_performance_comparison()

        print("\n🎉 All examples completed successfully!")

    except Exception as e:
        print(f"\n❌ Error running examples: {e}")
        import traceback
        traceback.print_exc()


if __name__ == "__main__":
    main()
