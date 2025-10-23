#!/usr/bin/env python3
"""
Create a subset BAM file with only rRNA chromosomes for testing.
"""


import pysam


def create_subset_bam(input_bam, output_bam, target_chroms):
    """Create a subset BAM file with only specified chromosomes."""

    print(f"📖 Reading BAM file: {input_bam}")

    # Open input BAM file
    with pysam.AlignmentFile(input_bam, "rb") as infile:
        # Create new header with only target chromosomes
        header_dict = infile.header.to_dict()
        new_sq = []

        # Add only target chromosomes to header
        for chrom in target_chroms:
            # Find the chromosome in original header
            for sq in infile.header['SQ']:
                if sq['SN'] == chrom:
                    new_sq.append(sq)
                    print(f"✅ Added chromosome: {sq['SN']} (length: {sq['LN']})")
                    break
            else:
                print(f"⚠️  Warning: Chromosome {chrom} not found in original BAM")

        # Update header with filtered chromosomes
        header_dict['SQ'] = new_sq
        new_header = pysam.AlignmentHeader.from_dict(header_dict)

        # Create output BAM file
        with pysam.AlignmentFile(output_bam, "wb", header=new_header) as outfile:
            total_reads = 0
            kept_reads = 0

            print("🔍 Processing reads...")

            for read in infile.fetch():
                total_reads += 1

                # Check if read is from target chromosome
                if read.reference_name in target_chroms:
                    outfile.write(read)
                    kept_reads += 1

                # Progress indicator
                if total_reads % 10000 == 0:
                    print(f"   Processed {total_reads:,} reads, kept {kept_reads:,}")

            print("✅ Completed!")
            print(f"   Total reads processed: {total_reads:,}")
            print(f"   Reads kept: {kept_reads:,}")
            print(f"   Kept ratio: {kept_reads/total_reads*100:.1f}%")

    # Create index for the new BAM file
    print("📝 Creating BAM index...")
    pysam.index(output_bam)
    print(f"✅ Index created: {output_bam}.bai")

if __name__ == "__main__":
    input_bam = "/home/yec/Desktop/S818.genes.bam"
    output_bam = "/home/yec/Desktop/test_rRNA.bam"

    target_chroms = [
        "rRNA-Hsa-nucleus_locus",
        "rRNA-Hsa-mitochrondria_locus",
        "rRNA-Hsa-5S"
    ]

    print("🧬 Creating rRNA subset BAM file")
    print("=" * 50)
    print(f"Input: {input_bam}")
    print(f"Output: {output_bam}")
    print(f"Target chromosomes: {', '.join(target_chroms)}")
    print()

    create_subset_bam(input_bam, output_bam, target_chroms)

    print(f"\n🎉 Subset BAM file created: {output_bam}")
    print(f"   Index file: {output_bam}.bai")
