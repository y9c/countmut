"""Tests for core mutation counting functionality."""

import tempfile
from pathlib import Path

import pysam
import pytest

from countmut.core import (
    count_mutations,
    determine_actual_strand,
    get_motif,
    read_fasta_index,
    reverse_complement,
)


@pytest.fixture
def test_data_dir():
    """Fixture to provide the path to test data."""
    # Assuming test.bam and test.fa are directly in ~/Desktop for now
    return Path("/home/yec/Desktop/")


class TestUtilityFunctions:
    """Test core utility functions."""

    def test_reverse_complement(self):
        """Test reverse complement function."""
        assert reverse_complement("ATGC") == "GCAT"
        assert reverse_complement("atgc") == "gcat"
        assert reverse_complement("ATGCN") == "NGCAT"
        assert reverse_complement("") == ""
        assert reverse_complement("AAAAA") == "TTTTT"
        assert reverse_complement("TTTTT") == "AAAAA"
        assert reverse_complement("GGGGG") == "CCCCC"
        assert reverse_complement("CCCCC") == "GGGGG"

    def test_get_motif(self):
        """Test motif extraction function."""
        seq = "ATGCGATCG"

        # Normal case - within bounds
        assert get_motif(seq, 2, 6) == "GCGA"

        # Start before sequence - should add N padding
        assert get_motif(seq, -2, 3) == "NNATG"

        # End after sequence - should add N padding
        assert get_motif(seq, 6, 12) == "TCGNNN"

        # Both before and after - N padding on both sides
        assert get_motif(seq, -1, 10) == "NATGCGATCGN"

        # Edge cases
        assert get_motif(seq, 0, 0) == ""
        assert get_motif("", 0, 5) == "NNNNN"

    def test_read_fasta_index(self, tmp_path):
        """Test reading FASTA index file."""
        # Create a mock .fai file
        fai_content = """chr1\t1000\t5\t50\t51
chr2\t2000\t1010\t50\t51
chr3\t500\t3015\t50\t51
"""
        fasta_path = tmp_path / "test.fa"
        fai_path = tmp_path / "test.fa.fai"

        # Write the .fai file
        fai_path.write_text(fai_content)

        # Read the index
        chrom_lengths = read_fasta_index(str(fasta_path))

        # Verify results
        assert chrom_lengths == {"chr1": 1000, "chr2": 2000, "chr3": 500}

    def test_read_fasta_index_missing(self, tmp_path):
        """Test reading missing FASTA index file."""
        fasta_path = tmp_path / "nonexistent.fa"

        with pytest.raises(FileNotFoundError) as exc_info:
            read_fasta_index(str(fasta_path))

        assert "FASTA index file not found" in str(exc_info.value)
        assert "samtools faidx" in str(exc_info.value)


class TestDetermineActualStrand:
    """Test strand determination logic."""

    def create_mock_read(self, is_paired=True, is_read1=True, is_reverse=False):
        """Create a mock read for testing."""
        read = pysam.AlignedSegment()
        read.is_paired = is_paired
        read.is_read1 = is_read1
        read.is_read2 = not is_read1 if is_paired else False
        read.is_reverse = is_reverse
        return read

    def test_paired_read1_forward(self):
        """Test read1 forward strand."""
        read = self.create_mock_read(is_paired=True, is_read1=True, is_reverse=False)
        assert determine_actual_strand(read) == "+"

    def test_paired_read1_reverse(self):
        """Test read1 reverse strand."""
        read = self.create_mock_read(is_paired=True, is_read1=True, is_reverse=True)
        assert determine_actual_strand(read) == "-"

    def test_paired_read2_forward(self):
        """Test read2 forward strand (reversed logic)."""
        read = self.create_mock_read(is_paired=True, is_read1=False, is_reverse=False)
        assert determine_actual_strand(read) == "-"

    def test_paired_read2_reverse(self):
        """Test read2 reverse strand (reversed logic)."""
        read = self.create_mock_read(is_paired=True, is_read1=False, is_reverse=True)
        assert determine_actual_strand(read) == "+"

    def test_single_end_forward(self):
        """Test single-end forward strand."""
        read = self.create_mock_read(is_paired=False, is_read1=False, is_reverse=False)
        assert determine_actual_strand(read) == "+"

    def test_single_end_reverse(self):
        """Test single-end reverse strand."""
        read = self.create_mock_read(is_paired=False, is_read1=False, is_reverse=True)
        assert determine_actual_strand(read) == "-"


class TestCountMutations:
    """Test main count_mutations function."""

    def test_count_mutations_invalid_base(self):
        """Test count_mutations with invalid base parameters."""
        with tempfile.TemporaryDirectory() as tmpdir:
            # Create dummy files
            bam_path = Path(tmpdir) / "test.bam"
            fa_path = Path(tmpdir) / "test.fa"
            bam_path.touch()
            fa_path.touch()

            # Test with invalid reference base
            result = count_mutations(
                str(bam_path), str(fa_path), ref_base="X", mut_base="G"
            )
            assert result is False

            # Test with invalid mutation base
            result = count_mutations(
                str(bam_path), str(fa_path), ref_base="A", mut_base="X"
            )
            assert result is False


class TestIntegration:
    """Integration tests requiring real BAM/FASTA files."""

    @pytest.mark.skipif(
        not Path("/home/yec/Desktop/test_rRNA.bam").exists(),
        reason="Test files not available",
    )
    def test_count_mutations_real_data(self, test_data_dir, tmp_path):
        # Paths to test files
        input_bam = test_data_dir / "test_rRNA.bam"  # Corrected BAM file name
        reference_fasta = test_data_dir / "genes.fa"  # Corrected FASTA file name
        output_tsv = tmp_path / "output.tsv"

        # Run mutation counting with default parameters
        stats = count_mutations(
            samfile=str(input_bam),
            reference=str(reference_fasta),
            output_file=str(output_tsv),
            output_bam=None,
            ref_base="A",
            mut_base="G",
            ref_base2=None,
            mut_base2=None,
            bin_size=1000000000,
            threads=1,
            save_rest=False,
            region=None,
            force=True,
            strand="both",
            pad=0,
            trim_start=0,
            trim_end=0,
            max_unc=100,
            min_con=0,
            max_sub=100,
            min_baseq=0,
            min_mapq=0,
            verbose=False,
        )

        # Assert that the function returned a dictionary with expected statistics
        assert isinstance(stats, dict)
        assert "total_processed_regions" in stats
        assert "total_mutations_found" in stats
        assert "total_raw_reads" in stats
        assert "total_reads_processed" in stats
        assert "total_reads_skipped" in stats
        assert stats["total_mutations_found"] > 0
        assert output_tsv.exists()
