"""Tests for utility functions."""

from countmut.utils import format_duration, get_output_headers, write_output


class TestFormatDuration:
    """Test duration formatting function."""

    def test_format_duration_seconds(self):
        """Test formatting durations in seconds."""
        assert format_duration(0.5) == "0.50s"
        assert format_duration(12.34) == "12.34s"
        assert format_duration(59.99) == "59.99s"

    def test_format_duration_minutes(self):
        """Test formatting durations in minutes."""
        assert format_duration(60) == "1m 0s"
        assert format_duration(90) == "1m 30s"
        assert format_duration(125) == "2m 5s"
        assert format_duration(3599) == "59m 59s"

    def test_format_duration_hours(self):
        """Test formatting durations in hours."""
        assert format_duration(3600) == "1h 0m 0s"
        assert format_duration(3661) == "1h 1m 1s"
        assert format_duration(7200) == "2h 0m 0s"
        assert format_duration(7323) == "2h 2m 3s"

    def test_format_duration_edge_cases(self):
        """Test edge cases."""
        assert format_duration(0) == "0.00s"
        assert format_duration(59.5) == "59.50s"
        assert format_duration(60.5) == "1m 0s"  # Rounds down


class TestGetOutputHeaders:
    """Test output header generation."""

    def test_get_output_headers_minimal(self):
        """Test minimal headers without additional stats."""
        headers = get_output_headers(save_rest=False)

        assert len(headers) == 10
        assert headers[0] == "chrom"
        assert headers[1] == "pos"
        assert headers[2] == "strand"
        assert headers[3] == "motif"
        assert "o0" not in headers
        assert "o1" not in headers
        assert "o2" not in headers

    def test_get_output_headers_with_rest(self):
        """Test headers with additional statistics."""
        headers = get_output_headers(save_rest=True)

        assert len(headers) == 13
        assert headers[0] == "chrom"
        assert "o0" in headers
        assert "o1" in headers
        assert "o2" in headers


class TestWriteOutput:
    """Test output writing function."""

    def test_write_output_to_file(self, tmp_path):
        """Test writing output to a file."""
        output_file = tmp_path / "output.tsv"
        results = [
            ["chr1", 100, "+", "ATGC", 1, 2, 3, 4, 5, 6],
            ["chr1", 200, "-", "GCTA", 7, 8, 9, 10, 11, 12],
        ]

        write_output(results, str(output_file), save_rest=False)

        assert output_file.exists()
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check header
        assert "chrom" in lines[0]
        assert "pos" in lines[0]

        # Check data rows
        assert len(lines) == 3  # Header + 2 data rows
        assert "chr1" in lines[1]
        assert "chr1" in lines[2]

    def test_write_output_with_rest(self, tmp_path):
        """Test writing output with additional statistics."""
        output_file = tmp_path / "output_rest.tsv"
        results = [
            ["chr1", 100, "+", "ATGC", 1, 2, 3, 4, 5, 6, 7, 8, 9],
        ]

        write_output(results, str(output_file), save_rest=True)

        assert output_file.exists()
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check that o columns are in header
        assert "o0" in lines[0]
        assert "o1" in lines[0]
        assert "o2" in lines[0]

    def test_write_output_creates_directory(self, tmp_path):
        """Test that output creates parent directories."""
        output_file = tmp_path / "subdir" / "nested" / "output.tsv"
        results = [["chr1", 100, "+", "ATGC", 1, 2, 3, 4, 5, 6]]

        write_output(results, str(output_file), save_rest=False)

        assert output_file.exists()
        assert output_file.parent.exists()

    def test_write_output_empty_results(self, tmp_path):
        """Test writing empty results."""
        output_file = tmp_path / "empty.tsv"
        results = []

        write_output(results, str(output_file), save_rest=False)

        assert output_file.exists()
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Should only have header
        assert len(lines) == 1
        assert "chrom" in lines[0]

    def test_write_output_to_stdout(self, capsys):
        """Test writing output to stdout."""
        results = [
            ["chr1", 100, "+", "ATGC", 1, 2, 3, 4, 5, 6],
        ]

        write_output(results, None, save_rest=False)

        captured = capsys.readouterr()
        assert "chrom" in captured.out
        assert "chr1" in captured.out
        assert "100" in captured.out
