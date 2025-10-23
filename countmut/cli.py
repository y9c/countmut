#!/usr/bin/env python3
"""
CountMut CLI - Beautiful command-line interface for mutation counting

This module provides a modern CLI interface for CountMut with rich output,
progress tracking, and comprehensive help.

Author: Ye Chang
Date: 2025-10-23
"""

import os
import sys
import time
from importlib import metadata as importlib_metadata

import pysam
import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .core import count_mutations
from .utils import format_duration

# Configure rich-click
click.rich_click.TEXT_MARKUP = "rich"
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.ERRORS_SUGGESTION = (
    "Try running the '--help' flag for more information."
)
click.rich_click.ERRORS_EPILOGUE = "To find out more, visit [link=https://github.com/y9c/countmut]https://github.com/y9c/countmut[/link]"
click.rich_click.TEXT_EMOJIS = True

console = Console()


def validate_bam_file(bam_file: str, threads: int = 1) -> bool:
    """Validate BAM file and create index if needed."""
    try:
        # Check if index exists and is newer than BAM file
        index_file = bam_file + ".bai"
        bam_mtime = os.path.getmtime(bam_file)

        if not os.path.exists(index_file):
            console.print("📝 Creating BAM index...")
            pysam.index(bam_file, "-@", str(threads))
        else:
            # Check if index is older than BAM file
            index_mtime = os.path.getmtime(index_file)
            if index_mtime < bam_mtime:
                console.print("🔄 BAM index is older than BAM file, rebuilding...")
                pysam.index(bam_file, "-@", str(threads))

        # Check if sorted by coordinate
        with pysam.AlignmentFile(bam_file, "rb") as f:
            header = f.header
            if "HD" in header and "SO" in header["HD"]:
                return header["HD"]["SO"] == "coordinate"
            return False
    except Exception as e:
        console.print(f"[red]Validation error: {e}[/red]")
        return False



@click.command(
    name="countmut",
    no_args_is_help=True,
    context_settings={"help_option_names": ["-h", "--help"]},
    epilog="""
Examples:

# Basic usage
countmut -i input.bam -r reference.fa

# Save to file with custom parameters
countmut -i input.bam -r reference.fa -o mutations.tsv --ref-base T --mut-base C

# Use more threads and smaller bins
countmut -i input.bam -r reference.fa -t 16 -b 5000

# Save additional statistics
countmut -i input.bam -r reference.fa -s

# Process specific region
countmut -i input.bam -r reference.fa --region chr1:1000000-2000000
    """,
)
@click.option(
    "-i",
    "--input",
    "samfile",
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Input BAM file (coordinate-sorted, indexed)",
)
@click.option(
    "-r",
    "--reference",
    "reffile",
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Reference FASTA file",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=str),
    help="[bold]Output file[/bold] for mutation counts (TSV format). If not specified, prints to stdout.",
)
@click.option(
    "--ref-base",
    default="A",
    show_default=True,
    help="[bold]Reference base[/bold] to count mutations from (A, T, G, or C)",
)
@click.option(
    "--mut-base",
    default="G",
    show_default=True,
    help="[bold]Mutation base[/bold] to count (A, T, G, or C)",
)
@click.option(
    "-b",
    "--bin-size",
    type=int,
    default=10_000,
    show_default=True,
    help="[bold]Genomic bin size[/bold] for parallel processing (in base pairs)",
)
@click.option(
    "-t",
    "--threads",
    type=int,
    default=None,
    help="[bold]Number of threads[/bold] for parallel processing (default: auto-detect)",
)
@click.option(
    "-s",
    "--save-rest",
    is_flag=True,
    help="[bold]Save additional statistics[/bold] including y0, y1, y2 columns",
)
@click.option(
    "--region",
    type=str,
    help="[bold]Genomic region[/bold] to process (e.g., 'chr1:1000000-2000000')",
)
@click.option(
    "-f",
    "--force",
    is_flag=True,
    help="[bold]Overwrite output file[/bold] without prompting",
)
@click.option(
    "--strand",
    type=click.Choice(["both", "forward", "reverse"], case_sensitive=False),
    default="both",
    show_default=True,
    help="[bold]Strand processing[/bold]: 'both' (default), 'forward' (+ only), or 'reverse' (- only)",
)
@click.option(
    "--pad",
    type=int,
    default=15,
    show_default=True,
    help="[bold]Motif half-window padding[/bold] around each site",
)
@click.option(
    "--trim-start",
    type=int,
    default=2,
    show_default=True,
    help="[bold]Trim bases[/bold] at read 5' end when counting",
)
@click.option(
    "--trim-end",
    type=int,
    default=2,
    show_default=True,
    help="[bold]Trim bases[/bold] at read 3' end when counting",
)
@click.option(
    "--max-unc",
    type=int,
    default=3,
    show_default=True,
    help="[bold]Max unconverted threshold[/bold] (Zf) to consider converted",
)
@click.option(
    "--min-con",
    type=int,
    default=1,
    show_default=True,
    help="[bold]Min converted threshold[/bold] (Yf) to consider converted",
)
@click.option(
    "--max-sub",
    type=int,
    default=1,
    show_default=True,
    help="[bold]Max substitutions[/bold] (NS) to consider mapped",
)
@click.version_option(importlib_metadata.version("countmut"), "--version", "-v", prog_name="countmut", message="%(prog)s %(version)s")
def main(
    samfile: str,
    reffile: str,
    output: str,
    ref_base: str,
    mut_base: str,
    bin_size: int,
    threads: int,
    save_rest: bool,
    region: str,
    force: bool,
    strand: str,
    pad: int,
    trim_start: int,
    trim_end: int,
    max_unc: int,
    min_con: int,
    max_sub: int,
):
    """
    [bold blue]🧬 CountMut - Fast Mutation Counting from BAM Pileup[/bold blue]

    A fast, parallel tool for counting mutations from BAM pileup data with
    bisulfite conversion analysis and genomic window processing.

    [bold]Key Features:[/bold]

     • [bold green]Parallel processing[/bold green]: Multi-threaded genomic window processing\n
     • [bold green]Bisulfite analysis[/bold green]: Built-in conversion detection\n
     • [bold green]Flexible output[/bold green]: TSV format with optional statistics\n
     • [bold green]Memory efficient[/bold green]: Streaming processing for large files\n
     • [bold green]Rich output[/bold green]: Progress bars and detailed statistics\n
     • [bold green]Region support[/bold green]: Process specific genomic regions

    [bold]Input Requirements:[/bold]

     • Input BAM file, coordinate-sorted (required), indexed with .bai file (created automatically if missing)\n
     • Reference FASTA file (required)

    [bold]Output Format:[/bold]

     The output TSV file contains the following columns:
     • chrom: Chromosome name
     • pos: Position (1-based)
     • strand: Strand (+ or -)
     • motif: Sequence motif around the position
     • u0, d0: Drop counts (unconverted, total)
     • u1, d1: Clean counts (unconverted, total)
     • u2, d2: Unconverted counts (unconverted, total)
     • y0, y1, y2: Additional counts (if --save-rest is used)

    """

    # Banner disabled per request

    # Validate input files
    if not os.path.exists(samfile):
        console.print(f"[red]❌ Input BAM file '{samfile}' does not exist![/red]")
        return

    if not os.path.exists(reffile):
        console.print(f"[red]❌ Reference file '{reffile}' does not exist![/red]")
        return

    # Validate base parameters
    valid_bases = {"A", "T", "G", "C"}
    if ref_base.upper() not in valid_bases:
        console.print(
            f"[red]❌ Invalid reference base '{ref_base}'. Must be one of: {', '.join(valid_bases)}[/red]"
        )
        return

    if mut_base.upper() not in valid_bases:
        console.print(
            f"[red]❌ Invalid mutation base '{mut_base}'. Must be one of: {', '.join(valid_bases)}[/red]"
        )
        return

    # Validate numeric parameters
    if bin_size <= 0:
        console.print(f"[red]❌ Bin size must be positive, got: {bin_size}[/red]")
        return

    if threads is not None and threads <= 0:
        console.print(f"[red]❌ Thread count must be positive, got: {threads}[/red]")
        return

    # Check output file
    if output and os.path.exists(output):
        if not force:
            response = console.input(
                f"[yellow]⚠️  Output file '{output}' already exists. Overwrite? (y/N): [/yellow]"
            )
            if response.lower() != "y":
                console.print("[yellow]Operation cancelled.[/yellow]")
                return
        else:
            console.print(f"[yellow]⚠️  Overwriting existing file: {output}[/yellow]")

    # Create output directory if needed
    if output:
        output_dir = os.path.dirname(output)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    # Display configuration
    info_table = Table(show_header=False, box=None, padding=(0, 1))
    info_table.add_column(style="bold blue", justify="left")
    info_table.add_column(style="white", justify="left")

    info_table.add_row("Input BAM:", samfile)
    info_table.add_row("Reference:", reffile)
    info_table.add_row("Output:", output or "stdout")
    info_table.add_row("Reference base:", ref_base.upper())
    info_table.add_row("Mutation base:", mut_base.upper())
    info_table.add_row("Bin size:", f"{bin_size:,}")
    info_table.add_row("Threads:", str(threads or "auto"))
    info_table.add_row("Save additional stats:", "Yes" if save_rest else "No")
    if region:
        info_table.add_row("Region:", region)
    info_table.add_row("Pad:", str(pad))
    info_table.add_row("Trim start:", str(trim_start))
    info_table.add_row("Trim end:", str(trim_end))
    info_table.add_row("Max unconverted (Zf):", str(max_unc))
    info_table.add_row("Min converted (Yf):", str(min_con))
    info_table.add_row("Max substitutions (NS):", str(max_sub))

    console.print(
        Panel(
            info_table,
            title="[bold green]Processing configuration[/bold green]",
            border_style="green",
            padding=(1, 2),
        )
    )

    # Validate BAM file
    console.print("🔍 Validating BAM file...")
    if not validate_bam_file(samfile, threads or 1):
        console.print("[red]❌ BAM file validation failed![/red]")
        return
    console.print("✅ BAM file validation passed")

    # Process the file
    console.print("🚀 Starting mutation counting...")

    try:
        _start_time = time.time()
        success = count_mutations(
            samfile=samfile,
            reffile=reffile,
            output_file=output,
            ref_base=ref_base.upper(),
            mut_base=mut_base.upper(),
            bin_size=bin_size,
            threads=threads,
            save_rest=save_rest,
            region=region,
            strand=strand,
            pad=pad,
            trim_start=trim_start,
            trim_end=trim_end,
            max_unc=max_unc,
            min_con=min_con,
            max_sub=max_sub,
        )

        if success:
            _duration = time.time() - _start_time
            console.print(f"✅ Mutation counting completed successfully! (⏱️ {format_duration(_duration)})")
            if output:
                console.print(f"📄 Results saved to: {output}")
        else:
            console.print("[red]❌ Mutation counting failed![/red]")
            sys.exit(1)

    except KeyboardInterrupt:
        console.print("\n[yellow]⚠️  Processing interrupted by user[/yellow]")
        sys.exit(1)
    except Exception as e:
        console.print(f"[red]❌ Unexpected error: {e}[/red]")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
