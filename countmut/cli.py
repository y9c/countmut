#!/usr/bin/env python3
"""
CountMut CLI - Beautiful command-line interface for mutation counting

This module provides a modern CLI interface for CountMut with rich output,
progress tracking, and comprehensive help. Options are organized into logical
groups for better readability and user experience.

Author: Ye Chang
Date: 2025-10-23
Version: 0.0.2
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

# Configure option groups for better organization
click.rich_click.OPTION_GROUPS = {
    "countmut": [
        {
            "name": "Input/Output Options",
            "options": ["--input", "--reference", "--output", "--force"],
        },
        {
            "name": "Mutation Analysis",
            "options": ["--ref-base", "--mut-base", "--strand", "--region"],
        },
        {
            "name": "Quality Filters",
            "options": [
                "--min-baseq",
                "--min-mapq",
                "--max-sub",
                "--trim-start",
                "--trim-end",
                "--max-unc",
                "--min-con",
            ],
        },
        {
            "name": "Output Records",
            "options": ["--pad", "--save-rest"],
        },
        {
            "name": "Performance Options",
            "options": ["--threads", "--bin-size"],
        },
        {
            "name": "Alternative Mutation Tagging",
            "options": ["--ref-base2", "--mut-base2", "--output-bam"],
        },
        {
            "name": "Help & Version",
            "options": ["--help", "--version"],
        },
    ],
}

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
    "--output-bam",
    type=click.Path(path_type=str),
    help="[bold]Output BAM file[/bold] with alternative mutation tags (Yc, Zc). If not specified, a temporary file is used and then deleted.",
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
    "--ref-base2",
    default=None,
    help="[bold]Alternative reference base[/bold] for tagging (e.g., 'C')",
)
@click.option(
    "--mut-base2",
    default=None,
    help="[bold]Alternative mutation base[/bold] for tagging (e.g., 'T')",
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
    help="[bold]Save other bases[/bold] statistics (o0, o1, o2 columns)",
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
    "-p",
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
    help="[bold]Trim N bases[/bold] from read 5' end (fragment orientation)",
)
@click.option(
    "--trim-end",
    type=int,
    default=2,
    show_default=True,
    help="[bold]Trim N bases[/bold] from read 3' end (fragment orientation)",
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
    help="[bold]Max substitutions[/bold] (NS) for high-quality alignment",
)
@click.option(
    "--min-baseq",
    type=int,
    default=20,
    show_default=True,
    help="[bold]Min base quality[/bold] (Phred score) to count bases",
)
@click.option(
    "--min-mapq",
    type=int,
    default=0,
    show_default=True,
    help="[bold]Min mapping quality[/bold] (MAPQ) to count reads",
)
@click.version_option(
    importlib_metadata.version("countmut"),
    "--version",
    "-v",
    prog_name="countmut",
    message="%(prog)s %(version)s",
)
def main(
    samfile: str,
    reffile: str,
    output: str,
    output_bam: str,
    ref_base: str,
    mut_base: str,
    ref_base2: str,
    mut_base2: str,
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
    min_baseq: int,
    min_mapq: int,
):
    """
    [bold blue]🧬 CountMut - Ultra-fast Strand-aware Mutation Counter[/bold blue]

    Fast, parallel mutation counting from bisulfite sequencing (BS-seq, CAM-seq,
    GLORI-seq, eTAM-seq) BAM files with quality-based overlap deduplication.

    [bold]Key Features:[/bold]

     • [bold green]Ultra-fast[/bold green]: Direct FASTA index reading, shared file handles, BGZF multi-threading\n
     • [bold green]Bisulfite support[/bold green]: NS, Zf, Yf tag filtering for conversion analysis\n
     • [bold green]Accurate[/bold green]: Quality-based mate overlap deduplication prevents double-counting\n
     • [bold green]Parallel[/bold green]: Multi-threaded genomic window processing\n
     • [bold green]Flexible[/bold green]: Configurable filtering, strand-specific processing, auto-indexing

    [bold]Input Requirements:[/bold]

     • BAM file (coordinate-sorted, auto-creates .bai index if missing)\n
     • Reference FASTA file (auto-creates .fai index if missing)\n
     • [yellow]Required tags[/yellow]: NS, Zf, Yf (essential for bisulfite analysis)

    [bold]Output Format:[/bold]

     TSV file with columns:
     • chrom, pos, strand, motif: Position and sequence context
     • u0, u1, u2: Unconverted (reference base) counts [drop, clean, unconverted]
     • m0, m1, m2: Mutation (mutation base only) counts [drop, clean, unconverted]
     • o0, o1, o2: Other bases counts [drop, clean, unconverted] (with --save-rest)

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
    if ref_base2 and mut_base2:
        info_table.add_row("Alt. Reference base:", ref_base2.upper())
        info_table.add_row("Alt. Mutation base:", mut_base2.upper())
        info_table.add_row("Output BAM:", output_bam or "Temporary")
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
    info_table.add_row("Min base quality (Q):", str(min_baseq))
    info_table.add_row("Min mapping quality (MAPQ):", str(min_mapq))

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
            output_bam=output_bam,
            ref_base=ref_base.upper(),
            mut_base=mut_base.upper(),
            ref_base2=ref_base2,
            mut_base2=mut_base2,
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
            min_baseq=min_baseq,
            min_mapq=min_mapq,
        )

        if success:
            _duration = time.time() - _start_time
            console.print(
                f"✅ Mutation counting completed successfully! (⏱️ {format_duration(_duration)})"
            )
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
