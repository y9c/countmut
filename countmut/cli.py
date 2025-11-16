#!/usr/bin/env python3
"""
CountMut CLI - Beautiful command-line interface for mutation counting

This module provides a modern CLI interface for CountMut with rich output,
progress tracking, and comprehensive help. Options are organized into logical
groups for better readability and user experience.

Author: Ye Chang
Date: 2025-10-23
"""

import os
from importlib import metadata as importlib_metadata

import rich.box
import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .core import count_mutations

__version__ = importlib_metadata.version("countmut")

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

# Define option groups for better CLI help formatting
click.rich_click.OPTION_GROUPS = {
    "countmut": [
        {
            "name": "Input/Output Options",
            "options": [
                "--input",
                "--reference",
                "--output",
                "--output-bam",
                "--force",
            ],
        },
        {
            "name": "Mutation Options",
            "options": ["--ref-base", "--mut-base", "--ref-base2", "--mut-base2"],
        },
        {
            "name": "Tuning Options",
            "options": [
                "--bin-size",
                "--threads",
                "--save-rest",
                "--region",
                "--strand",
                "--pad",
                "--max-reads-per-chunk",
            ],
        },
        {
            "name": "Filtering Options",
            "options": [
                "--trim-start",
                "--trim-end",
                "--max-unc",
                "--min-con",
                "--max-sub",
                "--min-baseq",
                "--min-mapq",
            ],
        },
        {
            "name": "Miscellaneous Options",
            "options": [
                "--verbose",
                "--version",
                "--help",
            ],
        },
    ]
}

console = Console()


@click.command(
    cls=click.RichCommand,
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
)
@click.version_option(__version__, "-v", "--version", prog_name="countmut")
@click.option(
    "-i",
    "--input",
    "samfile",  # Parameter name for main function
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Input BAM file (coordinate-sorted, indexed)",
)
@click.option(
    "-r",
    "--reference",
    "reference",  # Parameter name for main function
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Reference FASTA file",
)
@click.option(
    "-o",
    "--output",
    "output_file",  # Parameter name for main function
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
    help="[dim](base-level)[/dim] [bold]Trim N bases[/bold] from read 5' end (fragment orientation)",
)
@click.option(
    "--trim-end",
    type=int,
    default=2,
    show_default=True,
    help="[dim](base-level)[/dim] [bold]Trim N bases[/bold] from read 3' end (fragment orientation)",
)
@click.option(
    "--max-unc",
    type=int,
    default=3,
    show_default=True,
    help="[dim](read-level)[/dim] [bold]Max unconverted threshold[/bold] (Zf) to consider converted",
)
@click.option(
    "--min-con",
    type=int,
    default=1,
    show_default=True,
    help="[dim](read-level)[/dim] [bold]Min converted threshold[/bold] (Yf) to consider converted",
)
@click.option(
    "--max-sub",
    type=int,
    default=1,
    show_default=True,
    help="[dim](read-level)[/dim] [bold]Max substitutions[/bold] (NS) for high-quality alignment",
)
@click.option(
    "--min-mapq",
    type=int,
    default=0,
    show_default=True,
    help="[dim](read-level)[/dim] [bold]Min mapping quality[/bold] (MAPQ) to count reads",
)
@click.option(
    "--min-baseq",
    type=int,
    default=20,
    show_default=True,
    help="[dim](base-level)[/dim] [bold]Min base quality[/bold] (Phred score) to count bases",
)
@click.option(
    "--max-reads-per-chunk",
    type=int,
    default=100_000,
    show_default=True,
    help="[bold]Maximum reads per chunk[/bold] before splitting region for parallel processing (default: 100,000)",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Enable verbose logging output.",
)
def main(
    samfile: str,
    reference: str,
    output_file: str | None,
    output_bam: str | None,
    ref_base: str,
    mut_base: str,
    ref_base2: str | None,
    mut_base2: str | None,
    bin_size: int,
    threads: int,
    save_rest: bool,
    region: str | None,
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
    max_reads_per_chunk: int,
    verbose: bool,
):
    """
    [bold green]🚀 Ultra-fast strand-aware mutation counter[/bold green]

    `countmut` is a high-performance tool for counting specific mutations
    (e.g., A->G) from BAM files, designed for bisulfite / m6A sequencing analysis.
    """
    console = Console()

    # Get absolute paths for better logging and error messages
    input_bam_abs = os.path.abspath(samfile)
    reference_fasta_abs = os.path.abspath(reference)
    output_file_abs = os.path.abspath(output_file) if output_file else None
    output_bam_abs = os.path.abspath(output_bam) if output_bam else None

    # Check for output file overwrite
    if output_file_abs and os.path.exists(output_file_abs):
        if not force:
            response = console.input(
                f"[yellow]⚠️  Output file '{output_file_abs}' already exists. Overwrite? (y/N): [/yellow]"
            )
            if response.lower() != "y":
                console.print("[yellow]Operation cancelled.[/yellow]")
                return
        else:
            console.print(
                f"[yellow]⚠️  Overwriting existing file: {output_file_abs}[/yellow]"
            )

    # Create output directory if needed
    if output_file_abs:
        output_dir = os.path.dirname(output_file_abs)
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)

    try:
        # Display processing configuration in a rich panel
        config_table = Table(box=rich.box.MINIMAL, show_header=False)
        config_table.add_column("Setting", style="bold")
        config_table.add_column("Value", style="cyan")

        config_table.add_row("Input BAM:", input_bam_abs)
        config_table.add_row("Reference:", reference_fasta_abs)
        config_table.add_row(
            "Output:", output_file_abs if output_file_abs else "(stdout)"
        )
        if output_bam_abs:
            config_table.add_row("Output BAM:", output_bam_abs)
        if ref_base2 and mut_base2:
            config_table.add_row("Reference base:", ref_base)
            config_table.add_row("Mutation base:", mut_base)
            config_table.add_row("Alt. Reference base:", ref_base2)
            config_table.add_row("Alt. Mutation base:", mut_base2)
        else:
            config_table.add_row("Reference base:", ref_base)
            config_table.add_row("Mutation base:", mut_base)
        config_table.add_row("Bin size:", f"{bin_size:,}")
        config_table.add_row("Threads:", str(threads))
        config_table.add_row("Save additional stats:", "Yes" if save_rest else "No")
        if region:
            config_table.add_row("Region:", region)
        config_table.add_row("Strand processing:", strand)
        config_table.add_row("Pad:", str(pad))
        config_table.add_row("Trim start:", str(trim_start))
        config_table.add_row("Trim end:", str(trim_end))
        config_table.add_row("Max unconverted (Zf):", str(max_unc))
        config_table.add_row("Min converted (Yf):", str(min_con))
        config_table.add_row("Max substitutions (NS):", str(max_sub))
        config_table.add_row("Min base quality (Q):", str(min_baseq))
        config_table.add_row("Min mapping quality (MAPQ):", str(min_mapq))
        config_table.add_row("Max reads per chunk:", f"{max_reads_per_chunk:,}")
        config_table.add_row("Verbose logging:", "Yes" if verbose else "No")

        config_panel = Panel(
            config_table,
            title="[bold blue]Processing Configuration[/bold blue]",
            border_style="blue",
            expand=False,
        )
        console.print(config_panel)

        stats = count_mutations(
            samfile=input_bam_abs,
            reference=reference_fasta_abs,
            output_file=output_file_abs,
            output_bam=output_bam_abs,
            ref_base=ref_base,
            mut_base=mut_base,
            ref_base2=ref_base2,
            mut_base2=mut_base2,
            bin_size=bin_size,
            threads=threads,
            save_rest=save_rest,
            region=region,
            force=force,  # Added force keyword argument
            strand=strand,
            pad=pad,
            trim_start=trim_start,
            trim_end=trim_end,
            max_unc=max_unc,
            min_con=min_con,
            max_sub=max_sub,
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            max_reads_per_chunk=max_reads_per_chunk,
            verbose=verbose,
        )

        # Display final statistics in a rich panel
        if stats:
            stats_table = Table(
                box=rich.box.MINIMAL, show_header=False
            )  # Updated to use rich.box.MINIMAL
            stats_table.add_column("Metric", style="bold")
            stats_table.add_column("Value", style="cyan")
            stats_table.add_row(
                "Regions processed:", f"{stats['total_processed_regions']:,}"
            )
            stats_table.add_row(
                "Regions skipped (no reads):", f"{stats['total_skipped_regions']:,}"
            )
            stats_table.add_row("Total raw reads:", f"{stats['total_raw_reads']:,}")
            stats_table.add_row(
                "Total reads processed:", f"{stats['total_reads_processed']:,}"
            )
            stats_table.add_row(
                "Total reads skipped:", f"{stats['total_reads_skipped']:,}"
            )
            stats_table.add_row(
                "Total mutations found:", f"{stats['total_mutations_found']:,}"
            )

            if stats["total_reads_skipped"] > 0:
                # Display detailed skipped reads only if there were skipped reads
                stats_table.add_section()
                stats_table.add_row("Skipped details:", "")
                stats_table.add_row(
                    "  Wrong strand:",
                    f"{stats.get('total_skipped_wrong_strand_agg', 0):,}",
                )
                stats_table.add_row(
                    "  Unmapped:", f"{stats.get('total_skipped_unmapped_agg', 0):,}"
                )  # New: Display unmapped
                stats_table.add_row(
                    "  Duplicate:", f"{stats.get('total_skipped_duplicate_agg', 0):,}"
                )  # New: Display duplicate
                stats_table.add_row(
                    "  Secondary:", f"{stats.get('total_skipped_secondary_agg', 0):,}"
                )  # New: Display secondary
                stats_table.add_row(
                    "  Failed mismatch filter:",
                    f"{stats.get('total_skipped_mismatch_filter_agg', 0):,}",
                )  # New: Display mismatch filter skipped
                stats_table.add_row(
                    "  Failed mapping quality filter:",
                    f"{stats.get('total_skipped_mapq_filter_agg', 0):,}",
                )  # New: Display mapq filter skipped
                stats_table.add_row(
                    "  Failed conversion filter:",
                    f"{stats.get('total_skipped_conversion_filter_agg', 0):,}",
                )  # New: Display conversion filter skipped
                stats_table.add_row(
                    "  Missing tags:",
                    f"{stats.get('total_skipped_missing_tags_agg', 0):,}",
                )
                stats_table.add_row(
                    "  No sequence:",
                    f"{stats.get('total_skipped_no_sequence_agg', 0):,}",
                )

            final_panel = Panel(
                stats_table,
                title="[bold green]Processing Summary[/bold green]",
                border_style="green",
                expand=False,
            )
            console.print(final_panel)

    except Exception as e:
        console.print(f"[red]Error during mutation counting: {e}[/red]")
        console.print("[red]Please check the logs for more details.[/red]")


if __name__ == "__main__":
    main()
