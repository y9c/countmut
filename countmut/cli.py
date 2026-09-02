#!/usr/bin/env python3
"""
CountMut CLI -- unified strand-aware counter with a C backend.

There is no --mode / --ref-base / --mut-base: one counter, and the output is
your choice.
   countmut -i x.bam -r ref.fa -o out.tsv                     composition table
   countmut -i x.bam -r ref.fa -o out.tsv --output-format "{pos+1}\\t{ref}\\t{t}/({c}+{t})"
                                                               custom columns
   countmut -i x.bam -r ref.fa -o out.vcf --vcf                allele VCF

The heavy computation runs in the bundled C core (``backend/countmut_core``);
this wrapper drives it and renders a rich summary.
"""

import os
import sys
from importlib import metadata as importlib_metadata

import rich.box
import rich_click as click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from .backend import run_backend
from .model import EngineConfig, FilterConfig, StrandConfig

try:
    __version__ = importlib_metadata.version("countmut")
except importlib_metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.2.1"

click.rich_click.TEXT_MARKUP = "rich"
click.rich_click.SHOW_ARGUMENTS = True
click.rich_click.GROUP_ARGUMENTS_OPTIONS = True
click.rich_click.STYLE_ERRORS_SUGGESTION = "magenta italic"
click.rich_click.ERRORS_SUGGESTION = (
    "Try running the '--help' flag for more information."
)

click.rich_click.OPTION_GROUPS = {
    "countmut": [
        {"name": "Input/Output", "options": ["--input", "--reference", "--output"]},
        {
            "name": "Mode & Engine",
            "options": ["--engine", "--region", "--threads"],
        },
        {
            "name": "Base/Allele / Output Options",
            "options": [
                "--output-format",
                "--strandless",
                "--count-indels",
                "--max-depth",
                "--vcf",
            ],
        },
        {
            "name": "Filters & Trimming (expressions)",
            "options": ["--expression", "--pile-expression"],
        },
        {"name": "Misc", "options": ["--verbose", "--version", "--help"]},
    ]
}

console = Console()


@click.command(
    cls=click.RichCommand,
    context_settings={"help_option_names": ["-h", "--help"]},
    no_args_is_help=True,
    name="countmut",
)
@click.version_option(__version__, "-v", "--version", prog_name="countmut")
@click.option(
    "-i",
    "--input",
    "samfile",
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Input BAM file (coordinate-sorted; index auto-created)",
)
@click.option(
    "-r",
    "--reference",
    type=click.Path(exists=True, path_type=str),
    required=True,
    help="Reference FASTA (index auto-created)",
)
@click.option(
    "-o",
    "--output",
    type=click.Path(path_type=str),
    default=None,
    help="Output file (default: stdout)",
)
@click.option(
    "--engine",
    type=click.Choice(["auto", "read-walk", "pileup"], case_sensitive=False),
    default="auto",
    show_default=True,
    help="BAM-walk strategy (auto picks per mode)",
)
@click.option(
    "--region", type=str, default=None, help="Region, e.g. 'chr1:1000000-2000000'"
)
@click.option(
    "-t", "--threads", type=int, default=None, help="Worker threads (default: auto)"
)
@click.option(
    "--strandless",
    is_flag=True,
    default=False,
    show_default=True,
    help="Collapse '+'/'-' strands into one row (base/allele mode; default is per-strand)",
)
@click.option(
    "--count-indels",
    is_flag=True,
    help="Append ins/del/ref_skip/fail columns (base mode)",
)
@click.option(
    "--max-depth",
    type=int,
    default=0,
    show_default=True,
    help="Per-position depth cap (pileup engine; 0 = unlimited, counts all reads)",
)
@click.option("--vcf", is_flag=True, help="Emit VCF in allele mode")
@click.option(
    "-e",
    "--expression",
    "read_expr",
    default=None,
    help='Lua read filter (e.g. "mapq >= 20 and flags & UNMAP == 0")',
)
@click.option(
    "-p",
    "--pile-expression",
    "pile_expr",
    default=None,
    help="Lua site filter (e.g. \"ref == 'A' and depth >= 5 and g > 2\")",
)
@click.option(
    "--output-format",
    "output_format",
    default=None,
    help=(
        "Output-row template: literal text plus {expr} placeholders over the site "
        'values, e.g. "{pos+1}\\t{ref}\\t{a}/({a}+{t})" (helpers: round(), int(), …). '
        "Default output is the per-base composition; --vcf switches to an allele "
        "VCF.  For a template, add --fmt-header."
    ),
)
@click.option(
    "--fmt-header",
    default=None,
    help="Header line for a custom --output-format template (default: none)",
)
@click.option(
    "--verbose",
    is_flag=True,
    default=False,
    help="Real-time per-region progress on stderr (large BAMs)",
)
def main(
    samfile,
    reference,
    output,
    engine,
    region,
    threads,
    strandless,
    count_indels,
    max_depth,
    vcf,
    read_expr,
    pile_expr,
    output_format,
    fmt_header,
    verbose,
):
    """[bold green]countmut: one counter, output format is yours[/bold green]."""
    # One counting core.  Output: per-base composition by default, allele VCF
    # with --vcf, or any row template via --output-format.
    output_expr = output_format if output_format else None
    # Panels/logs go to stderr so stdout stays pure data.
    console = Console(stderr=True)

    # The C core keeps conservative defaults for read acceptance.
    fcfg = FilterConfig(
        max_depth=max_depth,
    )
    scfg = StrandConfig(process="both", strandless=strandless)
    ecfg = EngineConfig(
        engine=engine,
        threads=threads,
        region=region,
        count_indels=count_indels,
        strandless=strandless,
        vcf=vcf,
        read_expr=read_expr,
        pile_expr=pile_expr,
        output_expr=output_expr,
        fmt_header=fmt_header,
        verbose=verbose,
    )

    config_table = Table(box=rich.box.MINIMAL, show_header=False)
    config_table.add_column("Setting", style="bold")
    config_table.add_column("Value", style="cyan")
    config_table.add_row("Input BAM:", os.path.abspath(samfile))
    config_table.add_row("Reference:", os.path.abspath(reference))
    config_table.add_row("Output:", os.path.abspath(output) if output else "(stdout)")
    config_table.add_row(
        "Output:",
        "allele VCF" if vcf else ("custom template" if output_expr else "composition"),
    )
    config_table.add_row("Engine:", engine)
    if read_expr:
        config_table.add_row("Read filter:", read_expr)
    if pile_expr:
        config_table.add_row("Site filter:", pile_expr)
    config_table.add_row("Threads:", str(threads or "auto"))
    if region:
        config_table.add_row("Region:", region)
    console.print(
        Panel(
            config_table,
            title="[bold blue]Processing Configuration[/bold blue]",
            border_style="blue",
            expand=False,
        )
    )

    # Make sure output dir exists
    if output:
        d = os.path.dirname(os.path.abspath(output))
        if d:
            os.makedirs(d, exist_ok=True)

    stats = run_backend(samfile, reference, output, fcfg=fcfg, scfg=scfg, ecfg=ecfg)

    if not stats.get("success", False):
        console.print(f"[red]Error during counting: {stats.get('error')}[/red]")
        sys.exit(1)

    summary = Table(box=rich.box.MINIMAL, show_header=False)
    summary.add_column("Metric", style="bold")
    summary.add_column("Value", style="cyan")
    summary.add_row("Backend:", "C core")
    summary.add_row("Sites / depth rows:", f"{stats['total_sites']:,}")
    summary.add_row("Wall-clock:", f"{stats['elapsed']:.3f}s")
    console.print(
        Panel(
            summary,
            title="[bold green]Processing Summary[/bold green]",
            border_style="green",
            expand=False,
        )
    )


if __name__ == "__main__":
    main()
