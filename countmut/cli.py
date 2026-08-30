#!/usr/bin/env python3
"""
CountMut CLI -- unified strand-aware counter with a C backend.

Usage: countmut --mode mutation --ref-base A --mut-base G ...
       countmut --mode base ...
       countmut --mode allele --vcf ...

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
from .model import EngineConfig, FilterConfig, MutationConfig, StrandConfig

try:
    __version__ = importlib_metadata.version("countmut")
except importlib_metadata.PackageNotFoundError:  # pragma: no cover
    __version__ = "0.1.5"

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
            "options": ["--mode", "--engine", "--region", "--threads"],
        },
        {
            "name": "Mutation Options",
            "options": ["--ref-base", "--mut-base", "--pad", "--save-rest"],
        },
        {
            "name": "Base/Allele Options",
            "options": [
                "--split-strand",
                "--count-indels",
                "--min-depth",
                "--min-allele-support",
                "--vcf",
            ],
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
    "--mode",
    type=click.Choice(["mutation", "base", "allele"], case_sensitive=False),
    default="mutation",
    show_default=True,
    help="Counting mode",
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
    "--ref-base", default="A", show_default=True, help="Reference base to count from"
)
@click.option(
    "--mut-base", default="G", show_default=True, help="Mutation base to count"
)
@click.option(
    "--pad", type=int, default=15, show_default=True, help="Motif half-window"
)
@click.option(
    "-s", "--save-rest", is_flag=True, help="Also emit o0/o1/o2 (other bases)"
)
@click.option(
    "--split-strand",
    is_flag=True,
    default=True,
    show_default=True,
    help="Emit separate '+'/'-' rows (base/allele mode)",
)
@click.option(
    "--count-indels",
    is_flag=True,
    help="Append ins/del/ref_skip/fail columns (base mode)",
)
@click.option(
    "--min-depth",
    type=int,
    default=0,
    show_default=True,
    help="Minimum site depth to report (base/allele mode)",
)
@click.option(
    "--max-depth",
    type=int,
    default=0,
    show_default=True,
    help="Per-position depth cap (pileup engine; 0 = unlimited, counts all reads)",
)
@click.option(
    "--min-allele-support",
    type=int,
    default=1,
    show_default=True,
    help="Minimum alt allele support (allele mode)",
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
    "--verbose",
    is_flag=True,
    default=False,
    help="Real-time per-region progress on stderr (large BAMs)",
)
def main(
    samfile,
    reference,
    output,
    mode,
    engine,
    region,
    threads,
    ref_base,
    mut_base,
    pad,
    save_rest,
    split_strand,
    count_indels,
    min_depth,
    max_depth,
    min_allele_support,
    vcf,
    read_expr,
    pile_expr,
    verbose,
):
    """[bold green]countmut: unified ultra-fast strand-aware counter[/bold green]."""
    # Panels/logs go to stderr so stdout stays pure data.
    console = Console(stderr=True)

    # The C core keeps conservative defaults for read acceptance and mutation
    # categorisation.
    fcfg = FilterConfig(
        min_mapq=0,
        min_baseq=20,
        max_sub=1,
        max_unc=3,
        min_con=1,
        trim_start=2,
        trim_end=2,
        max_depth=max_depth,
    )
    mcfg = (
        MutationConfig(
            ref_base=ref_base, mut_base=mut_base, pad=pad, save_rest=save_rest
        )
        if mode == "mutation"
        else None
    )
    scfg = StrandConfig(process="both", split=split_strand)
    ecfg = EngineConfig(
        engine=engine,
        mode=mode,
        threads=threads,
        region=region,
        count_indels=count_indels,
        split_strand=split_strand,
        min_depth=min_depth,
        min_allele_support=min_allele_support,
        vcf=vcf,
        read_expr=read_expr,
        pile_expr=pile_expr,
        verbose=verbose,
    )

    config_table = Table(box=rich.box.MINIMAL, show_header=False)
    config_table.add_column("Setting", style="bold")
    config_table.add_column("Value", style="cyan")
    config_table.add_row("Input BAM:", os.path.abspath(samfile))
    config_table.add_row("Reference:", os.path.abspath(reference))
    config_table.add_row("Output:", os.path.abspath(output) if output else "(stdout)")
    config_table.add_row("Mode:", mode)
    config_table.add_row("Engine:", engine)
    if mode == "mutation":
        config_table.add_row("Substitution:", f"{ref_base} -> {mut_base}")
        config_table.add_row("Motif pad:", str(pad))
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

    stats = run_backend(
        samfile, reference, output, fcfg=fcfg, mcfg=mcfg, scfg=scfg, ecfg=ecfg
    )

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
