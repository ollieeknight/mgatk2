"""CLI options and decorators for mgatk2."""

import click


def common_options(f):
    """Common options for run, tenx, and test commands."""
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing (no output files created)",
    )(f)
    f = click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")(f)
    f = click.option(
        "--threads",
        "-t",
        "ncores",
        default=None,
        type=int,
        help="Number of threads for parallel processing",
    )(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name (e.g chrM, MT, or M)",
    )(f)
    f = click.option(
        "--barcode-tag",
        "-bt",
        default="CB",
        show_default=True,
        help="BAM tag for cell barcode",
    )(f)
    f = click.option(
        "--barcodes",
        "-b",
        "barcode_file",
        default=None,
        type=click.Path(exists=True),
        help="Barcode file (e.g., barcodes.tsv). If not provided, performs bulk calling on all reads.",
    )(f)
    f = click.option(
        "--output",
        "-o",
        "output_dir",
        default="mgatk2",
        type=click.Path(),
        help="Output directory for analysis results",
    )(f)
    return click.option(
        "--input",
        "-i",
        "bam_path",
        default=".",
        type=click.Path(exists=True),
        help=(
            "Input BAM file or 10x outs/ directory "
            "[default: current directory, auto-detects possorted_bam.bam]"
        ),
    )(f)
