"""CLI options and decorators for mgatk2."""

import click


def common_options(f):
    """Common options for run, tenx, and call commands."""
    # Apply options in reverse order (last decorator applied first)
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing (no output files created)",
    )(f)
    f = click.option(
        "--format",
        "-f",
        "output_format",
        type=click.Choice(["txt", "hdf5"], case_sensitive=False),
        default="hdf5",
        show_default=True,
        help="Output format: txt (text files) or hdf5 (fast binary)",
    )(f)
    f = click.option(
        "--deduplication",
        "-d",
        "dedup_mode",
        type=click.Choice(
            ["alignment_and_fragment_length", "alignment_start", "none"],
            case_sensitive=False,
        ),
        default="alignment_and_fragment_length",
        show_default=True,
        help="Deduplication strategy",
    )(f)
    f = click.option(
        "--min-distance-from-end",
        "-e",
        "min_distance_from_end",
        default=5,
        type=int,
        show_default=True,
        help="Minimum distance from read ends (bp)",
    )(f)
    f = click.option(
        "--max-strand-bias",
        "-s",
        "max_strand_bias",
        default=1.0,
        type=float,
        show_default=True,
        help="Maximum strand bias (0-1)",
    )(f)
    f = click.option(
        "--min-reads",
        "-c",
        "min_reads",
        default=1,
        type=int,
        show_default=True,
        help="Minimum deduplicated reads per cell to include in analysis",
    )(f)
    f = click.option(
        "--mapq",
        "min_mapq",
        default=30,
        type=int,
        show_default=True,
        help="Minimum alignment/mapping quality",
    )(f)
    f = click.option(
        "--quality",
        "-q",
        "base_qual",
        default=20,
        type=int,
        show_default=True,
        help="Minimum base quality (Phred score)",
    )(f)
    f = click.option(
        "--memory",
        "-m",
        "max_memory",
        default=128,
        type=float,
        show_default=True,
        help="Maximum memory usage in GB",
    )(f)
    f = click.option(
        "--batch-size",
        "batch_size",
        default=250,
        type=int,
        show_default=True,
        help="Number of cells to process per batch",
    )(f)
    f = click.option(
        "--verbose",
        "-v",
        is_flag=True,
        default=True,
        show_default="on",
        help="Enable verbose logging",
    )(f)
    f = click.option(
        "--threads",
        "-t",
        "ncores",
        default=16,
        type=int,
        show_default=True,
        help="Number of threads for parallel processing",
    )(f)
    f = click.option(
        "--output",
        "-o",
        "output_dir",
        default="mgatk2/",
        type=click.Path(),
        show_default=True,
        help="Output directory for analysis results",
    )(f)
    f = click.option(
        "--min-barcode-reads",
        default=10,
        type=int,
        show_default=True,
        help="Minimum reads per barcode when auto-detecting from BAM",
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
        help="Barcode file (singlecell.csv, barcodes.tsv/csv, or auto-detect from BAM)",
    )(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name (e.g chrM, MT, or M)",
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


def tenx_options(f):
    """Options for tenx command with 10x-specific defaults."""
    # Apply options in reverse order
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing",
    )(f)
    f = click.option(
        "--format",
        "-f",
        "output_format",
        type=click.Choice(["txt", "hdf5"], case_sensitive=False),
        default="txt",
        show_default=True,
        help="Output format",
    )(f)
    f = click.option(
        "--deduplication",
        "-d",
        "dedup_mode",
        type=click.Choice(
            ["alignment_and_fragment_length", "alignment_start", "none"],
            case_sensitive=False,
        ),
        default="alignment_start",
        show_default=True,
        help="Deduplication strategy",
    )(f)
    f = click.option(
        "--min-distance-from-end",
        "-e",
        "min_distance_from_end",
        default=0,
        type=int,
        show_default=True,
        help="Minimum distance from read ends (bp)",
    )(f)
    f = click.option(
        "--max-strand-bias",
        "-s",
        "max_strand_bias",
        default=1.0,
        type=float,
        show_default=True,
        help="Maximum strand bias (0-1)",
    )(f)
    f = click.option(
        "--min-reads",
        "-c",
        "min_reads",
        default=0,
        type=int,
        show_default=True,
        help="Minimum deduplicated reads per cell",
    )(f)
    f = click.option(
        "--mapq",
        "min_mapq",
        default=0,
        type=int,
        show_default=True,
        help="Minimum alignment/mapping quality",
    )(f)
    f = click.option(
        "--quality",
        "-q",
        "base_qual",
        default=0,
        type=int,
        show_default=True,
        help="Minimum base quality (Phred score)",
    )(f)
    f = click.option(
        "--memory",
        "-m",
        "max_memory",
        type=float,
        help="Maximum memory usage in GB",
    )(f)
    f = click.option(
        "--batch-size",
        "batch_size",
        default=250,
        type=int,
        show_default=True,
        help="Number of cells to process per batch",
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
        "--output",
        "-o",
        "output_dir",
        default="mgatk2",
        type=click.Path(),
        help="Output directory",
    )(f)
    f = click.option(
        "--min-barcode-reads",
        default=10,
        type=int,
        show_default=True,
        help="Minimum reads per barcode when auto-detecting from BAM",
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
        help="Barcode file (singlecell.csv, barcodes.tsv/csv, or auto-detect from BAM)",
    )(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name",
    )(f)
    return click.option(
        "--input",
        "-i",
        "bam_path",
        default=".",
        type=click.Path(exists=True),
        help="Input BAM file or 10x outs/ directory",
    )(f)


def call_options(f):
    """Options for call command (bulk analysis, one BAM per cell)."""
    # Apply options in reverse order
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing",
    )(f)
    f = click.option(
        "--format",
        "-f",
        "output_format",
        type=click.Choice(["txt", "hdf5"], case_sensitive=False),
        default="hdf5",
        show_default=True,
        help="Output format: txt (text files) or hdf5 (fast binary)",
    )(f)
    f = click.option(
        "--deduplication",
        "-d",
        "dedup_mode",
        type=click.Choice(
            ["alignment_and_fragment_length", "alignment_start", "none"],
            case_sensitive=False,
        ),
        default="alignment_and_fragment_length",
        show_default=True,
        help="Deduplication strategy",
    )(f)
    f = click.option(
        "--min-distance-from-end",
        "-e",
        "min_distance_from_end",
        default=5,
        type=int,
        show_default=True,
        help="Minimum distance from read ends (bp)",
    )(f)
    f = click.option(
        "--max-strand-bias",
        "-s",
        "max_strand_bias",
        default=1.0,
        type=float,
        show_default=True,
        help="Maximum strand bias (0-1)",
    )(f)
    f = click.option(
        "--mapq",
        "min_mapq",
        default=30,
        type=int,
        show_default=True,
        help="Minimum alignment/mapping quality",
    )(f)
    f = click.option(
        "--quality",
        "-q",
        "base_qual",
        default=20,
        type=int,
        show_default=True,
        help="Minimum base quality (Phred score)",
    )(f)
    f = click.option(
        "--memory",
        "-m",
        "max_memory",
        default=128,
        type=float,
        show_default=True,
        help="Maximum memory usage in GB",
    )(f)
    f = click.option(
        "--verbose",
        "-v",
        is_flag=True,
        default=True,
        show_default="on",
        help="Enable verbose logging",
    )(f)
    f = click.option(
        "--threads",
        "-t",
        "ncores",
        default=16,
        type=int,
        show_default=True,
        help="Number of threads for parallel processing",
    )(f)
    f = click.option(
        "--output",
        "-o",
        "output_dir",
        default="mgatk2/",
        type=click.Path(),
        show_default=True,
        help="Output directory for analysis results",
    )(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name (e.g chrM, MT, or M)",
    )(f)
    return click.option(
        "--input",
        "-i",
        "bam_path",
        type=click.Path(exists=True),
        required=True,
        help="Input BAM file for bulk analysis (one BAM per cell)",
    )(f)
