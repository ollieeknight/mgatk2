"""CLI options and decorators for mgatk2."""

import click


def _paired_quality_options(f):
    """Shared paired evidence and threshold options."""
    options = [
        click.option(
            "--max-strand-bias",
            "-s",
            default=0.9,
            type=float,
            show_default=True,
            help=(
                "Maximum tumour alternate strand imbalance, |forward - reverse| / total (0-1). "
                "Above this the allele is flagged STRAND_BIAS."
            ),
        ),
        click.option(
            "--min-distance-from-end",
            "-e",
            default=5,
            type=int,
            show_default=True,
            help="Discard base observations within this many bases of either read end",
        ),
        click.option(
            "--mapq",
            "min_mapq",
            default=20,
            type=int,
            show_default=True,
            help="Minimum alignment/mapping quality",
        ),
        click.option(
            "--quality",
            "-q",
            "base_qual",
            default=20,
            type=int,
            show_default=True,
            help="Minimum base quality (Phred score)",
        ),
    ]
    for option in options:
        f = option(f)
    return f


def paired_options(f):
    """Options for paired tumour/normal mitochondrial evidence analysis."""
    f = click.option(
        "--dry-run",
        is_flag=True,
        help=(
            "Validate the configuration, open both alignments, and check the "
            "mitochondrial contig and index are present, then exit"
        ),
    )(f)
    f = click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")(f)
    f = click.option(
        "--input-is-consensus",
        is_flag=True,
        help="Declare upstream UMI-consensus inputs (requires --deduplication none)",
    )(f)
    f = click.option(
        "--circular-edge-bases",
        default=500,
        type=int,
        show_default=True,
        help=(
            "Bases at each end of a linear mitochondrial reference to flag "
            "CIRCULAR_EDGE_UNRESOLVED and exclude from callable territory"
        ),
    )(f)
    f = click.option(
        "--autosomal-median-depth",
        default=None,
        type=float,
        help=(
            "Median autosomal depth of the tumour. Enables the POSSIBLE_NUMT filter, "
            "which flags alternate support a single-copy NuMT could account for."
        ),
    )(f)
    f = click.option(
        "--custom-blacklist",
        default=None,
        type=click.Path(exists=True, dir_okay=False),
        help="BED of chrM positions to flag BLACKLIST and exclude from callable territory",
    )(f)
    f = click.option(
        "--max-normal-af",
        default=0.01,
        type=float,
        show_default=True,
        help="Normal allele fraction above this is flagged HIGH_NORMAL_AF",
    )(f)
    f = click.option(
        "--min-tumor-af",
        default=0.005,
        type=float,
        show_default=True,
        help="Tumour allele fraction below this is flagged LOW_TUMOR_AF",
    )(f)
    f = click.option(
        "--min-alt-observations",
        default=3,
        type=int,
        show_default=True,
        help="Fewer tumour alternate observations than this is flagged LOW_ALT_OBSERVATIONS",
    )(f)
    f = click.option(
        "--min-normal-depth",
        default=5,
        type=int,
        show_default=True,
        help="Normal depth below this is flagged LOW_NORMAL_DEPTH and is not callable",
    )(f)
    f = click.option(
        "--min-tumor-depth",
        default=10,
        type=int,
        show_default=True,
        help="Tumour depth below this is flagged LOW_TUMOR_DEPTH and is not callable",
    )(f)
    f = click.option(
        "--deduplication",
        "-d",
        type=click.Choice(
            ["alignment_and_fragment_length", "alignment_start", "none"],
            case_sensitive=False,
        ),
        default="alignment_and_fragment_length",
        show_default=True,
        help=(
            "Coordinate deduplication fallback for unmarked input. "
            "Use none for already-deduplicated or UMI-consensus input."
        ),
    )(f)
    f = _paired_quality_options(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name (e.g. chrM, MT, or M)",
    )(f)
    f = click.option(
        "--sample-name",
        required=True,
        help="Filename prefix for the VCF, index, and callable BED",
    )(f)
    f = click.option(
        "--output",
        "-o",
        "output_dir",
        required=True,
        type=click.Path(),
        help="Output directory for the VCF, index, and callable BED",
    )(f)
    f = click.option(
        "--reference",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Indexed reference FASTA defining REF; also decodes CRAM input",
    )(f)
    f = click.option(
        "--normal",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Autologous normal (comparator) BAM or CRAM",
    )(f)
    return click.option(
        "--tumor",
        required=True,
        type=click.Path(exists=True, dir_okay=False),
        help="Tumour (query) BAM or CRAM",
    )(f)


# Every single-cell command counts bases the same way; only the presets differ.
# `run` targets filtered HDF5, `tenx` reproduces original mgatk behaviour, and
# `call` treats each BAM in a directory as one bulk sample.
_PRESETS = {
    "run": {
        "output_format": "hdf5",
        "dedup_mode": "alignment_and_fragment_length",
        "base_qual": 20,
        "min_mapq": 30,
        "min_reads": 1,
        "min_distance_from_end": 5,
    },
    "tenx": {
        "output_format": "txt",
        "dedup_mode": "alignment_start",
        "base_qual": 0,
        "min_mapq": 0,
        "min_reads": 0,
        "min_distance_from_end": 0,
    },
    "call": {
        "output_format": "hdf5",
        "dedup_mode": "alignment_and_fragment_length",
        "base_qual": 20,
        "min_mapq": 30,
        "min_reads": 1,
        "min_distance_from_end": 5,
    },
}


def singlecell_options(preset: str):
    """Build the shared single-cell option surface for one command preset.

    `run`, `tenx`, and `call` previously carried three verbatim copies of the
    same seventeen options, which is how they drifted apart on defaults, help
    text, and which filters were exposed at all.
    """
    defaults = _PRESETS[preset]
    is_bulk = preset == "call"

    options = [
        click.option(
            "--input",
            "-i",
            "bam_path",
            **(
                {
                    "required": True,
                    "help": (
                        "Directory containing one BAM file per cell, each treated as an "
                        "independent bulk sample"
                    ),
                }
                if is_bulk
                else {
                    "default": ".",
                    "help": (
                        "Input BAM file or 10x outs/ directory "
                        "[default: current directory, auto-detects possorted_bam.bam]"
                    ),
                }
            ),
            type=click.Path(exists=True),
        ),
        click.option(
            "--genome",
            "-g",
            "mito_genome",
            default="chrM",
            show_default=True,
            help="Mitochondrial chromosome name (e.g chrM, MT, or M)",
        ),
    ]

    if not is_bulk:
        options += [
            click.option(
                "--barcodes",
                "-b",
                "barcode_file",
                default=None,
                type=click.Path(exists=True),
                help="Barcode file (singlecell.csv, barcodes.tsv/csv, or auto-detect from BAM)",
            ),
            click.option(
                "--barcode-tag",
                "-bt",
                default="CB",
                show_default=True,
                help="BAM tag for cell barcode",
            ),
            click.option(
                "--min-barcode-reads",
                default=10,
                type=int,
                show_default=True,
                help=(
                    "Minimum reads per barcode when auto-detecting from BAM. "
                    "Ignored when --barcodes is supplied or a 10x barcode file is found."
                ),
            ),
        ]

    options += [
        click.option(
            "--output",
            "-o",
            "output_dir",
            default="mgatk2",
            type=click.Path(),
            show_default=True,
            help="Output directory for analysis results",
        ),
        click.option(
            "--threads",
            "-t",
            "ncores",
            default=None,
            type=int,
            help=(
                "Number of BAM files processed concurrently"
                if is_bulk
                else "Number of concurrent barcode shards"
            )
            + " (default: auto-detect from SLURM or system)",
        ),
        click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging"),
        click.option(
            "--memory",
            "-m",
            "max_memory",
            default=128,
            type=float,
            show_default=True,
            help=(
                # A bulk sample is one pseudo-cell, so it never shards and the
                # budget can never bind.
                "Accepted for compatibility; a bulk sample is a single cell, so "
                "the budget never binds"
                if is_bulk
                else "Memory budget in GB shared by concurrent shard workers; sets shard width"
            ),
        ),
        click.option(
            "--quality",
            "-q",
            "base_qual",
            default=defaults["base_qual"],
            type=int,
            show_default=True,
            help="Minimum base quality (Phred score)",
        ),
        click.option(
            "--mapq",
            "min_mapq",
            default=defaults["min_mapq"],
            type=int,
            show_default=True,
            help="Minimum alignment/mapping quality",
        ),
    ]

    if not is_bulk:
        options.append(
            click.option(
                "--min-reads",
                "-c",
                "min_reads",
                default=defaults["min_reads"],
                type=int,
                show_default=True,
                help=(
                    "Minimum deduplicated reads per cell to include in analysis. "
                    "Floored at 1, so 0 and 1 behave identically."
                ),
            )
        )

    options += [
        click.option(
            "--max-strand-bias",
            "-s",
            "max_strand_bias",
            default=1.0,
            type=float,
            show_default=True,
            help=(
                "Maximum per-base strand imbalance, |forward - reverse| / total (0-1). "
                "Bases above this are zeroed; 1.0 disables the filter."
            ),
        ),
        click.option(
            "--min-distance-from-end",
            "-e",
            "min_distance_from_end",
            default=defaults["min_distance_from_end"],
            type=int,
            show_default=True,
            help="Minimum distance from read ends (bp)",
        ),
        click.option(
            "--nh-max",
            "nh_max",
            default=0,
            type=int,
            show_default=True,
            help="Max NH tag (multi-mapper filter). 0=disabled; 1 matches mgatk default.",
        ),
        click.option(
            "--nm-max",
            "nm_max",
            default=0,
            type=int,
            show_default=True,
            help="Max NM/nM tag (mismatch filter). 0=disabled; 4 matches mgatk default.",
        ),
        click.option(
            "--deduplication",
            "-d",
            "dedup_mode",
            type=click.Choice(
                ["alignment_and_fragment_length", "alignment_start", "none"],
                case_sensitive=False,
            ),
            default=defaults["dedup_mode"],
            show_default=True,
            help="Deduplication strategy",
        ),
        click.option(
            "--format",
            "-f",
            "output_format",
            type=click.Choice(["txt", "hdf5"], case_sensitive=False),
            default=defaults["output_format"],
            show_default=True,
            help="Output format: txt (text files) or hdf5 (fast binary)",
        ),
        click.option(
            "--no-tn5",
            "compute_tn5",
            is_flag=True,
            default=True,
            flag_value=False,
            help="Skip Tn5 cut site tracking. Use for non-ATAC assays (RNA-seq, WGS).",
        ),
        click.option(
            "--dry-run",
            is_flag=True,
            help=(
                "Show the configuration, check the input alignment carries the "
                "mitochondrial contig and an index, then exit without processing"
            ),
        ),
    ]

    def decorator(f):
        for option in reversed(options):
            f = option(f)
        return f

    return decorator
