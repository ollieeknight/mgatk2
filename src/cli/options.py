"""CLI options and decorators for mgatk2."""

import click


def _paired_quality_options(f):
    """Shared paired/WES evidence and migration-threshold options."""
    options = [
        click.option("--max-strand-bias", "-s", default=0.9, type=float, show_default=True),
        click.option("--min-distance-from-end", "-e", default=5, type=int, show_default=True),
        click.option("--mapq", "min_mapq", default=20, type=int, show_default=True),
        click.option("--quality", "-q", "base_qual", default=20, type=int, show_default=True),
    ]
    for option in options:
        f = option(f)
    return f


def paired_options(f):
    """Options for generic paired mitochondrial evidence analysis."""
    f = click.option("--dry-run", is_flag=True, help="Validate and show configuration only")(f)
    f = click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")(f)
    f = click.option(
        "--input-is-consensus",
        is_flag=True,
        help="Declare upstream UMI-consensus inputs (requires --deduplication none)",
    )(f)
    f = click.option("--circular-edge-bases", default=500, type=int, show_default=True)(f)
    f = click.option("--custom-blacklist", type=click.Path(exists=True, dir_okay=False))(f)
    f = click.option("--min-query-baseline-ratio", default=3.0, type=float, show_default=True)(f)
    f = click.option("--max-baseline-af", default=0.01, type=float, show_default=True)(f)
    f = click.option("--min-query-af", default=0.005, type=float, show_default=True)(f)
    f = click.option("--min-alt-observations", default=3, type=int, show_default=True)(f)
    f = click.option("--min-baseline-depth", default=5, type=int, show_default=True)(f)
    f = click.option("--min-query-depth", default=10, type=int, show_default=True)(f)
    f = click.option(
        "--deduplication",
        type=click.Choice(
            ["alignment_and_fragment_length", "alignment_start", "none"],
            case_sensitive=False,
        ),
        default="alignment_and_fragment_length",
        show_default=True,
    )(f)
    f = _paired_quality_options(f)
    f = click.option("--genome", "mito_genome", default="chrM", show_default=True)(f)
    f = click.option("--sample-name", required=True)(f)
    f = click.option("--output", "output_dir", required=True, type=click.Path())(f)
    f = click.option("--reference", required=True, type=click.Path(exists=True, dir_okay=False))(f)
    f = click.option("--baseline", required=True, type=click.Path(exists=True, dir_okay=False))(f)
    return click.option("--query", required=True, type=click.Path(exists=True, dir_okay=False))(f)


def common_options(f):
    """Common options for run, tenx, and call commands."""
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing (no output files created)",
    )(f)
    f = click.option(
        "--no-tn5",
        "compute_tn5",
        is_flag=True,
        default=True,
        flag_value=False,
        help="Skip Tn5 cut site tracking. Use for non-ATAC assays (RNA-seq, WGS).",
    )(f)
    f = click.option(
        "--pileup-mode",
        "pileup_mode",
        type=click.Choice(["classic", "fast"], case_sensitive=False),
        default="fast",
        show_default=True,
        help="Pileup engine: fast (vectorized numpy, default) or classic (reference).",
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
        default=None,
        type=int,
        show_default=True,
        help="Worker batch size: cells per parallel processing batch (default: matches cores)",
    )(f)
    f = click.option(
        "--verbose",
        "-v",
        is_flag=True,
        help="Enable verbose logging",
    )(f)
    f = click.option(
        "--threads",
        "-t",
        "ncores",
        default=None,
        type=int,
        help="Number of threads (default: auto-detect from SLURM or system)",
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
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing",
    )(f)
    f = click.option(
        "--no-tn5",
        "compute_tn5",
        is_flag=True,
        default=True,
        flag_value=False,
        help="Skip Tn5 cut site tracking. Use for non-ATAC assays.",
    )(f)
    f = click.option(
        "--pileup-mode",
        "pileup_mode",
        type=click.Choice(["classic", "fast"], case_sensitive=False),
        default="fast",
        show_default=True,
        help="Pileup engine: fast (vectorized numpy, default) or classic (reference).",
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
        "--nh-max",
        "nh_max",
        default=0,
        type=int,
        show_default=True,
        help="Max NH tag (multi-mapper filter). 0=disabled; 1 matches mgatk default.",
    )(f)
    f = click.option(
        "--nm-max",
        "nm_max",
        default=0,
        type=int,
        show_default=True,
        help="Max NM/nM tag (mismatch filter). 0=disabled; 4 matches mgatk default.",
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
        default=None,
        type=int,
        show_default=True,
        help="Worker batch size: cells per parallel processing batch (default: matches cores)",
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


def wes_options(f):
    """Options for wes command (somatic mito calling from paired tumour/normal CRAMs)."""
    f = click.option(
        "--dry-run",
        is_flag=True,
        help="Show configuration and exit without processing",
    )(f)
    f = click.option(
        "--sample-name",
        "-n",
        "sample_name",
        default="sample",
        show_default=True,
        help="Sample name prefix for output files",
    )(f)
    f = click.option(
        "--blacklist",
        "blacklist_build",
        default="none",
        type=click.Choice(["hg38", "hg19", "mm10", "mm9", "none"], case_sensitive=False),
        show_default=True,
        help=(
            "Bundled NuMT blacklist build to apply to chrM positions. "
            "Note: bundled BEDs are nuclear-side NuMT positions and contain no chrM rows; "
            "primary NuMT defence is --mapq 20. Use 'none' (default) or supply "
            "--custom-blacklist for a chrM-side BED."
        ),
    )(f)
    f = click.option(
        "--custom-blacklist",
        "custom_blacklist",
        default=None,
        type=click.Path(exists=True),
        help="Path to a custom chrM-side blacklist BED (overrides --blacklist)",
    )(f)
    f = click.option(
        "--min-tn-ratio",
        "min_tn_ratio",
        default=3.0,
        type=float,
        show_default=True,
        help="Minimum tumour/normal VAF ratio to call somatic",
    )(f)
    f = click.option(
        "--max-normal-af",
        "max_normal_af",
        default=0.01,
        type=float,
        show_default=True,
        help="Maximum normal VAF before calling germline",
    )(f)
    f = click.option(
        "--min-af",
        "min_tumour_af",
        default=0.005,
        type=float,
        show_default=True,
        help="Minimum tumour variant allele fraction to report",
    )(f)
    f = click.option(
        "--min-alt-reads",
        "min_tumour_alt_reads",
        default=3,
        type=int,
        show_default=True,
        help="Minimum absolute alt read count in tumour",
    )(f)
    f = click.option(
        "--min-normal-depth",
        "min_normal_depth",
        default=5,
        type=int,
        show_default=True,
        help="Minimum read depth in normal at a called site",
    )(f)
    f = click.option(
        "--min-depth",
        "min_tumour_depth",
        default=10,
        type=int,
        show_default=True,
        help="Minimum read depth in tumour at a called site",
    )(f)
    f = _paired_quality_options(f)
    f = click.option(
        "--verbose",
        "-v",
        is_flag=True,
        help="Enable verbose logging",
    )(f)
    f = click.option(
        "--output",
        "-o",
        "output_dir",
        default="mgatk2_wes/",
        type=click.Path(),
        show_default=True,
        help="Output directory",
    )(f)
    f = click.option(
        "--genome",
        "-g",
        "mito_genome",
        default="chrM",
        show_default=True,
        help="Mitochondrial chromosome name in the CRAM (chrM / MT / M)",
    )(f)
    f = click.option(
        "--reference",
        "-r",
        "reference",
        type=click.Path(exists=True),
        required=True,
        help="Reference FASTA used to encode the CRAM files (required for decoding)",
    )(f)
    f = click.option(
        "--normal",
        "normal_cram",
        type=click.Path(exists=True),
        required=True,
        help="Normal sample CRAM (e.g. CD14+ monocytes from CLONK-WES)",
    )(f)
    return click.option(
        "--tumour",
        "tumour_cram",
        type=click.Path(exists=True),
        required=True,
        help="Tumour sample CRAM (e.g. NK cells from CLONK-WES)",
    )(f)


def call_options(f):
    """Options for call command (bulk analysis, one BAM per cell)."""
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
        help="Enable verbose logging",
    )(f)
    f = click.option(
        "--threads",
        "-t",
        "ncores",
        default=None,
        type=int,
        help="Number of threads (default: auto-detect from SLURM or system)",
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
        help="Directory containing one bulk BAM per sample",
    )(f)
