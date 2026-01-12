"""Test command"""

import logging
import os
from pathlib import Path

import click

from core.exceptions import InvalidInputError
from utils.utils import validate_bam_file, validate_barcode_file

from ..options import common_options
from ..utils import (
    auto_detect_10x_structure,
    get_10x_parent_directory_name,
    normalize_mito_chr,
    parse_reads_number,
    run_pipeline_command,
    subsample_bam,
)

logger = logging.getLogger(__name__)


@click.command()
@common_options
@click.option(
    "--input-reads",
    "--n-reads",
    "n_reads",
    default="1M",
    help="Number of reads to randomly sample (e.g., 1M, 10M, 500K) [default: 1M]",
)
@click.option(
    "--batch-size",
    "batch_size",
    default=250,
    type=int,
    help="Number of cells to process per batch",
)
@click.option(
    "--quality",
    "-q",
    "base_qual",
    default=20,
    type=int,
    show_default=True,
    help="Minimum base quality (Phred score)",
)
@click.option(
    "--mapq",
    "min_mapq",
    default=30,
    type=int,
    show_default=True,
    help="Minimum alignment/mapping quality",
)
@click.option(
    "--min-reads",
    "-c",
    "min_reads",
    default=1,
    type=int,
    show_default=True,
    help="Minimum deduplicated reads per cell",
)
@click.option(
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
)
@click.option(
    "--format",
    "-f",
    "output_format",
    type=click.Choice(["txt", "hdf5"], case_sensitive=False),
    default="hdf5",
    show_default=True,
    help="Output format",
)
def test(
    bam_path,
    output_dir,
    barcode_file,
    barcode_tag,
    mito_genome,
    ncores,
    verbose,
    dry_run,
    n_reads,
    batch_size,
    base_qual,
    min_mapq,
    min_reads,
    dedup_mode,
    output_format,
):
    """Test run with random subsample of reads.

    If no barcode file is provided, automatically extracts barcodes from the BAM file.
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    try:
        n_reads_int = parse_reads_number(n_reads)
    except ValueError:
        logger.error(f"Invalid --input-reads format: {n_reads}")
        raise SystemExit(1) from None

    logger.info(f"Sampling {n_reads_int:,} reads for testing")

    try:
        bam_path, barcode_file = auto_detect_10x_structure(bam_path, barcode_file)

        validate_bam_file(bam_path)
        if barcode_file:
            validate_barcode_file(barcode_file)

        mito_chr = normalize_mito_chr(mito_genome)

        if dry_run:
            _show_test_config(
                bam_path,
                barcode_file,
                output_dir,
                n_reads_int,
                mito_chr,
                barcode_tag,
                ncores,
                batch_size,
                base_qual,
                min_mapq,
                min_reads,
                dedup_mode,
                output_format,
            )
            return

        os.makedirs(output_dir, exist_ok=True)

        # Create subsampled BAM
        test_bam = Path(output_dir) / "test_subsample.bam"
        subsample_bam(bam_path, str(test_bam), mito_chr, n_reads_int)

        run_pipeline_command(
            bam_path=str(test_bam),
            output_dir=output_dir,
            barcode_file=barcode_file,
            barcode_tag=barcode_tag,
            mito_genome=mito_genome,
            ncores=ncores,
            verbose=verbose,
            batch_size=batch_size,
            max_memory=None,
            base_qual=base_qual,
            min_mapq=min_mapq,
            min_reads=min_reads,
            max_strand_bias=1.0,
            min_distance_from_end=5,
            dedup_mode=dedup_mode,
            output_format=output_format,
            sequential=False,
            dry_run=dry_run,
            original_bam_path=bam_path,
            report_title=get_10x_parent_directory_name(bam_path),
            report_subtitle="mgatk2 test analysis",
            working_directory=os.getcwd(),
        )

    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except InvalidInputError as e:
        logger.error(f"Test setup failed: {e}")
        raise SystemExit(1) from None
    except Exception as e:
        logger.error(f"Test execution failed: {e}")
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from None


def _show_test_config(
    bam_path,
    barcode_file,
    output_dir,
    n_reads,
    mito_chr,
    barcode_tag,
    ncores,
    batch_size,
    base_qual,
    min_mapq,
    min_reads,
    dedup_mode,
    output_format,
):
    """Show test configuration in dry run mode."""
    config = [
        ("Input BAM", bam_path),
        ("Barcodes", barcode_file or "auto-extract"),
        ("Output directory", output_dir),
        ("Sample reads", f"{n_reads:,}"),
        ("Mitochondrial chr", mito_chr),
        ("Barcode tag", barcode_tag),
        ("Cores", ncores or "auto-detect"),
        ("Batch size", batch_size),
        ("Min base quality", base_qual),
        ("Min mapping quality", min_mapq),
        ("Min reads per cell", min_reads),
        ("Deduplication", dedup_mode),
        ("Output format", output_format),
    ]

    for label, value in config:
        logger.info(f"{label:<22} {value}")
