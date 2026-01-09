"""Main run command"""

import logging
import os

import click

from core.exceptions import InvalidInputError, ProcessingError

from ..options import common_options
from ..utils import get_10x_parent_directory_name, run_pipeline_command

logger = logging.getLogger(__name__)


@click.command()
@common_options
@click.option(
    "--batch-size",
    "batch_size",
    default=250,
    type=int,
    help="Number of cells to process per batch",
)
@click.option(
    "--memory",
    "-m",
    "max_memory",
    type=float,
    help="Maximum memory usage in GB",
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
    help="Minimum deduplicated reads per cell to include in analysis",
)
@click.option(
    "--max-strand-bias",
    "-s",
    "max_strand_bias",
    default=1.0,
    type=float,
    show_default=True,
    help="Maximum strand bias (0-1). Default, 1.0, keeps all data",
)
@click.option(
    "--min-distance-from-end",
    "-e",
    "min_distance_from_end",
    default=5,
    type=int,
    show_default=True,
    help="Minimum distance from read ends (bp). Default: 5bp, 0 keeps all bases",
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
    help="Output format: txt (text files) or hdf5 (fast binary)",
)
@click.option(
    "--sequential",
    "sequential",
    is_flag=True,
    help="Process cells sequentially, disabling multiprocessing",
)
def run(
    bam_path,
    output_dir,
    barcode_file,
    barcode_tag,
    mito_genome,
    ncores,
    verbose,
    dry_run,
    batch_size,
    max_memory,
    base_qual,
    min_mapq,
    min_reads,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    sequential,
):
    """Run mgatk2 with quality filtering, deduplication, and HDF5 output.

    If no barcode file is provided with -b/--barcodes, performs bulk calling on all reads.
    """
    if verbose:
        logging.basicConfig(level=logging.DEBUG)

    try:
        run_pipeline_command(
            bam_path=bam_path,
            output_dir=output_dir,
            barcode_file=barcode_file,
            barcode_tag=barcode_tag,
            mito_genome=mito_genome,
            ncores=ncores,
            verbose=verbose,
            batch_size=batch_size,
            max_memory=max_memory,
            base_qual=base_qual,
            min_mapq=min_mapq,
            min_reads=min_reads,
            max_strand_bias=max_strand_bias,
            min_distance_from_end=min_distance_from_end,
            dedup_mode=dedup_mode,
            output_format=output_format,
            sequential=sequential,
            dry_run=dry_run,
            original_bam_path=bam_path,
            report_title=get_10x_parent_directory_name(bam_path),
            report_subtitle="mgatk2 output analysis",
            working_directory=os.getcwd(),
        )

    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except (InvalidInputError, ProcessingError) as e:
        logger.error(f"{type(e).__name__}: {e}")
        raise SystemExit(1) from None
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise SystemExit(1) from None
