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
def run(
    bam_path,
    mito_genome,
    barcode_file,
    barcode_tag,
    min_barcode_reads,
    output_dir,
    ncores,
    verbose,
    batch_size,
    max_memory,
    base_qual,
    min_mapq,
    min_reads,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    dry_run,
):
    """Run mgatk2 with optimised defaults"""
    # Determine if processing should be sequential (auto-enabled when threads=1)
    sequential = ncores == 1

    try:
        run_pipeline_command(
            bam_path=bam_path,
            output_dir=output_dir,
            barcode_file=barcode_file,
            barcode_tag=barcode_tag,
            min_barcode_reads=min_barcode_reads,
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
