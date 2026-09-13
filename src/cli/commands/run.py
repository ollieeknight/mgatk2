"""Main run command"""

import logging
import os

import click
from click.core import ParameterSource

from cli.base import CONTEXT_SETTINGS
from core.exceptions import InvalidInputError, ProcessingError

from ..options import apply_assay_preset, singlecell_options
from ..utils import get_10x_parent_directory_name, run_pipeline_command

logger = logging.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("run")
@click.pass_context
def run(
    ctx,
    bam_path,
    mito_genome,
    barcode_file,
    barcode_tag,
    min_barcode_reads,
    output_dir,
    ncores,
    verbose,
    max_memory,
    base_qual,
    min_mapq,
    min_reads,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    dry_run,
    compute_tn5,
    nh_max,
    nm_max,
    assay,
    panel_bed,
):
    """Run mgatk2 with optimised defaults"""
    try:
        # The assay preset only fills options the user left at their default,
        # so an explicit flag always wins (with a warning when they disagree).
        explicit = {
            name
            for name in ("dedup_mode", "compute_tn5", "min_distance_from_end", "barcode_tag")
            if ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE
        }
        resolved = apply_assay_preset(
            assay,
            {
                "dedup_mode": dedup_mode,
                "compute_tn5": compute_tn5,
                "min_distance_from_end": min_distance_from_end,
                "barcode_tag": barcode_tag,
            },
            explicit,
        )
        dedup_mode = resolved["dedup_mode"]
        compute_tn5 = resolved["compute_tn5"]
        min_distance_from_end = resolved["min_distance_from_end"]
        barcode_tag = resolved["barcode_tag"]
        status = run_pipeline_command(
            bam_path=bam_path,
            output_dir=output_dir,
            barcode_file=barcode_file,
            barcode_tag=barcode_tag,
            min_barcode_reads=min_barcode_reads,
            mito_genome=mito_genome,
            ncores=ncores,
            verbose=verbose,
            max_memory=max_memory,
            base_qual=base_qual,
            min_mapq=min_mapq,
            min_reads=min_reads,
            max_strand_bias=max_strand_bias,
            min_distance_from_end=min_distance_from_end,
            dedup_mode=dedup_mode,
            output_format=output_format,
            dry_run=dry_run,
            compute_tn5=compute_tn5,
            nh_max=nh_max,
            nm_max=nm_max,
            assay=assay,
            panel_bed=panel_bed,
            original_bam_path=bam_path,
            report_title=get_10x_parent_directory_name(bam_path),
            report_subtitle="mgatk2 output analysis",
            working_directory=os.getcwd(),
        )
        if status:
            raise SystemExit(status)

    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except (InvalidInputError, ProcessingError) as e:
        logger.error(f"{type(e).__name__}: {e}")
        raise SystemExit(1) from None
    except Exception as e:
        logger.error(f"Unexpected error: {e}")
        raise SystemExit(1) from None
