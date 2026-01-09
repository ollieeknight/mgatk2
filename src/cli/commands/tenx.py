"""Tenx command for mgatk2."""

import click

from ..options import common_options
from ..utils import run_pipeline_command


@click.command()
@common_options
@click.option(
    "--sequential",
    "sequential",
    is_flag=True,
    help="Process cells sequentially, disabling multiprocessing",
)
def tenx(
    bam_path,
    output_dir,
    barcode_file,
    barcode_tag,
    mito_genome,
    ncores,
    verbose,
    dry_run,
    sequential,
):
    """Run mgatk2 with original mgatk package behaviour.

    Equivalent to 'mgatk2 run' with tenx-compatible defaults:
    - No quality filtering (base_qual=0, mapq=0)
    - Alignment-based deduplication only (no fragment length)
    - Text output format
    - Batch size optimized for 10x data

    If no barcode file is provided with -b/--barcodes, performs bulk calling on all reads.
    """
    run_pipeline_command(
        bam_path=bam_path,
        output_dir=output_dir,
        barcode_file=barcode_file,
        barcode_tag=barcode_tag,
        mito_genome=mito_genome,
        ncores=ncores,
        verbose=verbose,
        batch_size=9,
        max_memory=None,
        base_qual=0,
        min_mapq=0,
        min_reads=0,
        max_strand_bias=1.0,
        min_distance_from_end=0,
        dedup_mode="alignment_start",
        output_format="txt",
        sequential=sequential,
        dry_run=dry_run,
    )
