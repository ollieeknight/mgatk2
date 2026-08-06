"""Tenx command for mgatk2."""

import click

from ..options import tenx_options
from ..utils import run_pipeline_command


@click.command()
@tenx_options
def tenx(
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
    nh_max,
    nm_max,
    compute_tn5,
):
    """Run mgatk2 with original mgatk package behaviour"""
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
        nh_max=nh_max,
        nm_max=nm_max,
        compute_tn5=compute_tn5,
    )
    if status:
        raise SystemExit(status)
