"""Deprecated ``mgatk2 wes`` compatibility adapter."""

import click

from core.config import PairedConfig

from ..options import wes_options
from ..utils import normalise_mito_chr
from .paired import execute_paired


@click.command(short_help="Deprecated alias for paired mitochondrial analysis")
@wes_options
def wes(
    tumour_cram,
    normal_cram,
    reference,
    mito_genome,
    output_dir,
    base_qual,
    min_mapq,
    min_distance_from_end,
    max_strand_bias,
    min_tumour_depth,
    min_normal_depth,
    min_tumour_af,
    max_normal_af,
    min_tn_ratio,
    min_tumour_alt_reads,
    blacklist_build,
    custom_blacklist,
    sample_name,
    verbose,
    dry_run,
):
    """Compatibility wrapper; use ``mgatk2 paired`` for new analyses."""
    click.echo(
        "Warning: 'mgatk2 wes' is deprecated; use 'mgatk2 paired' with "
        "--query/--baseline. A legacy TSV projection will be written.",
        err=True,
    )
    try:
        config = PairedConfig(
            query=tumour_cram,
            baseline=normal_cram,
            reference=reference,
            output=output_dir,
            sample_name=sample_name,
            mito_chr=normalise_mito_chr(mito_genome),
            min_baseq=base_qual,
            min_mapq=min_mapq,
            min_distance_from_end=min_distance_from_end,
            max_strand_bias=max_strand_bias,
            deduplication="none",
            min_query_depth=min_tumour_depth,
            min_baseline_depth=min_normal_depth,
            min_alt_observations=min_tumour_alt_reads,
            min_query_af=min_tumour_af,
            max_baseline_af=max_normal_af,
            min_query_baseline_ratio=min_tn_ratio,
            custom_blacklist=custom_blacklist,
            input_is_consensus=True,
            write_legacy_tsv=True,
        )
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc
    if blacklist_build.lower() != "none":
        click.echo(
            "Bundled NUMT BEDs are nuclear-side resources and are not applied as chrM filters.",
            err=True,
        )
    if dry_run:
        click.echo("Configuration valid; dry run complete.")
        return
    result = execute_paired(config, verbose)
    click.echo(f"Complete: {result.candidates} candidates; legacy projection written.")
