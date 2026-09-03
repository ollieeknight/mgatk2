"""Paired tumour/normal mitochondrial evidence command."""

from __future__ import annotations

import logging
from dataclasses import asdict
from pathlib import Path

import click

from cli.base import CONTEXT_SETTINGS
from core.config import PairedConfig
from core.exceptions import MgatkError
from processing.paired_pileup import PairedResult, run_paired_pipeline

from ..options import paired_options
from ..utils import check_alignment, normalise_mito_chr

logger = logging.getLogger(__name__)


def execute_paired(config: PairedConfig, verbose: bool = False) -> PairedResult:
    """Execute paired analysis with command-level error translation."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    try:
        return run_paired_pipeline(config)
    except KeyboardInterrupt:
        raise click.Abort() from None
    except (MgatkError, OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc


@click.command(
    short_help="Paired tumour/normal mitochondrial SNV evidence",
    context_settings=CONTEXT_SETTINGS,
)
@paired_options
def paired(
    tumor,
    normal,
    reference,
    output_dir,
    sample_name,
    mito_genome,
    base_qual,
    min_mapq,
    min_distance_from_end,
    max_strand_bias,
    deduplication,
    min_tumor_depth,
    min_normal_depth,
    min_alt_observations,
    min_tumor_af,
    max_normal_af,
    custom_blacklist,
    autosomal_median_depth,
    circular_edge_bases,
    input_is_consensus,
    verbose,
    dry_run,
):
    """Compare mitochondrial evidence in a TUMOR sample to an autologous NORMAL."""
    try:
        config = PairedConfig(
            tumor=tumor,
            normal=normal,
            reference=reference,
            output=output_dir,
            sample_name=sample_name,
            mito_chr=normalise_mito_chr(mito_genome),
            min_baseq=base_qual,
            min_mapq=min_mapq,
            min_distance_from_end=min_distance_from_end,
            max_strand_bias=max_strand_bias,
            deduplication=deduplication,
            min_tumor_depth=min_tumor_depth,
            min_normal_depth=min_normal_depth,
            min_alt_observations=min_alt_observations,
            min_tumor_af=min_tumor_af,
            max_normal_af=max_normal_af,
            custom_blacklist=custom_blacklist,
            autosomal_median_depth=autosomal_median_depth,
            circular_edge_bases=circular_edge_bases,
            input_is_consensus=input_is_consensus,
        )
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc

    logger.info("Effective paired configuration: %s", asdict(config))
    if dry_run:
        try:
            for alignment in (config.tumor, config.normal):
                check_alignment(alignment, config.mito_chr, reference_filename=config.reference)
        except MgatkError as exc:
            raise click.ClickException(str(exc)) from exc
        click.echo("Configuration valid; both alignments readable; dry run complete.")
        return
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    result = execute_paired(config, verbose)
    click.echo(
        f"Complete: {result.evidence_positions} positions, "
        f"{result.candidates} candidates, {result.pass_candidates} PASS"
    )
