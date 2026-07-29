"""Generic paired mitochondrial evidence command."""

from __future__ import annotations

import logging
from dataclasses import asdict
from pathlib import Path

import click

from core.config import PairedConfig
from core.exceptions import MgatkError
from processing.paired_pileup import PairedResult, run_paired_pipeline

from ..options import paired_options
from ..utils import normalise_mito_chr

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


@click.command(short_help="Paired query/baseline mitochondrial SNV evidence")
@paired_options
def paired(
    query,
    baseline,
    reference,
    output_dir,
    sample_name,
    mito_genome,
    base_qual,
    min_mapq,
    min_distance_from_end,
    max_strand_bias,
    deduplication,
    min_query_depth,
    min_baseline_depth,
    min_alt_observations,
    min_query_af,
    max_baseline_af,
    min_query_baseline_ratio,
    custom_blacklist,
    circular_edge_bases,
    input_is_consensus,
    verbose,
    dry_run,
):
    """Compare mitochondrial evidence in a QUERY population to an autologous BASELINE."""
    try:
        config = PairedConfig(
            query=query,
            baseline=baseline,
            reference=reference,
            output=output_dir,
            sample_name=sample_name,
            mito_chr=normalise_mito_chr(mito_genome),
            min_baseq=base_qual,
            min_mapq=min_mapq,
            min_distance_from_end=min_distance_from_end,
            max_strand_bias=max_strand_bias,
            deduplication=deduplication,
            min_query_depth=min_query_depth,
            min_baseline_depth=min_baseline_depth,
            min_alt_observations=min_alt_observations,
            min_query_af=min_query_af,
            max_baseline_af=max_baseline_af,
            min_query_baseline_ratio=min_query_baseline_ratio,
            custom_blacklist=custom_blacklist,
            circular_edge_bases=circular_edge_bases,
            input_is_consensus=input_is_consensus,
        )
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc

    logger.info("Effective paired configuration: %s", asdict(config))
    if dry_run:
        click.echo("Configuration valid; dry run complete.")
        return
    Path(output_dir).mkdir(parents=True, exist_ok=True)
    result = execute_paired(config, verbose)
    click.echo(
        f"Complete: {result.evidence_positions} positions, "
        f"{result.candidates} candidates, {result.pass_candidates} PASS"
    )
