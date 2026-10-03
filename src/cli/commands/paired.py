"""Paired tumour/normal mitochondrial evidence command."""

import logging
from dataclasses import asdict
from pathlib import Path

import click

from cli.base import CONTEXT_SETTINGS
from core.config import PairedConfig
from core.exceptions import MgatkError
from processing.paired_pileup import run_paired_pipeline

from ..options import paired_options
from ..utils import check_alignment

logger = logging.getLogger(__name__)


@click.command(
    short_help="Paired tumour/normal mitochondrial SNV evidence",
    context_settings=CONTEXT_SETTINGS,
)
@paired_options
def paired(verbose, dry_run, **options):
    """Compare mitochondrial evidence in a TUMOR sample to an autologous NORMAL."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    try:
        config = PairedConfig(**options)
    except ValueError as exc:
        raise click.BadParameter(str(exc)) from exc
    logger.info("Effective paired configuration: %s", asdict(config))

    try:
        if dry_run:
            for alignment in (config.tumor, config.normal):
                check_alignment(alignment, config.mito_chr, reference_filename=config.reference)
            click.echo("Configuration valid; both alignments readable; dry run complete.")
            return
        Path(config.output).mkdir(parents=True, exist_ok=True)
        result = run_paired_pipeline(config)
    except (MgatkError, OSError, ValueError) as exc:
        raise click.ClickException(str(exc)) from exc

    click.echo(
        f"Complete: {result.evidence_positions} positions, "
        f"{result.candidates} candidates, {result.pass_candidates} PASS"
    )
