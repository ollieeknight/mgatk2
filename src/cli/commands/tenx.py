"""Tenx command for mgatk2."""

import click

from cli.base import CONTEXT_SETTINGS

from ..options import singlecell_options
from ..utils import run_pipeline_command


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("tenx")
def tenx(**options):
    """Run mgatk2 with original mgatk package behaviour"""
    raise SystemExit(run_pipeline_command(**options))
