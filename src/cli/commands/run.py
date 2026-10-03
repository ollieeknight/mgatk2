"""Main run command"""

import click
from click.core import ParameterSource

from cli.base import CONTEXT_SETTINGS

from ..options import apply_assay_preset, singlecell_options
from ..utils import run_pipeline_command


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("run")
@click.pass_context
def run(ctx, assay, **options):
    """Run mgatk2 with optimised defaults"""
    explicit = {
        name for name in options if ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE
    }
    options = apply_assay_preset(assay, options, explicit)
    raise SystemExit(run_pipeline_command(assay=assay, **options))
