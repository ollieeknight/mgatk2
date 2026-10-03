"""Main run command"""

import click
from click.core import ParameterSource

from cli.base import CONTEXT_SETTINGS

from ..options import apply_assay_preset, singlecell_options
from ..utils import run_pipeline_command

# Options an --assay preset may fill.
_PRESET_OPTIONS = ("dedup_mode", "compute_tn5", "min_distance_from_end", "barcode_tag")


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("run")
@click.pass_context
def run(ctx, assay, **options):
    """Run mgatk2 with optimised defaults"""
    # The assay preset only fills options the user left at their default,
    # so an explicit flag always wins (with a warning when they disagree).
    explicit = {
        name
        for name in _PRESET_OPTIONS
        if ctx.get_parameter_source(name) is ParameterSource.COMMANDLINE
    }
    options.update(
        apply_assay_preset(assay, {name: options[name] for name in _PRESET_OPTIONS}, explicit)
    )
    raise SystemExit(run_pipeline_command(assay=assay, **options))
