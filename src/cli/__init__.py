"""CLI main module for mgatk2."""

import logging

from .base import cli, main
from .commands import call, hardmask_fasta, paired, run, tenx

logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

cli.add_command(run)
cli.add_command(tenx)
cli.add_command(call)
cli.add_command(paired)
cli.add_command(hardmask_fasta)

__all__ = ["cli", "main"]
