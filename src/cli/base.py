"""Base CLI setup for mgatk2."""

from importlib.metadata import version

import click


@click.group()
@click.version_option(version=version("mgatk2"))
def cli():
    """\b
    mgatk2: An improved mitochondrial genome analysis toolkit
    Inspired by Caleb Lareau's original mgatk for 10x single-cell ATAC-seq data"""


def main():
    """Entry point for mgatk2 CLI."""
    cli()
