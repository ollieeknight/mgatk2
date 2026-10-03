"""Hard-mask a reference FASTA with the bundled nuclear NUMT blacklists."""

import logging
from pathlib import Path

import click

from cli.base import CONTEXT_SETTINGS
from utils.masking import (
    detect_fasta_chr_prefix,
    get_blacklist_path,
    load_blacklist_regions,
    mask_fasta,
    normalise_bed_chromosomes,
)

logger = logging.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@click.option(
    "--input-fasta",
    "-i",
    type=click.Path(exists=True, dir_okay=False),
    required=True,
    help="Input reference genome FASTA (plain or .gz)",
)
@click.option(
    "--output-fasta",
    "-o",
    type=click.Path(dir_okay=False),
    required=True,
    help="Output hard-masked FASTA; written gzipped when the name ends in .gz",
)
@click.option(
    "--genome",
    "-g",
    required=True,
    help="Genome build (hg38, hg19, GRCh38, GRCh37, mm10, mm9, GRCm38, GRCm37)",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def hardmask_fasta(input_fasta, output_fasta, genome, verbose):
    """Hard-mask reference genome FASTA with blacklists"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    try:
        use_chr_prefix = detect_fasta_chr_prefix(Path(input_fasta))
        blacklist_path = get_blacklist_path(genome)
        numt_regions = normalise_bed_chromosomes(
            load_blacklist_regions(blacklist_path), use_chr_prefix
        )
        logger.info(
            "Masking %s NUMT regions from %s (%s) in %s",
            len(numt_regions),
            blacklist_path.name,
            "chr-prefixed" if use_chr_prefix else "unprefixed",
            input_fasta,
        )
        logger.debug("Blacklist %s; first regions: %s", blacklist_path, numt_regions[:3])
        stats = mask_fasta(Path(input_fasta), Path(output_fasta), numt_regions)
    except (ValueError, OSError) as e:
        raise click.ClickException(str(e)) from e

    logger.info(
        "Wrote %s: %s contigs read, %s regions and %s bases masked",
        output_fasta,
        stats["chromosomes_processed"],
        stats["regions_masked"],
        stats["bases_masked"],
    )
