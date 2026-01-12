"""Mask command"""

import logging
from pathlib import Path

import click

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--input-fasta",
    "-i",
    type=click.Path(exists=True),
    required=True,
    help="Input reference genome FASTA file",
)
@click.option(
    "--output-fasta",
    "-o",
    type=click.Path(),
    required=True,
    help="Output hard-masked FASTA file",
)
@click.option(
    "--genome",
    "-g",
    type=str,
    required=True,
    help="Genome build (hg38, hg19, GRCh38, GRCh37, mm10, mm9, GRCm38, GRCm37)",
)
@click.option(
    "--mt-chrom",
    default=None,
    help="Mitochondrial chromosome name (auto-detected from genome if not provided)",
)
@click.option(
    "--mask-numts/--no-mask-numts",
    default=True,
    show_default=True,
    help="Mask NUMT regions in nuclear chromosomes (recommended)",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
def hardmask_fasta(input_fasta, output_fasta, genome, mt_chrom, mask_numts, verbose):
    """Hard-mask reference genome FASTA with blacklists"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    try:
        from utils.genome_utils import (
            get_blacklist_path,
            get_mt_chrom_name,
            load_blacklist_regions,
            normalize_genome_name,
        )
        from utils.masking import mask_fasta

        # Normalise genome name (handles GRCh38 → hg38, etc.)
        genome_normalized = normalize_genome_name(genome)
        logger.info("Genome:                 %s → %s", genome, genome_normalized)

        # Get blacklist file
        blacklist_path = get_blacklist_path(genome_normalized)
        logger.info("Blacklist file:         %s", blacklist_path)

        # Get MT chromosome name
        if mt_chrom is None:
            mt_chrom = get_mt_chrom_name(genome_normalized)
        logger.info("MT chromosome:          %s", mt_chrom)
        logger.info("Input FASTA:            %s", Path(input_fasta).resolve())
        logger.info("Output FASTA:           %s", Path(output_fasta).resolve())

        # Load blacklist regions
        numt_regions = None
        if mask_numts:
            numt_regions = load_blacklist_regions(blacklist_path)

        # Perform masking
        logger.info("Masking genome...")
        mask_fasta(
            input_fasta=Path(input_fasta),
            output_fasta=Path(output_fasta),
            numt_regions=numt_regions,
            mt_chrom=mt_chrom,
            mt_blacklist_positions=None,
        )

    except ValueError as e:
        logger.error("Error: %s", e)
        raise SystemExit(1) from e
    except FileNotFoundError as e:
        logger.error("File not found: %s", e)
        raise SystemExit(1) from e
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from e
