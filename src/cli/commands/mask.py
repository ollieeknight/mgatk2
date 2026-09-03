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
    normalise_genome_name,
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
    type=str,
    required=True,
    help="Genome build (hg38, hg19, GRCh38, GRCh37, mm10, mm9, GRCm38, GRCm37)",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    help="Enable verbose logging",
)
def hardmask_fasta(input_fasta, output_fasta, genome, verbose):
    """Hard-mask reference genome FASTA with blacklists"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    try:
        genome_normalised = normalise_genome_name(genome)
        logger.info("Genome:                 %s -> %s", genome, genome_normalised)

        use_chr_prefix = detect_fasta_chr_prefix(Path(input_fasta))
        logger.info("Input FASTA:            %s", Path(input_fasta).resolve())
        logger.info("  Compressed:           %s", str(input_fasta).endswith(".gz"))
        logger.info("  Chr naming:           %s", "chr1, chr2..." if use_chr_prefix else "1, 2...")
        logger.info("Output FASTA:           %s", Path(output_fasta).resolve())

        blacklist_path = get_blacklist_path(genome_normalised)
        logger.info("Blacklist file:         %s", blacklist_path)
        numt_regions = normalise_bed_chromosomes(
            load_blacklist_regions(blacklist_path), use_chr_prefix
        )
        logger.info("Loaded %s NUMT regions to mask", len(numt_regions))
        logger.debug("Sample normalised regions: %s", numt_regions[:3])

        logger.info("Masking genome...")
        stats = mask_fasta(
            input_fasta=Path(input_fasta),
            output_fasta=Path(output_fasta),
            numt_regions=numt_regions,
        )
        logger.info(
            "Masking complete: %s contigs read, %s regions and %s bases masked",
            stats["chromosomes_processed"],
            stats["regions_masked"],
            stats["bases_masked"],
        )

    except (ValueError, FileNotFoundError) as e:
        logger.error("%s: %s", type(e).__name__, e)
        raise SystemExit(1) from e
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from e
