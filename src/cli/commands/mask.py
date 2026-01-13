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
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    default=False,
    show_default="off",
    help="Enable verbose logging",
)
def hardmask_fasta(input_fasta, output_fasta, genome, mt_chrom, mask_numts, verbose):
    """Hard-mask reference genome FASTA with blacklists"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    try:
        from utils.genome_utils import (
            detect_fasta_chr_prefix,
            get_blacklist_path,
            get_mt_chrom_name,
            load_blacklist_regions,
            normalise_bed_chromosomes,
            normalise_genome_name,
        )
        from utils.masking import mask_fasta

        # Normalise genome name (handles GRCh38 → hg38, etc.)
        genome_normalised = normalise_genome_name(genome)
        logger.info("Genome:                 %s → %s", genome, genome_normalised)

        # Get blacklist file
        blacklist_path = get_blacklist_path(genome_normalised)
        logger.info("Blacklist file:         %s", blacklist_path)

        # Detect FASTA chromosome naming convention
        use_chr_prefix = detect_fasta_chr_prefix(Path(input_fasta))
        logger.info("Input FASTA:            %s", Path(input_fasta).resolve())
        logger.info("  Compressed:           %s", str(input_fasta).endswith(".gz"))
        logger.info("  Chr naming:           %s", "chr1, chr2..." if use_chr_prefix else "1, 2...")

        # Get MT chromosome name and normalise it
        if mt_chrom is None:
            mt_chrom = get_mt_chrom_name(genome_normalised)
            # Normalise MT chromosome name to match FASTA
            if not use_chr_prefix and mt_chrom.startswith("chr"):
                mt_chrom = mt_chrom[3:]
            elif use_chr_prefix and not mt_chrom.startswith("chr"):
                mt_chrom = f"chr{mt_chrom}"
        logger.info("MT chromosome:          %s", mt_chrom)
        logger.info("Output FASTA:           %s", Path(output_fasta).resolve())

        # Load and normalise blacklist regions
        numt_regions = None
        if mask_numts:
            numt_regions = load_blacklist_regions(blacklist_path)
            numt_regions = normalise_bed_chromosomes(numt_regions, use_chr_prefix)
            logger.info("Loaded %s NUMT regions to mask", len(numt_regions))
            if verbose:
                logger.debug(
                    "Sample normalised regions: %s",
                    numt_regions[:3] if len(numt_regions) >= 3 else numt_regions,
                )

        # Perform masking
        logger.info("Masking genome...")
        mask_fasta(
            input_fasta=Path(input_fasta),
            output_fasta=Path(output_fasta),
            numt_regions=numt_regions,
            mt_chrom=mt_chrom,
            mt_blacklist_positions=None,
        )

        # Log results
        logger.info("✓ Masking complete!")

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
