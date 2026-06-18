"""WES somatic mitochondrial variant calling command."""

import logging
from pathlib import Path

import click

from core.config import WesConfig
from core.exceptions import InvalidInputError, ProcessingError

from ..options import wes_options
from ..utils import normalise_mito_chr, setup_file_logging

logger = logging.getLogger(__name__)


@click.command()
@wes_options
def wes(
    tumour_cram,
    normal_cram,
    reference,
    mito_genome,
    output_dir,
    base_qual,
    min_mapq,
    min_distance_from_end,
    max_strand_bias,
    min_tumour_depth,
    min_normal_depth,
    min_tumour_af,
    max_normal_af,
    min_tn_ratio,
    min_tumour_alt_reads,
    blacklist_build,
    custom_blacklist,
    sample_name,
    verbose,
    dry_run,
):
    """Somatic mitochondrial variant calling from paired tumour/normal CRAMs"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)

    mito_chr = normalise_mito_chr(mito_genome)

    config = WesConfig(
        min_baseq=base_qual,
        min_mapq=min_mapq,
        min_distance_from_end=min_distance_from_end,
        max_strand_bias=max_strand_bias,
        mito_chr=mito_chr,
        min_tumour_depth=min_tumour_depth,
        min_normal_depth=min_normal_depth,
        min_tumour_af=min_tumour_af,
        max_normal_af=max_normal_af,
        min_tn_ratio=min_tn_ratio,
        min_tumour_alt_reads=min_tumour_alt_reads,
        blacklist_build=blacklist_build,
        custom_blacklist=custom_blacklist,
    )

    logger.info("mgatk2 wes")
    logger.info("  Tumour CRAM:         %s", tumour_cram)
    logger.info("  Normal CRAM:         %s", normal_cram)
    logger.info("  Reference FASTA:     %s", reference)
    logger.info("  Mito chromosome:     %s", mito_chr)
    logger.info("  Output directory:    %s", output_dir)
    logger.info("  Sample name:         %s", sample_name)
    logger.info("  Min tumour depth:    %s", min_tumour_depth)
    logger.info("  Min normal depth:    %s", min_normal_depth)
    logger.info("  Min tumour AF:       %s", min_tumour_af)
    logger.info("  Max normal AF:       %s", max_normal_af)
    logger.info("  Min T/N ratio:       %s", min_tn_ratio)
    logger.info("  Min alt reads:       %s", min_tumour_alt_reads)
    logger.info("  Max strand bias:     %s", max_strand_bias)
    logger.info("  Min MAPQ:            %s", min_mapq)
    logger.info("  Min base quality:    %s", base_qual)
    logger.info("  Blacklist build:     %s", blacklist_build)
    if custom_blacklist:
        logger.info("  Custom blacklist:    %s", custom_blacklist)

    if dry_run:
        logger.info("Dry run — exiting without processing")
        return

    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    setup_file_logging(output_path / "wes.log")

    try:
        from processing.wes_pileup import run_wes_pipeline

        result = run_wes_pipeline(
            tumour_cram=tumour_cram,
            normal_cram=normal_cram,
            reference=reference,
            output_dir=output_dir,
            config=config,
            sample_name=sample_name,
        )

        logger.info(
            "Complete: %d PASS variants → %s",
            result["pass_variants"],
            result["output_tsv"],
        )

    except KeyboardInterrupt:
        raise SystemExit(130) from None
    except (InvalidInputError, ProcessingError) as e:
        logger.error("%s: %s", type(e).__name__, e)
        raise SystemExit(1) from None
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        if verbose:
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from None
