"""Call commands for mgatk2"""

import logging
from pathlib import Path

import click

from ..options import call_options
from ..utils import (
    normalise_mito_chr,
    run_pipeline_command,
)

logger = logging.getLogger(__name__)


@click.command()
@call_options
def call(
    bam_path,
    mito_genome,
    output_dir,
    ncores,
    verbose,
    max_memory,
    base_qual,
    min_mapq,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    dry_run,
):
    """Run mgatk2 and treat each bam file as a single cell"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    sequential = ncores == 1

    try:
        input_path = Path(bam_path)
        bam_files = sorted(input_path.glob("*.bam"))

        if not bam_files:
            logger.error(f"No BAM files (*.bam) found in directory: {bam_path}")
            logger.info("Try: ls *.bam to check for BAM files in the directory")
            raise SystemExit(1)

        mito_chr = normalise_mito_chr(mito_genome)

        logger.info("Auto-detected %s BAM files:", len(bam_files))
        for bam_file in bam_files:
            logger.info(f"  {bam_file.name} ({bam_file.stat().st_size / (1024**3):.2f} GB)")

        if dry_run:
            _show_call_configuration(bam_files, output_dir, mito_chr, ncores)
            return

        for i, bam_file in enumerate(bam_files, 1):
            sample_name = bam_file.stem  # filename without .bam
            sample_output = Path(output_dir) / sample_name

            logger.info(f"[{i}/{len(bam_files)}] Processing: {bam_file.name} → {sample_output}")

            sample_output.mkdir(parents=True, exist_ok=True)

            status = run_pipeline_command(
                bam_path=str(bam_file),
                output_dir=str(sample_output),
                barcode_file="bulk",
                barcode_tag="CB",
                min_barcode_reads=10,
                mito_genome=mito_chr,
                ncores=ncores,
                verbose=verbose,
                batch_size=1,
                max_memory=max_memory,
                base_qual=base_qual,
                min_mapq=min_mapq,
                min_reads=0,
                max_strand_bias=max_strand_bias,
                min_distance_from_end=min_distance_from_end,
                dedup_mode=dedup_mode,
                output_format=output_format,
                sequential=sequential,
                dry_run=False,
            )
            if status:
                raise SystemExit(status)
            logger.info("Completed: %s", sample_name)

        logger.info("Analysis completed for all %s BAM files", len(bam_files))

    except KeyboardInterrupt:
        logger.info("Bulk analysis interrupted by user")
        raise SystemExit(130) from None
    except Exception as e:
        logger.error(f"Bulk analysis failed: {e}")
        if logger.isEnabledFor(logging.DEBUG):
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from None


def _show_call_configuration(bam_files, output_dir, mito_chr, ncores):
    """Show call configuration in dry run mode."""
    logger.info("BAM files found:        %s", len(bam_files))
    for bam_file in bam_files:
        logger.info("  - %s", bam_file.name)
    logger.info("Output directory:       %s", output_dir)
    logger.info("Mitochondrial chr:      %s", mito_chr)
    logger.info("Cores:                  %s", ncores or "auto-detect")
