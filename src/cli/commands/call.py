"""Call commands for mgatk2"""

import logging
from pathlib import Path

import click

from ..utils import (
    normalize_mito_chr,
    run_pipeline_command,
    setup_file_logging,
    split_and_index_bam,
)

logger = logging.getLogger(__name__)


@click.command()
@click.option(
    "--input-dir",
    "-i",
    default=".",
    type=click.Path(exists=True),
    help=("Directory containing BAM files"),
)
@click.option(
    "--output-dir",
    "-o",
    default="mgatk2",
    type=click.Path(),
    help="Output directory",
)
@click.option(
    "--mito-chr",
    "-g",
    "mito_genome",
    default="chrM",
    show_default=True,
    help="Mitochondrial chromosome name",
)
@click.option(
    "--ncores",
    "-c",
    default=None,
    type=int,
    help="Number of cores",
)
@click.option(
    "--mode",
    type=click.Choice(["run", "tenx"], case_sensitive=False),
    default="run",
    show_default=True,
    help="Processing mode: run or tenx",
)
@click.option("--verbose", "-v", is_flag=True, help="Enable verbose logging")
@click.option(
    "--dry-run",
    is_flag=True,
    help="Show configuration and detected BAM files, then exit without processing",
)
@click.option(
    "--split-and-index",
    is_flag=True,
    help="Split BAM files by mitochondrial reads and create indices (preprocessing only)",
)
def call(input_dir, output_dir, mito_genome, ncores, mode, verbose, dry_run, split_and_index):
    """Bulk analysis of multiple BAM files with autodetection"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    logger.info("mgatk2 call: Bulk analysis with autodetection")

    try:
        # Find all BAM files in directory
        input_path = Path(input_dir)
        bam_files = sorted(input_path.glob("*.bam"))

        if not bam_files:
            logger.error(f"No BAM files (*.bam) found in directory: {input_dir}")
            logger.info("Try: ls *.bam to check for BAM files in the directory")
            raise SystemExit(1)

        # Normalize mito chr
        mito_chr = normalize_mito_chr(mito_genome)

        # Show detected files and configuration
        logger.info("Auto-detected %s BAM files:", len(bam_files))
        for bam_file in bam_files:
            logger.info(f"  {bam_file.name} ({bam_file.stat().st_size / (1024**3):.2f} GB)")

        if dry_run:
            _show_call_configuration(bam_files, output_dir, mito_chr, ncores, mode, split_and_index)
            return

        # Setup logging to file
        log_file = Path(output_dir) / "output.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        setup_file_logging(log_file)

        # Process each BAM
        for i, bam_file in enumerate(bam_files, 1):
            sample_name = bam_file.stem  # filename without .bam
            sample_output = Path(output_dir) / sample_name

            logger.info(f"[{i}/{len(bam_files)}] Processing: {bam_file.name} → {sample_output}")

            if split_and_index:
                # Only split and index
                split_and_index_bam(bam_file, sample_output, mito_chr)
                logger.info("Split and indexed: %s", sample_name)
                continue

            # Full analysis pipeline
            sample_output.mkdir(parents=True, exist_ok=True)

            # Run pipeline based on mode
            pipeline_args = {
                "bam_path": str(bam_file),
                "output_dir": str(sample_output),
                "barcode_file": None,  # Process all reads (bulk mode)
                "barcode_tag": "CB",
                "mito_genome": mito_chr,
                "ncores": ncores,
                "verbose": verbose,
                "batch_size": None,
                "max_memory": None,
                "max_strand_bias": 1.0,
                "sequential": False,
                "dry_run": False,
            }

            if mode.lower() == "tenx":
                pipeline_args.update(
                    {
                        "base_qual": 0,
                        "min_mapq": 0,
                        "min_reads": 0,
                        "min_distance_from_end": 0,
                        "dedup_mode": "alignment_start",
                        "output_format": "txt",
                    }
                )
            else:  # run mode
                pipeline_args.update(
                    {
                        "base_qual": 20,
                        "min_mapq": 30,
                        "min_reads": 10,
                        "min_distance_from_end": 5,
                        "dedup_mode": "alignment_and_fragment_length",
                        "output_format": "hdf5",
                    }
                )

            run_pipeline_command(**pipeline_args)
            logger.info("Completed: %s", sample_name)

        if split_and_index:
            logger.info(f"Split and indexing completed for all {len(bam_files)} BAM files")
        else:
            logger.info("Analysis completed for all %s BAM files!", len(bam_files))

    except KeyboardInterrupt:
        logger.info("Bulk analysis interrupted by user")
        raise SystemExit(130) from None
    except Exception as e:
        logger.error(f"Bulk analysis failed: {e}")
        if logger.isEnabledFor(logging.DEBUG):
            import traceback

            traceback.print_exc()
        raise SystemExit(1) from None


def _show_call_configuration(bam_files, output_dir, mito_chr, ncores, mode, split_and_index):
    """Show call configuration in dry run mode."""
    logger.info("BAM files found:        %s", len(bam_files))
    for bam_file in bam_files:
        logger.info("  - %s", bam_file.name)
    logger.info("Output directory:       %s", output_dir)
    logger.info("Mitochondrial chr:      %s", mito_chr)
    logger.info("Processing mode:        %s", mode)
    logger.info("Cores:                  %s", ncores or "auto-detect")
    logger.info("Split and index only:   %s", split_and_index)
