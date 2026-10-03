"""Call commands for mgatk2"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import click

from cli.base import CONTEXT_SETTINGS
from core.exceptions import MgatkError
from processing.processors import MP_CONTEXT

from ..options import singlecell_options
from ..utils import check_alignment, determine_cores, run_pipeline_command

logger = logging.getLogger(__name__)


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("call")
def call(bam_path, output_dir, ncores, dry_run, **options):
    """Run mgatk2 and treat each bam file as a single cell"""
    input_path = Path(bam_path)
    if input_path.is_file():
        logger.error(
            "'call' expects a directory of one-BAM-per-cell files, not a single BAM: %s. "
            "Use mgatk2 run -i <bam> for a single multi-cell BAM.",
            bam_path,
        )
        raise SystemExit(1)

    bam_files = sorted(input_path.glob("*.bam"))
    if not bam_files:
        logger.error("No BAM files (*.bam) found in directory: %s", bam_path)
        raise SystemExit(1)

    logger.info("Found %s BAM files", len(bam_files))

    if dry_run:
        try:
            for bam_file in bam_files:
                check_alignment(str(bam_file), options["mito_genome"])
        except MgatkError as e:
            logger.error("%s", e)
            raise SystemExit(1) from None
        return

    tasks = [
        {
            **options,
            "bam_path": str(bam_file),
            "output_dir": str(Path(output_dir) / bam_file.stem),
            "barcode_file": "bulk",
            # A bulk sample is one pseudo-cell and never shards, so the
            # thread budget goes on running samples side by side instead.
            "ncores": 1,
            "min_reads": 0,
        }
        for bam_file in bam_files
    ]

    workers = min(determine_cores(ncores), len(tasks))
    logger.info("Processing %s BAM files on %s worker(s)", len(tasks), workers)
    # One process per sample: each installs its own root-logger file handler.
    if workers <= 1:
        statuses = [run_pipeline_command(**task) for task in tasks]
    else:
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=mp.get_context(MP_CONTEXT)
        ) as pool:
            futures = [pool.submit(run_pipeline_command, **task) for task in tasks]
            statuses = [future.result() for future in futures]

    failed = [bam.name for bam, status in zip(bam_files, statuses, strict=True) if status]
    if failed:
        logger.error("%s of %s BAM files failed: %s", len(failed), len(tasks), ", ".join(failed))
        raise SystemExit(1)
    logger.info("Analysis completed for all %s BAM files", len(tasks))
