"""Call commands for mgatk2"""

import logging
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import click

from cli.base import CONTEXT_SETTINGS
from processing.processors import MP_CONTEXT

from ..options import singlecell_options
from ..utils import (
    check_alignment,
    determine_cores,
    normalise_mito_chr,
    run_pipeline_command,
)

logger = logging.getLogger(__name__)


def _run_one_sample(arguments: dict) -> tuple[str, int]:
    """Process one bulk BAM in its own process.

    Each sample owns a whole pipeline, including the root-logger file handler
    that `setup_file_logging` installs, so samples must not share a process.
    """
    sample_name = arguments.pop("sample_name")
    return sample_name, run_pipeline_command(**arguments)


@click.command(context_settings=CONTEXT_SETTINGS)
@singlecell_options("call")
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
    compute_tn5,
    nh_max,
    nm_max,
):
    """Run mgatk2 and treat each bam file as a single cell"""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
        logger.setLevel(logging.DEBUG)

    try:
        input_path = Path(bam_path)

        if input_path.is_file():
            logger.error(
                f"'call' expects a directory of one-BAM-per-cell files, not a single BAM: {bam_path}"
            )
            logger.info(
                "Try: mgatk2 run -i <bam> for a single multi-cell BAM (e.g. Tapestri/10x .cells.bam)"
            )
            raise SystemExit(1)

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
            for bam_file in bam_files:
                check_alignment(str(bam_file), mito_chr)
            _show_call_configuration(bam_files, output_dir, mito_chr, ncores)
            return

        tasks = []
        for bam_file in bam_files:
            sample_output = Path(output_dir) / bam_file.stem
            sample_output.mkdir(parents=True, exist_ok=True)
            tasks.append(
                {
                    "sample_name": bam_file.stem,
                    "bam_path": str(bam_file),
                    "output_dir": str(sample_output),
                    "barcode_file": "bulk",
                    "barcode_tag": "CB",
                    "min_barcode_reads": 10,
                    "mito_genome": mito_chr,
                    # Each sample is one bulk pseudo-cell, so a sample never
                    # shards; the thread budget buys concurrency across samples.
                    "ncores": 1,
                    "verbose": verbose,
                    "max_memory": max_memory,
                    "base_qual": base_qual,
                    "min_mapq": min_mapq,
                    "min_reads": 0,
                    "max_strand_bias": max_strand_bias,
                    "min_distance_from_end": min_distance_from_end,
                    "dedup_mode": dedup_mode,
                    "output_format": output_format,
                    "dry_run": False,
                    "compute_tn5": compute_tn5,
                    "nh_max": nh_max,
                    "nm_max": nm_max,
                }
            )

        workers = min(determine_cores(ncores), len(tasks))
        logger.info("Processing %s BAM files on %s worker(s)", len(tasks), workers)

        failures = 0
        if workers <= 1:
            for index, task in enumerate(tasks, 1):
                name = task["sample_name"]
                logger.info("[%s/%s] Processing: %s", index, len(tasks), name)
                failures += bool(_run_one_sample(dict(task))[1])
        else:
            with ProcessPoolExecutor(
                max_workers=workers, mp_context=mp.get_context(MP_CONTEXT)
            ) as pool:
                futures = [pool.submit(_run_one_sample, dict(task)) for task in tasks]
                for index, future in enumerate(as_completed(futures), 1):
                    name, status = future.result()
                    failures += bool(status)
                    logger.info(
                        "[%s/%s] %s: %s",
                        index,
                        len(tasks),
                        name,
                        "failed" if status else "complete",
                    )

        if failures:
            logger.error("%s of %s BAM files failed", failures, len(tasks))
            raise SystemExit(1)

        logger.info("Analysis completed for all %s BAM files", len(bam_files))

    except SystemExit:
        raise
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
