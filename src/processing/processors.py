"""Shard scheduling for mgatk2."""

import logging
import multiprocessing as mp
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from tqdm import tqdm

from processing.pileup import plan_shards, scan_shard

logger = logging.getLogger(__name__)

# fork on Linux: cheap, and the parent heap is small by design here.
# spawn on macOS: fork is unsafe with the Objective-C runtime (Python 3.12+).
MP_CONTEXT = "fork" if platform.system() == "Linux" else "spawn"


def build_tasks(bam_path, config, barcodes, reference_filename=None) -> list[tuple]:
    """Split the barcode list into contiguous shards that fit the memory budget."""
    per_shard = plan_shards(len(barcodes), config)
    return [
        (str(bam_path), config, barcodes[lo : lo + per_shard], lo, reference_filename)
        for lo in range(0, len(barcodes), per_shard)
    ]


def process_shards(bam_path, config, barcodes, writer, reference_filename=None) -> dict:
    """Scan chrM once per shard, writing each finished shard straight to disk."""
    tasks = build_tasks(bam_path, config, barcodes, reference_filename)
    n_cells = len(barcodes)
    workers = min(config.performance.n_cores, len(tasks))

    logger.info(
        "Counting %s cells in %s shard(s) of up to %s cells on %s worker(s)",
        f"{n_cells:,}",
        len(tasks),
        len(tasks[0][2]),
        workers,
    )

    totals = {"total_reads": 0, "duplicate_reads": 0, "kept_reads": 0, "cells_passed": 0}

    def absorb(result):
        writer.write_shard(result, barcodes)
        totals["total_reads"] = max(totals["total_reads"], result.total_reads)
        totals["duplicate_reads"] += result.duplicate_reads
        totals["kept_reads"] += int(result.n_reads.sum())
        totals["cells_passed"] += int(np.count_nonzero(result.kept))

    with tqdm(total=n_cells, desc="Counting cells", unit="cell") as progress:
        if workers <= 1:
            for task in tasks:
                absorb(scan_shard(task))
                progress.update(len(task[2]))
        else:
            with ProcessPoolExecutor(
                max_workers=workers, mp_context=mp.get_context(MP_CONTEXT)
            ) as pool:
                futures = {pool.submit(scan_shard, task): task for task in tasks}
                for future in as_completed(futures):
                    absorb(future.result())
                    progress.update(len(futures[future][2]))

    return totals
