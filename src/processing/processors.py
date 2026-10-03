"""Shard scheduling for mgatk2."""

import logging
import multiprocessing as mp
import platform
from collections import deque
from concurrent.futures import ProcessPoolExecutor

import numpy as np
from tqdm import tqdm

from processing.pileup import plan_shards, scan_shard

logger = logging.getLogger(__name__)

# fork on Linux: cheap, and the parent heap is small by design here.
# spawn on macOS: fork is unsafe with the Objective-C runtime (Python 3.12+).
MP_CONTEXT = "fork" if platform.system() == "Linux" else "spawn"


def process_shards(bam_path, config, barcodes, writer) -> dict:
    """Scan chrM once per shard, writing each finished shard straight to disk."""
    # Contiguous shards, so each writes one chunk-aligned block of HDF5 columns.
    per_shard = plan_shards(len(barcodes), config)
    tasks = [
        (str(bam_path), config, barcodes[lo : lo + per_shard], lo)
        for lo in range(0, len(barcodes), per_shard)
    ]
    n_cells = len(barcodes)
    workers = min(config.n_cores, len(tasks))

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
            return totals

        # Shards are absorbed in barcode order so every output is reproducible.
        # Submitting at most `workers` ahead keeps finished shards waiting on a
        # slow predecessor from piling up in memory.
        with ProcessPoolExecutor(
            max_workers=workers, mp_context=mp.get_context(MP_CONTEXT)
        ) as pool:
            pending = deque()
            for task in tasks:
                pending.append((pool.submit(scan_shard, task), len(task[2])))
                if len(pending) == workers:
                    future, n = pending.popleft()
                    absorb(future.result())
                    progress.update(n)
            for future, n in pending:
                absorb(future.result())
                progress.update(n)

    return totals
