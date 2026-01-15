"""Cell processing for mgatk2."""

import gc
import logging
import multiprocessing as mp
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from tqdm import tqdm

from processing.pileup import PileupGenerator

logger = logging.getLogger(__name__)

# Use spawn context on HPC/cluster environments for better stability
# Fork can have issues with certain file handles and threading on HPC
MP_CONTEXT = "spawn"


def process_barcode_worker(args):
    barcode, reads, config = args
    if not reads or len(reads) < config.min_reads_per_cell:
        return None

    try:
        pileup_gen = PileupGenerator(config)
        pileup = pileup_gen.generate_pileup(reads)
        pileup = pileup_gen.filter_strand_bias(pileup)

        if not pileup:
            return None

        n_reads = len(reads)
        n_paired = sum(1 for r in reads if r.is_paired)
        n_proper = sum(1 for r in reads if r.is_proper_pair)
        fwd = sum(1 for r in reads if not r.is_reverse)

        depths = np.array([c["depth"] for c in pileup.values()])
        mean_cov = depths.mean()

        return {
            "barcode": barcode,
            "pileup": pileup,
            "n_reads": n_reads,
            "qc": {
                "barcode": barcode,
                "total_reads": n_reads,
                "total_fragments": n_reads // 2 if n_paired > 0 else n_reads,
                "proper_pairs": n_proper,
                "proper_pair_pct": 100.0 * n_proper / n_reads if n_reads > 0 else 0,
                "mean_coverage": mean_cov,
                "median_coverage": np.median(depths),
                "coverage_std": depths.std(),
                "coverage_cv": depths.std() / mean_cov if mean_cov > 0 else 0,
                "positions_covered": len(depths),
                "positions_10x": (depths >= 10).sum(),
                "positions_50x": (depths >= 50).sum(),
                "mean_mapq": np.mean([r.mapping_quality for r in reads]),
                "forward_reads": fwd,
                "reverse_reads": n_reads - fwd,
                "strand_ratio": fwd / (n_reads - fwd) if fwd < n_reads else 0,
            },
        }
    except Exception as e:
        logger.error(f"Error processing {barcode}: {e}")
        return None


class CellProcessor:
    def __init__(self, config, output_dir):
        self.config = config
        self.output_dir = output_dir

    def process_cells_direct(self, reads_by_barcode, incremental_writer=None):
        barcodes = list(reads_by_barcode.keys())
        results = []
        failed = 0

        with tqdm(total=len(barcodes), desc="Processing cells", unit="cell") as pbar:
            for bc in barcodes:
                result = process_barcode_worker((bc, reads_by_barcode.pop(bc), self.config))

                if result:
                    if incremental_writer:
                        incremental_writer.write_cell(result)
                        results.append({"barcode": result["barcode"], "n_reads": result["n_reads"]})
                    else:
                        results.append(result)
                else:
                    failed += 1
                pbar.update(1)

        if failed > 0:
            logger.warning(f"{failed} cells failed")
        gc.collect()
        return results

    def process_cells_progressive(self, reads_by_barcode, incremental_writer=None):
        barcodes = list(reads_by_barcode.keys())
        n_cells = len(barcodes)
        total_reads = sum(len(r) for r in reads_by_barcode.values())
        avg = total_reads / n_cells if n_cells > 0 else 0

        logger.info(f"Processing {n_cells} cells, {avg:.0f} avg reads/cell")

        # Check if sequential processing is forced via config
        if self.config.performance.sequential:
            logger.info("Sequential processing enabled via config")
            return self.process_cells_direct(reads_by_barcode, incremental_writer)

        # For very large datasets, use sequential to avoid memory issues
        if total_reads > 50_000_000:
            logger.info("Large dataset (>50M reads), using sequential processing")
            return self.process_cells_direct(reads_by_barcode, incremental_writer)

        logger.info(f"Using parallel processing with {self.config.performance.n_cores} cores")
        return self._process_parallel(
            reads_by_barcode, barcodes, self.config.performance.batch_size, incremental_writer
        )

    def _process_parallel(self, reads_by_barcode, barcodes, batch_size, incremental_writer=None):
        n_cells = len(barcodes)
        results = []
        n_batches = (n_cells + batch_size - 1) // batch_size
        
        logger.info(f"Starting parallel processing: {n_batches} batches of ~{batch_size} cells each")

        with ProcessPoolExecutor(
            max_workers=self.config.performance.n_cores, mp_context=mp.get_context(MP_CONTEXT)
        ) as executor:
            with tqdm(total=n_cells, desc="Processing cells", unit="cell") as pbar:
                for i in range(n_batches):
                    start = i * batch_size
                    end = min(start + batch_size, n_cells)
                    batch_bcs = barcodes[start:end]
                    args = [(bc, reads_by_barcode.pop(bc), self.config) for bc in batch_bcs]

                    try:
                        futures = [executor.submit(process_barcode_worker, a) for a in args]
                        batch_results = [f.result() for f in as_completed(futures) if f.result()]
                        results.extend(batch_results)

                        if incremental_writer:
                            for r in batch_results:
                                incremental_writer.write_cell(r)
                    except Exception as e:
                        logger.error(f"Batch {i + 1} failed: {e}, falling back to sequential")
                        pbar.close()
                        return self.process_cells_direct(reads_by_barcode, incremental_writer)

                    pbar.update(len(args))
                    gc.collect()

        logger.info(f"Processed {len(results):,}/{n_cells:,} cells")
        return results
