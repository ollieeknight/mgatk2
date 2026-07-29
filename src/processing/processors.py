"""Cell processing for mgatk2."""

import gc
import logging
import multiprocessing as mp
import platform
from concurrent.futures import ProcessPoolExecutor, as_completed

import numpy as np
from tqdm import tqdm

from processing.pileup import PileupGenerator

logger = logging.getLogger(__name__)

# fork on Linux: copy-on-write, no re-import cost, no pickle roundtrip for reads.
# spawn on macOS: fork is unsafe with Objective-C runtime (Python 3.12+).
MP_CONTEXT = "fork" if platform.system() == "Linux" else "spawn"


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

        depths = np.array([c["depth"] for c in pileup.values()])
        mean_cov = depths.mean()
        positions_covered = len(depths)
        coverage_breadth = positions_covered / config.mito_length if config.mito_length > 0 else 0

        return {
            "barcode": barcode,
            "pileup": pileup,
            "n_reads": n_reads,
            "qc": {
                "barcode": barcode,
                "total_reads": n_reads,
                "total_fragments": n_reads - n_paired // 2,
                "mean_depth": mean_cov,
                "coverage_breadth": coverage_breadth,
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

        logger.info(f"Processing {n_cells} cells at an average of {avg:.0f} reads/cell")

        if self.config.performance.sequential:
            logger.info("Sequential processing enabled via config")
            return self.process_cells_direct(reads_by_barcode, incremental_writer)

        # Sequential fallback: threshold on avg bases/cell, not reads/cell.
        # Scales correctly for long reads (1000bp) vs short reads (100bp).
        sample_cells = [reads_by_barcode[b] for b in barcodes[:20] if reads_by_barcode[b]]
        if sample_cells:
            sample_reads = [r for cell in sample_cells for r in cell]
            avg_read_len = (
                sum(len(r.query_sequence) for r in sample_reads) / len(sample_reads)
                if sample_reads
                else 100
            )
        else:
            avg_read_len = 100
        avg_bases = avg * avg_read_len
        if avg_bases > 250_000:  # ~2500 reads × 100bp; scales for long reads
            logger.info(
                f"High bases/cell ({avg_bases:.0f} > 250k, ~{avg:.0f} reads × {avg_read_len:.0f}bp), "
                "using sequential processing"
            )
            return self.process_cells_direct(reads_by_barcode, incremental_writer)

        return self._process_parallel(
            reads_by_barcode,
            barcodes,
            self.config.performance.worker_batch_size,
            incremental_writer,
        )

    def _process_parallel(self, reads_by_barcode, barcodes, batch_size, incremental_writer=None):
        n_cells = len(barcodes)
        results = []
        n_batches = (n_cells + batch_size - 1) // batch_size

        with ProcessPoolExecutor(
            max_workers=self.config.performance.n_cores, mp_context=mp.get_context(MP_CONTEXT)
        ) as executor:
            with tqdm(total=n_cells, desc="Processing cells", unit="cell") as pbar:
                for i in range(n_batches):
                    start = i * batch_size
                    end = min(start + batch_size, n_cells)
                    batch_bcs = barcodes[start:end]
                    args = [(bc, reads_by_barcode.pop(bc), self.config) for bc in batch_bcs]

                    futures = {}
                    batch_results = []
                    for arg in args:
                        try:
                            futures[executor.submit(process_barcode_worker, arg)] = arg
                        except Exception as exc:
                            logger.error(
                                "Worker submission failed; processing cell directly: %s", exc
                            )
                            result = process_barcode_worker(arg)
                            if result:
                                batch_results.append(result)

                    for future in as_completed(futures):
                        try:
                            result = future.result()
                        except Exception as exc:
                            logger.error("Worker failed; processing cell directly: %s", exc)
                            result = process_barcode_worker(futures[future])
                        if result:
                            batch_results.append(result)

                    results.extend(batch_results)
                    if incremental_writer:
                        for result in batch_results:
                            incremental_writer.write_cell(result)

                    pbar.update(len(args))
                    gc.collect()

        logger.info(f"Processed {len(results):,}/{n_cells:,} cells")
        return results
