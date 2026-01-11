"""Quality control metrics"""

import datetime
import logging
from collections import defaultdict
from importlib.metadata import version

import numpy as np

from core.config import PipelineConfig

logger = logging.getLogger(__name__)


class QCCalculator:
    """Calculate quality control metrics for mtDNA sequencing"""

    def __init__(self, config: PipelineConfig):
        self.config = config

    def calculate_position_stats(self, cell_results: list[dict]) -> dict[int, dict]:
        # Collect depths for each position across all cells
        position_data: dict[int, dict] = defaultdict(lambda: {"depths": [], "n_cells": 0})

        # Collect data for each position
        for result in cell_results:
            pileup = result["pileup"]
            for pos, counts in pileup.items():
                pos_1based = pos + 1
                depth = counts["depth"]
                if depth > 0:
                    position_data[pos_1based]["depths"].append(depth)
                    position_data[pos_1based]["n_cells"] += 1

        # Calculate statistics for each position
        stats = {}
        for pos in range(1, self.config.mito_length + 1):
            if pos in position_data:
                depths = position_data[pos]["depths"]
                n_cells = position_data[pos]["n_cells"]
                n_cells_10x = sum(1 for d in depths if d >= 10)
                n_cells_50x = sum(1 for d in depths if d >= 50)
                mean_cov = np.mean(depths)
                median_cov = np.median(depths)
                std_cov = np.std(depths)
                cv = std_cov / mean_cov if mean_cov > 0 else 0

                stats[pos] = {
                    "n_cells_any": n_cells,
                    "n_cells_10x": n_cells_10x,
                    "n_cells_50x": n_cells_50x,
                    "mean_cov": mean_cov,
                    "median_cov": median_cov,
                    "cv": cv,
                }
            else:
                stats[pos] = {
                    "n_cells_any": 0,
                    "n_cells_10x": 0,
                    "n_cells_50x": 0,
                    "mean_cov": 0,
                    "median_cov": 0,
                    "cv": 0,
                }

        return stats

    def collect_run_metadata(
        self, bam_path: str, output_dir: str, n_cells_input: int, n_cells_passed: int
    ) -> dict:
        """Collect run metadata and parameters"""
        return {
            "mgatk_version": version("mgatk2"),
            "run_date": datetime.datetime.now().isoformat(),
            "input_bam": str(bam_path),
            "output_dir": str(output_dir),
            "reference": self.config.mito_chr,
            "reference_length": self.config.mito_length,
            "cells_total": n_cells_input,
            "cells_passed_qc": n_cells_passed,
            "cells_failed_qc": n_cells_input - n_cells_passed,
            "parameters": {
                "min_base_quality": self.config.quality.min_baseq,
                "min_mapping_quality": self.config.quality.min_mapq,
                "min_reads_per_cell": self.config.min_reads_per_cell,
                "max_strand_bias": self.config.quality.max_strand_bias,
                "skip_deduplication": self.config.dedup.skip,
                "use_fragment_length_dedup": self.config.dedup.use_fragment_length,
                "barcode_tag": self.config.barcode_tag,
                "mito_chr": self.config.mito_chr,
                "n_cores": self.config.performance.n_cores,
            },
        }
