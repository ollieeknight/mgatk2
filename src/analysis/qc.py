"""Quality control metrics"""

import datetime
import logging
from importlib.metadata import version

from core.config import PipelineConfig

logger = logging.getLogger(__name__)


class QCCalculator:
    """Calculate quality control metrics for mtDNA sequencing"""

    def __init__(self, config: PipelineConfig):
        self.config = config

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
