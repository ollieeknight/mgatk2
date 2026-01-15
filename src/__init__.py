"""mgatk2: Mitochondrial genome analysis toolkit.

For 10x scATAC- and scASAP-seq data.
"""

import logging

from core.config import (
    DeduplicationConfig,
    PerformanceConfig,
    PipelineConfig,
    QualityThresholds,
    SimpleRead,
)
from core.exceptions import (
    BAMReadError,
    InsufficientDataError,
    InvalidInputError,
    MgatkError,
    ProcessingError,
)
from core.pipeline import MtDNAPipeline, run_pipeline
from processing.pileup import PileupGenerator
from processing.processors import CellProcessor, process_barcode_worker
from processing.readers import BAMReader
from utils.utils import (
    load_singlecell_csv,
    validate_bam_file,
    validate_barcode_file,
)

__version__ = "1.0.0"

# Configure logging
logger = logging.getLogger(__name__)
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

# Export main classes and functions
__all__ = [
    # Configuration
    "SimpleRead",
    "PipelineConfig",
    "QualityThresholds",
    "DeduplicationConfig",
    "PerformanceConfig",
    # Exceptions
    "MgatkError",
    "InvalidInputError",
    "ProcessingError",
    "BAMReadError",
    "InsufficientDataError",
    # Core functionality
    "PileupGenerator",
    "MtDNAPipeline",
    "run_pipeline",
    "BAMReader",
    "CellProcessor",
    "process_barcode_worker",
    # Utilities
    "load_singlecell_csv",
    "validate_bam_file",
    "validate_barcode_file",
]
