"""Core pipeline components."""

from .config import (
    DeduplicationConfig,
    PerformanceConfig,
    PipelineConfig,
    QualityThresholds,
    SimpleRead,
)
from .exceptions import (
    BAMFormatError,
    BAMReadError,
    InvalidInputError,
    MgatkError,
    NoBarcodeTagsError,
    NoChrMReadsError,
    ProcessingError,
)
from .pipeline import MtDNAPipeline, run_pipeline

__all__ = [
    # From config
    "PipelineConfig",
    "QualityThresholds",
    "DeduplicationConfig",
    "PerformanceConfig",
    "SimpleRead",
    # From pipeline
    "run_pipeline",
    "MtDNAPipeline",
    # From exceptions
    "MgatkError",
    "InvalidInputError",
    "ProcessingError",
    "BAMReadError",
    "NoChrMReadsError",
    "NoBarcodeTagsError",
    "BAMFormatError",
]
