"""Core configuration and exceptions."""

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

__all__ = [
    "PipelineConfig",
    "QualityThresholds",
    "DeduplicationConfig",
    "PerformanceConfig",
    "SimpleRead",
    "MgatkError",
    "InvalidInputError",
    "ProcessingError",
    "BAMReadError",
    "NoChrMReadsError",
    "NoBarcodeTagsError",
    "BAMFormatError",
]
