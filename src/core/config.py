"""Configuration classes."""

from dataclasses import dataclass

import numpy


@dataclass
class QualityThresholds:
    """Quality filtering parameters."""

    min_baseq: int = 20
    min_mapq: int = 30
    max_strand_bias: float = 1.0
    min_distance_from_end: int = 5


@dataclass
class DeduplicationConfig:
    """Deduplication parameters."""

    skip: bool = False
    use_fragment_length: bool = True


@dataclass
class PerformanceConfig:
    """Resource management."""

    n_cores: int = 8
    batch_size: int = 8
    max_memory_gb: float = 128.0
    sequential: bool = False


@dataclass
class SimpleRead:
    """Lightweight BAM read"""

    reference_start: int
    is_reverse: bool
    mapping_quality: int
    query_sequence: bytes
    query_qualities: "numpy.ndarray"
    cigar: list[tuple[int, int]]
    is_proper_pair: bool = False
    is_paired: bool = False
    template_length: int = 0

    def get_aligned_pairs(self) -> list[tuple[int, int]]:
        """Get aligned (query_pos, ref_pos) pairs."""
        pairs = []
        ref_pos = self.reference_start
        query_pos = 0

        for op, length in self.cigar:
            if op == 0:
                for _ in range(length):
                    pairs.append((query_pos, ref_pos))
                    query_pos += 1
                    ref_pos += 1
            elif op == 1:
                query_pos += length
            elif op in [2, 3]:
                ref_pos += length
            elif op == 4:
                query_pos += length
            elif op in [7, 8]:
                for _ in range(length):
                    pairs.append((query_pos, ref_pos))
                    query_pos += 1
                    ref_pos += 1
        return pairs


class PipelineConfig:
    """Pipeline configuration."""

    def __init__(
        self,
        min_baseq: int = 20,
        min_mapq: int = 30,
        max_strand_bias: float = 0.9,
        skip_deduplication: bool = False,
        use_fragment_length_dedup: bool = True,
        n_cores: int = 8,
        batch_size: int | None = None,
        max_memory_gb: float = 128.0,
        sequential: bool = False,
        min_reads_per_cell: int = 1,
        barcode_tag: str = "CB",
        mito_chr: str = "chrM",
        mito_length: int = 16569,
        **kwargs,
    ):
        self.quality = QualityThresholds(
            min_baseq=min_baseq, min_mapq=min_mapq, max_strand_bias=max_strand_bias
        )
        self.dedup = DeduplicationConfig(
            skip=skip_deduplication, use_fragment_length=use_fragment_length_dedup
        )
        self.performance = PerformanceConfig(
            n_cores=n_cores,
            batch_size=batch_size or n_cores,
            max_memory_gb=max_memory_gb,
            sequential=sequential,
        )
        self.min_reads_per_cell = min_reads_per_cell
        self.barcode_tag = barcode_tag
        self.mito_chr = mito_chr
        self.mito_length = mito_length
