"""Configuration classes."""

from dataclasses import dataclass
from pathlib import Path

import numpy


@dataclass
class QualityThresholds:
    """Quality filtering parameters"""

    min_baseq: int = 20
    min_mapq: int = 30
    max_strand_bias: float = 1.0
    min_distance_from_end: int = 5
    nh_max: int = 0  # 0 = disabled; 1 matches mgatk tenx default
    nm_max: int = 0  # 0 = disabled; 4 matches mgatk tenx default


@dataclass
class DeduplicationConfig:
    """Deduplication parameters."""

    skip: bool = False
    use_fragment_length: bool = True


@dataclass
class PerformanceConfig:
    """Resource management."""

    n_cores: int = 8  # CPU cores for parallel processing
    worker_batch_size: int = 8  # Cells per parallel batch
    io_batch_size: int = 100  # Cells before HDF5 flush (dynamic: 10% of barcodes)
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
    query_name: str = ""
    reference_end: int = 0
    is_read1: bool = False
    is_read2: bool = False
    is_qcfail: bool = False
    is_duplicate: bool = False
    read_group: str | None = None

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


@dataclass
class PairedConfig:
    """Configuration for query/baseline mitochondrial evidence analysis."""

    query: str
    baseline: str
    reference: str
    output: str
    sample_name: str
    min_baseq: int = 20
    min_mapq: int = 20
    min_distance_from_end: int = 5
    mito_chr: str = "chrM"
    deduplication: str = "alignment_and_fragment_length"
    min_query_depth: int = 10
    min_baseline_depth: int = 5
    min_alt_observations: int = 3
    min_query_af: float = 0.005
    max_baseline_af: float = 0.01
    min_query_baseline_ratio: float = 3.0
    max_strand_bias: float = 0.9
    custom_blacklist: str | None = None
    input_is_consensus: bool = False
    shifted_reference_supplied: bool = False
    circular_edge_bases: int = 500
    evidence_schema_version: str = "1.0"
    candidate_schema_version: str = "1.0"
    qc_schema_version: str = "1.0"
    write_legacy_tsv: bool = True

    def __post_init__(self) -> None:
        for name in (
            "min_baseq",
            "min_mapq",
            "min_distance_from_end",
            "min_query_depth",
            "min_baseline_depth",
            "min_alt_observations",
            "circular_edge_bases",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")
        for name in ("min_query_af", "max_baseline_af", "max_strand_bias"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.min_query_baseline_ratio <= 0:
            raise ValueError("min_query_baseline_ratio must be positive")
        if self.deduplication not in {
            "alignment_and_fragment_length",
            "alignment_start",
            "none",
        }:
            raise ValueError(f"Unsupported deduplication mode: {self.deduplication}")
        if Path(self.query).resolve() == Path(self.baseline).resolve():
            raise ValueError("query and baseline must be different files")
        if not self.sample_name or any(c in self.sample_name for c in "/\\"):
            raise ValueError("sample_name must be a non-empty filename prefix")
        if self.input_is_consensus and self.deduplication != "none":
            raise ValueError("consensus inputs require --deduplication none")


class PipelineConfig:
    """Pipeline configuration."""

    def __init__(
        self,
        min_baseq: int = 20,
        min_mapq: int = 30,
        max_strand_bias: float = 0.9,
        min_distance_from_end: int = 5,
        skip_deduplication: bool = False,
        use_fragment_length_dedup: bool = True,
        n_cores: int = 8,
        worker_batch_size: int | None = None,
        io_batch_size: int | None = None,
        max_memory_gb: float = 128.0,
        sequential: bool = False,
        min_reads_per_cell: int = 1,
        barcode_tag: str = "CB",
        mito_chr: str = "chrM",
        mito_length: int = 16569,
        nh_max: int = 0,
        nm_max: int = 0,
        pileup_mode: str = "fast",
        compute_tn5: bool = True,
    ):
        self.quality = QualityThresholds(
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            max_strand_bias=max_strand_bias,
            min_distance_from_end=min_distance_from_end,
            nh_max=nh_max,
            nm_max=nm_max,
        )
        self.dedup = DeduplicationConfig(
            skip=skip_deduplication, use_fragment_length=use_fragment_length_dedup
        )
        self.performance = PerformanceConfig(
            n_cores=n_cores,
            worker_batch_size=worker_batch_size or n_cores,
            io_batch_size=io_batch_size or 100,
            max_memory_gb=max_memory_gb,
            sequential=sequential,
        )
        self.min_reads_per_cell = min_reads_per_cell
        self.barcode_tag = barcode_tag
        self.mito_chr = mito_chr
        self.mito_length = mito_length
        self.pileup_mode = pileup_mode  # "classic" | "fast"
        self.compute_tn5 = compute_tn5
