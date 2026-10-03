"""Configuration classes."""

from dataclasses import dataclass
from pathlib import Path


@dataclass(slots=True)
class SimpleRead:
    """Lightweight BAM read. Used by the paired/bulk fragment path only."""

    reference_start: int
    is_reverse: bool
    mapping_quality: int
    query_sequence: bytes
    query_qualities: bytes
    cigar: list[tuple[int, int]]
    is_proper_pair: bool = False
    is_paired: bool = False
    template_length: int = 0
    query_name: str = ""
    is_read1: bool = False
    is_read2: bool = False


@dataclass
class PairedConfig:
    """Configuration for tumour/normal mitochondrial evidence analysis."""

    tumor: str
    normal: str
    reference: str
    output: str
    sample_name: str
    min_baseq: int = 20
    min_mapq: int = 20
    min_distance_from_end: int = 5
    mito_chr: str = "chrM"
    deduplication: str = "alignment_and_fragment_length"
    min_tumor_depth: int = 10
    min_normal_depth: int = 5
    min_alt_observations: int = 3
    min_tumor_af: float = 0.005
    max_normal_af: float = 0.01
    max_strand_bias: float = 0.9
    custom_blacklist: str | None = None
    autosomal_median_depth: float | None = None
    input_is_consensus: bool = False
    circular_edge_bases: int = 500

    def __post_init__(self) -> None:
        for name in (
            "min_baseq",
            "min_mapq",
            "min_distance_from_end",
            "min_tumor_depth",
            "min_normal_depth",
            "min_alt_observations",
            "circular_edge_bases",
        ):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be non-negative")
        for name in ("min_tumor_af", "max_normal_af", "max_strand_bias"):
            if not 0 <= getattr(self, name) <= 1:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.deduplication not in {
            "alignment_and_fragment_length",
            "alignment_start",
            "none",
        }:
            raise ValueError(f"Unsupported deduplication mode: {self.deduplication}")
        if Path(self.tumor).resolve() == Path(self.normal).resolve():
            raise ValueError("tumor and normal must be different files")
        if not self.sample_name or any(c in self.sample_name for c in "/\\"):
            raise ValueError("sample_name must be a non-empty filename prefix")
        if self.input_is_consensus and self.deduplication != "none":
            raise ValueError("consensus inputs require --deduplication none")
        if self.autosomal_median_depth is not None and self.autosomal_median_depth < 0:
            raise ValueError("autosomal_median_depth must be non-negative")


@dataclass
class PipelineConfig:
    """Single-cell counting configuration."""

    min_baseq: int = 20
    min_mapq: int = 30
    max_strand_bias: float = 1.0
    min_distance_from_end: int = 5
    nh_max: int = 0  # 0 disables
    nm_max: int = 0  # 0 disables
    skip_deduplication: bool = False
    use_fragment_length_dedup: bool = True
    n_cores: int = 8
    max_memory_gb: float = 128.0
    min_reads_per_cell: int = 1
    barcode_tag: str = "CB"
    mito_chr: str = "chrM"
    mito_length: int = 16569
    compute_tn5: bool = True
    # 1-based amplicon panel positions; coverage breadth is scoped to these.
    panel_positions: frozenset[int] | None = None

    def bytes_per_cell(self) -> int:
        # uint32 base counts (4 bases x 2 strands) plus uint32 Tn5 cuts (2 strands).
        return self.mito_length * (4 * 2 + 2) * 4
