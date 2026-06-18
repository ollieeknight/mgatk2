"""Configuration classes."""

from dataclasses import dataclass

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
class WesConfig:
    """Configuration for mgatk2 wes somatic mitochondrial variant calling."""

    # Quality thresholds
    min_baseq: int = 20
    min_mapq: int = 20  # lower than sc — mito reads multi-map more in WES
    min_distance_from_end: int = 5

    # Mito genome
    mito_chr: str = "chrM"
    mito_length: int = 16569

    # Depth thresholds
    min_tumour_depth: int = 10
    min_normal_depth: int = 5

    # Variant AF thresholds
    min_tumour_af: float = 0.005  # 0.5% — low heteroplasmy detection
    max_normal_af: float = 0.10  # 10% — germline filter
    min_tn_ratio: float = 3.0  # tumour VAF / normal VAF

    # Strand bias for bulk (|fwd-rev| / (fwd+rev))
    max_strand_bias: float = 0.9

    # Minimum alt reads in tumour (absolute count guard)
    min_tumour_alt_reads: int = 3

    # Custom chrM-side blacklist BED (path to file). The bundled NuMT BEDs are
    # nuclear positions (used for reference masking) and contain no chrM rows.
    # Primary NuMT defence is min_mapq=20; supply a custom BED for chrM-side exclusions.
    blacklist_build: str = "none"  # hg38 | hg19 | mm10 | mm9 | none
    custom_blacklist: str | None = None  # path to custom chrM-side BED

    def to_pipeline_config(self) -> "PipelineConfig":
        """Create a PipelineConfig for use with PileupGenerator."""
        return PipelineConfig(
            min_baseq=self.min_baseq,
            min_mapq=self.min_mapq,
            min_distance_from_end=self.min_distance_from_end,
            skip_deduplication=True,  # CLONK-WES CRAMs are already fgbio-deduplicated
            mito_chr=self.mito_chr,
            mito_length=self.mito_length,
        )


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
        **kwargs,
    ):
        self.quality = QualityThresholds(
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            max_strand_bias=max_strand_bias,
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
