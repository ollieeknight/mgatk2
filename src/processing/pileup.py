"""Streaming per-cell base counting.

One shard = one contiguous slice of the barcode list. A shard worker opens the
BAM itself, streams chrM once, and accumulates directly into dense per-cell
count arrays. Reads are never materialised as Python objects, so worker memory
is O(cells in shard) rather than O(reads in BAM).
"""

from __future__ import annotations

import warnings
from dataclasses import dataclass

import numpy as np
import pysam

from core.config import PipelineConfig

# ASCII -> base index; 4 means "not ACGT" and is discarded.
BASE_LUT = np.full(256, 4, dtype=np.uint8)
for _char, _idx in (("A", 0), ("C", 1), ("G", 2), ("T", 3)):
    BASE_LUT[ord(_char)] = _idx
    BASE_LUT[ord(_char.lower())] = _idx

# Reused so the hot loop never allocates a position vector.
_POSITIONS = np.arange(1 << 16, dtype=np.int32)

CIGAR_MATCH = (0, 7, 8)
CIGAR_REF_ONLY = (2, 3)
CIGAR_QUERY_ONLY = (1, 4)


@dataclass
class ShardResult:
    """Finished counts for one contiguous block of barcodes."""

    offset: int  # first column index in the output matrices
    counts: np.ndarray  # (n_cells, mito_length, 4, 2) uint16
    tn5: np.ndarray  # (n_cells, mito_length, 2) uint16
    depth: np.ndarray  # (n_cells, mito_length) uint16
    base_totals: np.ndarray  # (mito_length, 4) int64, summed over cells
    n_reads: np.ndarray  # (n_cells,) int64, reads kept per cell
    n_paired: np.ndarray  # (n_cells,) int64, of which flagged as paired
    mean_depth: np.ndarray
    median_depth: np.ndarray
    max_depth: np.ndarray
    total_bases: np.ndarray
    coverage_breadth: np.ndarray
    kept: np.ndarray  # (n_cells,) bool, passed min_reads_per_cell
    total_reads: int  # chrM records seen (identical across shards)
    duplicate_reads: int  # duplicates dropped for this shard's cells


def plan_shards(n_cells: int, config: PipelineConfig) -> int:
    """Cells per shard, sized so all concurrent shards fit the memory budget."""
    n_cores = max(1, config.performance.n_cores)
    budget = config.performance.max_memory_gb * 1e9 * 0.6
    affordable = max(1, int(budget / (n_cores * config.bytes_per_cell())))
    even_split = max(1, -(-n_cells // n_cores))
    return min(even_split, affordable)


def scan_shard(task: tuple) -> ShardResult:
    """Stream chrM once and count bases for this shard's barcodes."""
    bam_path, config, barcodes, offset, reference_filename = task
    return _Shard(bam_path, config, barcodes, offset, reference_filename).run()


class _Shard:
    def __init__(self, bam_path, config: PipelineConfig, barcodes, offset, reference_filename):
        self.bam_path = str(bam_path)
        self.config = config
        self.offset = offset
        self.reference_filename = reference_filename
        self.index_of = {bc: i for i, bc in enumerate(barcodes)}
        self.n_cells = len(barcodes)
        self.is_bulk = list(barcodes) == ["bulk"]

        length = config.mito_length
        global _POSITIONS
        if length > _POSITIONS.size:
            _POSITIONS = np.arange(length, dtype=np.int32)
        self.counts = np.zeros((self.n_cells, length, 4, 2), dtype=np.uint32)
        self.tn5 = (
            np.zeros((self.n_cells, length, 2), dtype=np.uint32) if config.compute_tn5 else None
        )
        self.n_reads = np.zeros(self.n_cells, dtype=np.int64)
        self.n_paired = np.zeros(self.n_cells, dtype=np.int64)

    def _open(self):
        if self.bam_path.lower().endswith(".cram"):
            return pysam.AlignmentFile(
                self.bam_path, "rc", reference_filename=self.reference_filename
            )
        # BGZF decompression threads: the shard's wall time is decode-bound.
        return pysam.AlignmentFile(self.bam_path, "rb", threads=2)

    def run(self) -> ShardResult:
        total_reads, duplicates = self._accumulate()
        self._apply_strand_bias()
        return self._summarise(total_reads, duplicates)

    def _accumulate(self) -> tuple[int, int]:
        config = self.config
        quality = config.quality
        length = config.mito_length
        tag = config.barcode_tag
        nh_max = quality.nh_max
        nm_max = quality.nm_max
        min_mapq = quality.min_mapq
        min_baseq = quality.min_baseq
        min_dist = quality.min_distance_from_end
        dedup = not config.dedup.skip
        use_fragment_length = config.dedup.use_fragment_length
        index_of = self.index_of
        is_bulk = self.is_bulk
        counts = self.counts
        tn5 = self.tn5
        n_reads = self.n_reads
        n_paired = self.n_paired
        seen: list[set] = [set() for _ in range(self.n_cells)] if dedup else []

        total_reads = 0
        duplicates = 0
        default_quals = np.full(1024, 60, dtype=np.uint8)

        with self._open() as bam:
            for read in bam.fetch(config.mito_chr):
                total_reads += 1

                if read.is_unmapped or read.is_secondary or read.is_supplementary:
                    continue

                if is_bulk:
                    cell = 0
                else:
                    barcode = read.get_tag(tag) if read.has_tag(tag) else None
                    cell = index_of.get(barcode)
                    if cell is None:
                        continue

                # Filters kept in step with mgatk's filterClipBam.py.
                if nh_max > 0 and read.has_tag("NH") and read.get_tag("NH") > nh_max:
                    continue
                if nm_max > 0:
                    nm = next(
                        (read.get_tag(t) for t in ("NM", "nM") if read.has_tag(t)),
                        None,
                    )
                    if nm is not None and nm > nm_max:
                        continue

                if dedup:
                    key = (read.reference_start << 1) | int(read.is_reverse)
                    if use_fragment_length:
                        key |= abs(read.template_length or 0) << 20
                    cell_seen = seen[cell]
                    if key in cell_seen:
                        duplicates += 1
                        continue
                    cell_seen.add(key)

                if read.mapping_quality < min_mapq:
                    continue

                sequence = read.query_sequence
                if not sequence:
                    continue

                # Counted only once the read is certain to contribute bases, so
                # this stays equal to the Tn5 cut total.
                n_reads[cell] += 1
                n_paired[cell] += read.is_paired

                strand = 1 if read.is_reverse else 0

                if tn5 is not None:
                    # Tn5 inserts at the read's outermost aligned reference base.
                    cut = (read.reference_end - 1) if read.is_reverse else read.reference_start
                    if 0 <= cut < length:
                        tn5[cell, cut, strand] += 1

                read_length = len(sequence)
                bases = BASE_LUT[np.frombuffer(sequence.encode("ascii"), dtype=np.uint8)]
                raw_quals = read.query_qualities
                if raw_quals is None:
                    if read_length > default_quals.size:
                        default_quals = np.full(read_length, 60, dtype=np.uint8)
                    quals = default_quals[:read_length]
                else:
                    quals = np.frombuffer(raw_quals, dtype=np.uint8)

                self._count_read(
                    counts[cell],
                    bases,
                    quals,
                    read.reference_start,
                    read.cigartuples or [],
                    read_length,
                    strand,
                    min_baseq,
                    min_dist,
                    length,
                )

        return total_reads, duplicates

    @staticmethod
    def _count_read(
        cell_counts,
        bases,
        quals,
        reference_start,
        cigar,
        read_length,
        strand,
        min_baseq,
        min_dist,
        mito_length,
    ):
        """Add one read's usable bases into a (mito_length, 4, 2) cell block."""
        # Query positions outside this window are too close to a read end.
        window_start = min_dist
        window_end = read_length - min_dist

        ref_pos = reference_start
        query_pos = 0

        for op, op_length in cigar:
            if op in CIGAR_MATCH:
                q_start = max(query_pos, window_start)
                q_end = min(query_pos + op_length, window_end)
                # Clip the same amount off the reference side to stay in register.
                r_start = ref_pos + (q_start - query_pos)
                r_end = r_start + (q_end - q_start)
                if r_start < 0:
                    q_start -= r_start
                    r_start = 0
                if r_end > mito_length:
                    q_end -= r_end - mito_length
                    r_end = mito_length

                if q_start < q_end:
                    block_bases = bases[q_start:q_end]
                    usable = (quals[q_start:q_end] >= min_baseq) & (block_bases < 4)
                    # Reference positions within a block are unique, so a plain
                    # fancy-index add is exact (no np.add.at needed).
                    cell_counts[_POSITIONS[r_start:r_end][usable], block_bases[usable], strand] += 1

                query_pos += op_length
                ref_pos += op_length
            elif op in CIGAR_REF_ONLY:
                ref_pos += op_length
            elif op in CIGAR_QUERY_ONLY:
                query_pos += op_length

    def _apply_strand_bias(self):
        """Zero any base whose observations come too heavily from one strand."""
        max_bias = self.config.quality.max_strand_bias
        if max_bias >= 1.0:
            return  # a ratio can never exceed 1.0, so the filter is a no-op

        # Same metric `paired` applies: |forward - reverse| / total.
        totals = self.counts.sum(axis=3, keepdims=True)
        imbalance = np.abs(np.diff(self.counts.astype(np.int64), axis=3))
        biased = (imbalance > max_bias * totals) & (totals > 0)
        self.counts[np.broadcast_to(biased, self.counts.shape)] = 0

    def _summarise(self, total_reads: int, duplicates: int) -> ShardResult:
        length = self.config.mito_length
        kept = self.n_reads >= max(1, self.config.min_reads_per_cell)
        if not kept.all():
            self.counts[~kept] = 0
            if self.tn5 is not None:
                self.tn5[~kept] = 0

        depth = self.counts.sum(axis=(2, 3), dtype=np.int64)
        positions_covered = (depth > 0).sum(axis=1)
        total_bases = depth.sum(axis=1)

        # Depth statistics ignore uncovered positions, matching mgatk2 <= 1.2.
        covered_only = depth.astype(np.float32)
        covered_only[depth == 0] = np.nan
        with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
            # An entirely uncovered cell is expected, not exceptional.
            warnings.simplefilter("ignore", RuntimeWarning)
            mean_depth = np.where(positions_covered > 0, total_bases / positions_covered, 0.0)
            median_depth = np.nan_to_num(np.nanmedian(covered_only, axis=1))
        del covered_only

        base_totals = self.counts.sum(axis=(0, 3), dtype=np.int64)
        np.minimum(self.counts, 65535, out=self.counts)
        if self.tn5 is None:
            tn5 = np.zeros((self.n_cells, length, 2), dtype=np.uint16)
        else:
            np.minimum(self.tn5, 65535, out=self.tn5)
            tn5 = self.tn5.astype(np.uint16)

        return ShardResult(
            offset=self.offset,
            counts=self.counts.astype(np.uint16),
            tn5=tn5,
            depth=np.minimum(depth, 65535).astype(np.uint16),
            base_totals=base_totals,
            n_reads=self.n_reads,
            n_paired=self.n_paired,
            mean_depth=mean_depth,
            median_depth=median_depth,
            max_depth=np.minimum(depth.max(axis=1), 65535).astype(np.uint16),
            total_bases=total_bases,
            coverage_breadth=positions_covered / length,
            kept=kept,
            total_reads=total_reads,
            duplicate_reads=duplicates,
        )
