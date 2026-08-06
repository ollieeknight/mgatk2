"""Per-allele quality histograms and the statistics derived from them.

Storing a histogram per (position, allele) rather than a running sum is what
makes median and rank-sum statistics possible at all: a mean pooled over
reference and alternate observations, which is what mgatk2 reported before
schema 2.0, cannot separate a real allele from an artefact.
"""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

# Bin counts are chosen to cover the realistic range of each metric exactly.
BASEQ_BINS = 96  # Illumina tops out near 42; UMI-consensus qualities reach the 90s
MAPQ_BINS = 64  # bwa/bowtie2 cap at 60
DISTANCE_BINS = 64  # bases from the nearest read end, in DISTANCE_SCALE-wide bins
DISTANCE_SCALE = 2

BASE_INDEX = {"A": 0, "C": 1, "G": 2, "T": 3}


class QualityHistograms:
    """Per-position, per-allele counts and quality distributions for one sample."""

    def __init__(self, length: int):
        self.length = length
        self.counts = np.zeros((length, 4, 2), dtype=np.int64)  # allele x strand
        self.orientation = np.zeros((length, 4, 2), dtype=np.int64)  # allele x F1R2/F2R1
        self.baseq = np.zeros((length, 4, BASEQ_BINS), dtype=np.int32)
        self.mapq = np.zeros((length, 4, MAPQ_BINS), dtype=np.int32)
        self.distance = np.zeros((length, 4, DISTANCE_BINS), dtype=np.int32)
        self.clipped = np.zeros(length, dtype=np.int64)
        self.overlap_disagreements = np.zeros(length, dtype=np.int64)
        self._buffer: list[list[int]] = [[] for _ in range(7)]
        self._buffered = 0

    # Observations arrive one at a time from a Python loop; scalar NumPy
    # indexing there would dominate runtime, so buffer and apply in blocks.
    FLUSH_EVERY = 1_000_000

    def add(
        self,
        position: int,
        allele: int,
        strand: int,
        base_quality: int,
        mapping_quality: int,
        distance_from_end: int,
        orientation: int,
    ) -> None:
        buffer = self._buffer
        buffer[0].append(position)
        buffer[1].append(allele)
        buffer[2].append(strand)
        buffer[3].append(base_quality)
        buffer[4].append(mapping_quality)
        buffer[5].append(distance_from_end)
        buffer[6].append(orientation)
        self._buffered += 1
        if self._buffered >= self.FLUSH_EVERY:
            self.flush()

    def flush(self) -> None:
        if not self._buffered:
            return
        position, allele, strand, baseq, mapq, distance, orientation = (
            np.asarray(values, dtype=np.int64) for values in self._buffer
        )
        cell = position * 4 + allele

        np.add.at(self.counts.reshape(-1), cell * 2 + strand, 1)
        oriented = orientation >= 0
        if oriented.any():
            np.add.at(self.orientation.reshape(-1), cell[oriented] * 2 + orientation[oriented], 1)
        np.add.at(self.baseq.reshape(-1), cell * BASEQ_BINS + np.minimum(baseq, BASEQ_BINS - 1), 1)
        np.add.at(self.mapq.reshape(-1), cell * MAPQ_BINS + np.minimum(mapq, MAPQ_BINS - 1), 1)
        np.add.at(
            self.distance.reshape(-1),
            cell * DISTANCE_BINS + np.minimum(distance // DISTANCE_SCALE, DISTANCE_BINS - 1),
            1,
        )

        self._buffer = [[] for _ in range(7)]
        self._buffered = 0

    def depth(self) -> np.ndarray:
        return self.counts.sum(axis=(1, 2))

    def allele_counts(self) -> np.ndarray:
        """(length, 4) observations per allele, both strands."""
        return self.counts.sum(axis=2)


def histogram_median(histogram: np.ndarray, scale: int = 1) -> float:
    """Lower median of the distribution a histogram represents."""
    total = int(histogram.sum())
    if total == 0:
        return 0.0
    cumulative = np.cumsum(histogram)
    return float(int(np.searchsorted(cumulative, (total + 1) // 2)) * scale)


def rank_sum(alternate: np.ndarray, reference: np.ndarray) -> tuple[float, float]:
    """Mann-Whitney z-score and two-sided p-value for alternate versus reference.

    A large |z| means the alternate observations sit systematically higher or
    lower on this metric than the reference ones at the same position, which is
    the signature of an artefact rather than a real allele.
    """
    n_alternate = int(alternate.sum())
    n_reference = int(reference.sum())
    if n_alternate == 0 or n_reference == 0:
        return 0.0, 1.0

    combined = alternate.astype(np.float64) + reference
    below = np.cumsum(combined) - combined
    mid_ranks = below + (combined + 1) / 2

    rank_total = float((alternate * mid_ranks).sum())
    u_statistic = rank_total - n_alternate * (n_alternate + 1) / 2
    total = n_alternate + n_reference
    expected = n_alternate * n_reference / 2

    ties = float((combined**3 - combined).sum())
    if total < 2:
        return 0.0, 1.0
    variance = n_alternate * n_reference / 12 * ((total + 1) - ties / (total * (total - 1)))
    if variance <= 0:
        return 0.0, 1.0

    z_score = (u_statistic - expected) / np.sqrt(variance)
    return float(z_score), float(2 * norm.sf(abs(z_score)))
