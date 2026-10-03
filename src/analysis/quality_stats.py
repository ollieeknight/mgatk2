"""Per-allele quality histograms and the statistics derived from them.

Storing a histogram per (position, allele) rather than a running sum is what
makes median and rank-sum statistics possible at all: a mean pooled over
reference and alternate observations, which is what mgatk2 reported before
v1.3, cannot separate a real allele from an artefact.
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
        self.counts = np.zeros((length, 4, 2), dtype=np.int64)  # allele x strand
        self.orientation = np.zeros((length, 4, 2), dtype=np.int64)  # allele x F1R2/F2R1
        self.baseq = np.zeros((length, 4, BASEQ_BINS), dtype=np.int32)
        self.mapq = np.zeros((length, 4, MAPQ_BINS), dtype=np.int32)
        self.distance = np.zeros((length, 4, DISTANCE_BINS), dtype=np.int32)

    def add(self, observations: dict[str, np.ndarray]) -> None:
        """Count a batch of observations from `fragment_observations`."""
        cell = observations["position"] * 4 + observations["allele"]
        orientation = observations["orientation"]
        oriented = orientation >= 0
        _count(self.counts, cell * 2 + observations["strand"])
        _count(self.orientation, cell[oriented] * 2 + orientation[oriented])
        _count(
            self.baseq, cell * BASEQ_BINS + np.minimum(observations["base_quality"], BASEQ_BINS - 1)
        )
        _count(
            self.mapq, cell * MAPQ_BINS + np.minimum(observations["mapping_quality"], MAPQ_BINS - 1)
        )
        distance_bin = np.minimum(observations["distance"] // DISTANCE_SCALE, DISTANCE_BINS - 1)
        _count(self.distance, cell * DISTANCE_BINS + distance_bin)

    def depth(self) -> np.ndarray:
        return self.counts.sum(axis=(1, 2))

    def allele_counts(self) -> np.ndarray:
        """(length, 4) observations per allele, both strands."""
        return self.counts.sum(axis=2)


def _count(target: np.ndarray, flat_index: np.ndarray) -> None:
    flat = target.reshape(-1)
    flat += np.bincount(flat_index, minlength=flat.size).astype(flat.dtype)


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
    variance = n_alternate * n_reference / 12 * ((total + 1) - ties / (total * (total - 1)))
    if variance <= 0:
        return 0.0, 1.0

    z_score = (u_statistic - expected) / np.sqrt(variance)
    return float(z_score), float(2 * norm.sf(abs(z_score)))
