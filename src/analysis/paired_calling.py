"""Pure tumour/normal SNV candidate construction and classification."""

from __future__ import annotations

from math import isfinite, log10

import numpy as np
from scipy.stats import beta, binom, fisher_exact, hypergeom

from analysis.quality_stats import (
    BASE_INDEX,
    DISTANCE_SCALE,
    QualityHistograms,
    histogram_median,
    rank_sum,
)
from core.config import PairedConfig

# Below this, an estimated substitution error rate is indistinguishable from
# having seen no errors at all, and the binomial test stops being meaningful.
MIN_ERROR_RATE = 1e-6

# Shared by every per-allele artefact test; deliberately strict, because all of
# them are run at every observed allele across the genome.
ARTEFACT_P = 0.001


def _fisher_p(table: list[list[int]]) -> float:
    """Two-sided Fisher exact p-value; 1 when undefined."""
    value = float(fisher_exact(table).pvalue)
    return value if isfinite(value) else 1.0


def phred(probability: float) -> float:
    """Convert an adjusted p-value into a capped phred-scaled quality."""
    if probability <= 0:
        return 1000.0
    return min(1000.0, max(0.0, -10 * log10(probability)))


def binomial_ci(successes: np.ndarray, total: np.ndarray, alpha: float = 0.05):
    """Clopper-Pearson confidence intervals, elementwise."""
    # Substitute a valid shape where the bound is fixed, so ppf never sees zero.
    low = np.where(
        successes == 0, 0.0, beta.ppf(alpha / 2, np.maximum(successes, 1), total - successes + 1)
    )
    high = np.where(
        successes == total,
        1.0,
        beta.ppf(1 - alpha / 2, successes + 1, np.maximum(total - successes, 1)),
    )
    return low, high


def _enrichment_p(tumor_alt, tumor_depth, normal_alt, normal_depth) -> np.ndarray:
    """One-sided Fisher exact p-values for tumour alternate enrichment, elementwise.

    Same formula as scipy's fisher_exact(alternative="greater"), vectorised:
    an empty row or column gives 1.
    """
    tumor_other = tumor_depth - tumor_alt
    normal_other = normal_depth - normal_alt
    degenerate = ((tumor_depth == 0) | (normal_depth == 0) | (tumor_alt + normal_alt == 0)) | (
        tumor_other + normal_other == 0
    )
    with np.errstate(all="ignore"):
        p = hypergeom.cdf(
            tumor_other,
            tumor_depth + normal_depth,
            tumor_depth,
            tumor_other + normal_other,
        )
    p = np.where(degenerate | ~np.isfinite(p), 1.0, p)
    return np.minimum(p, 1.0)


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Return monotone Benjamini-Hochberg adjusted p-values."""
    order = sorted(range(len(p_values)), key=p_values.__getitem__)
    adjusted = [1.0] * len(p_values)
    running = 1.0
    for rank_index in range(len(order) - 1, -1, -1):
        original_index = order[rank_index]
        rank = rank_index + 1
        running = min(running, p_values[original_index] * len(order) / rank)
        adjusted[original_index] = min(1.0, running)
    return adjusted


def _allele_quality(
    histograms: QualityHistograms, index: int, allele: int | None
) -> dict[str, float]:
    """Median base quality, mapping quality, and distance from read end."""
    if allele is None:
        return {"mbq": 0.0, "mmq": 0.0, "mpos": 0.0}
    return {
        "mbq": histogram_median(histograms.baseq[index, allele]),
        "mmq": histogram_median(histograms.mapq[index, allele]),
        "mpos": histogram_median(histograms.distance[index, allele], DISTANCE_SCALE),
    }


def _rank_sums(
    histograms: QualityHistograms, index: int, alternate: int, reference: int | None
) -> dict[str, float]:
    """Alternate-versus-reference rank-sum z-scores and p-values."""
    results = {}
    for name, source in (
        ("rsbq", histograms.baseq),
        ("rsmq", histograms.mapq),
        ("rspos", histograms.distance),
    ):
        if reference is None:
            z_score, p_value = 0.0, 1.0
        else:
            z_score, p_value = rank_sum(source[index, alternate], source[index, reference])
        results[name] = z_score
        results[f"{name}_p"] = p_value
    return results


def _strand_bias(forward: int, reverse: int) -> float:
    total = forward + reverse
    return abs(forward - reverse) / total if total else 0.0


def _pair(array, index: int, allele: int | None) -> tuple[int, int]:
    """Forward/reverse or F1R2/F2R1 counts for one allele; zeros for an N reference."""
    if allele is None:
        return 0, 0
    return int(array[index, allele, 0]), int(array[index, allele, 1])


def unresolved_edge(position: int, length: int, config: PairedConfig) -> bool:
    """Within --circular-edge-bases of either end of a linear reference."""
    edge = config.circular_edge_bases
    return position <= edge or position > length - edge


def construct_candidates(
    chromosome: str,
    reference: str,
    tumor: QualityHistograms,
    normal: QualityHistograms,
    config: PairedConfig,
    blacklist: set[int],
    error_rates: dict[str, float],
) -> list[dict]:
    """Create one row for every observed non-reference SNV allele."""
    candidates = []
    tumor_depths = tumor.depth()
    normal_depths = normal.depth()

    for index, ref in enumerate(reference):
        reference_allele = BASE_INDEX.get(ref)
        for alt, alternate_allele in BASE_INDEX.items():
            if alt == ref:
                continue
            tumor_alt_fwd, tumor_alt_rev = _pair(tumor.counts, index, alternate_allele)
            normal_alt_fwd, normal_alt_rev = _pair(normal.counts, index, alternate_allele)
            tumor_alt = tumor_alt_fwd + tumor_alt_rev
            normal_alt = normal_alt_fwd + normal_alt_rev
            if tumor_alt == normal_alt == 0:
                continue

            tumor_ref_fwd, tumor_ref_rev = _pair(tumor.counts, index, reference_allele)
            normal_ref_fwd, normal_ref_rev = _pair(normal.counts, index, reference_allele)
            tumor_depth = int(tumor_depths[index])
            normal_depth = int(normal_depths[index])

            strand_p = _fisher_p([[tumor_alt_fwd, tumor_alt_rev], [tumor_ref_fwd, tumor_ref_rev]])

            tumor_alt_f1r2, tumor_alt_f2r1 = _pair(tumor.orientation, index, alternate_allele)
            tumor_ref_f1r2, tumor_ref_f2r1 = _pair(tumor.orientation, index, reference_allele)
            normal_alt_f1r2, normal_alt_f2r1 = _pair(normal.orientation, index, alternate_allele)
            normal_ref_f1r2, normal_ref_f2r1 = _pair(normal.orientation, index, reference_allele)
            # Only meaningful when both mates were present; single-end and
            # orphan input leaves every orientation count at zero.
            orientation_p = (
                _fisher_p([[tumor_alt_f1r2, tumor_alt_f2r1], [tumor_ref_f1r2, tumor_ref_f2r1]])
                if (tumor_alt_f1r2 + tumor_alt_f2r1) and (tumor_ref_f1r2 + tumor_ref_f2r1)
                else 1.0
            )

            row = {
                "chrom": chromosome,
                "pos": index + 1,
                "ref": ref,
                "alt": alt,
                "normal_dp": normal_depth,
                "normal_ref_count": normal_ref_fwd + normal_ref_rev,
                "normal_ac": normal_alt,
                "normal_af": normal_alt / normal_depth if normal_depth else 0.0,
                "tumor_dp": tumor_depth,
                "tumor_ref_count": tumor_ref_fwd + tumor_ref_rev,
                "tumor_ac": tumor_alt,
                "tumor_af": tumor_alt / tumor_depth if tumor_depth else 0.0,
                "normal_ref_fwd": normal_ref_fwd,
                "normal_ref_rev": normal_ref_rev,
                "normal_alt_fwd": normal_alt_fwd,
                "normal_alt_rev": normal_alt_rev,
                "tumor_ref_fwd": tumor_ref_fwd,
                "tumor_ref_rev": tumor_ref_rev,
                "tumor_alt_fwd": tumor_alt_fwd,
                "tumor_alt_rev": tumor_alt_rev,
                "normal_ref_f1r2": normal_ref_f1r2,
                "normal_ref_f2r1": normal_ref_f2r1,
                "normal_alt_f1r2": normal_alt_f1r2,
                "normal_alt_f2r1": normal_alt_f2r1,
                "tumor_ref_f1r2": tumor_ref_f1r2,
                "tumor_ref_f2r1": tumor_ref_f2r1,
                "tumor_alt_f1r2": tumor_alt_f1r2,
                "tumor_alt_f2r1": tumor_alt_f2r1,
                "error_rate": error_rates.get(f"{ref}>{alt}", MIN_ERROR_RATE),
                "strand_p": strand_p,
                "orient_p": orientation_p,
            }
            # Quality medians are tumour-derived only: nothing filters on the
            # normal-side values.
            for label, allele in (("ref", reference_allele), ("alt", alternate_allele)):
                for metric, value in _allele_quality(tumor, index, allele).items():
                    row[f"tumor_{label}_{metric}"] = value
            row.update(_rank_sums(tumor, index, alternate_allele, reference_allele))
            candidates.append(row)

    # The allele-fraction statistics are vectorised: scipy's per-call overhead
    # dominates at full depth, where nearly every position has some alternate.
    def column(name):
        return np.array([row[name] for row in candidates], dtype=np.int64)

    tumor_alt, tumor_depth = column("tumor_ac"), column("tumor_dp")
    normal_alt, normal_depth = column("normal_ac"), column("normal_dp")
    error_rate = np.array([row["error_rate"] for row in candidates], dtype=np.float64)
    tumor_ci = binomial_ci(tumor_alt, tumor_depth)
    normal_ci = binomial_ci(normal_alt, normal_depth)
    statistics = {
        "tumor_af_ci_low": tumor_ci[0],
        "tumor_af_ci_high": tumor_ci[1],
        "normal_af_ci_low": normal_ci[0],
        "normal_af_ci_high": normal_ci[1],
        "seq_p": np.where(tumor_alt > 0, binom.sf(tumor_alt - 1, tumor_depth, error_rate), 1.0),
        "enrich_p": _enrichment_p(tumor_alt, tumor_depth, normal_alt, normal_depth),
    }
    for name, values in statistics.items():
        for row, value in zip(candidates, values.tolist(), strict=True):
            row[name] = value

    q_values = benjamini_hochberg([row["enrich_p"] for row in candidates])
    for row, q_value in zip(candidates, q_values, strict=True):
        row["enrich_q"] = q_value
        row["qual"] = round(phred(q_value), 2)
        row["filter"] = ";".join(_filters(row, config, blacklist, len(reference))) or "PASS"
    return candidates


def _filters(row: dict, config: PairedConfig, blacklist: set[int], length: int) -> list[str]:
    """Technical and statistical flags; an empty list means PASS."""
    filters = []
    if row["tumor_dp"] < config.min_tumor_depth:
        filters.append("LOW_TUMOR_DEPTH")
    if row["normal_dp"] < config.min_normal_depth:
        filters.append("LOW_NORMAL_DEPTH")
    if row["tumor_ac"] < config.min_alt_observations:
        filters.append("LOW_ALT_OBSERVATIONS")
    if row["tumor_af"] < config.min_tumor_af:
        filters.append("LOW_TUMOR_AF")
    if row["normal_af"] > config.max_normal_af:
        filters.append("HIGH_NORMAL_AF")
    if row["pos"] in blacklist:
        filters.append("BLACKLIST")
    if unresolved_edge(row["pos"], length, config):
        filters.append("CIRCULAR_EDGE_UNRESOLVED")
    if (
        row["strand_p"] < ARTEFACT_P
        or _strand_bias(row["tumor_alt_fwd"], row["tumor_alt_rev"]) > config.max_strand_bias
    ):
        filters.append("STRAND_BIAS")
    if row["orient_p"] < ARTEFACT_P:
        filters.append("ORIENTATION_BIAS")
    if row["seq_p"] > ARTEFACT_P:
        filters.append("WEAK_EVIDENCE")
    if row["enrich_q"] > 0.05:
        filters.append("NOT_SIGNIFICANT")
    # A negative z means the alternate observations sit systematically lower on
    # the metric than the reference ones at the same position, which is the
    # artefact direction; a positive skew is not evidence against the allele.
    for name, flag in (("rsbq", "BASE_QUAL"), ("rsmq", "MAP_QUAL"), ("rspos", "POSITION")):
        if row[name] < 0 and row[f"{name}_p"] < ARTEFACT_P:
            filters.append(flag)
    # A NuMT carried at one autosomal copy contributes roughly the autosomal
    # median depth of reads, so alternate support at or below that is not
    # separable from nuclear leakage.
    if (
        config.autosomal_median_depth is not None
        and row["tumor_ac"] <= config.autosomal_median_depth
    ):
        filters.append("POSSIBLE_NUMT")
    return filters
