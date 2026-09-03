"""Pure tumour/normal SNV candidate construction and classification."""

from __future__ import annotations

from math import isfinite, log10

from scipy.stats import beta, binom, fisher_exact

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


def _fisher_p(table: list[list[int]], alternative: str = "two-sided") -> float:
    value = float(fisher_exact(table, alternative=alternative).pvalue)
    return value if isfinite(value) else 1.0


def phred(probability: float) -> float:
    """Convert an adjusted p-value into a capped phred-scaled quality."""
    if probability <= 0:
        return 1000.0
    return min(1000.0, max(0.0, -10 * log10(probability)))


def binomial_ci(successes: int, total: int, alpha: float = 0.05) -> tuple[float, float]:
    """Clopper-Pearson confidence interval."""
    if total == 0:
        return 0.0, 1.0
    low = 0.0 if successes == 0 else float(beta.ppf(alpha / 2, successes, total - successes + 1))
    high = (
        1.0
        if successes == total
        else float(beta.ppf(1 - alpha / 2, successes + 1, total - successes))
    )
    return low, high


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    """Return monotone Benjamini-Hochberg adjusted p-values."""
    if not p_values:
        return []
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
    if reference is None:
        return {
            "rsbq": 0.0,
            "rsbq_p": 1.0,
            "rsmq": 0.0,
            "rsmq_p": 1.0,
            "rspos": 0.0,
            "rspos_p": 1.0,
        }
    results = {}
    for name, source in (
        ("rsbq", histograms.baseq),
        ("rsmq", histograms.mapq),
        ("rspos", histograms.distance),
    ):
        z_score, p_value = rank_sum(source[index, alternate], source[index, reference])
        results[name] = z_score
        results[f"{name}_p"] = p_value
    return results


def _strand_bias(forward: int, reverse: int) -> float:
    total = forward + reverse
    return abs(forward - reverse) / total if total else 0.0


def _orientation(histograms: QualityHistograms, index: int, allele: int | None) -> tuple[int, int]:
    if allele is None:
        return 0, 0
    return (
        int(histograms.orientation[index, allele, 0]),
        int(histograms.orientation[index, allele, 1]),
    )


def construct_candidates(
    evidence_rows: list[dict],
    tumor: QualityHistograms,
    normal: QualityHistograms,
    config: PairedConfig,
    blacklist: set[int],
    error_rates: dict[str, float],
) -> list[dict]:
    """Create one row for every observed non-reference SNV allele."""
    candidates = []
    p_values = []
    length = len(evidence_rows)

    for evidence in evidence_rows:
        index = evidence["pos"] - 1
        ref = evidence["ref"]
        reference_allele = BASE_INDEX.get(ref)

        for alt, alternate_allele in BASE_INDEX.items():
            if alt == ref:
                continue
            tumor_alt_fwd = evidence[f"tumor_{alt.lower()}_fwd"]
            tumor_alt_rev = evidence[f"tumor_{alt.lower()}_rev"]
            normal_alt_fwd = evidence[f"normal_{alt.lower()}_fwd"]
            normal_alt_rev = evidence[f"normal_{alt.lower()}_rev"]
            tumor_alt = tumor_alt_fwd + tumor_alt_rev
            normal_alt = normal_alt_fwd + normal_alt_rev
            if tumor_alt == normal_alt == 0:
                continue

            known_reference = ref in BASE_INDEX
            reference_key = ref.lower()
            tumor_ref_fwd = evidence[f"tumor_{reference_key}_fwd"] if known_reference else 0
            tumor_ref_rev = evidence[f"tumor_{reference_key}_rev"] if known_reference else 0
            normal_ref_fwd = evidence[f"normal_{reference_key}_fwd"] if known_reference else 0
            normal_ref_rev = evidence[f"normal_{reference_key}_rev"] if known_reference else 0
            tumor_ref = tumor_ref_fwd + tumor_ref_rev
            normal_ref = normal_ref_fwd + normal_ref_rev
            tumor_depth = evidence["tumor_dp"]
            normal_depth = evidence["normal_dp"]
            tumor_af = tumor_alt / tumor_depth if tumor_depth else 0.0
            normal_af = normal_alt / normal_depth if normal_depth else 0.0

            enrichment_p = _fisher_p(
                [
                    [tumor_alt, max(0, tumor_depth - tumor_alt)],
                    [normal_alt, max(0, normal_depth - normal_alt)],
                ],
                alternative="greater",
            )

            error_rate = error_rates.get(f"{ref}>{alt}", MIN_ERROR_RATE)
            sequencing_error_p = (
                float(binom.sf(tumor_alt - 1, tumor_depth, error_rate))
                if tumor_depth and tumor_alt
                else 1.0
            )

            strand_p = _fisher_p([[tumor_alt_fwd, tumor_alt_rev], [tumor_ref_fwd, tumor_ref_rev]])
            tumor_alt_f1r2, tumor_alt_f2r1 = _orientation(tumor, index, alternate_allele)
            tumor_ref_f1r2, tumor_ref_f2r1 = _orientation(tumor, index, reference_allele)
            normal_alt_f1r2, normal_alt_f2r1 = _orientation(normal, index, alternate_allele)
            normal_ref_f1r2, normal_ref_f2r1 = _orientation(normal, index, reference_allele)
            # Only meaningful when both mates were present; single-end and
            # orphan input leaves every orientation count at zero.
            orientation_p = (
                _fisher_p([[tumor_alt_f1r2, tumor_alt_f2r1], [tumor_ref_f1r2, tumor_ref_f2r1]])
                if (tumor_alt_f1r2 + tumor_alt_f2r1) and (tumor_ref_f1r2 + tumor_ref_f2r1)
                else 1.0
            )

            tumor_ci = binomial_ci(tumor_alt, tumor_depth)
            normal_ci = binomial_ci(normal_alt, normal_depth)

            row = {
                "chrom": evidence["chrom"],
                "pos": evidence["pos"],
                "ref": ref,
                "alt": alt,
                "normal_dp": normal_depth,
                "normal_ref_count": normal_ref,
                "normal_ac": normal_alt,
                "normal_af": normal_af,
                "tumor_dp": tumor_depth,
                "tumor_ref_count": tumor_ref,
                "tumor_ac": tumor_alt,
                "tumor_af": tumor_af,
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
                "normal_af_ci_low": normal_ci[0],
                "normal_af_ci_high": normal_ci[1],
                "tumor_af_ci_low": tumor_ci[0],
                "tumor_af_ci_high": tumor_ci[1],
                "error_rate": error_rate,
                "seq_p": sequencing_error_p,
                "enrich_p": enrichment_p,
                "enrich_q": 1.0,
                "strand_p": strand_p,
                "orient_p": orientation_p,
            }

            # Quality medians are tumour-derived only: nothing filters on the
            # normal-side values and they doubled the record for no reader.
            for label, allele in (("ref", reference_allele), ("alt", alternate_allele)):
                for metric, value in _allele_quality(tumor, index, allele).items():
                    row[f"tumor_{label}_{metric}"] = value
            row.update(_rank_sums(tumor, index, alternate_allele, reference_allele))

            candidates.append(row)
            p_values.append(enrichment_p)

    for row, q_value in zip(candidates, benjamini_hochberg(p_values), strict=True):
        row["enrich_q"] = q_value
        row["qual"] = round(phred(q_value), 2)
        row["filter"] = ";".join(_filters(row, config, blacklist, length)) or "PASS"
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
    if not config.shifted_reference_supplied and (
        row["pos"] <= config.circular_edge_bases or row["pos"] > length - config.circular_edge_bases
    ):
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
