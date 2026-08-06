"""Pure query/baseline SNV candidate construction and classification."""

from __future__ import annotations

from math import isfinite

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

# Shared by the strand and orientation tests; deliberately strict, because both
# are run at every observed allele across the genome.
ARTEFACT_P = 0.001


def _fisher_p(table: list[list[int]], alternative: str = "two-sided") -> float:
    value = float(fisher_exact(table, alternative=alternative).pvalue)
    return value if isfinite(value) else 1.0


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
        return {"MBQ": 0.0, "MMQ": 0.0, "MPOS": 0.0}
    return {
        "MBQ": histogram_median(histograms.baseq[index, allele]),
        "MMQ": histogram_median(histograms.mapq[index, allele]),
        "MPOS": histogram_median(histograms.distance[index, allele], DISTANCE_SCALE),
    }


def _rank_sums(
    histograms: QualityHistograms, index: int, alternate: int, reference: int | None
) -> dict[str, float]:
    """Alternate-versus-reference rank-sum z-scores and p-values."""
    if reference is None:
        return dict.fromkeys(("RSBQ", "RSBQ_P", "RSMQ", "RSMQ_P", "RSPOS", "RSPOS_P"), 0.0) | {
            "RSBQ_P": 1.0,
            "RSMQ_P": 1.0,
            "RSPOS_P": 1.0,
        }
    results = {}
    for name, source in (
        ("RSBQ", histograms.baseq),
        ("RSMQ", histograms.mapq),
        ("RSPOS", histograms.distance),
    ):
        z_score, p_value = rank_sum(source[index, alternate], source[index, reference])
        results[name] = z_score
        results[f"{name}_P"] = p_value
    return results


def _strand_bias(forward: int, reverse: int) -> float:
    total = forward + reverse
    return abs(forward - reverse) / total if total else 0.0


def construct_candidates(
    evidence_rows: list[dict],
    query: QualityHistograms,
    baseline: QualityHistograms,
    config: PairedConfig,
    blacklist: set[int],
    error_rates: dict[str, float],
) -> list[dict]:
    """Create one row for every observed non-reference SNV allele."""
    candidates = []
    p_values = []
    length = len(evidence_rows)

    for evidence in evidence_rows:
        index = evidence["POS"] - 1
        ref = evidence["REF"]
        reference_allele = BASE_INDEX.get(ref)

        for alt, alternate_allele in BASE_INDEX.items():
            if alt == ref:
                continue
            query_alt_fwd = evidence[f"QUERY_{alt}_FWD"]
            query_alt_rev = evidence[f"QUERY_{alt}_REV"]
            baseline_alt_fwd = evidence[f"BASELINE_{alt}_FWD"]
            baseline_alt_rev = evidence[f"BASELINE_{alt}_REV"]
            query_alt = query_alt_fwd + query_alt_rev
            baseline_alt = baseline_alt_fwd + baseline_alt_rev
            if query_alt == baseline_alt == 0:
                continue

            known_reference = ref in BASE_INDEX
            query_ref_fwd = evidence[f"QUERY_{ref}_FWD"] if known_reference else 0
            query_ref_rev = evidence[f"QUERY_{ref}_REV"] if known_reference else 0
            baseline_ref_fwd = evidence[f"BASELINE_{ref}_FWD"] if known_reference else 0
            baseline_ref_rev = evidence[f"BASELINE_{ref}_REV"] if known_reference else 0
            query_ref = query_ref_fwd + query_ref_rev
            baseline_ref = baseline_ref_fwd + baseline_ref_rev
            query_depth = evidence["QUERY_DEPTH"]
            baseline_depth = evidence["BASELINE_DEPTH"]
            query_af = query_alt / query_depth if query_depth else 0.0
            baseline_af = baseline_alt / baseline_depth if baseline_depth else 0.0

            # Smoothed so a zero-count baseline yields a finite, ordered ratio.
            ratio = ((query_alt + 0.5) / (query_depth + 1)) / (
                (baseline_alt + 0.5) / (baseline_depth + 1)
            )
            enrichment_p = _fisher_p(
                [
                    [query_alt, max(0, query_depth - query_alt)],
                    [baseline_alt, max(0, baseline_depth - baseline_alt)],
                ],
                alternative="greater",
            )

            error_rate = error_rates.get(f"{ref}>{alt}", MIN_ERROR_RATE)
            sequencing_error_p = (
                float(binom.sf(query_alt - 1, query_depth, error_rate))
                if query_depth and query_alt
                else 1.0
            )

            strand_p = _fisher_p([[query_alt_fwd, query_alt_rev], [query_ref_fwd, query_ref_rev]])
            alt_f1r2, alt_f2r1 = (
                int(query.orientation[index, alternate_allele, 0]),
                int(query.orientation[index, alternate_allele, 1]),
            )
            ref_f1r2, ref_f2r1 = (
                (
                    int(query.orientation[index, reference_allele, 0]),
                    int(query.orientation[index, reference_allele, 1]),
                )
                if reference_allele is not None
                else (0, 0)
            )
            # Only meaningful when both mates were present; single-end and
            # orphan input leaves every orientation count at zero.
            orientation_p = (
                _fisher_p([[alt_f1r2, alt_f2r1], [ref_f1r2, ref_f2r1]])
                if (alt_f1r2 + alt_f2r1) and (ref_f1r2 + ref_f2r1)
                else 1.0
            )

            query_ci = binomial_ci(query_alt, query_depth)
            baseline_ci = binomial_ci(baseline_alt, baseline_depth)

            row = {
                "CHROM": evidence["CHROM"],
                "POS": evidence["POS"],
                "REF": ref,
                "ALT": alt,
                "QUERY_DEPTH": query_depth,
                "QUERY_REF_COUNT": query_ref,
                "QUERY_REF_FWD": query_ref_fwd,
                "QUERY_REF_REV": query_ref_rev,
                "QUERY_ALT_COUNT": query_alt,
                "QUERY_ALT_FWD": query_alt_fwd,
                "QUERY_ALT_REV": query_alt_rev,
                "QUERY_AF": query_af,
                "QUERY_AF_CI_LOW": query_ci[0],
                "QUERY_AF_CI_HIGH": query_ci[1],
                "BASELINE_DEPTH": baseline_depth,
                "BASELINE_REF_COUNT": baseline_ref,
                "BASELINE_REF_FWD": baseline_ref_fwd,
                "BASELINE_REF_REV": baseline_ref_rev,
                "BASELINE_ALT_COUNT": baseline_alt,
                "BASELINE_ALT_FWD": baseline_alt_fwd,
                "BASELINE_ALT_REV": baseline_alt_rev,
                "BASELINE_AF": baseline_af,
                "BASELINE_AF_CI_LOW": baseline_ci[0],
                "BASELINE_AF_CI_HIGH": baseline_ci[1],
                "QUERY_BASELINE_RATIO": ratio,
                "ERROR_RATE": error_rate,
                "SEQUENCING_ERROR_P": sequencing_error_p,
                "ENRICHMENT_P": enrichment_p,
                "ENRICHMENT_Q": 1.0,
                "STRAND_P": strand_p,
                "ORIENTATION_P": orientation_p,
            }

            for name, histograms, alt_pair, ref_pair in (
                ("QUERY", query, (alt_f1r2, alt_f2r1), (ref_f1r2, ref_f2r1)),
                (
                    "BASELINE",
                    baseline,
                    (
                        int(baseline.orientation[index, alternate_allele, 0]),
                        int(baseline.orientation[index, alternate_allele, 1]),
                    ),
                    (
                        (
                            int(baseline.orientation[index, reference_allele, 0]),
                            int(baseline.orientation[index, reference_allele, 1]),
                        )
                        if reference_allele is not None
                        else (0, 0)
                    ),
                ),
            ):
                row[f"{name}_REF_F1R2"], row[f"{name}_REF_F2R1"] = ref_pair
                row[f"{name}_ALT_F1R2"], row[f"{name}_ALT_F2R1"] = alt_pair
                for label, allele in (("REF", reference_allele), ("ALT", alternate_allele)):
                    for metric, value in _allele_quality(histograms, index, allele).items():
                        row[f"{name}_{label}_{metric}"] = value

            row.update(
                {
                    f"QUERY_{key}": value
                    for key, value in _rank_sums(
                        query, index, alternate_allele, reference_allele
                    ).items()
                }
            )

            candidates.append(row)
            p_values.append(enrichment_p)

    for row, q_value in zip(candidates, benjamini_hochberg(p_values), strict=True):
        row["ENRICHMENT_Q"] = q_value
        row["FILTER"] = ";".join(_filters(row, config, blacklist, length)) or "PASS"
    return candidates


def _filters(row: dict, config: PairedConfig, blacklist: set[int], length: int) -> list[str]:
    """Technical and statistical flags; an empty list means PASS."""
    filters = []
    if row["QUERY_DEPTH"] < config.min_query_depth:
        filters.append("LOW_QUERY_DEPTH")
    if row["BASELINE_DEPTH"] < config.min_baseline_depth:
        filters.append("LOW_BASELINE_DEPTH")
    if row["QUERY_ALT_COUNT"] < config.min_alt_observations:
        filters.append("LOW_ALT_OBSERVATIONS")
    if row["QUERY_AF"] < config.min_query_af:
        filters.append("LOW_QUERY_AF")
    if row["BASELINE_AF"] > config.max_baseline_af:
        filters.append("HIGH_BASELINE_AF")
    if row["POS"] in blacklist:
        filters.append("BLACKLIST")
    if not config.shifted_reference_supplied and (
        row["POS"] <= config.circular_edge_bases or row["POS"] > length - config.circular_edge_bases
    ):
        filters.append("CIRCULAR_EDGE_UNRESOLVED")
    if (
        row["STRAND_P"] < ARTEFACT_P
        or _strand_bias(row["QUERY_ALT_FWD"], row["QUERY_ALT_REV"]) > config.max_strand_bias
    ):
        filters.append("STRAND_BIAS")
    if row["ORIENTATION_P"] < ARTEFACT_P:
        filters.append("ORIENTATION_BIAS")
    if row["SEQUENCING_ERROR_P"] > ARTEFACT_P:
        filters.append("WEAK_EVIDENCE")
    if row["ENRICHMENT_Q"] > 0.05:
        filters.append("NOT_SIGNIFICANT")
    # A NuMT carried at one autosomal copy contributes roughly the autosomal
    # median depth of reads, so alternate support at or below that is not
    # separable from nuclear leakage.
    if (
        config.autosomal_median_depth is not None
        and row["QUERY_ALT_COUNT"] <= config.autosomal_median_depth
    ):
        filters.append("POSSIBLE_NUMT")
    return filters
