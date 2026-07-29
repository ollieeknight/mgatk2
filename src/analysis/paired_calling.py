"""Pure query/baseline SNV candidate construction and classification."""

from __future__ import annotations

from math import isfinite

from scipy.stats import beta, fisher_exact

from core.config import PairedConfig


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


def _legacy_filters(row: dict, config: PairedConfig, blacklisted: bool) -> list[str]:
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
    raw_ratio = row["QUERY_AF"] / row["BASELINE_AF"] if row["BASELINE_AF"] > 0 else float("inf")
    if raw_ratio < config.min_query_baseline_ratio:
        filters.append("LOW_QUERY_BASELINE_RATIO")
    alt_total = row["QUERY_FWD"] + row["QUERY_REV"]
    strand_bias = abs(row["QUERY_FWD"] - row["QUERY_REV"]) / alt_total if alt_total else 0
    if strand_bias > config.max_strand_bias:
        filters.append("STRAND_BIAS")
    if blacklisted:
        filters.append("BLACKLIST")
    return filters


def construct_candidates(
    evidence_rows: list[dict],
    config: PairedConfig,
    blacklist: set[int],
) -> list[dict]:
    """Create one row for every observed non-reference SNV allele."""
    candidates = []
    p_values = []
    for evidence in evidence_rows:
        ref = evidence["REF"]
        for alt in "ACGT":
            if alt == ref:
                continue
            query_alt = evidence[f"QUERY_{alt}_FWD"] + evidence[f"QUERY_{alt}_REV"]
            baseline_alt = evidence[f"BASELINE_{alt}_FWD"] + evidence[f"BASELINE_{alt}_REV"]
            if query_alt == baseline_alt == 0:
                continue
            query_ref = (
                evidence[f"QUERY_{ref}_FWD"] + evidence[f"QUERY_{ref}_REV"] if ref in "ACGT" else 0
            )
            baseline_ref = (
                evidence[f"BASELINE_{ref}_FWD"] + evidence[f"BASELINE_{ref}_REV"]
                if ref in "ACGT"
                else 0
            )
            query_depth = evidence["QUERY_DEPTH"]
            baseline_depth = evidence["BASELINE_DEPTH"]
            query_af = query_alt / query_depth if query_depth else 0.0
            baseline_af = baseline_alt / baseline_depth if baseline_depth else 0.0
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
            query_alt_fwd = evidence[f"QUERY_{alt}_FWD"]
            query_alt_rev = evidence[f"QUERY_{alt}_REV"]
            query_ref_fwd = evidence.get(f"QUERY_{ref}_FWD", 0)
            query_ref_rev = evidence.get(f"QUERY_{ref}_REV", 0)
            strand_p = _fisher_p([[query_alt_fwd, query_alt_rev], [query_ref_fwd, query_ref_rev]])
            query_ci = binomial_ci(query_alt, query_depth)
            baseline_ci = binomial_ci(baseline_alt, baseline_depth)
            row = {
                "CHROM": evidence["CHROM"],
                "POS": evidence["POS"],
                "REF": ref,
                "ALT": alt,
                "QUERY_DEPTH": query_depth,
                "BASELINE_DEPTH": baseline_depth,
                "QUERY_REF_COUNT": query_ref,
                "QUERY_ALT_COUNT": query_alt,
                "BASELINE_REF_COUNT": baseline_ref,
                "BASELINE_ALT_COUNT": baseline_alt,
                "QUERY_FWD": query_alt_fwd,
                "QUERY_REV": query_alt_rev,
                "BASELINE_FWD": evidence[f"BASELINE_{alt}_FWD"],
                "BASELINE_REV": evidence[f"BASELINE_{alt}_REV"],
                "QUERY_AF": query_af,
                "BASELINE_AF": baseline_af,
                "QUERY_BASELINE_RATIO": ratio,
                "QUERY_AF_CI_LOW": query_ci[0],
                "QUERY_AF_CI_HIGH": query_ci[1],
                "BASELINE_AF_CI_LOW": baseline_ci[0],
                "BASELINE_AF_CI_HIGH": baseline_ci[1],
                "ENRICHMENT_P": enrichment_p,
                "ENRICHMENT_Q": 1.0,
                "STRAND_P": strand_p,
            }
            for sample in ("QUERY", "BASELINE"):
                for allele_name, allele in (("REF", ref), ("ALT", alt)):
                    for orientation in ("F1R2", "F2R1"):
                        row[f"{sample}_{allele_name}_{orientation}"] = evidence.get(
                            f"{sample}_{allele}_{orientation}", 0
                        )
            legacy = _legacy_filters(row, config, evidence["POS"] in blacklist)
            row["LEGACY_FILTER"] = "PASS" if not legacy else "|".join(legacy)
            candidates.append(row)
            p_values.append(enrichment_p)

    for row, q_value in zip(candidates, benjamini_hochberg(p_values), strict=True):
        row["ENRICHMENT_Q"] = q_value
        filters = []
        if row["QUERY_DEPTH"] < config.min_query_depth:
            filters.append("LOW_QUERY_DEPTH")
        if row["BASELINE_DEPTH"] < config.min_baseline_depth:
            filters.append("LOW_BASELINE_DEPTH")
        if row["POS"] in blacklist:
            filters.append("BLACKLIST")
        if not config.shifted_reference_supplied and (
            row["POS"] <= config.circular_edge_bases
            or row["POS"] > len(evidence_rows) - config.circular_edge_bases
        ):
            filters.append("CIRCULAR_EDGE_UNRESOLVED")
        alt_total = row["QUERY_FWD"] + row["QUERY_REV"]
        strand_bias = abs(row["QUERY_FWD"] - row["QUERY_REV"]) / alt_total if alt_total else 0
        if row["STRAND_P"] < 0.001 or strand_bias > config.max_strand_bias:
            filters.append("STRAND_BIAS")
        if q_value > 0.05:
            filters.append("NOT_SIGNIFICANT")
        row["FILTER"] = "PASS" if not filters else "|".join(filters)
    return candidates
