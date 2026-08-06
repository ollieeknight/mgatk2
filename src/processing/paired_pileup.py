"""Evidence-first query/baseline mitochondrial SNV pipeline."""

from __future__ import annotations

import hashlib
import logging
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

import numpy as np
import pysam

from analysis.paired_calling import MIN_ERROR_RATE, construct_candidates
from analysis.quality_stats import BASE_INDEX, QualityHistograms
from core.config import PairedConfig, PipelineConfig
from core.exceptions import InvalidInputError
from data.blacklists import load_blacklist_positions
from file_io.paired_writers import write_paired_outputs
from processing.fragments import resolve_fragment_observations
from processing.readers import BAMReader

logger = logging.getLogger(__name__)


@dataclass
class PairedResult:
    """Completed output paths and summary counts."""

    outputs: dict[str, str]
    evidence_positions: int
    candidates: int
    pass_candidates: int
    callable_positions: int


def load_fasta_reference(reference_path: str, requested_chromosome: str) -> tuple[str, str, str]:
    """Load the FASTA-defined mitochondrial reference and its checksum."""
    path = Path(reference_path)
    if not path.exists():
        raise InvalidInputError(f"Reference FASTA not found: {path}")
    if not Path(f"{path}.fai").exists():
        raise InvalidInputError(f"Reference FASTA index not found: {path}.fai")
    try:
        with pysam.FastaFile(str(path)) as fasta:
            chromosome = next(
                (
                    name
                    for name in (requested_chromosome, "chrM", "MT", "M", "chrMT")
                    if name in fasta.references
                ),
                None,
            )
            if chromosome is None:
                raise InvalidInputError(
                    f"No mitochondrial chromosome found in {path}; available: "
                    f"{', '.join(fasta.references[:10])}"
                )
            sequence = fasta.fetch(chromosome).upper()
    except InvalidInputError:
        raise
    except Exception as exc:
        raise InvalidInputError(f"Cannot read reference FASTA {path}: {exc}") from exc
    if not sequence or any(base not in "ACGTN" for base in sequence):
        raise InvalidInputError("Mitochondrial FASTA sequence is empty or contains invalid bases")
    return chromosome, sequence, hashlib.sha256(sequence.encode("ascii")).hexdigest()


def collect_sample_evidence(
    alignment_path: str, config: PairedConfig, reference_length: int
) -> tuple[QualityHistograms, dict]:
    """Collect fragment-collapsed, strand-specific evidence for one sample."""
    pipeline_config = PipelineConfig(
        min_baseq=config.min_baseq,
        min_mapq=config.min_mapq,
        min_distance_from_end=config.min_distance_from_end,
        mito_chr=config.mito_chr,
        mito_length=reference_length,
        compute_tn5=False,
    )
    fragments, stats = BAMReader(
        alignment_path,
        pipeline_config,
        reference_filename=config.reference,
    ).collect_bulk_reads(config.deduplication)
    if stats["reference_length"] != reference_length:
        raise InvalidInputError(
            f"{alignment_path} header length for {config.mito_chr} is "
            f"{stats['reference_length']}, but FASTA length is {reference_length}"
        )

    histograms = QualityHistograms(reference_length)
    overlap_totals = {
        "overlap_positions": 0,
        "overlap_agreements": 0,
        "overlap_disagreements": 0,
    }
    orientation_index = {"F1R2": 0, "F2R1": 1}
    for fragment in fragments:
        observations, overlap = resolve_fragment_observations(
            fragment, config.min_baseq, config.min_distance_from_end
        )
        for key in overlap_totals:
            overlap_totals[key] += int(overlap[key])
        for position in overlap["disagreement_positions"]:
            if 0 <= position < reference_length:
                histograms.overlap_disagreements[position] += 1
        for position, observation in observations.items():
            if not 0 <= position < reference_length:
                continue
            histograms.add(
                position,
                BASE_INDEX[observation.base],
                int(observation.is_reverse),
                observation.base_quality,
                observation.mapping_quality,
                observation.distance_from_end,
                orientation_index.get(observation.orientation, -1),
            )
            if observation.clipped:
                histograms.clipped[position] += 1
    histograms.flush()
    stats.update(overlap_totals)
    stats["counted_observations"] = int(histograms.depth().sum())
    return histograms, stats


def estimate_error_rates(
    baseline: QualityHistograms, reference: str, max_real_allele_fraction: float = 0.01
) -> dict[str, float]:
    """Per-substitution sequencing error rate, learned from the baseline sample.

    A plain query-versus-baseline Fisher test assumes both samples share an
    error rate, so unequal depth alone can look significant. Estimating the
    rate for each REF>ALT substitution gives the caller an absolute noise floor
    to test against as well.
    """
    per_allele = baseline.allele_counts()
    depth = per_allele.sum(axis=1)
    reference_index = np.array([BASE_INDEX.get(base, -1) for base in reference], dtype=np.int64)
    with np.errstate(invalid="ignore", divide="ignore"):
        fraction = np.where(depth[:, None] > 0, per_allele / depth[:, None], 0.0)

    usable = (depth > 0) & (reference_index >= 0)
    rates = {}
    for reference_base, reference_position in BASE_INDEX.items():
        at_reference = usable & (reference_index == reference_position)
        for alternate_base, alternate_position in BASE_INDEX.items():
            if alternate_base == reference_base:
                continue
            # Exclude sites carrying a plausible real allele, or the estimate
            # absorbs the very heteroplasmy it is meant to distinguish.
            noise_only = at_reference & (
                fraction[:, alternate_position] <= max_real_allele_fraction
            )
            observed = int(per_allele[noise_only, alternate_position].sum())
            total = int(depth[noise_only].sum())
            rates[f"{reference_base}>{alternate_base}"] = max(
                observed / total if total else 0.0, MIN_ERROR_RATE
            )
    return rates


def build_position_evidence(
    chromosome: str,
    reference: str,
    query: QualityHistograms,
    baseline: QualityHistograms,
    config: PairedConfig,
    blacklist: set[int],
) -> list[dict]:
    """Build the all-position public evidence table."""
    rows = []
    length = len(reference)
    samples = (("QUERY", query), ("BASELINE", baseline))
    depths = {name: histograms.depth() for name, histograms in samples}

    for position, ref in enumerate(reference, start=1):
        index = position - 1
        query_depth = int(depths["QUERY"][index])
        baseline_depth = int(depths["BASELINE"][index])
        query_callable = query_depth >= config.min_query_depth
        baseline_callable = baseline_depth >= config.min_baseline_depth
        unresolved_edge = not config.shifted_reference_supplied and (
            position <= config.circular_edge_bases or position > length - config.circular_edge_bases
        )
        row = {
            "CHROM": chromosome,
            "POS": position,
            "REF": ref,
            "QUERY_DEPTH": query_depth,
            "BASELINE_DEPTH": baseline_depth,
            "QUERY_CLIPPED_OBSERVATIONS": int(query.clipped[index]),
            "BASELINE_CLIPPED_OBSERVATIONS": int(baseline.clipped[index]),
            "QUERY_OVERLAP_DISAGREEMENTS": int(query.overlap_disagreements[index]),
            "BASELINE_OVERLAP_DISAGREEMENTS": int(baseline.overlap_disagreements[index]),
            "QUERY_CALLABLE": query_callable,
            "BASELINE_CALLABLE": baseline_callable,
            "_JOINT_CALLABLE": (
                query_callable
                and baseline_callable
                and position not in blacklist
                and not unresolved_edge
            ),
        }
        for name, histograms in samples:
            for base, base_index in BASE_INDEX.items():
                row[f"{name}_{base}_FWD"] = int(histograms.counts[index, base_index, 0])
                row[f"{name}_{base}_REV"] = int(histograms.counts[index, base_index, 1])
        rows.append(row)
    return rows


def _depth_summary(histograms: QualityHistograms) -> dict:
    depths = histograms.depth().astype(float)
    return {
        "breadth": float(np.count_nonzero(depths) / len(depths)),
        "depth_quantiles": {
            str(quantile): float(np.quantile(depths, quantile))
            for quantile in (0, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
        },
    }


def _package_version() -> str:
    try:
        return version("mgatk2")
    except PackageNotFoundError:
        return "unknown"


def _git_commit() -> str | None:
    try:
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=Path(__file__).parents[2],
            capture_output=True,
            text=True,
            check=False,
        )
    except OSError:
        return None
    return result.stdout.strip() if result.returncode == 0 else None


def run_paired_pipeline(config: PairedConfig) -> PairedResult:
    """Run the single public paired-analysis orchestration seam."""
    started = time.monotonic()
    chromosome, reference, checksum = load_fasta_reference(config.reference, config.mito_chr)
    config.mito_chr = chromosome
    blacklist = load_blacklist_positions(
        build="none", custom_bed=config.custom_blacklist, mito_chr=chromosome
    )
    query, query_stats = collect_sample_evidence(config.query, config, len(reference))
    baseline, baseline_stats = collect_sample_evidence(config.baseline, config, len(reference))
    evidence = build_position_evidence(chromosome, reference, query, baseline, config, blacklist)
    error_rates = estimate_error_rates(baseline, reference)
    candidates = construct_candidates(evidence, query, baseline, config, blacklist, error_rates)
    callable_positions = sum(row["_JOINT_CALLABLE"] for row in evidence)
    qc = {
        "schema_version": config.qc_schema_version,
        "evidence_schema_version": config.evidence_schema_version,
        "candidate_schema_version": config.candidate_schema_version,
        "mgatk2_version": _package_version(),
        "git_commit": _git_commit(),
        "command_line": sys.argv,
        "parameters": asdict(config),
        "inputs": {
            "query": {
                "path": str(Path(config.query).resolve()),
                "type": Path(config.query).suffix.lower().lstrip("."),
                "declared_upstream_consensus": config.input_is_consensus,
                "statistics": query_stats,
                **_depth_summary(query),
            },
            "baseline": {
                "path": str(Path(config.baseline).resolve()),
                "type": Path(config.baseline).suffix.lower().lstrip("."),
                "declared_upstream_consensus": config.input_is_consensus,
                "statistics": baseline_stats,
                **_depth_summary(baseline),
            },
        },
        "reference": {
            "path": str(Path(config.reference).resolve()),
            "chromosome": chromosome,
            "length": len(reference),
            "sha256": checksum,
            "standard_reference": True,
            "shifted_reference_supplied": config.shifted_reference_supplied,
        },
        "deduplication": config.deduplication,
        "snv_only": True,
        "indels_called": False,
        "blacklist_numt_strategy": (
            "user_chrM_blacklist_and_MAPQ"
            if config.custom_blacklist
            else "MAPQ_only_no_chrM_blacklist"
        ),
        "circular_edge": {
            "bases": config.circular_edge_bases,
            "status": (
                "SHIFTED_REFERENCE_SUPPLIED"
                if config.shifted_reference_supplied
                else "CIRCULAR_EDGE_UNRESOLVED"
            ),
        },
        "counts": {
            "evidence_positions": len(evidence),
            "callable_positions": callable_positions,
            "candidates": len(candidates),
            "pass_candidates": sum(row["FILTER"] == "PASS" for row in candidates),
        },
        "substitution_error_rates": error_rates,
        "numt_strategy": (
            "autosomal_median_depth_and_MAPQ"
            if config.autosomal_median_depth is not None
            else "MAPQ_only"
        ),
        "elapsed_seconds": 0.0,
    }
    qc["elapsed_seconds"] = time.monotonic() - started
    outputs = write_paired_outputs(
        Path(config.output),
        config.sample_name,
        evidence,
        candidates,
        qc,
    )
    logger.info(
        "Paired analysis complete: %d positions, %d candidates, %d PASS",
        len(evidence),
        len(candidates),
        qc["counts"]["pass_candidates"],
    )
    return PairedResult(
        outputs=outputs,
        evidence_positions=len(evidence),
        candidates=len(candidates),
        pass_candidates=qc["counts"]["pass_candidates"],
        callable_positions=callable_positions,
    )
