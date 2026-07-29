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

from analysis.paired_calling import construct_candidates
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


def _empty_accumulator(length: int) -> list[dict]:
    return [
        {
            **{f"{base}_{strand}": 0 for base in "ACGT" for strand in ("FWD", "REV")},
            **{f"{base}_{orientation}": 0 for base in "ACGT" for orientation in ("F1R2", "F2R1")},
            "depth": 0,
            "mapq_sum": 0,
            "baseq_sum": 0,
            "read_position_sum": 0,
            "clipped_observations": 0,
            "overlap_disagreements": 0,
        }
        for _position in range(length)
    ]


def collect_sample_evidence(
    alignment_path: str, config: PairedConfig, reference_length: int
) -> tuple[list[dict], dict]:
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

    evidence = _empty_accumulator(reference_length)
    overlap_totals = {
        "overlap_positions": 0,
        "overlap_agreements": 0,
        "overlap_disagreements": 0,
    }
    for fragment in fragments:
        observations, overlap = resolve_fragment_observations(
            fragment, config.min_baseq, config.min_distance_from_end
        )
        for key in overlap_totals:
            overlap_totals[key] += int(overlap[key])
        for position in overlap["disagreement_positions"]:
            if 0 <= position < reference_length:
                evidence[position]["overlap_disagreements"] += 1
        for position, observation in observations.items():
            if not 0 <= position < reference_length:
                continue
            row = evidence[position]
            strand = "REV" if observation.is_reverse else "FWD"
            row[f"{observation.base}_{strand}"] += 1
            if observation.orientation:
                row[f"{observation.base}_{observation.orientation}"] += 1
            row["depth"] += 1
            row["mapq_sum"] += observation.mapping_quality
            row["baseq_sum"] += observation.base_quality
            row["read_position_sum"] += observation.read_position
            row["clipped_observations"] += int(observation.clipped)
    stats.update(overlap_totals)
    stats["counted_observations"] = sum(row["depth"] for row in evidence)
    return evidence, stats


def _mean(total: int, count: int) -> float:
    return total / count if count else 0.0


def build_position_evidence(
    chromosome: str,
    reference: str,
    query: list[dict],
    baseline: list[dict],
    config: PairedConfig,
    blacklist: set[int],
) -> list[dict]:
    """Build the all-position public evidence table."""
    rows = []
    length = len(reference)
    for position, ref in enumerate(reference, start=1):
        query_row = query[position - 1]
        baseline_row = baseline[position - 1]
        query_callable = query_row["depth"] >= config.min_query_depth
        baseline_callable = baseline_row["depth"] >= config.min_baseline_depth
        unresolved_edge = not config.shifted_reference_supplied and (
            position <= config.circular_edge_bases or position > length - config.circular_edge_bases
        )
        row = {
            "CHROM": chromosome,
            "POS": position,
            "REF": ref,
            "QUERY_DEPTH": query_row["depth"],
            "BASELINE_DEPTH": baseline_row["depth"],
            "QUERY_COUNTED_FRAGMENTS": query_row["depth"],
            "BASELINE_COUNTED_FRAGMENTS": baseline_row["depth"],
            "QUERY_MEAN_MAPQ": _mean(query_row["mapq_sum"], query_row["depth"]),
            "BASELINE_MEAN_MAPQ": _mean(baseline_row["mapq_sum"], baseline_row["depth"]),
            "QUERY_MEAN_BASEQ": _mean(query_row["baseq_sum"], query_row["depth"]),
            "BASELINE_MEAN_BASEQ": _mean(baseline_row["baseq_sum"], baseline_row["depth"]),
            "QUERY_MEAN_READ_POSITION": _mean(query_row["read_position_sum"], query_row["depth"]),
            "BASELINE_MEAN_READ_POSITION": _mean(
                baseline_row["read_position_sum"], baseline_row["depth"]
            ),
            "QUERY_CLIPPED_OBSERVATIONS": query_row["clipped_observations"],
            "BASELINE_CLIPPED_OBSERVATIONS": baseline_row["clipped_observations"],
            "QUERY_OVERLAP_DISAGREEMENTS": query_row["overlap_disagreements"],
            "BASELINE_OVERLAP_DISAGREEMENTS": baseline_row["overlap_disagreements"],
            "QUERY_CALLABLE": query_callable,
            "BASELINE_CALLABLE": baseline_callable,
            "_JOINT_CALLABLE": (
                query_callable
                and baseline_callable
                and position not in blacklist
                and not unresolved_edge
            ),
        }
        for sample, source in (("QUERY", query_row), ("BASELINE", baseline_row)):
            for base in "ACGT":
                for strand in ("FWD", "REV"):
                    row[f"{sample}_{base}_{strand}"] = source[f"{base}_{strand}"]
                for orientation in ("F1R2", "F2R1"):
                    row[f"{sample}_{base}_{orientation}"] = source[f"{base}_{orientation}"]
        rows.append(row)
    return rows


def _depth_summary(evidence: list[dict]) -> dict:
    depths = np.array([row["depth"] for row in evidence], dtype=float)
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
    candidates = construct_candidates(evidence, config, blacklist)
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
            "legacy_pass_candidates": sum(row["LEGACY_FILTER"] == "PASS" for row in candidates),
        },
        "legacy_projection_written": config.write_legacy_tsv,
        "elapsed_seconds": 0.0,
    }
    qc["elapsed_seconds"] = time.monotonic() - started
    outputs = write_paired_outputs(
        Path(config.output),
        config.sample_name,
        evidence,
        candidates,
        qc,
        config.write_legacy_tsv,
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
