"""Evidence-first tumour/normal mitochondrial SNV pipeline."""

from __future__ import annotations

import hashlib
import logging
import sys
from dataclasses import asdict, dataclass
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pysam

from analysis.paired_calling import MIN_ERROR_RATE, construct_candidates, unresolved_edge
from analysis.quality_stats import BASE_INDEX, QualityHistograms
from core.config import PairedConfig, PipelineConfig
from core.exceptions import InvalidInputError
from data.blacklists import load_bed_positions
from file_io.paired_writers import write_paired_outputs
from processing.fragments import fragment_observations
from processing.readers import BAMReader, resolve_mito_contig

logger = logging.getLogger(__name__)

FRAGMENTS_PER_CHUNK = 100_000


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
            chromosome = resolve_mito_contig(fasta.references, requested_chromosome, path)
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
    overlap_totals = {"overlap_positions": 0, "overlap_agreements": 0, "overlap_disagreements": 0}
    # Chunked so the per-base arrays stay bounded at any depth.
    for start in range(0, len(fragments), FRAGMENTS_PER_CHUNK):
        observations, overlap = fragment_observations(
            fragments[start : start + FRAGMENTS_PER_CHUNK],
            config.min_baseq,
            config.min_distance_from_end,
        )
        in_range = observations["position"] < reference_length
        histograms.add({name: values[in_range] for name, values in observations.items()})
        for key in overlap_totals:
            overlap_totals[key] += overlap[key]
    stats.update(overlap_totals)
    stats["counted_observations"] = int(histograms.depth().sum())
    return histograms, stats


def estimate_error_rates(
    normal: QualityHistograms, reference: str, max_real_allele_fraction: float
) -> dict[str, float]:
    """Per-substitution sequencing error rate, learned from the normal sample.

    A plain tumour-versus-normal Fisher test assumes both samples share an
    error rate, so unequal depth alone can look significant. Estimating the
    rate for each REF>ALT substitution gives the caller an absolute noise floor
    to test against as well.
    """
    per_allele = normal.allele_counts()
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


def callable_mask(
    reference: str,
    tumor: QualityHistograms,
    normal: QualityHistograms,
    config: PairedConfig,
    blacklist: set[int],
) -> np.ndarray:
    """Positions deep enough in both samples, off the blacklist, and off the edges."""
    length = len(reference)
    positions = np.arange(1, length + 1)
    callable_ = (tumor.depth() >= config.min_tumor_depth) & (
        normal.depth() >= config.min_normal_depth
    )
    callable_ &= ~np.isin(positions, list(blacklist))
    callable_ &= ~np.array([unresolved_edge(p, length, config) for p in positions], dtype=bool)
    return callable_


def _depth_summary(histograms: QualityHistograms) -> dict:
    depths = histograms.depth().astype(float)
    return {
        "breadth": float(np.count_nonzero(depths) / len(depths)),
        "depth_quantiles": {
            str(quantile): float(np.quantile(depths, quantile))
            for quantile in (0, 0.25, 0.5, 0.75, 0.9, 0.99, 1)
        },
    }


def run_paired_pipeline(config: PairedConfig) -> PairedResult:
    """Run the paired analysis and write the VCF, index, and callable BED."""
    chromosome, reference, checksum = load_fasta_reference(config.reference, config.mito_chr)
    config.mito_chr = chromosome
    blacklist = (
        load_bed_positions(config.custom_blacklist, chromosome)
        if config.custom_blacklist
        else set()
    )
    tumor, tumor_stats = collect_sample_evidence(config.tumor, config, len(reference))
    normal, normal_stats = collect_sample_evidence(config.normal, config, len(reference))
    callable_ = callable_mask(reference, tumor, normal, config, blacklist)
    error_rates = estimate_error_rates(normal, reference, config.max_normal_af)
    candidates = construct_candidates(
        chromosome, reference, tumor, normal, config, blacklist, error_rates
    )
    callable_positions = int(callable_.sum())
    qc = {
        "mgatk2_version": version("mgatk2"),
        "command_line": sys.argv,
        "parameters": asdict(config),
        "inputs": {
            "normal": {
                "path": str(Path(config.normal).resolve()),
                "type": Path(config.normal).suffix.lower().lstrip("."),
                "declared_upstream_consensus": config.input_is_consensus,
                "statistics": normal_stats,
                **_depth_summary(normal),
            },
            "tumor": {
                "path": str(Path(config.tumor).resolve()),
                "type": Path(config.tumor).suffix.lower().lstrip("."),
                "declared_upstream_consensus": config.input_is_consensus,
                "statistics": tumor_stats,
                **_depth_summary(tumor),
            },
        },
        "reference": {
            "path": str(Path(config.reference).resolve()),
            "chromosome": chromosome,
            "length": len(reference),
            "sha256": checksum,
        },
        "deduplication": config.deduplication,
        "snv_only": True,
        "indels_called": False,
        "blacklist_numt_strategy": (
            "user_chrM_blacklist_and_MAPQ"
            if config.custom_blacklist
            else "MAPQ_only_no_chrM_blacklist"
        ),
        "circular_edge_bases": config.circular_edge_bases,
        "counts": {
            "evidence_positions": len(reference),
            "callable_positions": callable_positions,
            "candidates": len(candidates),
            "pass_candidates": sum(row["filter"] == "PASS" for row in candidates),
        },
        "substitution_error_rates": error_rates,
        "error_rate_real_allele_exclusion": config.max_normal_af,
        "numt_strategy": (
            "autosomal_median_depth_and_MAPQ"
            if config.autosomal_median_depth is not None
            else "MAPQ_only"
        ),
    }
    outputs = write_paired_outputs(
        Path(config.output), config.sample_name, chromosome, callable_, candidates, qc
    )
    logger.info(
        "Paired analysis complete: %d positions, %d candidates, %d PASS",
        len(reference),
        len(candidates),
        qc["counts"]["pass_candidates"],
    )
    return PairedResult(
        outputs=outputs,
        evidence_positions=len(reference),
        candidates=len(candidates),
        pass_candidates=qc["counts"]["pass_candidates"],
        callable_positions=callable_positions,
    )
