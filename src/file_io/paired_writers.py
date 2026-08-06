"""Atomic writers for paired mitochondrial evidence outputs."""

from __future__ import annotations

import csv
import gzip
import json
import math
import os
import tempfile
from pathlib import Path

import pysam

EVIDENCE_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "QUERY_DEPTH",
    "BASELINE_DEPTH",
    *[
        f"{sample}_{base}_{strand}"
        for sample in ("QUERY", "BASELINE")
        for base in "ACGT"
        for strand in ("FWD", "REV")
    ],
    "QUERY_CLIPPED_OBSERVATIONS",
    "BASELINE_CLIPPED_OBSERVATIONS",
    "QUERY_OVERLAP_DISAGREEMENTS",
    "BASELINE_OVERLAP_DISAGREEMENTS",
    "QUERY_CALLABLE",
    "BASELINE_CALLABLE",
]

# Grouped by sample throughout, so the table reads in the same order as the
# VCF sample columns.
_SAMPLE_CANDIDATE_COLUMNS = [
    "{sample}_DEPTH",
    "{sample}_REF_COUNT",
    "{sample}_REF_FWD",
    "{sample}_REF_REV",
    "{sample}_ALT_COUNT",
    "{sample}_ALT_FWD",
    "{sample}_ALT_REV",
    "{sample}_AF",
    "{sample}_AF_CI_LOW",
    "{sample}_AF_CI_HIGH",
    "{sample}_REF_F1R2",
    "{sample}_REF_F2R1",
    "{sample}_ALT_F1R2",
    "{sample}_ALT_F2R1",
    "{sample}_REF_MBQ",
    "{sample}_REF_MMQ",
    "{sample}_REF_MPOS",
    "{sample}_ALT_MBQ",
    "{sample}_ALT_MMQ",
    "{sample}_ALT_MPOS",
]

CANDIDATE_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    *[column.format(sample="QUERY") for column in _SAMPLE_CANDIDATE_COLUMNS],
    *[column.format(sample="BASELINE") for column in _SAMPLE_CANDIDATE_COLUMNS],
    "QUERY_BASELINE_RATIO",
    "ERROR_RATE",
    "SEQUENCING_ERROR_P",
    "ENRICHMENT_P",
    "ENRICHMENT_Q",
    "STRAND_P",
    "ORIENTATION_P",
    "QUERY_RSBQ",
    "QUERY_RSBQ_P",
    "QUERY_RSMQ",
    "QUERY_RSMQ_P",
    "QUERY_RSPOS",
    "QUERY_RSPOS_P",
    "FILTER",
]

VCF_FILTERS = {
    "LOW_QUERY_DEPTH": "Query depth is below the configured reportability threshold",
    "LOW_BASELINE_DEPTH": "Baseline depth is below the configured reportability threshold",
    "LOW_ALT_OBSERVATIONS": "Fewer query alternate observations than required",
    "LOW_QUERY_AF": "Query allele fraction is below the configured threshold",
    "HIGH_BASELINE_AF": "Baseline allele fraction is above the configured threshold",
    "BLACKLIST": "Position overlaps the user-supplied mitochondrial blacklist",
    "CIRCULAR_EDGE_UNRESOLVED": "Linear-reference edge lacks shifted-reference evidence",
    "STRAND_BIAS": "Alternate evidence is associated with read strand",
    "ORIENTATION_BIAS": "Alternate evidence is associated with read-pair orientation",
    "WEAK_EVIDENCE": "Alternate support is consistent with the estimated substitution error rate",
    "NOT_SIGNIFICANT": "One-sided query enrichment BH q-value exceeds 0.05",
    "POSSIBLE_NUMT": "Alternate support is within the depth a single-copy NuMT would contribute",
}

VCF_INFO = (
    ("DP", 1, "Integer", "Total counted fragment observations across both samples"),
    ("QBRATIO", 1, "Float", "Smoothed finite query/baseline allele fraction ratio"),
    ("ERR", 1, "Float", "Estimated per-substitution sequencing error rate"),
    ("SEQP", 1, "Float", "Binomial p-value for query alternate support against ERR"),
    ("EP", 1, "Float", "One-sided Fisher exact query enrichment p-value"),
    ("EQ", 1, "Float", "Benjamini-Hochberg adjusted enrichment p-value"),
    ("SP", 1, "Float", "Fisher exact strand-bias p-value in the query"),
    ("OP", 1, "Float", "Fisher exact read-orientation-bias p-value in the query"),
    ("RSBQ", 1, "Float", "Alt-vs-ref base quality rank-sum z-score in the query"),
    ("RSBQ_P", 1, "Float", "Base quality rank-sum p-value in the query"),
    ("RSMQ", 1, "Float", "Alt-vs-ref mapping quality rank-sum z-score in the query"),
    ("RSMQ_P", 1, "Float", "Mapping quality rank-sum p-value in the query"),
    ("RSPOS", 1, "Float", "Alt-vs-ref read-end distance rank-sum z-score in the query"),
    ("RSPOS_P", 1, "Float", "Read-end distance rank-sum p-value in the query"),
)

VCF_FORMAT = (
    ("GT", 1, "String", "Genotype"),
    ("DP", 1, "Integer", "Counted A/C/G/T fragment observations"),
    ("AD", "R", "Integer", "Reference and alternate fragment depths"),
    ("AF", "A", "Float", "Alternate allele fraction"),
    (
        "SB",
        4,
        "Integer",
        "Reference forward, reference reverse, alternate forward, alternate reverse",
    ),
    ("F1R2", "R", "Integer", "Reference and alternate F1R2 observations"),
    ("F2R1", "R", "Integer", "Reference and alternate F2R1 observations"),
    ("MBQ", "R", "Integer", "Median base quality of reference and alternate observations"),
    ("MMQ", "R", "Integer", "Median mapping quality of reference and alternate observations"),
    ("MPOS", "R", "Integer", "Median distance from read end for reference and alternate"),
    ("AFCI", 2, "Float", "Clopper-Pearson 95% interval for the alternate allele fraction"),
)


def _temporary_path(destination: Path, suffix: str = "") -> Path:
    handle, name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=suffix, dir=destination.parent
    )
    os.close(handle)
    return Path(name)


def _write_gzip_tsv(path: Path, columns: list[str], rows: list[dict]) -> None:
    temporary = _temporary_path(path)
    try:
        with gzip.open(temporary, "wt", newline="") as handle:
            writer = csv.DictWriter(
                handle, fieldnames=columns, delimiter="\t", extrasaction="ignore"
            )
            writer.writeheader()
            writer.writerows(rows)
        with gzip.open(temporary, "rt") as handle:
            if sum(1 for _line in handle) - 1 != len(rows):
                raise OSError(f"Row-count validation failed for {path}")
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_json(path: Path, value: dict) -> None:
    temporary = _temporary_path(path)
    try:
        with temporary.open("w") as handle:
            json.dump(value, handle, indent=2, sort_keys=True, allow_nan=False)
            handle.write("\n")
        json.loads(temporary.read_text())
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def _write_callable_bed(path: Path, rows: list[dict]) -> None:
    source = _temporary_path(path, ".bed")
    compressed = Path(f"{source}.gz")
    try:
        with source.open("w") as handle:
            for row in rows:
                if row["_JOINT_CALLABLE"]:
                    handle.write(f"{row['CHROM']}\t{row['POS'] - 1}\t{row['POS']}\n")
        pysam.tabix_compress(str(source), str(compressed), force=True)
        with gzip.open(compressed, "rt") as handle:
            list(handle)
        os.replace(compressed, path)
    finally:
        source.unlink(missing_ok=True)
        compressed.unlink(missing_ok=True)


def _add_vcf_headers(header: pysam.VariantHeader, qc: dict) -> None:
    # pysam already emits ##fileformat; a second one would violate VCF 4.3 §1.1.
    header.add_meta("source", value=f"mgatk2-{qc['mgatk2_version']}")
    header.add_meta("reference", value=f"file://{qc['reference']['path']}")
    header.add_meta("query_sample", value="QUERY")
    header.add_meta("baseline_sample", value="BASELINE")
    header.add_meta("mgatk2_query", value=qc["inputs"]["query"]["path"])
    header.add_meta("mgatk2_baseline", value=qc["inputs"]["baseline"]["path"])
    header.contigs.add(qc["reference"]["chromosome"], length=qc["reference"]["length"])
    for identifier, description in VCF_FILTERS.items():
        header.filters.add(identifier, None, None, description)
    for identifier, number, value_type, description in VCF_INFO:
        header.info.add(identifier, number, value_type, description)
    for identifier, number, value_type, description in VCF_FORMAT:
        header.formats.add(identifier, number, value_type, description)
    header.add_sample("QUERY")
    header.add_sample("BASELINE")


def _phred(probability: float) -> float:
    """Convert an adjusted p-value into a capped phred-scaled QUAL."""
    if probability <= 0:
        return 1000.0
    return min(1000.0, max(0.0, -10 * math.log10(probability)))


def _write_vcf(path: Path, candidates: list[dict], qc: dict) -> None:
    source = _temporary_path(path, ".vcf")
    compressed = Path(f"{source}.gz")
    try:
        header = pysam.VariantHeader()
        _add_vcf_headers(header, qc)
        with pysam.VariantFile(source, "w", header=header) as output:
            for row in candidates:
                record = output.new_record(
                    contig=row["CHROM"],
                    start=row["POS"] - 1,
                    stop=row["POS"],
                    alleles=(row["REF"], row["ALT"]),
                )
                record.qual = _phred(row["ENRICHMENT_Q"])
                for flag in row["FILTER"].split(";"):
                    record.filter.add(flag)

                record.info["DP"] = row["QUERY_DEPTH"] + row["BASELINE_DEPTH"]
                record.info["QBRATIO"] = row["QUERY_BASELINE_RATIO"]
                record.info["ERR"] = row["ERROR_RATE"]
                record.info["SEQP"] = row["SEQUENCING_ERROR_P"]
                record.info["EP"] = row["ENRICHMENT_P"]
                record.info["EQ"] = row["ENRICHMENT_Q"]
                record.info["SP"] = row["STRAND_P"]
                record.info["OP"] = row["ORIENTATION_P"]
                for key in ("RSBQ", "RSBQ_P", "RSMQ", "RSMQ_P", "RSPOS", "RSPOS_P"):
                    record.info[key] = row[f"QUERY_{key}"]

                for sample in ("QUERY", "BASELINE"):
                    reference_count = row[f"{sample}_REF_COUNT"]
                    alternate_count = row[f"{sample}_ALT_COUNT"]
                    call = record.samples[sample]
                    # Somatic convention: the query carries the allele, the
                    # baseline is the reference-state comparator.
                    call["GT"] = (0, 1) if sample == "QUERY" else (0, 0)
                    call["DP"] = row[f"{sample}_DEPTH"]
                    call["AD"] = (reference_count, alternate_count)
                    call["AF"] = (row[f"{sample}_AF"],)
                    call["SB"] = (
                        row[f"{sample}_REF_FWD"],
                        row[f"{sample}_REF_REV"],
                        row[f"{sample}_ALT_FWD"],
                        row[f"{sample}_ALT_REV"],
                    )
                    call["F1R2"] = (row[f"{sample}_REF_F1R2"], row[f"{sample}_ALT_F1R2"])
                    call["F2R1"] = (row[f"{sample}_REF_F2R1"], row[f"{sample}_ALT_F2R1"])
                    for metric in ("MBQ", "MMQ", "MPOS"):
                        call[metric] = (
                            int(row[f"{sample}_REF_{metric}"]),
                            int(row[f"{sample}_ALT_{metric}"]),
                        )
                    call["AFCI"] = (row[f"{sample}_AF_CI_LOW"], row[f"{sample}_AF_CI_HIGH"])
                output.write(record)
        pysam.tabix_compress(str(source), str(compressed), force=True)
        pysam.tabix_index(str(compressed), preset="vcf", force=True)
        with pysam.VariantFile(compressed) as check:
            if sum(1 for _record in check) != len(candidates):
                raise OSError(f"Record-count validation failed for {path}")
        os.replace(compressed, path)
        os.replace(Path(f"{compressed}.tbi"), Path(f"{path}.tbi"))
    finally:
        source.unlink(missing_ok=True)
        compressed.unlink(missing_ok=True)
        Path(f"{compressed}.tbi").unlink(missing_ok=True)


def write_paired_outputs(
    output_dir: Path,
    sample_name: str,
    evidence: list[dict],
    candidates: list[dict],
    qc: dict,
) -> dict[str, str]:
    """Write and validate the complete paired output contract."""
    output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "evidence": output_dir / f"{sample_name}.mt_evidence.tsv.gz",
        "candidates": output_dir / f"{sample_name}.mt_candidates.tsv.gz",
        "vcf": output_dir / f"{sample_name}.mt_variants.vcf.gz",
        "vcf_index": output_dir / f"{sample_name}.mt_variants.vcf.gz.tbi",
        "callable_bed": output_dir / f"{sample_name}.mt_callable.bed.gz",
        "qc": output_dir / f"{sample_name}.mt_qc.json",
        "log": output_dir / f"{sample_name}.paired.log",
    }
    _write_gzip_tsv(paths["evidence"], EVIDENCE_COLUMNS, evidence)
    _write_gzip_tsv(paths["candidates"], CANDIDATE_COLUMNS, candidates)
    _write_vcf(paths["vcf"], candidates, qc)
    _write_callable_bed(paths["callable_bed"], evidence)
    with paths["log"].open("w") as handle:
        handle.write(
            f"mgatk2 paired: {len(evidence)} evidence positions, {len(candidates)} candidates\n"
        )
    qc["outputs"] = {key: str(value) for key, value in paths.items()}
    _write_json(paths["qc"], qc)
    return {key: str(value) for key, value in paths.items()}
