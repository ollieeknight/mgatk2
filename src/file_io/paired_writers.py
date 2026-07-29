"""Atomic writers for paired mitochondrial evidence outputs."""

from __future__ import annotations

import csv
import gzip
import json
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
    "QUERY_COUNTED_FRAGMENTS",
    "BASELINE_COUNTED_FRAGMENTS",
    "QUERY_MEAN_MAPQ",
    "BASELINE_MEAN_MAPQ",
    "QUERY_MEAN_BASEQ",
    "BASELINE_MEAN_BASEQ",
    "QUERY_MEAN_READ_POSITION",
    "BASELINE_MEAN_READ_POSITION",
    "QUERY_CLIPPED_OBSERVATIONS",
    "BASELINE_CLIPPED_OBSERVATIONS",
    "QUERY_OVERLAP_DISAGREEMENTS",
    "BASELINE_OVERLAP_DISAGREEMENTS",
    "QUERY_CALLABLE",
    "BASELINE_CALLABLE",
]

CANDIDATE_COLUMNS = [
    "CHROM",
    "POS",
    "REF",
    "ALT",
    "QUERY_DEPTH",
    "BASELINE_DEPTH",
    "QUERY_REF_COUNT",
    "QUERY_ALT_COUNT",
    "BASELINE_REF_COUNT",
    "BASELINE_ALT_COUNT",
    "QUERY_FWD",
    "QUERY_REV",
    "BASELINE_FWD",
    "BASELINE_REV",
    "QUERY_AF",
    "BASELINE_AF",
    "QUERY_BASELINE_RATIO",
    "QUERY_AF_CI_LOW",
    "QUERY_AF_CI_HIGH",
    "BASELINE_AF_CI_LOW",
    "BASELINE_AF_CI_HIGH",
    "ENRICHMENT_P",
    "ENRICHMENT_Q",
    "STRAND_P",
    "LEGACY_FILTER",
    "FILTER",
]


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
    header.add_meta("fileformat", value="VCFv4.3")
    header.add_meta("source", value=f"mgatk2-{qc['mgatk2_version']}")
    header.add_meta("mgatk2_query", value=qc["inputs"]["query"]["path"])
    header.add_meta("mgatk2_baseline", value=qc["inputs"]["baseline"]["path"])
    header.contigs.add(qc["reference"]["chromosome"], length=qc["reference"]["length"])
    for identifier, description in {
        "LOW_QUERY_DEPTH": "Query depth is below the configured reportability threshold",
        "LOW_BASELINE_DEPTH": "Baseline depth is below the configured reportability threshold",
        "BLACKLIST": "Position overlaps the user-supplied mitochondrial blacklist",
        "CIRCULAR_EDGE_UNRESOLVED": "Linear-reference edge lacks shifted-reference evidence",
        "STRAND_BIAS": "Alternate evidence is associated with read strand",
        "NOT_SIGNIFICANT": "One-sided query enrichment BH q-value exceeds 0.05",
    }.items():
        header.filters.add(identifier, None, None, description)
    for identifier, number, value_type, description in (
        ("QUERY_BASELINE_RATIO", 1, "Float", "Smoothed finite query/baseline AF ratio"),
        ("ENRICHMENT_P", 1, "Float", "One-sided Fisher exact enrichment p-value"),
        ("ENRICHMENT_Q", 1, "Float", "Benjamini-Hochberg adjusted enrichment p-value"),
        ("LEGACY_FILTER", 1, "String", "Migration-era deterministic filter result"),
    ):
        header.info.add(identifier, number, value_type, description)
    header.formats.add("DP", 1, "Integer", "Counted A/C/G/T fragment observations")
    header.formats.add("AD", "R", "Integer", "Reference and alternate fragment depths")
    header.formats.add("AF", "A", "Float", "Alternate allele fraction")
    header.formats.add("F1R2", "R", "Integer", "Reference and alternate F1R2 observations")
    header.formats.add("F2R1", "R", "Integer", "Reference and alternate F2R1 observations")
    header.add_sample("QUERY")
    header.add_sample("BASELINE")


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
                if row["FILTER"] == "PASS":
                    record.filter.add("PASS")
                else:
                    for flag in row["FILTER"].split("|"):
                        record.filter.add(flag)
                record.info["QUERY_BASELINE_RATIO"] = row["QUERY_BASELINE_RATIO"]
                record.info["ENRICHMENT_P"] = row["ENRICHMENT_P"]
                record.info["ENRICHMENT_Q"] = row["ENRICHMENT_Q"]
                record.info["LEGACY_FILTER"] = row["LEGACY_FILTER"]
                for sample in ("QUERY", "BASELINE"):
                    ref_count = row[f"{sample}_REF_COUNT"]
                    alt_count = row[f"{sample}_ALT_COUNT"]
                    record.samples[sample]["DP"] = row[f"{sample}_DEPTH"]
                    record.samples[sample]["AD"] = (ref_count, alt_count)
                    record.samples[sample]["AF"] = (row[f"{sample}_AF"],)
                    record.samples[sample]["F1R2"] = (
                        row.get(f"{sample}_REF_F1R2", 0),
                        row.get(f"{sample}_ALT_F1R2", 0),
                    )
                    record.samples[sample]["F2R1"] = (
                        row.get(f"{sample}_REF_F2R1", 0),
                        row.get(f"{sample}_ALT_F2R1", 0),
                    )
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


def _write_legacy(path: Path, candidates: list[dict]) -> None:
    columns = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "NORMAL_DEPTH",
        "TUMOUR_DEPTH",
        "NORMAL_REF_COUNT",
        "NORMAL_ALT_COUNT",
        "NORMAL_FWD",
        "NORMAL_REV",
        "TUMOUR_REF_COUNT",
        "TUMOUR_ALT_COUNT",
        "TUMOUR_FWD",
        "TUMOUR_REV",
        "NORMAL_VAF",
        "TUMOUR_VAF",
        "TN_RATIO",
        "STRAND_BIAS",
        "FILTER",
    ]
    rows = []
    for row in candidates:
        total = row["QUERY_FWD"] + row["QUERY_REV"]
        rows.append(
            {
                "CHROM": row["CHROM"],
                "POS": row["POS"],
                "REF": row["REF"],
                "ALT": row["ALT"],
                "NORMAL_DEPTH": row["BASELINE_DEPTH"],
                "TUMOUR_DEPTH": row["QUERY_DEPTH"],
                "NORMAL_REF_COUNT": row["BASELINE_REF_COUNT"],
                "NORMAL_ALT_COUNT": row["BASELINE_ALT_COUNT"],
                "NORMAL_FWD": row["BASELINE_FWD"],
                "NORMAL_REV": row["BASELINE_REV"],
                "TUMOUR_REF_COUNT": row["QUERY_REF_COUNT"],
                "TUMOUR_ALT_COUNT": row["QUERY_ALT_COUNT"],
                "TUMOUR_FWD": row["QUERY_FWD"],
                "TUMOUR_REV": row["QUERY_REV"],
                "NORMAL_VAF": row["BASELINE_AF"],
                "TUMOUR_VAF": row["QUERY_AF"],
                "TN_RATIO": row["QUERY_BASELINE_RATIO"],
                "STRAND_BIAS": abs(row["QUERY_FWD"] - row["QUERY_REV"]) / total if total else 0,
                "FILTER": row["LEGACY_FILTER"],
            }
        )
    temporary = _temporary_path(path)
    try:
        with temporary.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary, path)
    finally:
        temporary.unlink(missing_ok=True)


def write_paired_outputs(
    output_dir: Path,
    sample_name: str,
    evidence: list[dict],
    candidates: list[dict],
    qc: dict,
    write_legacy: bool = False,
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
    if write_legacy:
        paths["legacy_tsv"] = output_dir / f"{sample_name}.mito_somatic.tsv"
        _write_legacy(paths["legacy_tsv"], candidates)
    with paths["log"].open("w") as handle:
        handle.write(
            f"mgatk2 paired: {len(evidence)} evidence positions, {len(candidates)} candidates\n"
        )
        if write_legacy:
            handle.write("Legacy compatibility projection written; not a primary callset.\n")
    qc["outputs"] = {key: str(value) for key, value in paths.items()}
    _write_json(paths["qc"], qc)
    return {key: str(value) for key, value in paths.items()}
