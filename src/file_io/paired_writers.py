"""Atomic writers for paired mitochondrial evidence outputs."""

from __future__ import annotations

import gzip
import json
import os
import tempfile
from pathlib import Path

import pysam

# Normal first, then tumour: the sample order every somatic VCF consumer
# expects, and the order the tables read in.
SAMPLES = ("normal", "tumor")

VCF_FILTERS = {
    "LOW_TUMOR_DEPTH": "Tumour depth is below the configured reportability threshold",
    "LOW_NORMAL_DEPTH": "Normal depth is below the configured reportability threshold",
    "LOW_ALT_OBSERVATIONS": "Fewer tumour alternate observations than required",
    "LOW_TUMOR_AF": "Tumour allele fraction is below the configured threshold",
    "HIGH_NORMAL_AF": "Normal allele fraction is above the configured threshold",
    "BLACKLIST": "Position overlaps the user-supplied mitochondrial blacklist",
    "CIRCULAR_EDGE_UNRESOLVED": "Linear-reference edge lacks shifted-reference evidence",
    "STRAND_BIAS": "Alternate evidence is associated with read strand",
    "ORIENTATION_BIAS": "Alternate evidence is associated with read-pair orientation",
    "WEAK_EVIDENCE": "Alternate support is consistent with the estimated substitution error rate",
    "NOT_SIGNIFICANT": "One-sided tumour enrichment BH q-value exceeds 0.05",
    "BASE_QUAL": "Alternate observations carry systematically lower base quality than reference",
    "MAP_QUAL": "Alternate observations carry systematically lower mapping quality than reference",
    "POSITION": "Alternate observations sit systematically closer to read ends than reference",
    "POSSIBLE_NUMT": "Alternate support is within the depth a single-copy NuMT would contribute",
}

VCF_INFO = (
    ("DP", 1, "Integer", "Total counted fragment observations across both samples"),
    ("ERR", 1, "Float", "Estimated per-substitution sequencing error rate"),
    ("SEQP", 1, "Float", "Binomial p-value for tumour alternate support against ERR"),
    ("EP", 1, "Float", "One-sided Fisher exact tumour enrichment p-value"),
    ("EQ", 1, "Float", "Benjamini-Hochberg adjusted enrichment p-value"),
    ("SP", 1, "Float", "Fisher exact strand-bias p-value in the tumour"),
    ("OP", 1, "Float", "Fisher exact read-orientation-bias p-value in the tumour"),
    ("MBQ", "R", "Integer", "Median tumour base quality of reference and alternate observations"),
    ("MMQ", "R", "Integer", "Median tumour mapping quality of reference and alternate"),
    ("MPOS", "R", "Integer", "Median tumour distance from read end for reference and alternate"),
    ("RSBQ", 1, "Float", "Alt-vs-ref base quality rank-sum z-score in the tumour"),
    ("RSBQ_P", 1, "Float", "Base quality rank-sum p-value in the tumour"),
    ("RSMQ", 1, "Float", "Alt-vs-ref mapping quality rank-sum z-score in the tumour"),
    ("RSMQ_P", 1, "Float", "Mapping quality rank-sum p-value in the tumour"),
    ("RSPOS", 1, "Float", "Alt-vs-ref read-end distance rank-sum z-score in the tumour"),
    ("RSPOS_P", 1, "Float", "Read-end distance rank-sum p-value in the tumour"),
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
    ("AFCI", 2, "Float", "Clopper-Pearson 95% interval for the alternate allele fraction"),
)


def _temporary_path(destination: Path, suffix: str = "") -> Path:
    handle, name = tempfile.mkstemp(
        prefix=f".{destination.name}.", suffix=suffix, dir=destination.parent
    )
    os.close(handle)
    return Path(name)


def _write_callable_bed(path: Path, evidence: list[dict]) -> None:
    """Positions reportable in both samples, as a bgzipped BED.

    The one output the VCF cannot stand in for: it carries variant sites only,
    so callable territory is unrecoverable after the run. Same companion file
    GATK CallableLoci and Strelka emit alongside a somatic VCF, and the
    denominator for any per-callable-base mutation burden.
    """
    source = _temporary_path(path, ".bed")
    compressed = Path(f"{source}.gz")
    try:
        with source.open("w") as handle:
            for row in evidence:
                if row["_joint_callable"]:
                    handle.write(f"{row['chrom']}\t{row['pos'] - 1}\t{row['pos']}\n")
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
    header.add_meta("tumor_sample", value="TUMOR")
    header.add_meta("normal_sample", value="NORMAL")
    header.add_meta("mgatk2_tumor", value=qc["inputs"]["tumor"]["path"])
    header.add_meta("mgatk2_normal", value=qc["inputs"]["normal"]["path"])
    # The whole QC record, on one unstructured meta line, because the VCF is
    # now the only artefact and provenance cannot live in a file that no longer
    # exists. Read it back with:
    #   bcftools view -h out.vcf.gz | sed -n 's/^##mgatk2_qc=//p' | jq
    header.add_meta("mgatk2_qc", value=json.dumps(qc, sort_keys=True, allow_nan=False))
    header.contigs.add(qc["reference"]["chromosome"], length=qc["reference"]["length"])
    for identifier, description in VCF_FILTERS.items():
        header.filters.add(identifier, None, None, description)
    for identifier, number, value_type, description in VCF_INFO:
        header.info.add(identifier, number, value_type, description)
    for identifier, number, value_type, description in VCF_FORMAT:
        header.formats.add(identifier, number, value_type, description)
    for sample in SAMPLES:
        header.add_sample(sample.upper())


def _write_vcf(path: Path, candidates: list[dict], qc: dict) -> None:
    source = _temporary_path(path, ".vcf")
    compressed = Path(f"{source}.gz")
    try:
        header = pysam.VariantHeader()
        _add_vcf_headers(header, qc)
        with pysam.VariantFile(source, "w", header=header) as output:
            for row in candidates:
                record = output.new_record(
                    contig=row["chrom"],
                    start=row["pos"] - 1,
                    stop=row["pos"],
                    alleles=(row["ref"], row["alt"]),
                )
                record.qual = row["qual"]
                for flag in row["filter"].split(";"):
                    record.filter.add(flag)

                record.info["DP"] = row["tumor_dp"] + row["normal_dp"]
                record.info["ERR"] = row["error_rate"]
                record.info["SEQP"] = row["seq_p"]
                record.info["EP"] = row["enrich_p"]
                record.info["EQ"] = row["enrich_q"]
                record.info["SP"] = row["strand_p"]
                record.info["OP"] = row["orient_p"]
                for metric in ("mbq", "mmq", "mpos"):
                    record.info[metric.upper()] = (
                        int(row[f"tumor_ref_{metric}"]),
                        int(row[f"tumor_alt_{metric}"]),
                    )
                for key in ("rsbq", "rsbq_p", "rsmq", "rsmq_p", "rspos", "rspos_p"):
                    record.info[key.upper()] = row[key]

                for sample in SAMPLES:
                    call = record.samples[sample.upper()]
                    # Somatic convention: the tumour carries the allele, the
                    # normal is the reference-state comparator.
                    call["GT"] = (0, 1) if sample == "tumor" else (0, 0)
                    call["DP"] = row[f"{sample}_dp"]
                    call["AD"] = (row[f"{sample}_ref_count"], row[f"{sample}_ac"])
                    call["AF"] = (row[f"{sample}_af"],)
                    call["SB"] = (
                        row[f"{sample}_ref_fwd"],
                        row[f"{sample}_ref_rev"],
                        row[f"{sample}_alt_fwd"],
                        row[f"{sample}_alt_rev"],
                    )
                    call["F1R2"] = (row[f"{sample}_ref_f1r2"], row[f"{sample}_alt_f1r2"])
                    call["F2R1"] = (row[f"{sample}_ref_f2r1"], row[f"{sample}_alt_f2r1"])
                    call["AFCI"] = (row[f"{sample}_af_ci_low"], row[f"{sample}_af_ci_high"])
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
        "vcf": output_dir / f"{sample_name}.mt_variants.vcf.gz",
        "vcf_index": output_dir / f"{sample_name}.mt_variants.vcf.gz.tbi",
        "callable_bed": output_dir / f"{sample_name}.mt_callable.bed.gz",
    }
    qc["outputs"] = {key: str(value) for key, value in paths.items()}
    _write_vcf(paths["vcf"], candidates, qc)
    _write_callable_bed(paths["callable_bed"], evidence)
    return {key: str(value) for key, value in paths.items()}
