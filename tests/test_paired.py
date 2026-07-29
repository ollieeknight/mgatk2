import gzip
import json

import pysam
import pytest
from click.testing import CliRunner

from analysis.paired_calling import benjamini_hochberg, construct_candidates
from cli import cli
from core.config import PairedConfig
from core.exceptions import InvalidInputError
from processing.paired_pileup import _git_commit, run_paired_pipeline


def _config(paired_files, output, query="query_bam", baseline="baseline_bam", **kwargs):
    values = {
        "query": str(paired_files[query]),
        "baseline": str(paired_files[baseline]),
        "reference": str(paired_files["reference"]),
        "output": str(output),
        "sample_name": "pair",
        "min_distance_from_end": 0,
        "min_query_depth": 1,
        "min_baseline_depth": 1,
        "circular_edge_bases": 2,
    }
    values.update(kwargs)
    return PairedConfig(**values)


def _evidence(query_alt, query_depth, baseline_alt, baseline_depth):
    row = {
        "CHROM": "chrM",
        "POS": 1,
        "REF": "A",
        "QUERY_DEPTH": query_depth,
        "BASELINE_DEPTH": baseline_depth,
    }
    for sample in ("QUERY", "BASELINE"):
        for base in "ACGT":
            for strand in ("FWD", "REV"):
                row[f"{sample}_{base}_{strand}"] = 0
            for orientation in ("F1R2", "F2R1"):
                row[f"{sample}_{base}_{orientation}"] = 0
    row["QUERY_C_FWD"] = query_alt
    row["QUERY_A_FWD"] = query_depth - query_alt
    row["BASELINE_C_FWD"] = baseline_alt
    row["BASELINE_A_FWD"] = baseline_depth - baseline_alt
    return row


def _paired_args(paired_files, output):
    return [
        "paired",
        "--query",
        str(paired_files["query_bam"]),
        "--baseline",
        str(paired_files["baseline_bam"]),
        "--reference",
        str(paired_files["reference"]),
        "--output",
        str(output),
        "--sample-name",
        "pair",
    ]


def test_benjamini_hochberg_correction():
    assert benjamini_hochberg([0.01, 0.04, 0.03]) == [0.03, 0.04, 0.04]


def test_candidate_counts_and_uncertainty(paired_files, tmp_path):
    config = _config(paired_files, tmp_path)
    shallow = construct_candidates([_evidence(3, 10, 0, 5)], config, set())[0]
    deep = construct_candidates([_evidence(3, 10, 0, 500)], config, set())[0]

    assert shallow["ENRICHMENT_P"] > deep["ENRICHMENT_P"]
    assert shallow["BASELINE_AF_CI_HIGH"] > deep["BASELINE_AF_CI_HIGH"]

    row = _evidence(3, 10, 0, 10)
    row["QUERY_G_FWD"] = 2
    row["QUERY_A_FWD"] = 5

    candidate = construct_candidates([row], config, set())[0]

    assert candidate["QUERY_REF_COUNT"] == 5
    assert candidate["QUERY_REF_COUNT"] != candidate["QUERY_DEPTH"] - candidate["QUERY_ALT_COUNT"]


def test_legacy_filters_have_a_stable_order(paired_files, tmp_path):
    config = _config(
        paired_files,
        tmp_path,
        min_query_depth=10,
        min_baseline_depth=5,
    )
    row = construct_candidates([_evidence(1, 2, 1, 2)], config, {1})[0]

    assert row["LEGACY_FILTER"].split("|") == [
        "LOW_QUERY_DEPTH",
        "LOW_BASELINE_DEPTH",
        "LOW_ALT_OBSERVATIONS",
        "HIGH_BASELINE_AF",
        "LOW_QUERY_BASELINE_RATIO",
        "STRAND_BIAS",
        "BLACKLIST",
    ]


def test_paired_dry_run(paired_files, tmp_path):
    result = CliRunner().invoke(cli, [*_paired_args(paired_files, tmp_path), "--dry-run"])

    assert result.exit_code == 0
    assert "dry run complete" in result.output
    assert not (tmp_path / "pair.mt_qc.json").exists()


def test_paired_rejects_the_same_input_twice(paired_files, tmp_path):
    arguments = _paired_args(paired_files, tmp_path)
    arguments[4] = str(paired_files["query_bam"])

    result = CliRunner().invoke(cli, [*arguments, "--dry-run"])

    assert result.exit_code != 0
    assert "must be different" in result.output


def test_fasta_sets_ref_and_all_positions_are_written(paired_files, tmp_path):
    result = run_paired_pipeline(_config(paired_files, tmp_path))

    assert result.evidence_positions == 40
    assert result.candidates == 1


def test_bam_and_cram_give_the_same_counts(paired_files, tmp_path):
    bam = run_paired_pipeline(_config(paired_files, tmp_path / "bam"))
    cram = run_paired_pipeline(
        _config(paired_files, tmp_path / "cram", query="query_cram", baseline="baseline_cram")
    )

    assert (bam.evidence_positions, bam.candidates, bam.callable_positions) == (
        cram.evidence_positions,
        cram.candidates,
        cram.callable_positions,
    )
    with gzip.open(bam.outputs["evidence"], "rt") as bam_evidence:
        with gzip.open(cram.outputs["evidence"], "rt") as cram_evidence:
            assert bam_evidence.read() == cram_evidence.read()


def test_reference_length_must_match_the_alignment(paired_files, tmp_path):
    reference = tmp_path / "wrong.fa"
    reference.write_text(">chrM\n" + "A" * 41 + "\n")
    pysam.faidx(str(reference))
    config = _config(paired_files, tmp_path / "out")
    config.reference = str(reference)

    with pytest.raises(InvalidInputError, match="FASTA length"):
        run_paired_pipeline(config)


def test_excluded_reads_are_reported(paired_files, alignment_factory, tmp_path):
    query = alignment_factory(
        tmp_path / "filtered.bam",
        paired_files["reference"],
        [
            {"name": "kept", "start": 1},
            {"name": "duplicate", "start": 2, "flag": 1024},
            {"name": "qcfail", "start": 3, "flag": 512},
            {"name": "secondary", "start": 4, "flag": 256},
            {"name": "supplementary", "start": 5, "flag": 2048},
            {"name": "lowmapq", "start": 6, "mapq": 10},
            {"name": "missingqual", "start": 7, "qualities": None},
        ],
    )
    config = _config(paired_files, tmp_path / "policy")
    config.query = str(query)

    run_paired_pipeline(config)
    stats = json.loads((tmp_path / "policy/pair.mt_qc.json").read_text())["inputs"]["query"][
        "statistics"
    ]

    assert stats["preexisting_duplicate_reads"] == 1
    assert stats["qc_failed_reads"] == 1
    assert stats["secondary_reads"] == 1
    assert stats["supplementary_reads"] == 1
    assert stats["low_mapq_reads"] == 1
    assert stats["missing_quality_reads"] == 1
    assert stats["retained_reads"] == 1


def test_identical_inputs_have_no_pass_candidates(paired_files, alignment_factory, tmp_path):
    baseline = alignment_factory(
        tmp_path / "query-copy.bam",
        paired_files["reference"],
        [
            {"name": f"query{index}", "start": start, "sequence": sequence}
            for index, (start, sequence) in enumerate(
                (
                    (3, "AAAAAAACAAAAAAA"),
                    (4, "AAAAAACAAAAAAAA"),
                    (5, "AAAAACAAAAAAAAA"),
                )
            )
        ],
    )
    config = _config(paired_files, tmp_path / "negative")
    config.baseline = str(baseline)

    assert run_paired_pipeline(config).pass_candidates == 0


def test_provenance_tolerates_missing_git(monkeypatch):
    def missing_git(*args, **kwargs):
        raise FileNotFoundError

    monkeypatch.setattr("processing.paired_pileup.subprocess.run", missing_git)

    assert _git_commit() is None


def test_outputs_are_valid_and_repeatable(paired_files, tmp_path):
    config = _config(paired_files, tmp_path, write_legacy_tsv=True)
    result = run_paired_pipeline(config)

    assert set(result.outputs) == {
        "evidence",
        "candidates",
        "vcf",
        "vcf_index",
        "callable_bed",
        "qc",
        "log",
        "legacy_tsv",
    }
    with gzip.open(result.outputs["evidence"], "rt") as evidence:
        assert sum(1 for _line in evidence) == 41
    with pysam.VariantFile(result.outputs["vcf"]) as vcf:
        assert list(vcf)

    qc = json.loads((tmp_path / "pair.mt_qc.json").read_text())
    assert qc["schema_version"] == "1.0"
    assert qc["reference"]["sha256"]
    assert qc["snv_only"] is True

    rerun = run_paired_pipeline(config)
    assert (tmp_path / "pair.paired.log").read_text().count("mgatk2 paired:") == 1
    assert rerun.outputs == result.outputs


def test_wes_adapter_warns_and_writes_legacy_output(paired_files, tmp_path):
    result = CliRunner().invoke(
        cli,
        [
            "wes",
            "--tumour",
            str(paired_files["query_bam"]),
            "--normal",
            str(paired_files["baseline_bam"]),
            "--reference",
            str(paired_files["reference"]),
            "--output",
            str(tmp_path),
            "--sample-name",
            "legacy",
            "--min-distance-from-end",
            "0",
        ],
    )

    assert result.exit_code == 0, result.output
    assert result.output.count("deprecated") == 1
    assert (tmp_path / "legacy.mito_somatic.tsv").exists()
    assert (tmp_path / "legacy.mt_evidence.tsv.gz").exists()
