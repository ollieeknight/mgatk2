import gzip
import json

import numpy as np
import pysam
import pytest
from click.testing import CliRunner

from analysis.paired_calling import benjamini_hochberg, construct_candidates
from analysis.quality_stats import QualityHistograms, histogram_median, rank_sum
from cli import cli
from core.config import PairedConfig
from core.exceptions import InvalidInputError
from processing.paired_pileup import _git_commit, run_paired_pipeline


def _config(paired_files, output, tumor="tumor_bam", normal="normal_bam", **kwargs):
    values = {
        "tumor": str(paired_files[tumor]),
        "normal": str(paired_files[normal]),
        "reference": str(paired_files["reference"]),
        "output": str(output),
        "sample_name": "pair",
        "min_distance_from_end": 0,
        "min_tumor_depth": 1,
        "min_normal_depth": 1,
        "circular_edge_bases": 2,
    }
    values.update(kwargs)
    return PairedConfig(**values)


def _evidence(tumor_alt, tumor_depth, normal_alt, normal_depth):
    row = {
        "chrom": "chrM",
        "pos": 1,
        "ref": "A",
        "tumor_dp": tumor_depth,
        "normal_dp": normal_depth,
    }
    for sample in ("tumor", "normal"):
        for base in "acgt":
            for strand in ("fwd", "rev"):
                row[f"{sample}_{base}_{strand}"] = 0
    row["tumor_c_fwd"] = tumor_alt
    row["tumor_a_fwd"] = tumor_depth - tumor_alt
    row["normal_c_fwd"] = normal_alt
    row["normal_a_fwd"] = normal_depth - normal_alt
    return row


def _candidates(rows, config, blacklist=frozenset(), error_rates=None):
    """Call construct_candidates with empty histograms; counts come from rows."""
    length = max(row["pos"] for row in rows)
    return construct_candidates(
        rows,
        QualityHistograms(length),
        QualityHistograms(length),
        config,
        set(blacklist),
        error_rates or {},
    )


def _qc(vcf_path):
    """Read the QC record back out of the VCF header, the only place it lives."""
    with pysam.VariantFile(vcf_path) as vcf:
        for record in vcf.header.records:
            if record.key == "mgatk2_qc":
                return json.loads(record.value)
    raise AssertionError("no mgatk2_qc header record")


def _paired_args(paired_files, output):
    return [
        "paired",
        "--tumor",
        str(paired_files["tumor_bam"]),
        "--normal",
        str(paired_files["normal_bam"]),
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
    shallow = _candidates([_evidence(3, 10, 0, 5)], config)[0]
    deep = _candidates([_evidence(3, 10, 0, 500)], config)[0]

    assert shallow["enrich_p"] > deep["enrich_p"]
    assert shallow["normal_af_ci_high"] > deep["normal_af_ci_high"]

    row = _evidence(3, 10, 0, 10)
    row["tumor_g_fwd"] = 2
    row["tumor_a_fwd"] = 5

    candidate = _candidates([row], config)[0]

    assert candidate["tumor_ref_count"] == 5
    assert candidate["tumor_ref_count"] != candidate["tumor_dp"] - candidate["tumor_ac"]


def test_filters_have_a_stable_order(paired_files, tmp_path):
    config = _config(
        paired_files,
        tmp_path,
        min_tumor_depth=10,
        min_normal_depth=5,
        circular_edge_bases=0,
    )
    row = _candidates([_evidence(1, 2, 1, 2)], config, blacklist={1})[0]

    assert row["filter"].split(";") == [
        "LOW_TUMOR_DEPTH",
        "LOW_NORMAL_DEPTH",
        "LOW_ALT_OBSERVATIONS",
        "HIGH_NORMAL_AF",
        "BLACKLIST",
        "STRAND_BIAS",
        "NOT_SIGNIFICANT",
    ]


def test_sequencing_error_rate_gates_weak_alternate_support(paired_files, tmp_path):
    config = _config(paired_files, tmp_path, min_tumor_af=0.0, min_alt_observations=1)
    # 5 alt reads in 1000 is 0.5%: noise at a 1% error rate, signal at 1e-6.
    row = _evidence(5, 1000, 0, 1000)

    noisy = _candidates([row], config, error_rates={"A>C": 0.01})[0]
    clean = _candidates([row], config, error_rates={"A>C": 1e-6})[0]

    assert noisy["seq_p"] > clean["seq_p"]
    assert "WEAK_EVIDENCE" in noisy["filter"]
    assert "WEAK_EVIDENCE" not in clean["filter"]


def test_numt_filter_needs_autosomal_depth(paired_files, tmp_path):
    row = _evidence(20, 1000, 0, 1000)

    without = _candidates([row], _config(paired_files, tmp_path))[0]
    with_depth = _candidates(
        [row], _config(paired_files, tmp_path / "numt", autosomal_median_depth=30.0)
    )[0]

    assert "POSSIBLE_NUMT" not in without["filter"]
    assert "POSSIBLE_NUMT" in with_depth["filter"]


def test_histogram_median_and_rank_sum_separate_distributions():
    low = np.zeros(96, dtype=np.int32)
    low[10] = 40
    high = np.zeros(96, dtype=np.int32)
    high[40] = 40

    assert histogram_median(low) == 10
    assert histogram_median(low, scale=2) == 20

    z_score, p_value = rank_sum(low, high)
    assert z_score < 0 and p_value < 0.001
    assert rank_sum(low, low) == (0.0, 1.0)
    assert rank_sum(low, np.zeros(96, dtype=np.int32)) == (0.0, 1.0)


def test_paired_dry_run(paired_files, tmp_path):
    result = CliRunner().invoke(cli, [*_paired_args(paired_files, tmp_path), "--dry-run"])

    assert result.exit_code == 0
    assert "dry run complete" in result.output
    assert not (tmp_path / "pair.mt_variants.vcf.gz").exists()


def test_paired_dry_run_rejects_a_missing_contig(paired_files, tmp_path):
    result = CliRunner().invoke(
        cli, [*_paired_args(paired_files, tmp_path), "--dry-run", "-g", "chrM_absent"]
    )

    assert result.exit_code != 0
    assert "is absent from" in result.output


def test_paired_rejects_the_same_input_twice(paired_files, tmp_path):
    arguments = _paired_args(paired_files, tmp_path)
    arguments[4] = str(paired_files["tumor_bam"])

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
        _config(paired_files, tmp_path / "cram", tumor="tumor_cram", normal="normal_cram")
    )

    assert (bam.evidence_positions, bam.candidates, bam.callable_positions) == (
        cram.evidence_positions,
        cram.candidates,
        cram.callable_positions,
    )
    with pysam.VariantFile(bam.outputs["vcf"]) as bam_vcf:
        with pysam.VariantFile(cram.outputs["vcf"]) as cram_vcf:
            assert [str(record) for record in bam_vcf] == [str(record) for record in cram_vcf]


def test_reference_length_must_match_the_alignment(paired_files, tmp_path):
    reference = tmp_path / "wrong.fa"
    reference.write_text(">chrM\n" + "A" * 41 + "\n")
    pysam.faidx(str(reference))
    config = _config(paired_files, tmp_path / "out")
    config.reference = str(reference)

    with pytest.raises(InvalidInputError, match="FASTA length"):
        run_paired_pipeline(config)


def test_excluded_reads_are_reported(paired_files, alignment_factory, tmp_path):
    tumor = alignment_factory(
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
    config.tumor = str(tumor)

    result = run_paired_pipeline(config)
    stats = _qc(result.outputs["vcf"])["inputs"]["tumor"]["statistics"]

    assert stats["preexisting_duplicate_reads"] == 1
    assert stats["qc_failed_reads"] == 1
    assert stats["secondary_reads"] == 1
    assert stats["supplementary_reads"] == 1
    assert stats["low_mapq_reads"] == 1
    assert stats["missing_quality_reads"] == 1
    assert stats["retained_reads"] == 1


def test_identical_inputs_have_no_pass_candidates(paired_files, alignment_factory, tmp_path):
    normal = alignment_factory(
        tmp_path / "tumor-copy.bam",
        paired_files["reference"],
        [
            {"name": f"tumor{index}", "start": start, "sequence": sequence}
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
    config.normal = str(normal)

    assert run_paired_pipeline(config).pass_candidates == 0


def test_provenance_tolerates_missing_git(monkeypatch):
    def missing_git(*args, **kwargs):
        raise FileNotFoundError

    monkeypatch.setattr("processing.paired_pileup.subprocess.run", missing_git)

    assert _git_commit() is None


def test_outputs_are_valid_and_repeatable(paired_files, tmp_path):
    output = tmp_path / "out"
    config = _config(paired_files, output)
    result = run_paired_pipeline(config)

    assert set(result.outputs) == {"vcf", "vcf_index", "callable_bed"}
    assert sorted(path.name for path in output.iterdir()) == [
        "pair.mt_callable.bed.gz",
        "pair.mt_variants.vcf.gz",
        "pair.mt_variants.vcf.gz.tbi",
    ]
    with pysam.VariantFile(result.outputs["vcf"]) as vcf:
        assert list(vcf)

    qc = _qc(result.outputs["vcf"])
    assert qc["schema_version"] == "3.0"
    assert qc["reference"]["sha256"]
    assert qc["snv_only"] is True
    assert qc["counts"]["evidence_positions"] == 40
    assert qc["counts"]["callable_positions"] == result.callable_positions
    with gzip.open(result.outputs["callable_bed"], "rt") as bed:
        intervals = [line.split() for line in bed]
    assert len(intervals) == result.callable_positions
    assert all(start == str(int(end) - 1) for _chrom, start, end in intervals)

    rerun = run_paired_pipeline(config)
    assert rerun.outputs == result.outputs


def test_rank_sum_filters_fire_only_on_degraded_alternates(paired_files, tmp_path):
    """A low-quality alternate is flagged; a high-quality one is not."""
    config = _config(paired_files, tmp_path, min_tumor_af=0.0, min_alt_observations=1)
    row = _evidence(20, 100, 0, 100)

    def _candidate(alternate_bin):
        tumor = QualityHistograms(1)
        for source in (tumor.baseq, tumor.mapq, tumor.distance):
            source[0, 0, 30] = 80  # reference allele A
            source[0, 1, alternate_bin] = 20  # alternate allele C
        return construct_candidates(
            [row], tumor, QualityHistograms(1), config, set(), {"A>C": 1e-6}
        )[0]

    degraded = _candidate(5)
    healthy = _candidate(50)

    assert degraded["rsbq"] < 0
    for flag in ("BASE_QUAL", "MAP_QUAL", "POSITION"):
        assert flag in degraded["filter"]
        assert flag not in healthy["filter"]


def test_every_paired_option_is_accepted_by_the_command(paired_files, tmp_path):
    """The full option surface, exercised through the CLI rather than the config."""
    blacklist = tmp_path / "chrM.bed"
    blacklist.write_text("chrM\t0\t1\n")
    output = tmp_path / "cli-out"

    result = CliRunner().invoke(
        cli,
        [
            "paired",
            "--tumor",
            str(paired_files["tumor_bam"]),
            "--normal",
            str(paired_files["normal_bam"]),
            "--reference",
            str(paired_files["reference"]),
            "--output",
            str(output),
            "--sample-name",
            "pair",
            "--genome",
            "M",
            "--quality",
            "0",
            "--mapq",
            "0",
            "--min-distance-from-end",
            "0",
            "--max-strand-bias",
            "1.0",
            "--deduplication",
            "none",
            "--min-tumor-depth",
            "1",
            "--min-normal-depth",
            "1",
            "--min-alt-observations",
            "1",
            "--min-tumor-af",
            "0.0",
            "--max-normal-af",
            "0.5",
            "--custom-blacklist",
            str(blacklist),
            "--autosomal-median-depth",
            "0.5",
            "--circular-edge-bases",
            "2",
            "--input-is-consensus",
            "--verbose",
        ],
    )

    assert result.exit_code == 0, result.output
    assert (output / "pair.mt_variants.vcf.gz").exists()
    assert (output / "pair.mt_variants.vcf.gz.tbi").exists()
    assert (output / "pair.mt_callable.bed.gz").exists()

    qc = _qc(str(output / "pair.mt_variants.vcf.gz"))
    # --genome M must reach the reference as the FASTA's own contig name.
    assert qc["reference"]["chromosome"] == "chrM"
    assert qc["deduplication"] == "none"
    assert qc["numt_strategy"] == "autosomal_median_depth_and_MAPQ"
    assert qc["blacklist_numt_strategy"] == "user_chrM_blacklist_and_MAPQ"
    # The error-rate exclusion follows --max-normal-af rather than a constant.
    assert qc["error_rate_real_allele_exclusion"] == 0.5


def test_consensus_input_requires_deduplication_none(paired_files, tmp_path):
    result = CliRunner().invoke(
        cli,
        [
            "paired",
            "--tumor",
            str(paired_files["tumor_bam"]),
            "--normal",
            str(paired_files["normal_bam"]),
            "--reference",
            str(paired_files["reference"]),
            "--output",
            str(tmp_path / "out"),
            "--sample-name",
            "pair",
            "--input-is-consensus",
        ],
    )

    assert result.exit_code != 0
    assert "deduplication none" in result.output


def test_every_declared_vcf_field_is_written(paired_files, tmp_path):
    """Anything declared in the header but never populated is a silent gap."""
    from file_io.paired_writers import VCF_FORMAT, VCF_INFO

    result = run_paired_pipeline(_config(paired_files, tmp_path / "out", min_alt_observations=1))

    with pysam.VariantFile(result.outputs["vcf"]) as vcf:
        records = list(vcf)
    assert records

    for record in records:
        for identifier, *_ in VCF_INFO:
            assert identifier in record.info, identifier
        for sample in record.samples.values():
            for identifier, *_ in VCF_FORMAT:
                assert sample.get(identifier) is not None, identifier
