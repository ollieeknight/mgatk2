import gzip
import json
from importlib.metadata import version
from pathlib import Path

import numpy as np
import pysam
import pytest
from click.testing import CliRunner

from analysis.paired_calling import benjamini_hochberg, construct_candidates
from analysis.quality_stats import QualityHistograms, histogram_median, rank_sum
from cli import cli
from core.config import PairedConfig
from core.exceptions import InvalidInputError
from processing.paired_pileup import run_paired_pipeline


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


def _histograms(alt, depth):
    """One-position histograms with `alt` forward C reads and the rest forward A."""
    histograms = QualityHistograms(1)
    histograms.counts[0, 1, 0] = alt
    histograms.counts[0, 0, 0] = depth - alt
    return histograms


def _candidate(config, tumor, normal, blacklist=(), error_rates=None):
    """The A>C candidate at position 1 of a one-base reference."""
    return construct_candidates(
        "chrM", "A", tumor, normal, config, set(blacklist), error_rates or {}
    )[0]


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
    shallow = _candidate(config, _histograms(3, 10), _histograms(0, 5))
    deep = _candidate(config, _histograms(3, 10), _histograms(0, 500))

    assert shallow["enrich_p"] > deep["enrich_p"]
    assert shallow["normal_af_ci_high"] > deep["normal_af_ci_high"]

    # A third allele must not be counted as reference.
    tumor = _histograms(3, 8)
    tumor.counts[0, 2, 0] = 2
    candidate = _candidate(config, tumor, _histograms(0, 10))

    assert candidate["tumor_ref_count"] == 5
    assert candidate["tumor_dp"] == 10


def test_filters_have_a_stable_order(paired_files, tmp_path):
    config = _config(
        paired_files,
        tmp_path,
        min_tumor_depth=10,
        min_normal_depth=5,
        circular_edge_bases=0,
    )
    row = _candidate(config, _histograms(1, 2), _histograms(1, 2), blacklist={1})

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
    tumor, normal = _histograms(5, 1000), _histograms(0, 1000)

    noisy = _candidate(config, tumor, normal, error_rates={"A>C": 0.01})
    clean = _candidate(config, tumor, normal, error_rates={"A>C": 1e-6})

    assert noisy["seq_p"] > clean["seq_p"]
    assert "WEAK_EVIDENCE" in noisy["filter"]
    assert "WEAK_EVIDENCE" not in clean["filter"]


def test_numt_filter_needs_autosomal_depth(paired_files, tmp_path):
    tumor, normal = _histograms(20, 1000), _histograms(0, 1000)

    without = _candidate(_config(paired_files, tmp_path), tumor, normal)
    with_depth = _candidate(
        _config(paired_files, tmp_path, autosomal_median_depth=30.0), tumor, normal
    )

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
    args = [*_paired_args(paired_files, tmp_path / "out"), "--dry-run"]

    result = CliRunner().invoke(cli, args)
    assert result.exit_code == 0
    assert not (tmp_path / "out").exists()

    result = CliRunner().invoke(cli, [*args, "-g", "chrM_absent"])
    assert result.exit_code != 0
    assert "is absent from" in result.output


def test_bam_and_cram_give_the_same_counts(paired_files, tmp_path):
    bam = run_paired_pipeline(_config(paired_files, tmp_path / "bam"))
    cram = run_paired_pipeline(
        _config(paired_files, tmp_path / "cram", tumor="tumor_cram", normal="normal_cram")
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


def test_outputs_are_valid_and_repeatable(paired_files, tmp_path):
    output = tmp_path / "out"
    config = _config(paired_files, output)
    result = run_paired_pipeline(config)

    assert sorted(path.name for path in output.iterdir()) == [
        "pair.mt_callable.bed.gz",
        "pair.mt_variants.vcf.gz",
        "pair.mt_variants.vcf.gz.tbi",
    ]
    with pysam.VariantFile(result.outputs["vcf"]) as vcf:
        assert list(vcf)

    qc = _qc(result.outputs["vcf"])
    assert qc["mgatk2_version"] == version("mgatk2")
    assert qc["reference"]["sha256"]
    assert qc["snv_only"] is True
    assert qc["circular_edge_bases"] == 2
    assert "shifted_reference_supplied" not in qc["parameters"]
    assert qc["counts"]["evidence_positions"] == 40
    assert qc["counts"]["callable_positions"] == result.callable_positions
    with gzip.open(result.outputs["callable_bed"], "rt") as bed:
        intervals = [line.split() for line in bed]
    assert len(intervals) == result.callable_positions
    assert all(start == str(int(end) - 1) for _chrom, start, end in intervals)

    first = {key: Path(path).read_bytes() for key, path in result.outputs.items()}
    rerun = run_paired_pipeline(config)
    assert {key: Path(path).read_bytes() for key, path in rerun.outputs.items()} == first


def test_rank_sum_filters_fire_only_on_degraded_alternates(paired_files, tmp_path):
    """A low-quality alternate is flagged; a high-quality one is not."""
    config = _config(paired_files, tmp_path, min_tumor_af=0.0, min_alt_observations=1)

    def with_alternate_at(alternate_bin):
        tumor = _histograms(20, 100)
        for source in (tumor.baseq, tumor.mapq, tumor.distance):
            source[0, 0, 30] = 80  # reference allele A
            source[0, 1, alternate_bin] = 20  # alternate allele C
        return _candidate(config, tumor, _histograms(0, 100), error_rates={"A>C": 1e-6})

    degraded = with_alternate_at(5)
    healthy = with_alternate_at(50)

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
