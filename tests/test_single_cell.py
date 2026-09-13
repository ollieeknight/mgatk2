import importlib
import logging
import subprocess
import sys

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from cli.options import ASSAY_PRESETS, apply_assay_preset
from cli.utils import auto_detect_10x_structure, load_panel_positions, setup_file_logging
from core.config import PipelineConfig
from core.exceptions import InvalidInputError, NoBarcodeTagsError
from core.pipeline import run_pipeline
from file_io import IncrementalHDF5Writer
from processing.pileup import plan_shards, scan_shard
from processing.readers import BAMReader
from utils.utils import load_barcode_csv


@pytest.mark.parametrize("module", ["analysis", "cli", "core", "file_io", "processing", "utils"])
def test_packages_import_cleanly(module):
    result = subprocess.run(
        [sys.executable, "-c", f"import {module}"],
        capture_output=True,
        text=True,
    )

    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize("command", ["run", "tenx", "call", "paired", "hardmask-fasta"])
def test_every_command_help_renders(command):
    """Every registered command must at least parse its own option surface."""
    from cli import cli

    result = CliRunner().invoke(cli, [command, "--help"])

    assert result.exit_code == 0, result.output
    assert "--help" in result.output


def test_every_command_option_documents_itself():
    """An option with no help text is undiscoverable from the CLI."""
    import click

    from cli import cli

    context = click.Context(cli)
    undocumented = [
        f"{name} {parameter.opts}"
        for name, command in cli.commands.items()
        for parameter in command.get_params(context)
        if isinstance(parameter, click.Option) and not parameter.help
    ]

    assert undocumented == []


def test_read_end_threshold_reaches_the_pipeline():
    config = PipelineConfig(min_distance_from_end=0)

    assert config.quality.min_distance_from_end == 0


def test_umi_dedup_flag_reaches_the_pipeline():
    config = PipelineConfig(use_umi_dedup=True)

    assert config.dedup.use_umi is True
    assert config.dedup.umi_tag == "UB"


def test_umi_dedup_defaults_off():
    config = PipelineConfig()

    assert config.dedup.use_umi is False


def test_bulk_mode_skips_barcode_discovery(monkeypatch, tmp_path):
    captured = {}

    class FakePipeline:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def run(self):
            return {"cells_processed": 1}

    monkeypatch.setattr("core.pipeline.MtDNAPipeline", FakePipeline)
    result = run_pipeline(
        bam_path=str(tmp_path / "sample.bam"),
        barcode_file="bulk",
        output_dir=str(tmp_path / "output"),
    )

    assert result == {"cells_processed": 1}
    assert captured["barcodes"] == ["bulk"]


def test_barcode_less_bam_is_rejected(paired_files):
    with pytest.raises(NoBarcodeTagsError):
        BAMReader(
            str(paired_files["tumor_bam"]),
            PipelineConfig(mito_length=40),
            {"cell-1"},
        )


A, C, _, T = 0, 1, 2, 3
FWD, REV = 0, 1


def counting_config(**overrides):
    """Config that counts every base, so tests assert on the kernel alone."""
    defaults = {
        "mito_length": 40,
        "min_mapq": 0,
        "min_baseq": 0,
        "min_distance_from_end": 0,
        "max_strand_bias": 1.0,
    }
    return PipelineConfig(**{**defaults, **overrides})


def test_shard_counts_by_strand_and_deduplicates(barcoded_bam):
    result = scan_shard((str(barcoded_bam), counting_config(), ["cell-1", "cell-2"], 0, None))

    # r2 duplicates r1 exactly; the untagged barcode is never counted.
    assert result.duplicate_reads == 1
    assert result.n_reads.tolist() == [2, 1]

    # cell-1: r1 (ACGT x5 at 0) and r3 (ACGT x5 at 4) both put an A at position 4.
    assert result.counts[0, 0, A, FWD] == 1
    assert result.counts[0, 4, A, FWD] == 2
    assert result.counts[0, 1, C, FWD] == 1

    # cell-2 is a 20bp reverse-strand run of C.
    assert result.counts[1, :20, C, REV].tolist() == [1] * 20
    assert result.counts[1, :, :, FWD].sum() == 0

    # Tn5 cut sites: forward reads at their start, reverse reads at their end.
    assert result.tn5[0, 0, FWD] == 1
    assert result.tn5[0, 4, FWD] == 1
    assert result.tn5[1, 19, REV] == 1


def test_umi_dedup_collapses_reads_sharing_a_umi_at_different_starts(tmp_path, alignment_factory):
    """Same molecule, fragmented differently: positional dedup would miss this."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "umi_same.bam",
        reference,
        [
            {"name": "r1", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1", "UB": "AAAACCCCGGGG"}},
            {"name": "r2", "start": 3, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1", "UB": "AAAACCCCGGGG"}},
        ],
    )

    result = scan_shard((str(bam), counting_config(use_umi_dedup=True), ["cell-1"], 0, None))

    assert result.duplicate_reads == 1
    assert result.n_reads.tolist() == [1]


def test_umi_dedup_keeps_reads_sharing_a_start_with_different_umis(tmp_path, alignment_factory):
    """Two independent molecules that happen to start at the same base.

    `alignment_start` dedup would wrongly collapse these to one; this is the
    2.49x overcount measured against real scRNA data.
    """
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "umi_diff.bam",
        reference,
        [
            {"name": "r1", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1", "UB": "AAAACCCCGGGG"}},
            {"name": "r2", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1", "UB": "TTTTGGGGCCCC"}},
        ],
    )

    result = scan_shard((str(bam), counting_config(use_umi_dedup=True), ["cell-1"], 0, None))

    assert result.duplicate_reads == 0
    assert result.n_reads.tolist() == [2]


def test_umi_dedup_falls_back_to_position_when_tag_is_missing(tmp_path, alignment_factory):
    """A read without a UB tag must still be counted, not silently dropped."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "umi_missing.bam",
        reference,
        [{"name": "r1", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}}],
    )

    result = scan_shard((str(bam), counting_config(use_umi_dedup=True), ["cell-1"], 0, None))

    assert result.n_reads.tolist() == [1]
    assert result.duplicate_reads == 0


def test_insertion_keeps_query_and_reference_in_register(tmp_path, alignment_factory):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "insertion.bam",
        reference,
        [
            {
                "name": "ins",
                "start": 0,
                "sequence": "AAAAA" + "TTT" + "CCCCC",
                "cigar": ((0, 5), (1, 3), (0, 5)),
                "tags": {"CB": "cell-1"},
            }
        ],
    )

    result = scan_shard((str(bam), counting_config(), ["cell-1"], 0, None))

    assert result.counts[0, :5, A, FWD].tolist() == [1] * 5
    assert result.counts[0, 5:10, C, FWD].tolist() == [1] * 5
    assert result.counts[0, :, T, :].sum() == 0


def test_read_without_cigar_is_skipped(tmp_path, alignment_factory):
    """A placed read with no CIGAR has no aligned span, so it must not be counted."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "nocigar.bam",
        reference,
        [
            {"name": "aligned", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}},
            # Mapped, but CIGAR is absent: reference_end is None.
            {
                "name": "nocigar",
                "start": 10,
                "sequence": "A" * 20,
                "cigar": None,
                "flag": 16,
                "tags": {"CB": "cell-1"},
            },
        ],
    )

    result = scan_shard((str(bam), counting_config(), ["cell-1"], 0, None))

    # Only the genuinely aligned read contributes a read, bases, and a cut site.
    assert result.n_reads.tolist() == [1]
    assert result.counts[0].sum() == 20
    assert result.tn5[0, :, REV].sum() == 0


def test_min_distance_from_end_trims_both_read_ends(barcoded_bam):
    result = scan_shard(
        (str(barcoded_bam), counting_config(min_distance_from_end=2), ["cell-2"], 0, None)
    )

    # 20bp read, 2bp clipped at each end: only reference positions 2..17 survive.
    assert result.counts[0, :2, C, REV].sum() == 0
    assert result.counts[0, 2:18, C, REV].tolist() == [1] * 16
    assert result.counts[0, 18:, C, REV].sum() == 0


def test_strand_bias_filter_drops_single_stranded_bases(barcoded_bam):
    result = scan_shard(
        (str(barcoded_bam), counting_config(max_strand_bias=0.9), ["cell-2"], 0, None)
    )

    # Every observation is on one strand, so a 0.9 ceiling removes all of them.
    assert result.counts.sum() == 0
    assert result.depth.sum() == 0


def test_min_reads_per_cell_zeroes_failing_cells(barcoded_bam):
    result = scan_shard(
        (str(barcoded_bam), counting_config(min_reads_per_cell=2), ["cell-1", "cell-2"], 0, None)
    )

    assert result.kept.tolist() == [True, False]
    assert result.counts[1].sum() == 0
    assert result.mean_depth[1] == 0


def test_hdf5_output_matches_shard_counts(barcoded_bam, tmp_path):
    config = counting_config()
    barcodes = ["cell-1", "cell-2"]
    result = scan_shard((str(barcoded_bam), config, barcodes, 0, None))

    writer = IncrementalHDF5Writer(tmp_path, config, barcodes)
    writer.write_shard(result, barcodes)
    writer.finalize(tmp_path / "qc")

    with h5py.File(tmp_path / "output" / "counts.h5") as handle:
        assert [b.decode() for b in handle["barcode"][:]] == barcodes
        np.testing.assert_array_equal(handle["A_fwd"][:], result.counts[:, :, A, FWD].T)
        np.testing.assert_array_equal(handle["C_rev"][:], result.counts[:, :, C, REV].T)
    with h5py.File(tmp_path / "output" / "metadata.h5") as handle:
        np.testing.assert_array_equal(handle["coverage"][:], result.depth.T)
        assert handle["reference"][0] == b"A"


def test_both_reports_render_from_one_hdf5_run(barcoded_bam, tmp_path):
    """The two report flavours share a loader; neither may drift from the HDF5 layout."""
    import json

    from analysis.report import generate_html_report, generate_scrna_html_report

    config = counting_config()
    barcodes = ["cell-1", "cell-2"]
    writer = IncrementalHDF5Writer(tmp_path, config, barcodes)
    writer.write_shard(scan_shard((str(barcoded_bam), config, barcodes, 0, None)), barcodes)
    writer.finalize(tmp_path / "qc")
    (tmp_path / "qc" / "run_config.json").write_text(
        json.dumps({"mgatk_version": "test", "parameters": {"min_base_quality": 20}})
    )

    for generate in (generate_html_report, generate_scrna_html_report):
        report = generate(tmp_path, "sample")
        assert report.exists()
        page = report.read_text()
        assert "data:image/png;base64," in page
        assert "min_base_quality" in page


def test_plan_shards_caps_cells_by_memory_budget():
    config = PipelineConfig(n_cores=4, max_memory_gb=1.0)
    per_shard = plan_shards(100_000, config)

    assert per_shard >= 1
    assert per_shard * 4 * config.bytes_per_cell() <= 1.0e9
    # A small run still splits evenly across the cores rather than over-sharding.
    assert plan_shards(40, PipelineConfig(n_cores=4, max_memory_gb=128.0)) == 10


def test_file_log_captures_module_messages(tmp_path):
    log_path = tmp_path / "output.log"
    module_logger = logging.getLogger("core.pipeline")
    previous_level = module_logger.level
    module_logger.setLevel(logging.INFO)
    setup_file_logging(log_path)

    module_logger.info("pipeline message")
    for handler in logging.getLogger().handlers:
        handler.flush()
    module_logger.setLevel(previous_level)

    assert "pipeline message" in log_path.read_text()


def test_call_uses_bulk_mode(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.call")
    (tmp_path / "sample.bam").touch()
    calls = []
    monkeypatch.setattr(
        command_module,
        "run_pipeline_command",
        lambda **kwargs: calls.append(kwargs) or 0,
    )

    result = CliRunner().invoke(
        command_module.call,
        ["--input", str(tmp_path), "--output", str(tmp_path / "output")],
    )

    assert result.exit_code == 0, result.output
    assert calls[0]["barcode_file"] == "bulk"


def test_call_rejects_single_bam_file(caplog, tmp_path):
    command_module = importlib.import_module("cli.commands.call")
    bam = tmp_path / "one.bam"
    bam.touch()

    with caplog.at_level(logging.INFO):
        result = CliRunner().invoke(
            command_module.call,
            ["--input", str(bam), "--output", str(tmp_path / "output")],
        )

    assert result.exit_code == 1
    assert "mgatk2 run" in caplog.text


def test_auto_detect_finds_10x_multi_single_sample(tmp_path):
    count_dir = tmp_path / "outs" / "per_sample_outs" / "sampleA" / "count"
    count_dir.mkdir(parents=True)
    (count_dir / "sample_alignments.bam").touch()
    (count_dir / "sample_filtered_barcodes.csv").write_text("ref,AAAA-1\n")

    bam_path, barcode_file = auto_detect_10x_structure(str(tmp_path))

    assert bam_path.endswith("sample_alignments.bam")
    assert barcode_file.endswith("sample_filtered_barcodes.csv")


def test_auto_detect_rejects_multiple_10x_multi_samples(tmp_path):
    for sample in ("sampleA", "sampleB"):
        count_dir = tmp_path / "outs" / "per_sample_outs" / sample / "count"
        count_dir.mkdir(parents=True)
        (count_dir / "sample_alignments.bam").touch()

    with pytest.raises(InvalidInputError) as excinfo:
        auto_detect_10x_structure(str(tmp_path))

    assert "sampleA" in str(excinfo.value)
    assert "sampleB" in str(excinfo.value)


def test_load_barcode_csv_reads_atac_singlecell_schema(tmp_path):
    csv_file = tmp_path / "singlecell.csv"
    csv_file.write_text("barcode,is__cell_barcode\nAAAA-1,1\nCCCC-1,0\n")

    barcodes, metadata = load_barcode_csv(str(csv_file))

    assert barcodes == ["AAAA-1"]
    assert metadata is not None


def test_load_barcode_csv_reads_10x_multi_schema(tmp_path):
    csv_file = tmp_path / "sample_filtered_barcodes.csv"
    csv_file.write_text("GRCh38,AAAA-1\nGRCh38,CCCC-1\n")

    barcodes, metadata = load_barcode_csv(str(csv_file))

    assert barcodes == ["AAAA-1", "CCCC-1"]
    assert metadata is None


def test_load_barcode_csv_rejects_empty_file(tmp_path):
    csv_file = tmp_path / "empty.csv"
    csv_file.write_text("")

    with pytest.raises(InvalidInputError):
        load_barcode_csv(str(csv_file))


def test_run_returns_pipeline_failures(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.run")
    bam = tmp_path / "sample.bam"
    bam.touch()
    monkeypatch.setattr(command_module, "run_pipeline_command", lambda **kwargs: 1)

    result = CliRunner().invoke(command_module.run, ["--input", str(bam)])

    assert result.exit_code == 1


def test_low_mapq_reads_are_not_counted_as_retained(barcoded_bam):
    """Discarded reads must not inflate n_reads, or Tn5 totals stop matching."""
    result = scan_shard(
        (str(barcoded_bam), counting_config(min_mapq=61), ["cell-1", "cell-2"], 0, None)
    )

    assert result.n_reads.sum() == 0
    assert result.tn5.sum() == 0
    assert result.kept.sum() == 0


def test_tn5_cut_total_equals_retained_read_count(barcoded_bam):
    for min_mapq in (0, 61):
        result = scan_shard(
            (
                str(barcoded_bam),
                counting_config(min_mapq=min_mapq),
                ["cell-1", "cell-2"],
                0,
                None,
            )
        )

        assert result.tn5.sum() == result.n_reads.sum()


def test_strand_bias_uses_the_same_metric_as_paired(tmp_path, alignment_factory):
    """|forward - reverse| / total, matching analysis.paired_calling._strand_bias."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "bias.bam",
        reference,
        [
            {"name": "f", "start": 0, "sequence": "A" * 4, "tags": {"CB": "cell-1"}},
            {"name": "r", "start": 0, "sequence": "A" * 4, "flag": 16, "tags": {"CB": "cell-1"}},
        ],
    )

    def counts_at_zero(max_strand_bias):
        result = scan_shard(
            (str(bam), counting_config(max_strand_bias=max_strand_bias), ["cell-1"], 0, None)
        )
        return int(result.counts[0, 0, A, :].sum())

    # One forward and one reverse observation: |1 - 1| / 2 == 0, so even the
    # strictest ceiling keeps it. The old max/total metric never fell below 0.5
    # and would have zeroed this position for any ceiling under 0.5.
    assert counts_at_zero(0.0) == 2
    assert counts_at_zero(0.4) == 2


@pytest.mark.parametrize("column", ["is__cell_barcode", "is_cell_barcode", "is_cell"])
def test_load_barcode_csv_accepts_every_cell_flag_spelling(tmp_path, column):
    csv_file = tmp_path / "singlecell.csv"
    csv_file.write_text(f"barcode,{column}\nAAAA-1,1\nCCCC-1,0\n")

    barcodes, metadata = load_barcode_csv(str(csv_file))

    assert barcodes == ["AAAA-1"]
    assert metadata is not None


@pytest.mark.parametrize("command", ["run", "tenx", "call"])
def test_umi_deduplication_choice_is_accepted(command, monkeypatch, tmp_path):
    """--deduplication umi must parse on every single-cell command."""
    command_module = importlib.import_module(f"cli.commands.{command}")
    monkeypatch.setattr(command_module, "run_pipeline_command", lambda **kwargs: 0)
    if command == "call":
        (tmp_path / "sample.bam").touch()

    result = CliRunner().invoke(
        getattr(command_module, command),
        [
            "--input",
            str(tmp_path),
            "--output",
            str(tmp_path / "output"),
            "--deduplication",
            "umi",
        ],
    )

    assert result.exit_code == 0, result.output


def test_deduplication_umi_forwards_use_umi_dedup(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.run")
    (tmp_path / "sample.bam").touch()
    calls = []
    monkeypatch.setattr(
        command_module,
        "run_pipeline_command",
        lambda **kwargs: calls.append(kwargs) or 0,
    )

    result = CliRunner().invoke(
        command_module.run,
        [
            "--input",
            str(tmp_path / "sample.bam"),
            "--output",
            str(tmp_path / "output"),
            "--deduplication",
            "umi",
        ],
    )

    assert result.exit_code == 0, result.output
    assert calls[0]["dedup_mode"] == "umi"


def test_call_forwards_no_tn5(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.call")
    (tmp_path / "sample.bam").touch()
    calls = []
    monkeypatch.setattr(
        command_module,
        "run_pipeline_command",
        lambda **kwargs: calls.append(kwargs) or 0,
    )

    result = CliRunner().invoke(
        command_module.call,
        ["--input", str(tmp_path), "--output", str(tmp_path / "output"), "--no-tn5"],
    )

    assert result.exit_code == 0, result.output
    assert calls[0]["compute_tn5"] is False


# The presets `run`, `tenx`, and `call` are the whole reason three copies of the
# single-cell option surface existed. Pinning them here is what makes one
# shared builder safe.
EXPECTED_DEFAULTS = {
    "run": {
        "bam_path": ".",
        "output_format": "hdf5",
        "dedup_mode": "alignment_and_fragment_length",
        "base_qual": 20,
        "min_mapq": 30,
        "min_reads": 1,
        "min_distance_from_end": 5,
        "max_strand_bias": 1.0,
        "max_memory": 128.0,
        "barcode_tag": "CB",
        "min_barcode_reads": 10,
        "mito_genome": "chrM",
        "compute_tn5": True,
        "nh_max": 0,
        "nm_max": 0,
        "output_dir": "mgatk2",
    },
    "tenx": {
        "bam_path": ".",
        "output_format": "txt",
        "dedup_mode": "alignment_start",
        "base_qual": 0,
        "min_mapq": 0,
        "min_reads": 0,
        "min_distance_from_end": 0,
        "max_strand_bias": 1.0,
        "max_memory": 128.0,
        "barcode_tag": "CB",
        "min_barcode_reads": 10,
        "mito_genome": "chrM",
        "compute_tn5": True,
        "nh_max": 0,
        "nm_max": 0,
        "output_dir": "mgatk2",
    },
    "call": {
        "output_format": "hdf5",
        "dedup_mode": "alignment_and_fragment_length",
        "base_qual": 20,
        "min_mapq": 30,
        "min_distance_from_end": 5,
        "max_strand_bias": 1.0,
        "max_memory": 128.0,
        "mito_genome": "chrM",
        "compute_tn5": True,
        "nh_max": 0,
        "nm_max": 0,
        "output_dir": "mgatk2",
    },
    "paired": {
        "base_qual": 20,
        "min_mapq": 20,
        "min_distance_from_end": 5,
        "max_strand_bias": 0.9,
        "deduplication": "alignment_and_fragment_length",
        "min_tumor_depth": 10,
        "min_normal_depth": 5,
        "min_alt_observations": 3,
        "min_tumor_af": 0.005,
        "max_normal_af": 0.01,
        "circular_edge_bases": 500,
        "mito_genome": "chrM",
        "autosomal_median_depth": None,
        "custom_blacklist": None,
        "input_is_consensus": False,
    },
}


@pytest.mark.parametrize("command", sorted(EXPECTED_DEFAULTS))
def test_command_defaults_are_pinned(command):
    import click

    from cli import cli

    context = click.Context(cli)
    # get_default resolves Click 8.5's UNSET sentinel to the value the command
    # actually receives; parameter.default does not.
    actual = {
        parameter.name: parameter.get_default(context)
        for parameter in cli.commands[command].get_params(context)
    }

    for name, value in EXPECTED_DEFAULTS[command].items():
        assert actual[name] == value, name


def test_call_has_no_barcode_options():
    """A bulk BAM has no barcodes; offering the options would be a lie."""
    import click

    from cli import cli

    context = click.Context(cli)
    names = {p.name for p in cli.commands["call"].get_params(context)}

    assert not names & {"barcode_file", "barcode_tag", "min_barcode_reads", "min_reads"}


@pytest.mark.parametrize("command", ["run", "tenx", "call", "paired", "hardmask-fasta"])
def test_short_help_flag_is_accepted(command):
    from cli import cli

    assert CliRunner().invoke(cli, [command, "-h"]).exit_code == 0


@pytest.mark.parametrize("command", ["run", "tenx", "call"])
def test_multimapper_filters_are_exposed_on_every_singlecell_command(command):
    """The counting kernel supports NH/NM everywhere; only tenx used to say so."""
    import click

    from cli import cli

    context = click.Context(cli)
    names = {p.name for p in cli.commands[command].get_params(context)}

    assert {"nh_max", "nm_max"} <= names


def test_bulk_runs_skip_the_per_cell_html_report(monkeypatch, tmp_path, barcoded_bam):
    """One pseudo-cell cannot populate a per-cell QC report."""
    from core.pipeline import run_pipeline

    called = []
    import analysis.report as report_module

    monkeypatch.setattr(
        report_module, "generate_scrna_html_report", lambda *a, **_: called.append(a)
    )
    monkeypatch.setattr(report_module, "generate_html_report", lambda *a, **_: called.append(a))

    run_pipeline(
        bam_path=str(barcoded_bam),
        barcode_file="bulk",
        output_dir=str(tmp_path / "bulk"),
        mito_chr="chrM",
        output_format="hdf5",
        n_cores=1,
        min_mapq=0,
        min_baseq=0,
    )

    assert called == []
    assert not (tmp_path / "bulk" / "mgatk2_report.html").exists()


def test_call_processes_bam_files_concurrently(monkeypatch, tmp_path):
    """--threads on `call` used to be inert: one bulk cell never shards."""
    command_module = importlib.import_module("cli.commands.call")
    for name in ("a", "b", "c"):
        (tmp_path / f"{name}.bam").touch()

    submitted = []

    def fake_run(**kwargs):
        submitted.append(kwargs["bam_path"])
        return 0

    monkeypatch.setattr(command_module, "run_pipeline_command", fake_run)
    # Exercise the serial branch, which shares the task construction.
    result = CliRunner().invoke(
        command_module.call,
        ["--input", str(tmp_path), "--output", str(tmp_path / "out"), "--threads", "1"],
    )

    assert result.exit_code == 0, result.output
    assert len(submitted) == 3
    # Each sample is one bulk pseudo-cell, so its own pipeline never shards.
    assert all(name.endswith(".bam") for name in submitted)


def test_call_reports_a_failing_sample(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.call")
    (tmp_path / "a.bam").touch()
    monkeypatch.setattr(command_module, "run_pipeline_command", lambda **kwargs: 1)

    result = CliRunner().invoke(
        command_module.call,
        ["--input", str(tmp_path), "--output", str(tmp_path / "out"), "--threads", "1"],
    )

    assert result.exit_code == 1


def test_tapestri_preset_disables_coordinate_deduplication():
    """Amplicon reads share start coordinates, so coordinate dedup must be off."""
    assert ASSAY_PRESETS["tapestri"]["dedup_mode"] == "none"


def test_scrna_preset_uses_umi_deduplication():
    """Measured against real data: alignment_start overcounts scRNA depth 2.49x."""
    assert ASSAY_PRESETS["scrna"]["dedup_mode"] == "umi"


def test_assay_preset_fills_options_the_user_did_not_set():
    resolved = apply_assay_preset(
        "tapestri",
        {"dedup_mode": "alignment_start", "barcode_tag": "CB", "compute_tn5": True},
        explicit=set(),
    )

    assert resolved["dedup_mode"] == "none"
    assert resolved["barcode_tag"] == "RG"
    assert resolved["compute_tn5"] is False


def test_explicit_flag_overrides_assay_preset_with_a_warning(caplog):
    with caplog.at_level(logging.WARNING):
        resolved = apply_assay_preset(
            "tapestri",
            {"dedup_mode": "alignment_start"},
            explicit={"dedup_mode"},
        )

    assert resolved["dedup_mode"] == "alignment_start"
    assert "tapestri" in caplog.text


def test_no_assay_leaves_every_option_untouched():
    values = {"dedup_mode": "alignment_start", "barcode_tag": "CB"}

    assert apply_assay_preset(None, values, explicit=set()) == values


def test_panel_bed_scopes_coverage_breadth_to_targeted_bases(tmp_path, alignment_factory):
    """A panel only targets part of chrM, so untargeted bases are not misses."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "panel.bam",
        reference,
        [{"name": "amp", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}}],
    )
    # The read covers positions 1-20; the panel targets exactly those.
    panel = tmp_path / "panel.bed"
    panel.write_text("chrM\t0\t20\n")

    config = counting_config(panel_positions=load_panel_positions(str(panel), "chrM"))
    result = scan_shard((str(bam), config, ["cell-1"], 0, None))

    # 20 of 20 targeted bases covered, rather than 20 of the full 40.
    assert result.coverage_breadth[0] == pytest.approx(1.0)


def test_coverage_breadth_uses_whole_contig_without_a_panel(barcoded_bam):
    result = scan_shard((str(barcoded_bam), counting_config(), ["cell-2"], 0, None))

    # cell-2 covers 20 bases of the 40bp test contig.
    assert result.coverage_breadth[0] == pytest.approx(0.5)


def test_panel_positions_are_one_based_inclusive(tmp_path):
    bed = tmp_path / "p.bed"
    bed.write_text("chrM\t0\t3\nchrOther\t0\t99\n")

    # BED is 0-based half-open [0,3); only chrM rows count.
    assert load_panel_positions(str(bed), "chrM") == frozenset({1, 2, 3})


@pytest.fixture
def amplicon_bam(tmp_path, alignment_factory):
    """One Tapestri-style amplicon: 30 molecules sharing a start coordinate."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    reads = [
        {
            "name": f"mol{i}",
            "start": 5,
            "sequence": "ACGT" * 5,
            "template_length": 20,
            "tags": {"RG": "cell-1"},
        }
        for i in range(30)
    ]
    return alignment_factory(tmp_path / "amplicon.bam", reference, reads)


def test_tapestri_preset_keeps_every_amplicon_molecule(amplicon_bam):
    """Coordinate dedup would collapse this amplicon to a single read."""
    resolved = apply_assay_preset(
        "tapestri", {"dedup_mode": "alignment_and_fragment_length"}, explicit=set()
    )
    config = counting_config(
        barcode_tag="RG",
        skip_deduplication=resolved["dedup_mode"] == "none",
    )

    result = scan_shard((str(amplicon_bam), config, ["cell-1"], 0, None))

    assert int(result.n_reads[0]) == 30
    assert result.duplicate_reads == 0


def test_coordinate_dedup_would_destroy_amplicon_data(amplicon_bam):
    """Guards the reason the tapestri preset exists."""
    config = counting_config(barcode_tag="RG", skip_deduplication=False)

    result = scan_shard((str(amplicon_bam), config, ["cell-1"], 0, None))

    assert int(result.n_reads[0]) == 1
    assert result.duplicate_reads == 29


def test_barcodes_are_discovered_from_the_read_group_tag(amplicon_bam):
    from file_io.barcode_extraction import extract_barcodes_from_bam

    barcodes = extract_barcodes_from_bam(
        str(amplicon_bam), barcode_tag="RG", mito_chr="chrM", min_reads=1
    )

    assert barcodes == ["cell-1"]


def test_run_exposes_assay_and_panel_bed():
    from cli import cli

    result = CliRunner().invoke(cli, ["run", "--help"])

    assert "--assay" in result.output
    assert "--panel-bed" in result.output
    assert "tapestri" in result.output


@pytest.mark.parametrize("command", ["tenx", "call"])
def test_assay_is_scoped_to_run(command):
    """tenx keeps mgatk-compatible defaults; call is a bulk per-BAM path."""
    from cli import cli

    result = CliRunner().invoke(cli, [command, "--help"])

    assert "--assay" not in result.output


def test_run_passes_assay_resolved_options_to_the_pipeline(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.run")
    captured = {}
    monkeypatch.setattr(
        command_module, "run_pipeline_command", lambda **kwargs: captured.update(kwargs) or 0
    )
    bam = tmp_path / "in.bam"
    bam.touch()

    result = CliRunner().invoke(
        command_module.run,
        ["--input", str(bam), "--output", str(tmp_path / "out"), "--assay", "tapestri"],
    )

    assert result.exit_code == 0, result.output
    assert captured["dedup_mode"] == "none"
    assert captured["barcode_tag"] == "RG"
    assert captured["compute_tn5"] is False
    assert captured["min_distance_from_end"] == 0


def test_assay_run_reaches_the_pipeline_without_mocking_the_command(tmp_path, amplicon_bam):
    """Guards the wiring a mocked run_pipeline_command cannot see."""
    from cli import cli

    result = CliRunner().invoke(
        cli,
        [
            "run",
            "--input",
            str(amplicon_bam),
            "--output",
            str(tmp_path / "out"),
            "--assay",
            "tapestri",
            "--threads",
            "1",
            "--min-barcode-reads",
            "1",
        ],
    )

    assert result.exit_code == 0, result.output


def test_assay_selects_the_report_instead_of_inferring_it():
    """A tapestri run has no transposition to plot, even beside a singlecell.csv."""
    from core.pipeline import MtDNAPipeline

    # Metadata present would otherwise select the Tn5 report.
    assert MtDNAPipeline.wants_tn5_report(barcode_metadata={"a": 1}, assay="tapestri") is False
    assert MtDNAPipeline.wants_tn5_report(barcode_metadata={"a": 1}, assay="scrna") is False
    assert MtDNAPipeline.wants_tn5_report(barcode_metadata={"a": 1}, assay="scatac") is True


def test_report_selection_falls_back_to_inference_without_an_assay():
    """Unchanged behaviour when --assay is not given."""
    from core.pipeline import MtDNAPipeline

    assert MtDNAPipeline.wants_tn5_report(barcode_metadata={"a": 1}, assay=None) is True
    assert MtDNAPipeline.wants_tn5_report(barcode_metadata=None, assay=None) is False


def test_run_config_records_the_assay_and_panel():
    """A run must be reproducible from its recorded parameters."""
    from analysis.qc import QCCalculator

    config = PipelineConfig(panel_positions=frozenset({1, 2, 3}))
    metadata = QCCalculator(config, assay="tapestri").collect_run_metadata("in.bam", "out", 1, 1)

    assert metadata["parameters"]["assay"] == "tapestri"
    assert metadata["parameters"]["panel_positions"] == 3


def test_conflict_warning_names_the_real_flag(caplog):
    """A warning naming a flag that does not exist is worse than none."""
    with caplog.at_level(logging.WARNING):
        apply_assay_preset("tapestri", {"dedup_mode": "alignment_start"}, explicit={"dedup_mode"})

    assert "--deduplication" in caplog.text
    assert "--dedup-mode" not in caplog.text
