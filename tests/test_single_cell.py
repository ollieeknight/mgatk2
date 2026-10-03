import importlib
import logging
from pathlib import Path

import click
import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from cli import cli
from cli.options import apply_assay_preset
from cli.utils import auto_detect_10x_structure, load_panel_positions
from core.config import PipelineConfig
from core.exceptions import InvalidInputError
from core.pipeline import run_pipeline, wants_tn5_report
from file_io.writers import IncrementalHDF5Writer
from processing.pileup import plan_shards, scan_shard
from processing.readers import BAMReader
from utils.utils import load_barcode_csv


def _options(command):
    """A command's options, resolved in its own context.

    A bare Context(cli) would make Click cache the help option without -h.
    """
    context = click.Context(command, **command.context_settings)
    return context, command.get_params(context)


def test_every_command_option_documents_itself():
    """An option with no help text is undiscoverable from the CLI."""
    undocumented = [
        f"{name} {parameter.opts}"
        for name, command in cli.commands.items()
        for parameter in _options(command)[1]
        if isinstance(parameter, click.Option) and not parameter.help
    ]

    assert undocumented == []


def test_barcode_less_bam_is_rejected(paired_files):
    with pytest.raises(InvalidInputError, match="barcode tag"):
        BAMReader(
            str(paired_files["tumor_bam"]),
            PipelineConfig(mito_length=40),
            check_barcode_tag=True,
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

    from analysis.report import generate_html_report

    config = counting_config()
    barcodes = ["cell-1", "cell-2"]
    writer = IncrementalHDF5Writer(tmp_path, config, barcodes)
    writer.write_shard(scan_shard((str(barcoded_bam), config, barcodes, 0, None)), barcodes)
    writer.finalize(tmp_path / "qc")
    (tmp_path / "qc" / "run_config.json").write_text(
        json.dumps({"mgatk_version": "test", "parameters": {"min_base_quality": 20}})
    )

    for tn5 in (True, False):
        report = generate_html_report(tmp_path, "sample", tn5=tn5)
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


def test_load_barcode_csv_reads_10x_multi_schema(tmp_path):
    csv_file = tmp_path / "sample_filtered_barcodes.csv"
    csv_file.write_text("GRCh38,AAAA-1\nGRCh38,CCCC-1\n")

    barcodes, metadata = load_barcode_csv(str(csv_file))

    assert barcodes == ["AAAA-1", "CCCC-1"]
    assert metadata is None


def test_tn5_cut_total_equals_retained_read_count(barcoded_bam):
    """Reads dropped by a filter must not count as retained, or the totals diverge."""
    for min_mapq, retained in ((0, 3), (61, 0)):
        config = counting_config(min_mapq=min_mapq)
        result = scan_shard((str(barcoded_bam), config, ["cell-1", "cell-2"], 0, None))

        assert result.n_reads.sum() == retained
        assert result.tn5.sum() == retained


def test_strand_bias_is_forward_minus_reverse_over_total(tmp_path, alignment_factory):
    """Same metric as paired: |forward - reverse| / total."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    bam = alignment_factory(
        tmp_path / "bias.bam",
        reference,
        [
            {"name": "f", "start": 0, "sequence": "A" * 4, "tags": {"CB": "cell-1"}},
            {"name": "r", "start": 0, "sequence": "A" * 4, "flag": 16, "tags": {"CB": "cell-1"}},
            {"name": "f-only", "start": 10, "sequence": "A" * 4, "tags": {"CB": "cell-1"}},
        ],
    )

    config = counting_config(max_strand_bias=0.0)
    counts = scan_shard((str(bam), config, ["cell-1"], 0, None)).counts[0, :, A, :].sum(axis=1)

    # Balanced position 0 has bias 0 and survives even a zero ceiling; the
    # single-stranded position 10 has bias 1 and is removed.
    assert counts[0] == 2
    assert counts[10] == 0


@pytest.mark.parametrize("column", ["is__cell_barcode", "is_cell_barcode", "is_cell"])
def test_load_barcode_csv_accepts_every_cell_flag_spelling(tmp_path, column):
    csv_file = tmp_path / "singlecell.csv"
    csv_file.write_text(f"barcode,{column}\nAAAA-1,1\nCCCC-1,0\n")

    barcodes, metadata = load_barcode_csv(str(csv_file))

    assert barcodes == ["AAAA-1"]
    assert metadata is not None


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
    context, parameters = _options(cli.commands[command])
    # get_default resolves Click 8.5's UNSET sentinel to the value the command
    # actually receives; parameter.default does not.
    actual = {parameter.name: parameter.get_default(context) for parameter in parameters}

    for name, value in EXPECTED_DEFAULTS[command].items():
        assert actual[name] == value, name


@pytest.mark.parametrize("command", ["run", "tenx", "call", "paired", "hardmask-fasta"])
def test_short_help_flag_is_accepted(command):
    assert CliRunner().invoke(cli, [command, "-h"]).exit_code == 0


def test_bulk_runs_skip_the_per_cell_html_report(tmp_path, barcoded_bam):
    """One pseudo-cell cannot populate a per-cell QC report."""
    output = tmp_path / "bulk"
    config = PipelineConfig(n_cores=1, min_mapq=0, min_baseq=0)
    run_pipeline(str(barcoded_bam), str(output), config, barcode_file="bulk")

    assert (output / "output" / "counts.h5").exists()
    assert not (output / "mgatk2_report.html").exists()


@pytest.mark.parametrize("status", [0, 1])
def test_call_runs_each_bam_as_a_bulk_sample(monkeypatch, tmp_path, status):
    command_module = importlib.import_module("cli.commands.call")
    for name in ("a", "b"):
        (tmp_path / f"{name}.bam").touch()
    calls = []
    monkeypatch.setattr(
        command_module, "run_pipeline_command", lambda **kwargs: calls.append(kwargs) or status
    )

    result = CliRunner().invoke(
        command_module.call,
        ["--input", str(tmp_path), "--output", str(tmp_path / "out"), "--threads", "1", "--no-tn5"],
    )

    assert result.exit_code == status, result.output
    assert [Path(call["output_dir"]).name for call in calls] == ["a", "b"]
    assert all(call["barcode_file"] == "bulk" and call["ncores"] == 1 for call in calls)
    assert all(call["compute_tn5"] is False for call in calls)


def test_assay_preset_fills_options_the_user_did_not_set():
    resolved = apply_assay_preset(
        "tapestri",
        {"dedup_mode": "alignment_start", "barcode_tag": "CB", "compute_tn5": True},
        explicit=set(),
    )

    assert resolved["dedup_mode"] == "none"
    assert resolved["barcode_tag"] == "RG"
    assert resolved["compute_tn5"] is False


def test_panel_bed_scopes_coverage_breadth_to_targeted_bases(tmp_path, barcoded_bam):
    """A panel only targets part of chrM, so untargeted bases are not misses."""
    panel = tmp_path / "panel.bed"
    panel.write_text("chrM\t0\t20\n")  # cell-2's read covers exactly positions 1-20

    def breadth(config):
        return scan_shard((str(barcoded_bam), config, ["cell-2"], 0, None)).coverage_breadth[0]

    assert breadth(counting_config()) == pytest.approx(0.5)
    panel_positions = load_panel_positions(str(panel), "chrM")
    assert breadth(counting_config(panel_positions=panel_positions)) == pytest.approx(1.0)


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


def test_coordinate_dedup_would_destroy_amplicon_data(amplicon_bam):
    """Why the tapestri preset turns deduplication off."""
    config = counting_config(barcode_tag="RG", skip_deduplication=False)

    result = scan_shard((str(amplicon_bam), config, ["cell-1"], 0, None))

    assert int(result.n_reads[0]) == 1
    assert result.duplicate_reads == 29


def test_tapestri_assay_keeps_every_amplicon_molecule(tmp_path, amplicon_bam):
    """--assay tapestri reaches the kernel: RG barcodes, no deduplication."""
    output = tmp_path / "out"
    result = CliRunner().invoke(
        cli,
        ["run", "-i", str(amplicon_bam), "-o", str(output), "--assay", "tapestri", "-t", "1"]
        + ["--min-barcode-reads", "1"],
    )

    assert result.exit_code == 0, result.output
    rows = (output / "qc" / "cell_stats.csv").read_text().splitlines()
    assert rows[1].split(",")[0] == "cell-1"
    assert rows[1].split(",")[-1] == "30"


def test_assay_selects_the_report_before_metadata_inference():
    # singlecell.csv metadata implies ATAC only when no assay is declared.
    assert wants_tn5_report(barcode_metadata={"a": 1}, assay="tapestri") is False
    assert wants_tn5_report(barcode_metadata={"a": 1}, assay="scrna") is False
    assert wants_tn5_report(barcode_metadata=None, assay="scatac") is True
    assert wants_tn5_report(barcode_metadata={"a": 1}, assay=None) is True
    assert wants_tn5_report(barcode_metadata=None, assay=None) is False


def test_explicit_flag_beats_the_preset_and_warns_with_the_real_flag(caplog):
    with caplog.at_level(logging.WARNING):
        resolved = apply_assay_preset(
            "tapestri", {"dedup_mode": "alignment_start"}, explicit={"dedup_mode"}
        )

    assert resolved["dedup_mode"] == "alignment_start"
    assert "--deduplication" in caplog.text


def test_mito_alias_resolves_before_barcode_discovery(tmp_path, barcoded_bam):
    """-g MT on a chrM BAM must resolve the contig before anything fetches it."""
    run_pipeline(
        str(barcoded_bam),
        str(tmp_path / "out"),
        PipelineConfig(mito_chr="MT", mito_length=40, n_cores=1, min_mapq=0, min_baseq=0),
        min_barcode_reads=1,
    )

    assert (tmp_path / "out" / "output" / "counts.h5").exists()


def test_dry_run_creates_no_files(tmp_path, barcoded_bam):
    output = tmp_path / "out"
    result = CliRunner().invoke(
        cli, ["run", "--input", str(barcoded_bam), "--output", str(output), "--dry-run"]
    )

    assert result.exit_code == 0, result.output
    assert not output.exists()
