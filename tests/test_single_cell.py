import importlib
import logging
import subprocess
import sys

import h5py
import numpy as np
import pytest
from click.testing import CliRunner

from cli.utils import auto_detect_10x_structure, setup_file_logging
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


def test_read_end_threshold_reaches_the_pipeline():
    config = PipelineConfig(min_distance_from_end=0)

    assert config.quality.min_distance_from_end == 0


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


A, C, G, T = 0, 1, 2, 3
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
