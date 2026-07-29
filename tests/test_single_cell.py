import importlib
import logging
import subprocess
import sys

import numpy as np
import pytest
from click.testing import CliRunner

from cli.utils import setup_file_logging
from core.config import PipelineConfig, SimpleRead
from core.exceptions import NoBarcodeTagsError
from core.pipeline import run_pipeline
from processing.processors import process_barcode_worker
from processing.readers import BAMReader


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
            str(paired_files["query_bam"]),
            PipelineConfig(mito_length=40),
            {"cell-1"},
        )


def test_paired_and_orphan_reads_count_as_fragments(monkeypatch):
    class FakePileupGenerator:
        def __init__(self, config):
            pass

        def generate_pileup(self, reads):
            return {0: {"depth": len(reads)}}

        def filter_strand_bias(self, pileup):
            return pileup

    monkeypatch.setattr("processing.processors.PileupGenerator", FakePileupGenerator)
    reads = [
        SimpleRead(0, False, 60, b"A", np.array([40]), [(0, 1)], is_paired=True),
        SimpleRead(0, True, 60, b"A", np.array([40]), [(0, 1)], is_paired=True),
        SimpleRead(0, False, 60, b"A", np.array([40]), [(0, 1)]),
    ]

    result = process_barcode_worker(("cell", reads, PipelineConfig(mito_length=1)))

    assert result["qc"]["total_fragments"] == 2


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


def test_run_returns_pipeline_failures(monkeypatch, tmp_path):
    command_module = importlib.import_module("cli.commands.run")
    bam = tmp_path / "sample.bam"
    bam.touch()
    monkeypatch.setattr(command_module, "run_pipeline_command", lambda **kwargs: 1)

    result = CliRunner().invoke(command_module.run, ["--input", str(bam)])

    assert result.exit_code == 1
