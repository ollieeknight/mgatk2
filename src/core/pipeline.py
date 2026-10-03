"""Single-cell orchestration: load barcodes, count shards, write outputs."""

import datetime
import gzip
import logging
import time
from importlib.metadata import version
from pathlib import Path

from analysis.report import generate_html_report
from core.config import PipelineConfig
from core.exceptions import InvalidInputError, ProcessingError
from file_io.barcode_extraction import extract_barcodes_from_bam
from file_io.formats import write_run_config, write_run_summary
from file_io.writers import IncrementalHDF5Writer, IncrementalTextWriter
from processing.processors import process_shards
from processing.readers import BAMReader
from utils.utils import load_barcode_csv

logger = logging.getLogger(__name__)


def load_barcodes(
    barcode_file: str | None, bam_path: str, config: PipelineConfig, min_barcode_reads: int
) -> tuple[list[str], dict | None]:
    """Barcodes to count, plus per-barcode metadata when a singlecell.csv supplies it."""
    if barcode_file == "bulk":
        return ["bulk"], None

    metadata = None
    if barcode_file is None:
        barcodes = extract_barcodes_from_bam(
            bam_path, config.barcode_tag, config.mito_chr, min_barcode_reads
        )
        if not barcodes:
            raise InvalidInputError(
                f"No barcodes found in BAM file with tag '{config.barcode_tag}' "
                f"and minimum {min_barcode_reads} reads"
            )
    elif barcode_file.endswith(".csv"):
        barcodes, metadata = load_barcode_csv(barcode_file)
    else:
        opener = gzip.open if barcode_file.endswith(".gz") else open
        with opener(barcode_file, "rt") as f:
            barcodes = [line.strip() for line in f if line.strip()]

    # A repeated barcode would map every read to its last column and leave the rest empty.
    barcodes = list(dict.fromkeys(barcodes))
    logger.info("Loaded %s barcodes", f"{len(barcodes):,}")
    return barcodes, metadata


def wants_tn5_report(barcode_metadata, assay) -> bool:
    """Whether the QC report should centre on Tn5 transposition.

    A declared assay decides directly. Without one, fall back to the
    historical inference: barcode metadata only comes from a 10x scATAC
    singlecell.csv, so its presence implies ATAC.
    """
    if assay is not None:
        return assay == "scatac"
    return barcode_metadata is not None


def run_metadata(
    config: PipelineConfig,
    assay: str | None,
    bam_path: str,
    output_dir: str,
    n_cells_input: int,
    n_cells_passed: int,
) -> dict:
    """The run record written to qc/run_config.json and qc/summary.txt."""
    return {
        "mgatk_version": version("mgatk2"),
        "run_date": datetime.datetime.now().isoformat(),
        "input_bam": str(bam_path),
        "output_dir": str(output_dir),
        "reference": config.mito_chr,
        "reference_length": config.mito_length,
        "cells_total": n_cells_input,
        "cells_passed_qc": n_cells_passed,
        "cells_failed_qc": n_cells_input - n_cells_passed,
        "parameters": {
            "min_base_quality": config.min_baseq,
            "min_mapping_quality": config.min_mapq,
            "min_reads_per_cell": config.min_reads_per_cell,
            "max_strand_bias": config.max_strand_bias,
            "min_distance_from_end": config.min_distance_from_end,
            "nh_max": config.nh_max,
            "nm_max": config.nm_max,
            "compute_tn5": config.compute_tn5,
            "skip_deduplication": config.skip_deduplication,
            "use_fragment_length_dedup": config.use_fragment_length_dedup,
            "barcode_tag": config.barcode_tag,
            "assay": assay,
            # Size only: the full position set is recoverable from the BED.
            "panel_positions": len(config.panel_positions) if config.panel_positions else 0,
            "mito_chr": config.mito_chr,
            "n_cores": config.n_cores,
            "max_memory_gb": config.max_memory_gb,
        },
    }


def run_pipeline(
    bam_path: str,
    output_dir: str,
    config: PipelineConfig,
    barcode_file: str | None = None,
    min_barcode_reads: int = 10,
    output_format: str = "hdf5",
    assay: str | None = None,
    title: str = "mgatk2",
) -> None:
    """Count every barcode in one BAM and write the matrices, QC tables, and report."""
    started = time.monotonic()
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolves the mitochondrial contig name and checks the barcode tag before
    # anything else reads the BAM.
    BAMReader(bam_path, config, check_barcode_tag=barcode_file != "bulk")
    barcodes, barcode_metadata = load_barcodes(barcode_file, bam_path, config, min_barcode_reads)

    if output_format == "hdf5":
        writer = IncrementalHDF5Writer(output_dir, config, barcodes, barcode_metadata)
    else:
        writer = IncrementalTextWriter(output_dir, config, barcodes)

    totals = process_shards(bam_path, config, barcodes, writer)
    cells_passed = totals["cells_passed"]
    if not cells_passed:
        raise ProcessingError("No cells passed quality filters")

    if not config.skip_deduplication and totals["total_reads"]:
        logger.info(
            "%s duplicate reads removed (%.1f%%)",
            f"{totals['duplicate_reads']:,}",
            totals["duplicate_reads"] / totals["total_reads"] * 100,
        )
    logger.info(
        "Kept %s reads from %s cells at an average of %.0f reads/cell",
        f"{totals['kept_reads']:,}",
        f"{cells_passed:,}",
        totals["kept_reads"] / cells_passed,
    )

    qc_dir = output_dir / "qc"
    writer.finalize(qc_dir)
    metadata = run_metadata(config, assay, bam_path, output_dir, len(barcodes), cells_passed)
    write_run_summary(metadata, qc_dir / "summary.txt")
    write_run_config(metadata, qc_dir / "run_config.json")

    # The report is per-cell QC; a bulk sample is one pseudo-cell.
    if output_format == "hdf5" and barcodes != ["bulk"]:
        logger.info("Generating HTML QC report...")
        try:
            generate_html_report(output_dir, title, tn5=wants_tn5_report(barcode_metadata, assay))
        except Exception as e:
            logger.warning("Failed to generate HTML report: %s", e)

    elapsed = datetime.timedelta(seconds=round(time.monotonic() - started))
    logger.info("Pipeline complete in %s", elapsed)
