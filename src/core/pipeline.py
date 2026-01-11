"""Main pipeline orchestration for mtDNA genotyping."""

import gc
import logging
import time
from pathlib import Path
from typing import Any

import pysam

from analysis.qc import QCCalculator
from core.config import (
    DeduplicationConfig,
    PerformanceConfig,
    PipelineConfig,
    QualityThresholds,
)
from core.exceptions import InvalidInputError
from file_io import IncrementalHDF5Writer, IncrementalTextWriter
from processing.processors import CellProcessor
from processing.readers import BAMReader

logger = logging.getLogger(__name__)


class MtDNAPipeline:
    """Single-pass mtDNA genotyping pipeline"""

    def __init__(
        self,
        bam_path: str,
        barcodes: list[str],
        output_dir: Path,
        config: PipelineConfig | None = None,
        output_format: str = "standard",
        barcode_metadata=None,
        sample_name: str = "mgatk2",
        report_title: str | None = None,
        report_subtitle: str | None = None,
        working_directory: str | None = None,
    ):
        self.bam_path = Path(bam_path)
        self.barcodes = set(barcodes)
        self.barcode_list = barcodes
        self.output_dir = Path(output_dir)
        self.config = config or PipelineConfig()
        self.output_format = output_format.lower()
        self.barcode_metadata = barcode_metadata
        self.sample_name = sample_name
        self.report_title = report_title or sample_name
        self.report_subtitle = report_subtitle or "mgatk2 output analysis"
        self.working_directory = working_directory

        if not self.bam_path.exists():
            raise InvalidInputError(f"BAM file not found: {bam_path}")

        index_path = Path(str(self.bam_path) + ".bai")
        if not index_path.exists():
            logger.warning("BAM index not found, creating: %s", index_path)
            pysam.index(str(self.bam_path))

        with pysam.AlignmentFile(str(self.bam_path), "rb") as bam:
            if self.config.mito_chr not in bam.references:
                # Try common mitochondrial chromosome names
                for alt_name in ["chrM", "MT", "M", "chrMT"]:
                    if alt_name in bam.references:
                        logger.warning(f"Using '{alt_name}' instead of '{self.config.mito_chr}'")
                        self.config.mito_chr = alt_name
                        break
                else:
                    available = ", ".join(bam.references[:10])
                    raise InvalidInputError(
                        f"Mitochondrial chromosome '{self.config.mito_chr}' not found. "
                        f"Available: {available}"
                    )

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self) -> dict[str, Any]:
        start_time = time.time()

        logger.info("")
        logger.info("Collecting reads from BAM by barcode...")
        reader = BAMReader(str(self.bam_path), self.config, self.barcodes)
        reads_by_barcode, stats = reader.collect_reads_by_barcode()

        if not reads_by_barcode:
            logger.error("No reads found for any barcodes!")
            return {}

        n_cells_input = len(reads_by_barcode)

        incremental_writer = None
        if self.output_format == "hdf5":
            incremental_writer = IncrementalHDF5Writer(
                self.output_dir,
                self.config,
                self.barcode_list,
                barcode_metadata=self.barcode_metadata,
            )
        else:
            incremental_writer = IncrementalTextWriter(
                self.output_dir, self.config, self.barcode_list
            )

        processor = CellProcessor(self.config, self.output_dir)
        cell_results = processor.process_cells_progressive(reads_by_barcode, incremental_writer)

        if not cell_results:
            logger.error("No cells passed quality filters")
            return {}

        logger.info("Cleaning up...")
        qc_calc = QCCalculator(self.config)
        qc_dir = self.output_dir / "qc"

        incremental_writer.finalize(qc_dir)

        run_metadata = qc_calc.collect_run_metadata(
            str(self.bam_path),
            str(self.output_dir),
            n_cells_input,
            len(cell_results),
        )
        from file_io import write_run_summary

        write_run_summary(run_metadata, qc_dir / "summary.txt")

        gc.collect()

        if self.barcode_metadata is not None and self.output_format == "hdf5":
            logger.info("Generating HTML QC report...")
            try:
                from analysis.report import generate_html_report

                generate_html_report(
                    self.output_dir,
                    self.sample_name,
                    title=self.report_title,
                    subtitle=self.report_subtitle,
                    working_directory=self.working_directory,
                    input_dir=str(self.bam_path.parent),
                )
            except ImportError:
                logger.warning("matplotlib not installed, skipping HTML report generation")
            except Exception as e:
                logger.warning("Failed to generate HTML report: %s", e)
        else:
            logger.info("Skipping HTML report generation (requires singlecell.csv metadata")

        elapsed = time.time() - start_time
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)
        seconds = int(elapsed % 60)

        if hours > 0:
            logger.info(f"Pipeline complete (Elapsed time: {hours}h {minutes}m {seconds}s)")
        elif minutes > 0:
            logger.info(f"Pipeline complete (Elapsed time: {minutes}m {seconds}s)")
        else:
            logger.info(f"Pipeline complete (Elapsed time: {seconds}s)")

        return {
            "cells_processed": n_cells_input,
            "cells_passed_qc": len(cell_results),
            "mean_reads": (
                sum(r["n_reads"] for r in cell_results) / len(cell_results) if cell_results else 0
            ),
        }


def run_pipeline(
    bam_path: str,
    barcode_file: str | None = None,
    output_dir: str = "",
    sample_name: str = "mgatk",
    min_baseq: int = 20,
    min_mapq: int = 30,
    min_reads_per_cell: int = 10,
    max_strand_bias: float = 1.0,
    min_distance_from_end: int = 5,
    skip_deduplication: bool = False,
    use_fragment_length_dedup: bool = True,
    write_cell_bams: bool = False,
    barcode_tag: str = "CB",
    mito_chr: str = "chrM",
    n_cores: int = 8,
    batch_size: int | None = None,
    max_memory_gb: float = 128.0,
    output_format: str = "standard",
    sequential: bool = False,
    report_title: str | None = None,
    report_subtitle: str | None = None,
    working_directory: str | None = None,
) -> dict[str, Any]:
    """Run the pipeline with individual parameters"""
    barcode_metadata = None

    if barcode_file is None:
        # Bulk calling mode - process all reads without barcode filtering
        logger.info("Running in bulk calling mode (no barcode filtering)")
        barcodes = ["bulk"]  # Single pseudo-barcode for bulk analysis
    elif barcode_file.endswith(".csv"):
        from utils.utils import load_singlecell_csv

        barcodes, barcode_metadata = load_singlecell_csv(barcode_file)
    else:
        with open(barcode_file) as f:
            barcodes = [line.strip() for line in f if line.strip()]

    if batch_size is None:
        batch_size = n_cores

    config = PipelineConfig(
        quality=QualityThresholds(
            min_baseq=min_baseq,
            min_mapq=min_mapq,
            max_strand_bias=max_strand_bias,
            min_distance_from_end=min_distance_from_end,
        ),
        dedup=DeduplicationConfig(
            skip=skip_deduplication, use_fragment_length=use_fragment_length_dedup
        ),
        performance=PerformanceConfig(
            n_cores=n_cores,
            batch_size=batch_size,
            max_memory_gb=max_memory_gb,
            sequential=sequential,
        ),
        min_reads_per_cell=min_reads_per_cell,
        barcode_tag=barcode_tag,
        mito_chr=mito_chr,
        write_cell_bams=write_cell_bams,
    )

    pipeline = MtDNAPipeline(
        bam_path=bam_path,
        barcodes=barcodes,
        output_dir=Path(output_dir),
        config=config,
        output_format=output_format,
        barcode_metadata=barcode_metadata,
        sample_name=sample_name,
        report_title=report_title,
        report_subtitle=report_subtitle,
        working_directory=working_directory,
    )

    return pipeline.run()
