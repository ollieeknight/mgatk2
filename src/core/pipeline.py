"""Main pipeline orchestration for mtDNA genotyping."""

import logging
import time
from pathlib import Path
from typing import Any

import pysam

from analysis.qc import QCCalculator
from core.config import (
    PipelineConfig,
)
from core.exceptions import InvalidInputError
from file_io import IncrementalHDF5Writer, IncrementalTextWriter
from processing.processors import process_shards
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

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run(self) -> dict[str, Any]:
        start_time = time.time()

        # Validates the chromosome name and barcode tag before any heavy work.
        BAMReader(str(self.bam_path), self.config, self.barcodes)

        n_cells_input = len(self.barcode_list)
        if self.output_format == "hdf5":
            writer = IncrementalHDF5Writer(
                self.output_dir,
                self.config,
                self.barcode_list,
                barcode_metadata=self.barcode_metadata,
            )
        else:
            writer = IncrementalTextWriter(self.output_dir, self.config, self.barcode_list)

        totals = process_shards(self.bam_path, self.config, self.barcode_list, writer)

        cells_passed = totals["cells_passed"]
        if not cells_passed:
            logger.error("No cells passed quality filters")
            return {}

        if not self.config.dedup.skip and totals["total_reads"]:
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

        qc_dir = self.output_dir / "qc"
        writer.finalize(qc_dir)

        run_metadata = QCCalculator(self.config).collect_run_metadata(
            str(self.bam_path),
            str(self.output_dir),
            n_cells_input,
            cells_passed,
        )
        from file_io import write_run_summary

        write_run_summary(run_metadata, qc_dir / "summary.txt")

        if self.output_format == "hdf5":
            logger.info("Generating HTML QC report...")
            try:
                if self.barcode_metadata is not None:
                    # scATAC-seq: use ATAC report with Tn5 transposition plot
                    from analysis.report import generate_html_report

                    generate_html_report(
                        self.output_dir,
                        self.sample_name,
                        title=self.report_title,
                        subtitle=self.report_subtitle,
                        working_directory=self.working_directory,
                        input_dir=str(self.bam_path.parent),
                    )
                else:
                    # scRNA-seq: use RNA report with read start sites plot
                    from analysis.report import generate_scrna_html_report

                    generate_scrna_html_report(
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

        elapsed = time.time() - start_time
        hours = int(elapsed // 3600)
        minutes = int((elapsed % 3600) // 60)
        seconds = int(elapsed % 60)

        if hours > 0:
            logger.info("Pipeline complete")
            logger.info(f"Elapsed time: {hours}h {minutes}m {seconds}s")
        elif minutes > 0:
            logger.info("Pipeline complete")
            logger.info(f"Elapsed time: {minutes}m {seconds}s")
        else:
            logger.info("Pipeline complete")
            logger.info(f"Elapsed time: {seconds}s")

        return {
            "cells_processed": n_cells_input,
            "cells_passed_qc": cells_passed,
            "mean_reads": totals["kept_reads"] / cells_passed,
        }


def run_pipeline(
    bam_path: str,
    barcode_file: str | None = None,
    output_dir: str = "",
    sample_name: str = "mgatk",
    min_baseq: int = 20,
    min_mapq: int = 30,
    min_reads_per_cell: int = 1,
    max_strand_bias: float = 1.0,
    min_distance_from_end: int = 5,
    skip_deduplication: bool = False,
    use_fragment_length_dedup: bool = True,
    nh_max: int = 0,
    nm_max: int = 0,
    compute_tn5: bool = True,
    barcode_tag: str = "CB",
    min_barcode_reads: int = 1,
    mito_chr: str = "chrM",
    n_cores: int = 16,
    max_memory_gb: float = 128.0,
    output_format: str = "standard",
    report_title: str | None = None,
    report_subtitle: str | None = None,
    working_directory: str | None = None,
) -> dict[str, Any]:
    """Run the pipeline with individual parameters"""
    barcode_metadata = None

    if barcode_file == "bulk":
        barcodes = ["bulk"]
    elif barcode_file is None:
        from file_io.barcode_extraction import extract_barcodes_from_bam

        logger.info("No barcode file provided - extracting barcodes from BAM")
        barcodes = extract_barcodes_from_bam(
            bam_path, barcode_tag=barcode_tag, mito_chr=mito_chr, min_reads=min_barcode_reads
        )
        if not barcodes:
            raise InvalidInputError(
                f"No barcodes found in BAM file with tag '{barcode_tag}' "
                f"and minimum {min_barcode_reads} reads"
            )
    elif barcode_file.endswith(".csv"):
        from utils.utils import load_barcode_csv

        barcodes, barcode_metadata = load_barcode_csv(barcode_file)
    else:
        with open(barcode_file) as f:
            barcodes = [line.strip() for line in f if line.strip()]

    config = PipelineConfig(
        min_baseq=min_baseq,
        min_mapq=min_mapq,
        max_strand_bias=max_strand_bias,
        min_distance_from_end=min_distance_from_end,
        skip_deduplication=skip_deduplication,
        use_fragment_length_dedup=use_fragment_length_dedup,
        nh_max=nh_max,
        nm_max=nm_max,
        compute_tn5=compute_tn5,
        n_cores=n_cores,
        max_memory_gb=max_memory_gb,
        min_reads_per_cell=min_reads_per_cell,
        barcode_tag=barcode_tag,
        mito_chr=mito_chr,
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
