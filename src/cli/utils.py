"""Utility functions for CLI operations."""

import csv
import logging
import os
import sys
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

from core.exceptions import InvalidInputError, ProcessingError
from utils.utils import validate_bam_file, validate_barcode_file

logger = logging.getLogger(__name__)


def auto_detect_10x_structure(
    bam_path: str, barcode_file: str | None = None
) -> tuple[str, str | None]:
    """Auto-detect 10x Genomics output structure"""
    path = Path(bam_path)

    if path.is_dir():
        candidates = [path / "possorted_bam.bam", path / "outs" / "possorted_bam.bam"]
        bam_file = next((candidate for candidate in candidates if candidate.exists()), None)
        if bam_file:
            bam_path = str(bam_file)
            if not barcode_file:
                barcode_file = _find_barcode_file(bam_file.parent)
        else:
            logger.warning(f"No possorted_bam.bam found in {path}")

    elif path.is_file() and not barcode_file:
        if path.parent.name == "outs":
            logger.info("Detected 10x BAM in outs directory")
            barcode_file = _find_barcode_file(path.parent)

    return str(Path(bam_path).resolve()), barcode_file


def _find_barcode_file(directory: Path) -> str | None:
    """Find barcode file in 10x directory."""
    singlecell = directory / "singlecell.csv"
    if singlecell.exists():
        return str(singlecell)

    for pattern in [
        "filtered_peak_bc_matrix/barcodes.tsv",
        "filtered_tf_bc_matrix/barcodes.tsv.gz",
    ]:
        bc_file = directory / pattern
        if bc_file.exists():
            return str(bc_file)

    logger.warning("No barcode file found")
    return None


def normalise_mito_chr(mito_genome: str) -> str:
    """Normalise mitochondrial chromosome name"""
    if mito_genome.upper() in ["M", "MT"]:
        return "chrM"
    if mito_genome in ["chrM", "chrMT"]:
        return mito_genome
    logger.warning("Unusual mitochondrial chromosome name: %s", mito_genome)
    return mito_genome


def get_10x_parent_directory_name(bam_path: str) -> str:
    """Extract the parent directory name when processing 10x data"""
    bam_path_obj = Path(bam_path)

    if bam_path_obj.parent.name == "outs":
        return bam_path_obj.parent.parent.name
    return bam_path_obj.parent.name if bam_path_obj.parent.name != "." else "mgatk2"


def setup_file_logging(log_file_path):
    """Write application logs to one deterministic file."""
    root_logger = logging.getLogger()
    for handler in list(root_logger.handlers):
        if getattr(handler, "_mgatk2_file_handler", False):
            root_logger.removeHandler(handler)
            handler.close()

    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setLevel(logging.INFO)
    file_handler._mgatk2_file_handler = True
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)
    root_logger.addHandler(file_handler)


def run_pipeline_command(
    bam_path,
    output_dir,
    barcode_file,
    barcode_tag,
    min_barcode_reads,
    mito_genome,
    ncores,
    verbose,
    batch_size,
    max_memory,
    base_qual,
    min_mapq,
    min_reads,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    sequential,
    dry_run=False,
    nh_max=0,
    nm_max=0,
    pileup_mode="classic",
    compute_tn5=True,
    original_bam_path=None,
    report_title=None,
    report_subtitle=None,
    working_directory=None,
):
    """Common pipeline execution logic"""

    if verbose:
        for logger_name in ["mgatk", __name__]:
            logging.getLogger(logger_name).setLevel(logging.DEBUG)

    __version__ = version("mgatk2")
    logger.info("mgatk2 version %s", __version__)

    if barcode_file != "bulk":
        bam_path, barcode_file = auto_detect_10x_structure(bam_path, barcode_file)

    name = "output_"

    skip_dedup = dedup_mode.lower() == "none"
    use_fragment_length_dedup = dedup_mode.lower() in [
        "alignment_and_fragment_length",
        "fragment-length",
        "hybrid",
    ]

    if not dry_run:
        log_file = Path(output_dir) / "output.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        setup_file_logging(log_file)

        cmd_args = sys.argv
        cmd_path = os.path.realpath(cmd_args[0]) if cmd_args else "mgatk2"
        full_command = f"{cmd_path} {' '.join(cmd_args[1:])}"
        logger.info("Command executed: %s", full_command)
        logger.info("Working directory: %s", os.getcwd())
        logger.info("Execution time: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    try:
        validate_bam_file(bam_path)
        if barcode_file and barcode_file != "bulk":
            validate_barcode_file(barcode_file)

        os.makedirs(output_dir, exist_ok=True)

        mito_chr = normalise_mito_chr(mito_genome)

        if barcode_file is None:
            logger.info("No barcode file provided - will extract barcodes from BAM")
        _log_configuration(
            bam_path,
            barcode_file,
            output_dir,
            mito_chr,
            ncores,
            barcode_tag,
            min_barcode_reads,
            base_qual,
            min_mapq,
            min_reads,
            max_strand_bias,
            min_distance_from_end,
            dedup_mode,
            output_format,
            sequential,
            batch_size,
            max_memory,
            nh_max=nh_max,
            nm_max=nm_max,
            pileup_mode=pileup_mode,
            compute_tn5=compute_tn5,
        )

        if dry_run:
            return 0

        if report_title is None:
            source_path = original_bam_path if original_bam_path else bam_path
            report_title = get_10x_parent_directory_name(source_path)
        if report_subtitle is None:
            report_subtitle = "mgatk2 output analysis"

        actual_cores = _determine_cores(ncores)

        worker_batch = batch_size if batch_size is not None else actual_cores
        logger.info(f"  Worker batch size:      {worker_batch} cells")

        run_args = {
            "bam_path": bam_path,
            "barcode_file": barcode_file,
            "output_dir": output_dir,
            "sample_name": name,
            "min_baseq": base_qual,
            "min_mapq": min_mapq,
            "min_reads_per_cell": min_reads,
            "output_format": output_format.lower(),
            "max_strand_bias": max_strand_bias,
            "min_distance_from_end": min_distance_from_end,
            "barcode_tag": barcode_tag,
            "min_barcode_reads": min_barcode_reads,
            "mito_chr": mito_chr,
            "n_cores": actual_cores,
            "worker_batch_size": worker_batch,
            "io_batch_size": None,  # Will be determined dynamically in pipeline
            "skip_deduplication": skip_dedup,
            "use_fragment_length_dedup": use_fragment_length_dedup,
            "nh_max": nh_max,
            "nm_max": nm_max,
            "pileup_mode": pileup_mode,
            "compute_tn5": compute_tn5,
            "sequential": sequential,
            "report_title": report_title,
            "report_subtitle": report_subtitle,
            "working_directory": working_directory,
        }

        if max_memory is not None:
            run_args["max_memory_gb"] = max_memory

        from core.pipeline import run_pipeline

        run_pipeline(**run_args)

        return 0

    except InvalidInputError as e:
        logger.error("Input validation failed: %s", e)
        return 1
    except ProcessingError as e:
        logger.error("Processing failed: %s", e)
        return 1
    except Exception as e:
        logger.error("Unexpected error: %s", e)
        if verbose:
            import traceback

            traceback.print_exc()
        return 1


def _determine_cores(ncores):
    """Determine number of cores to use."""
    import multiprocessing

    if ncores is None:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        slurm_ntasks = os.environ.get("SLURM_NTASKS")

        if slurm_cpus:
            try:
                actual_cores = int(slurm_cpus)
            except ValueError:
                actual_cores = max(1, multiprocessing.cpu_count())
        elif slurm_ntasks:
            try:
                actual_cores = int(slurm_ntasks)
            except ValueError:
                actual_cores = max(1, multiprocessing.cpu_count())
        else:
            actual_cores = max(1, multiprocessing.cpu_count())
    else:
        actual_cores = ncores

    return actual_cores


def _log_configuration(
    bam_path,
    barcode_file,
    output_dir,
    mito_chr,
    ncores,
    barcode_tag,
    min_barcode_reads,
    base_qual,
    min_mapq,
    min_reads,
    max_strand_bias,
    min_distance_from_end,
    dedup_mode,
    output_format,
    sequential,
    batch_size,
    max_memory,
    nh_max=0,
    nm_max=0,
    pileup_mode="classic",
    compute_tn5=True,
):
    """Log the pipeline configuration."""

    logger.info("  Input BAM:              %s", os.path.realpath(bam_path))
    logger.info(
        "  Input barcodes:         %s",
        (
            "bulk (all reads)"
            if barcode_file == "bulk"
            else os.path.realpath(barcode_file)
            if barcode_file
            else "None (auto-detect from BAM)"
        ),
    )
    logger.info("  Output directory:       %s", os.path.realpath(output_dir))
    logger.info("  BAM prefix:             %s", mito_chr)

    actual_cores = _determine_cores(ncores)
    if ncores is None:
        slurm_cpus = os.environ.get("SLURM_CPUS_PER_TASK")
        slurm_ntasks = os.environ.get("SLURM_NTASKS")

        if slurm_cpus:
            cores_msg = f"{actual_cores} (SLURM CPUS_PER_TASK)"
        elif slurm_ntasks:
            cores_msg = f"{actual_cores} (SLURM NTASKS)"
        else:
            cores_msg = f"{actual_cores} (all available)"
    else:
        cores_msg = str(actual_cores)
    logger.info("  Cores:                  %s", cores_msg)

    if batch_size is None:
        batch_size = actual_cores
        batch_msg = f"{batch_size} (matches cores)"
    else:
        batch_msg = str(batch_size)

    logger.info("  Barcode tag:            %s", barcode_tag)
    if barcode_file is None:
        logger.info("  Min barcode reads:      %s", min_barcode_reads)
    logger.info("  Min base quality:       %s", base_qual)
    logger.info("  Min mapping quality:    %s", min_mapq)
    logger.info("  Min reads per cell:     %s", min_reads)

    if dedup_mode.lower() == "none":
        dedup_display = "disabled"
    elif dedup_mode.lower() in ["alignment_and_fragment_length", "hybrid", "fragment"]:
        dedup_display = "alignment + strand + fragment length"
    else:
        dedup_display = "alignment + strand only"

    logger.info("  Max strand bias:        %s", max_strand_bias)
    logger.info("  Min dist from end:      %sbp", min_distance_from_end)
    logger.info("  NH max (multi-mapper):  %s", nh_max if nh_max > 0 else "disabled")
    logger.info("  NM max (mismatches):    %s", nm_max if nm_max > 0 else "disabled")
    logger.info("  Pileup mode:            %s", pileup_mode)
    logger.info("  Compute Tn5 cuts:       %s", compute_tn5)
    logger.info("  Deduplication:          %s", dedup_display)
    logger.info("  Worker batch size:      %s cells", batch_msg)

    format_display = "text files (.txt.gz)" if output_format == "txt" else "HDF5 (.h5)"
    logger.info("  Output format:          %s", format_display)
    logger.info("  Sequential processing:  %s", sequential)
    logger.info("  Using mitochondrial chromosome: %s", mito_chr)

    if barcode_file == "bulk":
        logger.info("  Barcodes:               bulk (all reads)")
        return
    if barcode_file and barcode_file.endswith(".csv"):
        n_barcodes = 0
        total_rows = 0
        column_found = None

        with open(barcode_file) as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames

            if headers is None:
                logger.warning("CSV file has no headers")
                return

            if "is__cell_barcode" in headers:
                column_found = "is__cell_barcode"
            elif "is_cell_barcode" in headers:
                column_found = "is_cell_barcode"
            elif "is_cell" in headers:
                column_found = "is_cell"

            for row in reader:
                total_rows += 1
                if column_found:
                    is_cell = row.get(column_found, "0")
                    if is_cell in ["1", "1.0", "True", "true", "TRUE"]:
                        n_barcodes += 1
    elif barcode_file:
        with open(barcode_file) as f:
            n_barcodes = sum(1 for line in f if line.strip())
    else:
        logger.info("  Barcodes:               bulk (all reads)")
        return
    logger.info("  Barcodes:               %s", n_barcodes)

    if max_memory:
        logger.info("  Max memory limit:       %sGB", max_memory)
