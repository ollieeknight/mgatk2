"""Utility functions for CLI operations."""

import csv
import logging
import os
import sys
from datetime import datetime
from importlib.metadata import version
from pathlib import Path

from core.exceptions import InvalidInputError, ProcessingError
from core.pipeline import run_pipeline
from utils.utils import validate_bam_file, validate_barcode_file

logger = logging.getLogger(__name__)


def auto_detect_10x_structure(
    bam_path: str, barcode_file: str | None = None
) -> tuple[str, str | None]:
    """Auto-detect 10x Genomics output structure"""
    path = Path(bam_path)

    # Directory input - find BAM file
    if path.is_dir():
        # Check for outs directory
        bam_file: Path | None = None
        if path.name == "outs" or "outs" in str(path):
            bam_file = path / "possorted_bam.bam"
            if not bam_file.exists():
                bam_file = path.parent / "outs" / "possorted_bam.bam"
        else:
            # Look for outs subdirectory
            outs_dir = path / "outs"
            bam_file = outs_dir / "possorted_bam.bam" if outs_dir.exists() else None

        if bam_file and bam_file.exists():
            logger.info(f"Found 10x BAM: {bam_file}")
            bam_path = str(bam_file)
            if not barcode_file:
                barcode_file = _find_barcode_file(bam_file.parent)
        else:
            logger.warning(f"No possorted_bam.bam found in {path}")

    # File input - check if in 10x structure
    elif path.is_file() and not barcode_file:
        if path.parent.name == "outs":
            logger.info("Detected 10x BAM in outs directory")
            barcode_file = _find_barcode_file(path.parent)

    return str(Path(bam_path).resolve()), barcode_file


def _find_barcode_file(directory: Path) -> str | None:
    """Find barcode file in 10x directory."""
    # Try singlecell.csv first
    singlecell = directory / "singlecell.csv"
    if singlecell.exists():
        logger.info(f"Found {singlecell}")
        return str(singlecell)

    # Try barcode files
    for pattern in [
        "filtered_peak_bc_matrix/barcodes.tsv",
        "filtered_tf_bc_matrix/barcodes.tsv.gz",
    ]:
        bc_file = directory / pattern
        if bc_file.exists():
            logger.info(f"Found {bc_file}")
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

    # If BAM is in an 'outs' directory, get the parent of 'outs'
    if bam_path_obj.parent.name == "outs":
        # Get the parent of 'outs' directory
        project_dir = bam_path_obj.parent.parent
        return project_dir.name

    # If we auto-detected from a directory containing 'outs', use that directory name
    if "outs" in str(bam_path_obj.parent):
        # Find the directory that contains 'outs'
        current = bam_path_obj.parent
        while current.parent != current:  # Not root
            if current.name == "outs":
                return current.parent.name
            current = current.parent

    # Fallback: use the directory name containing the BAM
    return bam_path_obj.parent.name if bam_path_obj.parent.name != "." else "mgatk2"


def setup_file_logging(log_file_path):
    """Setup file logging"""
    logger = logging.getLogger("mgatk")

    # Create file handler
    file_handler = logging.FileHandler(log_file_path, mode="w")
    file_handler.setLevel(logging.INFO)

    # Create formatter
    formatter = logging.Formatter("%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    file_handler.setFormatter(formatter)

    # Add handler to logger
    logger.addHandler(file_handler)


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
    original_bam_path=None,
    report_title=None,
    report_subtitle=None,
    working_directory=None,
):
    """Common pipeline execution logic"""

    if verbose:
        # Only enable DEBUG for mgatk2 loggers, not third-party libraries
        for logger_name in ['mgatk', __name__]:
            logging.getLogger(logger_name).setLevel(logging.DEBUG)

    __version__ = version("mgatk2")
    logger.info("mgatk2 version %s", __version__)

    # Auto-detect 10x Genomics structure only if barcode_file is not explicitly None
    if barcode_file != "bulk":
        bam_path, barcode_file = auto_detect_10x_structure(bam_path, barcode_file)

    name = "output_"

    # Convert dedup_mode
    skip_dedup = dedup_mode.lower() == "none"
    use_fragment_length_dedup = dedup_mode.lower() in [
        "alignment_and_fragment_length",
        "fragment-length",
        "hybrid",
    ]

    # Setup file logging if not in dry run mode
    if not dry_run:
        log_file = Path(output_dir) / "output.log"
        log_file.parent.mkdir(parents=True, exist_ok=True)
        setup_file_logging(log_file)

        # Log the command that was run
        cmd_args = sys.argv
        cmd_path = os.path.realpath(cmd_args[0]) if cmd_args else "mgatk2"
        full_command = f"{cmd_path} {' '.join(cmd_args[1:])}"
        logger.info("Command executed: %s", full_command)
        logger.info("Working directory: %s", os.getcwd())
        logger.info("Execution time: %s", datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

    try:
        validate_bam_file(bam_path)
        if barcode_file:
            validate_barcode_file(barcode_file)

        os.makedirs(output_dir, exist_ok=True)

        mito_chr = normalise_mito_chr(mito_genome)

        # Barcode file will be auto-extracted in run_pipeline if None
        if barcode_file is None:
            logger.info("No barcode file provided - will extract barcodes from BAM")
        else:
            logger.info("Using barcode file: %s", barcode_file)
        barcode_file_to_cleanup = None

        # Log configuration
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
        )

        if dry_run:
            return 0

        # Determine title for HTML report (use passed values or generate from original BAM path)
        if report_title is None:
            # Use original BAM path if provided, otherwise current bam_path
            source_path = original_bam_path if original_bam_path else bam_path
            report_title = get_10x_parent_directory_name(source_path)
        if report_subtitle is None:
            report_subtitle = "mgatk2 output analysis"

        # Run pipeline
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
            "n_cores": _determine_cores(ncores),
            "batch_size": batch_size or _determine_cores(ncores),
            "skip_deduplication": skip_dedup,
            "use_fragment_length_dedup": use_fragment_length_dedup,
            "sequential": sequential,
            "report_title": report_title,
            "report_subtitle": report_subtitle,
            "working_directory": working_directory,
        }

        if max_memory is not None:
            run_args["max_memory_gb"] = max_memory

        run_pipeline(**run_args)

        # Clean up auto-detected barcode file if it was created
        if barcode_file_to_cleanup:
            try:
                os.remove(barcode_file_to_cleanup)
                logger.info(f"Cleaned up auto-detected barcode file: {barcode_file_to_cleanup}")
            except Exception as e:
                logger.warning("Could not remove auto-detected barcode file: %s", e)

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
        # Check for SLURM environment variables
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
):
    """Log the pipeline configuration."""

    logger.info("  Input BAM:              %s", os.path.realpath(bam_path))
    logger.info(
        "  Input barcodes:         %s",
        os.path.realpath(barcode_file) if barcode_file else "None (auto-detect from BAM)",
    )
    logger.info("  Output directory:       %s", os.path.realpath(output_dir))
    logger.info("  BAM prefix:             %s", mito_chr)

    actual_cores = _determine_cores(ncores)
    if ncores is None:
        # Check for SLURM environment variables
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
        batch_size = 250

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
    logger.info("  Deduplication:          %s", dedup_display)
    logger.info("  Batch size:             %s cells", batch_size)

    format_display = "text files (.txt.gz)" if output_format == "txt" else "HDF5 (.h5)"
    logger.info("  Output format:          %s", format_display)
    logger.info("  Sequential processing:  %s", sequential)
    logger.info("  Using mitochondrial chromosome: %s", mito_chr)

    # Count only cell barcodes (is_cell_barcode = 1) if it's a CSV file
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

            # Determine which column to use
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
        # For plain text barcode files, count all non-empty lines
        with open(barcode_file) as f:
            n_barcodes = sum(1 for line in f if line.strip())
    else:
        # For bulk calling
        logger.info("  Barcodes:               bulk (all reads)")
        return
    logger.info("  Barcodes:               %s", n_barcodes)

    if max_memory:
        logger.info("  Max memory limit:       %sGB", max_memory)
