"""Shared plumbing for the single-cell commands."""

import logging
import os
import sys
from dataclasses import asdict
from importlib.metadata import version
from pathlib import Path

import pysam

from core.config import PipelineConfig
from core.exceptions import InvalidInputError, MgatkError
from core.pipeline import run_pipeline
from data.blacklists import load_bed_positions
from processing.readers import resolve_mito_contig
from utils.utils import has_alignment_index, validate_bam_file, validate_barcode_file

logger = logging.getLogger(__name__)


def auto_detect_10x_structure(
    bam_path: str, barcode_file: str | None = None
) -> tuple[str, str | None]:
    """Auto-detect 10x Genomics output structure (scATAC or 10x Multi)"""
    path = Path(bam_path)

    if path.is_dir():
        candidates = [path / "possorted_bam.bam", path / "outs" / "possorted_bam.bam"]
        bam_file = next((candidate for candidate in candidates if candidate.exists()), None)
        if bam_file is None:
            bam_file = _find_10x_multi_bam(path)
        if bam_file is None:
            logger.warning(f"No possorted_bam.bam found in {path}")
        else:
            bam_path = str(bam_file)
            if not barcode_file:
                barcode_file = _find_barcode_file(bam_file.parent)

    elif path.is_file() and not barcode_file and path.parent.name in ("outs", "count"):
        barcode_file = _find_barcode_file(path.parent)

    return str(Path(bam_path).resolve()), barcode_file


def _find_10x_multi_bam(root: Path) -> Path | None:
    """Find a 10x Multi per-sample BAM under <root>/[outs/]per_sample_outs/*/count/."""
    for per_sample_outs in [root / "per_sample_outs", root / "outs" / "per_sample_outs"]:
        matches = sorted(per_sample_outs.glob("*/count/sample_alignments.bam"))
        if len(matches) == 1:
            logger.info(f"Detected 10x Multi sample: {matches[0].parent.parent.name}")
            return matches[0]
        if len(matches) > 1:
            samples = ", ".join(m.parent.parent.name for m in matches)
            raise InvalidInputError(
                f"Multiple 10x Multi samples found under {per_sample_outs}: {samples}. "
                "Pass --input pointing at the specific sample's sample_alignments.bam."
            )
    return None


def _find_barcode_file(directory: Path) -> str | None:
    """Find barcode file in 10x directory."""
    for name in [
        "singlecell.csv",
        "sample_filtered_barcodes.csv",
        "filtered_peak_bc_matrix/barcodes.tsv",
        "filtered_tf_bc_matrix/barcodes.tsv.gz",
    ]:
        if (directory / name).exists():
            return str(directory / name)

    logger.warning("No barcode file found")
    return None


def load_panel_positions(panel_bed: str, mito_chr: str) -> frozenset[int]:
    """1-based targeted positions from an amplicon panel BED."""
    positions = load_bed_positions(panel_bed, mito_chr)
    if not positions:
        raise InvalidInputError(
            f"Panel BED {panel_bed} defines no positions on {mito_chr}; "
            "check the contig name in column 1"
        )
    return frozenset(positions)


def check_alignment(path: str, mito_chr: str, reference_filename: str | None = None) -> None:
    """Dry-run check: the contig and an index exist. Creates no files, indexes included."""
    try:
        with pysam.AlignmentFile(path, reference_filename=reference_filename) as alignment:
            references = alignment.references
    except Exception as exc:
        raise InvalidInputError(f"Cannot read alignment {path}: {exc}") from exc

    resolve_mito_contig(references, mito_chr, path)
    if not has_alignment_index(path):
        raise InvalidInputError(f"No index beside {path}; the run would have to build one")


def report_title(bam_path: str) -> str:
    """Name a run after its 10x run or sample directory, else the BAM's directory."""
    parent = Path(bam_path).resolve().parent
    # <run>/outs/possorted_bam.bam, or 10x Multi <sample>/count/sample_alignments.bam
    if parent.name in ("outs", "count"):
        return parent.parent.name
    return parent.name


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


def determine_cores(ncores: int | None) -> int:
    """--threads if given, else the SLURM allocation, else every CPU."""
    if ncores:
        return ncores
    for variable in ("SLURM_CPUS_PER_TASK", "SLURM_NTASKS"):
        value = os.environ.get(variable, "")
        if value.isdigit():
            return int(value)
    return os.cpu_count() or 1


def run_pipeline_command(
    bam_path,
    output_dir,
    mito_genome="chrM",
    barcode_file=None,
    barcode_tag="CB",
    min_barcode_reads=10,
    ncores=None,
    verbose=False,
    max_memory=128.0,
    base_qual=20,
    min_mapq=30,
    min_reads=1,
    max_strand_bias=1.0,
    min_distance_from_end=5,
    dedup_mode="alignment_and_fragment_length",
    output_format="hdf5",
    dry_run=False,
    nh_max=0,
    nm_max=0,
    compute_tn5=True,
    assay=None,
    panel_bed=None,
) -> int:
    """Run one sample from command-line options and return its exit status."""
    if verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    logger.info("mgatk2 version %s", version("mgatk2"))

    try:
        if barcode_file != "bulk":
            bam_path, barcode_file = auto_detect_10x_structure(bam_path, barcode_file)
            if barcode_file:
                validate_barcode_file(barcode_file)

        dedup_mode = dedup_mode.lower()
        config = PipelineConfig(
            min_baseq=base_qual,
            min_mapq=min_mapq,
            max_strand_bias=max_strand_bias,
            min_distance_from_end=min_distance_from_end,
            nh_max=nh_max,
            nm_max=nm_max,
            skip_deduplication=dedup_mode == "none",
            use_fragment_length_dedup=dedup_mode == "alignment_and_fragment_length",
            n_cores=determine_cores(ncores),
            max_memory_gb=max_memory,
            min_reads_per_cell=min_reads,
            barcode_tag=barcode_tag,
            mito_chr=mito_genome,
            compute_tn5=compute_tn5,
            panel_positions=load_panel_positions(panel_bed, mito_genome) if panel_bed else None,
        )

        if dry_run:
            check_alignment(bam_path, config.mito_chr)
            _log_configuration(bam_path, barcode_file, output_dir, output_format, config)
            return 0

        Path(output_dir).mkdir(parents=True, exist_ok=True)
        setup_file_logging(Path(output_dir) / "output.log")
        logger.info("Command executed: %s", " ".join(sys.argv))
        logger.info("Working directory: %s", os.getcwd())
        validate_bam_file(bam_path)
        _log_configuration(bam_path, barcode_file, output_dir, output_format, config)

        run_pipeline(
            bam_path,
            output_dir,
            config,
            barcode_file=barcode_file,
            min_barcode_reads=min_barcode_reads,
            output_format=output_format.lower(),
            assay=assay,
            title=report_title(bam_path),
        )
        return 0

    except KeyboardInterrupt:
        logger.info("Interrupted by user")
        return 130
    except MgatkError as e:
        logger.error("%s", e)
        return 1
    except Exception:
        logger.exception("Unexpected error")
        return 1


def _log_configuration(bam_path, barcode_file, output_dir, output_format, config):
    logger.info("  Input BAM:        %s", os.path.realpath(bam_path))
    logger.info("  Barcodes:         %s", barcode_file or "auto-detect from BAM")
    logger.info("  Output:           %s (%s)", os.path.realpath(output_dir), output_format)
    for name, value in asdict(config).items():
        if name == "panel_positions":
            value = len(value) if value else "none"
        logger.info("  %-17s %s", name + ":", value)
