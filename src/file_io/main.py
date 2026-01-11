"""Main save function for mgatk2 outputs."""

import logging
from collections import defaultdict
from pathlib import Path

import numpy as np

from core.config import PipelineConfig

from .formats import (
    write_cell_stats,
    write_position_stats,
    write_run_summary,
)
from .writers import write_parameters_json

logger = logging.getLogger(__name__)


def save_mgatk_outputs(
    output_dir: Path,
    cell_results: list[dict],
    config: PipelineConfig,
    cell_stats: list[dict] | None = None,
    position_stats: dict | None = None,
    run_metadata: dict | None = None,
    output_format: str = "hdf5",
):
    """Save output files"""

    # Create output directories
    output_data_dir = output_dir / "output"
    qc_dir = output_dir / "qc"

    output_data_dir.mkdir(exist_ok=True)
    qc_dir.mkdir(exist_ok=True)

    # Initialise data structures
    bases = ["A", "C", "G", "T"]
    strands = ["fwd", "rev"]
    base_data: dict[str, list] = {f"{base}_{strand}": [] for base in bases for strand in strands}
    coverage_data = []

    # Track reference allele by counting most common base at each position
    position_base_counts: dict[int, dict[str, int]] = defaultdict(lambda: {"A": 0, "C": 0, "G": 0, "T": 0})

    # Track depth per cell
    cell_depths = {}

    # Process each cell's pileup data
    for result in cell_results:
        barcode = result["barcode"]
        pileup = result["pileup"]

        # Calculate mean depth for this cell
        depths = [counts["depth"] for counts in pileup.values()]
        mean_depth = np.mean(depths) if depths else 0
        cell_depths[barcode] = mean_depth

        for pos, counts in pileup.items():
            # Convert 0-based BAM position to 1-based
            pos_1based = pos + 1

            # Coverage for this position
            total_cov = counts["depth"]
            if total_cov > 0:
                coverage_data.append((pos_1based, barcode, total_cov))

            # Base counts for each nucleotide
            for base in bases:
                fwd_count = counts.get(f"{base}_fwd", 0)
                rev_count = counts.get(f"{base}_rev", 0)
                total_count = counts.get(base, 0)

                # Save strand-specific counts separately
                if fwd_count > 0:
                    base_data[f"{base}_fwd"].append((pos_1based, barcode, fwd_count))
                if rev_count > 0:
                    base_data[f"{base}_rev"].append((pos_1based, barcode, rev_count))

                # Track base counts for reference determination
                if total_count > 0:
                    position_base_counts[pos_1based][base] += total_count

    logger.info("Writing output files...")

    # Determine reference allele at each position (most common base)
    ref_alleles = {}
    for pos in range(1, config.mito_length + 1):
        if pos in position_base_counts:
            counts = position_base_counts[pos]
            # Get base with maximum count
            max_base = max(bases, key=lambda b: counts[b])
            if counts[max_base] > 0:
                ref_alleles[pos] = max_base
            else:
                ref_alleles[pos] = "N"
        else:
            ref_alleles[pos] = "N"

    if output_format == "txt":
        # Original mgatk format
        _save_original_format(
            output_data_dir,
            base_data,
            coverage_data,
            cell_depths,
            ref_alleles,
            config,
            cell_stats,
            position_stats,
            run_metadata,
            qc_dir,
        )
    else:
        # HDF5 format (default)
        _save_hdf5_format(
            output_data_dir,
            base_data,
            coverage_data,
            cell_depths,
            ref_alleles,
            config,
            cell_stats,
            position_stats,
            run_metadata,
            qc_dir,
        )


def _save_hdf5_format(
    output_dir,
    base_data,
    coverage_data,
    cell_depths,
    ref_alleles,
    config,
    cell_stats,
    position_stats,
    run_metadata,
    qc_dir,
):
    """Save outputs in HDF5 format."""

    # Create HDF5 files

    # Write QC and metadata files
    _write_qc_and_metadata(qc_dir, cell_stats, position_stats, run_metadata, config)


def _save_original_format(
    output_dir,
    base_data,
    coverage_data,
    cell_depths,
    ref_alleles,
    config,
    cell_stats,
    position_stats,
    run_metadata,
    qc_dir,
):
    """Save outputs in original mgatk format."""
    import gzip

    bases = ["A", "C", "G", "T"]
    logger.info("Writing original mgatk format files...")

    # Write base count files (A.txt.gz, C.txt.gz, etc.)
    for base in bases:
        base_file = output_dir / f"{base}.txt.gz"
        logger.info("  Writing %s", base_file)

        with gzip.open(base_file, "wt") as f:
            f.write("pos,cell,fwd,rev\n")  # Header

            # Combine forward and reverse counts for each position/cell
            pos_cell_counts = {}

            # Collect forward counts
            for pos, cell, count in base_data[f"{base}_fwd"]:
                key = (pos, cell)
                if key not in pos_cell_counts:
                    pos_cell_counts[key] = [0, 0]  # [fwd, rev]
                pos_cell_counts[key][0] = count

            # Collect reverse counts
            for pos, cell, count in base_data[f"{base}_rev"]:
                key = (pos, cell)
                if key not in pos_cell_counts:
                    pos_cell_counts[key] = [0, 0]
                pos_cell_counts[key][1] = count

            # Write combined data
            for (pos, cell), (fwd, rev) in sorted(pos_cell_counts.items()):
                f.write(f"{pos},{cell},{fwd},{rev}\n")

    # Write coverage file
    coverage_file = output_dir / "coverage.txt.gz"
    logger.info("  Writing %s", coverage_file)

    with gzip.open(coverage_file, "wt") as f:
        f.write("pos,cell,cov\n")
        for pos, cell, cov in sorted(coverage_data):
            f.write(f"{pos},{cell},{cov}\n")

    # Write depth table
    depth_file = output_dir / "depthTable.txt"
    logger.info("  Writing %s", depth_file)

    with open(depth_file, "w") as f:
        f.write("cell\tmean_depth\n")
        for cell, depth in sorted(cell_depths.items()):
            f.write(f"{cell}\t{depth:.2f}\n")

    # Write reference allele file
    ref_file = output_dir / f"{config.mito_chr}_refAllele.txt"
    logger.info("  Writing %s", ref_file)

    with open(ref_file, "w") as f:
        f.write("pos\tref\n")
        for pos in range(1, config.mito_length + 1):
            ref_base = ref_alleles.get(pos, "N")
            f.write(f"{pos}\t{ref_base}\n")

    # Write QC and metadata files
    _write_qc_and_metadata(qc_dir, cell_stats, position_stats, run_metadata, config)


def _write_qc_and_metadata(qc_dir, cell_stats, position_stats, run_metadata, config):
    """Write QC and metadata files."""
    # Write QC files
    if cell_stats:
        write_cell_stats(cell_stats, qc_dir / "cell_stats.csv")

    if position_stats:
        write_position_stats(position_stats, qc_dir / "position_stats.csv", config.mito_length)

    # Write metadata
    metadata_dir = qc_dir.parent / "metadata"
    metadata_dir.mkdir(exist_ok=True)

    if run_metadata:
        write_run_summary(run_metadata, metadata_dir / "run_summary.txt")

        # Write parameters as JSON
        if "parameters" in run_metadata:
            write_parameters_json(run_metadata["parameters"], metadata_dir / "parameters.json")
