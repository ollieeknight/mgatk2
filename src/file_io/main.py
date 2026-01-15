"""Main save function for mgatk2 outputs."""

import logging

from .formats import write_cell_stats, write_position_stats, write_run_summary
from .writers import write_parameters_json

logger = logging.getLogger(__name__)


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
    _write_qc_and_metadata(
        qc_dir, cell_stats, position_stats, run_metadata, config
    )


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
    _write_qc_and_metadata(
        qc_dir, cell_stats, position_stats, run_metadata, config
    )


def _write_qc_and_metadata(
    qc_dir, cell_stats, position_stats, run_metadata, config
):
    """Write QC and metadata files."""
    # Write QC files
    if cell_stats:
        write_cell_stats(cell_stats, qc_dir / "cell_stats.csv")

    if position_stats:
        write_position_stats(
            position_stats,
            qc_dir / "position_stats.csv",
            config.mito_length,
        )

    # Write metadata
    metadata_dir = qc_dir.parent / "metadata"
    metadata_dir.mkdir(exist_ok=True)

    if run_metadata:
        write_run_summary(run_metadata, metadata_dir / "run_summary.txt")

        # Write parameters as JSON
        if "parameters" in run_metadata:
            write_parameters_json(
                run_metadata["parameters"],
                metadata_dir / "parameters.json",
            )
