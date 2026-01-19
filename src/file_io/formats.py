"""Format utilities and writers for CSV/TSV/JSON output"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def write_cell_stats(cell_stats: list[dict], output_path: Path):
    """Write per-cell QC statistics to CSV file"""
    if not cell_stats:
        return

    # Define the columns we want to output in order
    output_columns = ["barcode", "mean_depth", "coverage_breadth", "total_fragments", "total_reads"]

    with open(output_path, "w") as f:
        # Write header
        f.write(",".join(output_columns) + "\n")

        # Write data
        for stats in cell_stats:
            values = [str(stats.get(key, "NA")) for key in output_columns]
            f.write(",".join(values) + "\n")


def write_position_stats(position_stats: dict, output_path: Path, mito_length: int):
    """Write per-position statistics to CSV file."""
    with open(output_path, "w") as f:
        # Write header
        f.write(
            "position,n_cells_any_coverage,n_cells_10x,n_cells_50x,"
            "mean_depth,median_depth,coverage_cv\n"
        )

        # Write data for each position
        for pos in range(1, mito_length + 1):
            if pos in position_stats:
                stats = position_stats[pos]
                f.write(
                    f"{pos},{stats['n_cells_any']},{stats['n_cells_10x']},"
                    f"{stats['n_cells_50x']},{stats['mean_cov']:.2f},"
                    f"{stats['median_cov']:.2f},{stats['cv']:.4f}\n"
                )
            else:
                f.write(f"{pos},0,0,0,0.00,0.00,NA\n")


def write_run_summary(run_metadata: dict, output_path: Path):
    """Write run summary to text file."""
    with open(output_path, "w") as f:
        f.write("mgatk2 Run Summary\n")
        f.write("=" * 20 + "\n")

        for key, value in run_metadata.items():
            if key != "parameters":
                f.write(f"{key}: {value}\n")
