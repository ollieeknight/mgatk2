"""Format utilities and writers for CSV/TSV/JSON output"""

import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def write_cell_stats(cell_stats: list[dict], output_path: Path):
    """Write per-cell QC statistics to CSV file"""
    if not cell_stats:
        return

    output_columns = ["barcode", "mean_depth", "coverage_breadth", "total_fragments", "total_reads"]

    with open(output_path, "w") as f:
        f.write(",".join(output_columns) + "\n")

        for stats in cell_stats:
            values = [str(stats.get(key, "NA")) for key in output_columns]
            f.write(",".join(values) + "\n")


def write_run_summary(run_metadata: dict, output_path: Path):
    """Write run summary to text file."""
    with open(output_path, "w") as f:
        f.write("mgatk2 Run Summary\n")
        f.write("=" * 20 + "\n")

        for key, value in run_metadata.items():
            if key != "parameters":
                f.write(f"{key}: {value}\n")
