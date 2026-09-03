"""Per-cell CSV statistics, the run configuration, and the plaintext summary."""

import json
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


def write_run_config(run_metadata: dict, output_path: Path):
    """Write the run metadata as JSON.

    summary.txt is for reading; this is the same record for machines, so a
    downstream script never has to parse the summary's colon-separated lines.
    """
    with open(output_path, "w") as f:
        json.dump(run_metadata, f, indent=2, sort_keys=True, default=str)
        f.write("\n")


def write_run_summary(run_metadata: dict, output_path: Path):
    """Write run summary to text file."""
    parameters = run_metadata.get("parameters") or {}
    with open(output_path, "w") as f:
        f.write("mgatk2 Run Summary\n")
        f.write("=" * 20 + "\n")

        for key, value in run_metadata.items():
            if key != "parameters":
                f.write(f"{key}: {value}\n")

        # QCCalculator collects these; dropping them left summary.txt with no
        # record of what the run was actually configured with.
        if parameters:
            f.write("\nParameters\n")
            f.write("-" * 20 + "\n")
            for key, value in parameters.items():
                f.write(f"{key}: {value}\n")
