"""Per-cell CSV statistics, the run configuration, and the plaintext summary."""

import json
from pathlib import Path


def write_cell_stats(cell_stats: list[dict], output_path: Path):
    """Write per-cell QC statistics to CSV file"""
    if not cell_stats:
        return
    output_columns = ["barcode", "mean_depth", "coverage_breadth", "total_fragments", "total_reads"]

    with open(output_path, "w") as f:
        f.write(",".join(output_columns) + "\n")

        for stats in cell_stats:
            values = [str(stats[key]) for key in output_columns]
            f.write(",".join(values) + "\n")


def write_run_config(run_metadata: dict, output_path: Path):
    """The run metadata as JSON: the machine-readable twin of summary.txt."""
    with open(output_path, "w") as f:
        json.dump(run_metadata, f, indent=2, sort_keys=True, default=str)
        f.write("\n")


def write_run_summary(run_metadata: dict, output_path: Path):
    """Write run summary to text file."""
    parameters = run_metadata["parameters"]
    with open(output_path, "w") as f:
        f.write("mgatk2 Run Summary\n")
        f.write("=" * 20 + "\n")

        for key, value in run_metadata.items():
            if key != "parameters":
                f.write(f"{key}: {value}\n")

        f.write("\nParameters\n")
        f.write("-" * 20 + "\n")
        for key, value in parameters.items():
            f.write(f"{key}: {value}\n")
