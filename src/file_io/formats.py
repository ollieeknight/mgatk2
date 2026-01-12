"""Format utilities and writers for CSV/TSV/JSON output"""

import gzip
import json
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def write_cell_stats(cell_stats: list[dict], output_path: Path):
    """Write per-cell QC statistics to CSV file"""
    if not cell_stats:
        return

    # Get all unique keys from all cells
    all_keys: set[str] = set()
    for stats in cell_stats:
        all_keys.update(stats.keys())

    # Sort keys for consistent output
    sorted_keys = sorted(all_keys)

    with open(output_path, "w") as f:
        # Write header
        f.write(",".join(sorted_keys) + "\n")

        # Write data
        for stats in cell_stats:
            values = [str(stats.get(key, "NA")) for key in sorted_keys]
            f.write(",".join(values) + "\n")

    logger.info("Wrote cell stats to: %s", output_path)


def write_position_stats(position_stats: dict, output_path: Path, mito_length: int):
    """Write per-position statistics to CSV file."""
    with open(output_path, "w") as f:
        # Write header
        f.write(
            "position,n_cells_any_coverage,n_cells_10x,n_cells_50x,"
            "mean_coverage,median_coverage,coverage_cv\n"
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

    logger.info("Wrote position stats to: %s", output_path)


def write_run_summary(run_metadata: dict, output_path: Path):
    """Write run summary to text file."""
    with open(output_path, "w") as f:
        f.write("mgatk2 Run Summary\n")
        f.write("=" * 20 + "\n")

        for key, value in run_metadata.items():
            if key != "parameters":
                f.write(f"{key}: {value}\n")

    logger.info("Wrote run summary to: %s", output_path)


def write_parameters_json(parameters: dict, output_path: Path):
    """Write run parameters to JSON file."""
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        json.dump(parameters, f, indent=2, default=str)

    logger.info("Wrote parameters to: %s", output_path)


def write_reference_fasta(ref_alleles: dict, output_path: Path, chr_name: str, length: int):
    """Write reference sequence in FASTA format"""
    with open(output_path, "w") as f:
        f.write(f">{chr_name}\n")

        # Build sequence from ref_alleles
        sequence = "".join([ref_alleles.get(pos, "N") for pos in range(1, length + 1)])

        # Write sequence in 60-character lines (FASTA standard)
        for i in range(0, len(sequence), 60):
            f.write(sequence[i : i + 60] + "\n")

    logger.info("Wrote reference FASTA to: %s", output_path)


def write_tsv_file(data: list[dict], output_path: Path, delimiter="\t"):
    """Write data to TSV/CSV file"""
    if not data:
        return

    # Get all unique keys
    all_keys: set[str] = set()
    for row in data:
        all_keys.update(row.keys())
    sorted_keys = sorted(all_keys)

    # Handle compression if file ends with .gz
    if output_path.suffix == ".gz":
        f = gzip.open(output_path, "wt")
    else:
        f = open(output_path, "w")

    with f:
        # Write header
        f.write(delimiter.join(sorted_keys) + "\n")

        # Write data
        for row in data:
            values = [str(row.get(key, "")) for key in sorted_keys]
            f.write(delimiter.join(values) + "\n")

    logger.info("Wrote data to: %s", output_path)


def validate_file_format(file_path: Path, expected_format: str):
    """Validate file format based on extension"""
    file_path = Path(file_path)

    format_extensions = {
        "bam": [".bam"],
        "csv": [".csv", ".csv.gz"],
        "tsv": [".tsv", ".tsv.gz", ".txt"],
        "json": [".json"],
        "fasta": [".fasta", ".fa", ".fas"],
        "h5": [".h5", ".hdf5"],
    }

    valid_extensions = format_extensions.get(expected_format.lower(), [])

    if file_path.suffix == ".gz":
        stem_suffix = file_path.stem.split(".")[-1] if "." in file_path.stem else ""
        actual_extension = f".{stem_suffix}.gz"
    else:
        actual_extension = file_path.suffix

    if actual_extension not in valid_extensions:
        raise ValueError(
            f"File {file_path} has extension '{actual_extension}' "
            f"but expected format '{expected_format}' requires one of: {valid_extensions}"
        )

    return True


def detect_delimiter(file_path: Path) -> str:
    """Detect delimiter in CSV/TSV file"""
    import csv

    with open(file_path) as f:
        # Read first few lines
        sample = f.read(1024)
        sniffer = csv.Sniffer()
        return sniffer.sniff(sample).delimiter


def get_compression_type(file_path: Path) -> str:
    """Detect compression type of file"""
    if file_path.suffix == ".gz":
        return "gzip"
    if file_path.suffix in [".bz2"]:
        return "bzip2"
    return "none"
