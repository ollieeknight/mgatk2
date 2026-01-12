"""Utility functions for mgatk2."""

import csv
import logging
import os

import pysam

from core.exceptions import InvalidInputError

logger = logging.getLogger(__name__)


def load_barcodes(barcode_file: str) -> list[str] | None:
    """Load barcodes from file.

    Supports:
    - singlecell.csv (10x scATAC format with is__cell_barcode column)
    - barcodes.tsv or .csv (simple list, one per line)
    """
    if not barcode_file:
        return None

    # Check if this is a singlecell.csv file
    if barcode_file.endswith("singlecell.csv"):
        logger.info("Detected singlecell.csv format")
        barcodes, _ = load_singlecell_csv(barcode_file)
        return barcodes

    # Load as simple barcode list (TSV or CSV)
    try:
        with open(barcode_file) as f:
            barcodes = [line.strip() for line in f if line.strip()]
        logger.info(f"Loaded {len(barcodes)} barcodes from file")
        return barcodes
    except Exception as e:
        logger.error(f"Error loading barcodes: {e}")
        raise


def load_singlecell_csv(csv_file: str) -> tuple[list[str] | None, dict[str, list] | None]:
    """Load barcodes and metadata from cellranger-atac singlecell.csv file."""
    if csv_file is None:
        return None, None

    try:
        barcodes = []
        metadata: dict[str, list] = {}
        total_rows = 0

        with open(csv_file) as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames

            if headers is None:
                raise InvalidInputError("CSV file has no headers")

            if "is__cell_barcode" not in headers:
                raise InvalidInputError("singlecell.csv missing 'is__cell_barcode' column")

            for header in headers:
                metadata[header] = []

            for row in reader:
                total_rows += 1

                if row.get("is__cell_barcode") != "1":
                    continue

                barcode = row["barcode"]
                barcodes.append(barcode)

                for header in headers:
                    value = row[header]
                    if header not in ["barcode", "excluded_reason"]:
                        try:
                            value = int(value) if value else 0
                        except ValueError:
                            try:
                                value = float(value) if value else 0.0
                            except ValueError:
                                pass  # Keep as string
                    metadata[header].append(value)

        if len(barcodes) == 0:
            raise InvalidInputError(f"No cells found with is__cell_barcode == 1 in {csv_file}")

        logger.info(f"Performing analysis on {len(barcodes)} barcodes from singlecell.csv")

        return barcodes, metadata

    except FileNotFoundError as e:
        raise InvalidInputError(f"singlecell.csv file not found: {csv_file}") from e
    except Exception as e:
        logger.error("Error loading singlecell.csv: %s", e)
        raise


def validate_bam_file(bam_path: str) -> None:
    """Validate input BAM file exists and is properly formatted."""
    if not os.path.exists(bam_path):
        raise InvalidInputError(f'BAM file not found: "{bam_path}"')

    if not bam_path.endswith(".bam"):
        raise InvalidInputError(f"Input file must have .bam extension: {bam_path}")

    if not os.path.exists(bam_path + ".bai"):
        logger.info("Creating BAM index for: %s", bam_path)
        try:
            pysam.index(bam_path)
        except Exception as e:
            raise InvalidInputError(f"Failed to create BAM index: {e}") from e

        if not os.path.exists(bam_path + ".bai"):
            raise InvalidInputError(f"Cannot create index for BAM: {bam_path}")


def validate_barcode_file(barcode_file: str) -> None:
    """Validate barcode file exists and is readable."""
    if not barcode_file:
        return

    if not os.path.exists(barcode_file):
        raise InvalidInputError(f'Barcode file not found: "{barcode_file}"')

    try:
        with open(barcode_file) as f:
            if not f.readline().strip():
                raise InvalidInputError(f"Barcode file is empty: {barcode_file}")
    except Exception as e:
        raise InvalidInputError(f"Cannot read barcode file: {e}") from e
