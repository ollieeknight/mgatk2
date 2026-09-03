"""Utility functions for mgatk2."""

import csv
import logging
import os

import pysam

from core.exceptions import InvalidInputError

logger = logging.getLogger(__name__)

# CellRanger has shipped all three spellings; treat them as one schema.
CELL_FLAG_COLUMNS = ("is__cell_barcode", "is_cell_barcode", "is_cell")
TRUE_VALUES = frozenset({"1", "1.0", "True", "true", "TRUE"})


def load_singlecell_csv(
    csv_file: str,
) -> tuple[list[str] | None, dict[str, list] | None]:
    """Load barcodes and metadata from cellranger-atac singlecell.csv file."""
    if csv_file is None:
        return None, None

    try:
        barcodes = []
        metadata: dict[str, list] = {}

        with open(csv_file) as f:
            reader = csv.DictReader(f)
            headers = reader.fieldnames

            if headers is None:
                raise InvalidInputError("CSV file has no headers")

            cell_column = next(
                (name for name in CELL_FLAG_COLUMNS if name in headers),
                None,
            )
            if cell_column is None:
                raise InvalidInputError(
                    f"singlecell.csv missing a cell-flag column ({', '.join(CELL_FLAG_COLUMNS)})"
                )

            for header in headers:
                metadata[header] = []

            for row in reader:
                if row.get(cell_column, "0") not in TRUE_VALUES:
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
            raise InvalidInputError(f"No cells found with {cell_column} set in {csv_file}")

        return barcodes, metadata

    except FileNotFoundError as e:
        raise InvalidInputError(f"singlecell.csv file not found: {csv_file}") from e
    except Exception as e:
        logger.error("Error loading singlecell.csv: %s", e)
        raise


def load_barcode_csv(
    csv_file: str,
) -> tuple[list[str] | None, dict[str, list] | None]:
    """Load barcodes from a .csv file, detecting ATAC vs 10x Multi schema."""
    if csv_file is None:
        return None, None

    with open(csv_file) as f:
        header = f.readline()

    if not header.strip():
        raise InvalidInputError(f"Barcode file is empty: {csv_file}")

    fields = next(csv.reader([header]))
    if any(col in fields for col in CELL_FLAG_COLUMNS):
        return load_singlecell_csv(csv_file)

    # 10x Multi sample_filtered_barcodes.csv: headerless "<reference>,<barcode>" rows
    barcodes = []
    with open(csv_file) as f:
        for row in csv.reader(f):
            if row:
                barcodes.append(row[-1])

    if not barcodes:
        raise InvalidInputError(f"No barcodes found in {csv_file}")

    return barcodes, None


def validate_bam_file(bam_path: str) -> None:
    """Validate input BAM file exists and is properly formatted."""
    if not os.path.exists(bam_path):
        raise InvalidInputError(f'BAM file not found: "{bam_path}"')

    if not bam_path.endswith(".bam"):
        raise InvalidInputError(f"Input file must have .bam extension: {bam_path}")

    ensure_alignment_index(bam_path)


def has_alignment_index(alignment_path: str) -> bool:
    """True when any samtools index flavour is present beside the alignment."""
    stem = os.path.splitext(alignment_path)[0]
    return any(
        os.path.exists(candidate)
        for candidate in (
            alignment_path + ".bai",
            alignment_path + ".csi",
            alignment_path + ".crai",
            stem + ".bai",
            stem + ".csi",
            stem + ".crai",
        )
    )


def ensure_alignment_index(alignment_path: str) -> None:
    """Create an index beside the alignment when none of the flavours exist."""
    if has_alignment_index(alignment_path):
        return
    logger.info("Creating alignment index for: %s", alignment_path)
    try:
        pysam.index(alignment_path)
    except Exception as e:
        raise InvalidInputError(f"Failed to create alignment index: {e}") from e
    if not has_alignment_index(alignment_path):
        raise InvalidInputError(f"Cannot create index for: {alignment_path}")


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
