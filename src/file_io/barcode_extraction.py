"""Extract barcodes from BAM files"""

from __future__ import annotations

import logging

import pysam

logger = logging.getLogger(__name__)


def extract_barcodes_from_bam(
    bam_path: str, barcode_tag: str = "CB", mito_chr: str = "chrM", min_reads: int = 10
) -> list[str]:
    """Extract unique barcodes from a BAM file."""
    logger.info("Extracting barcodes from BAM file...")
    logger.info("  Looking for tag '%s' on chromosome '%s'", barcode_tag, mito_chr)

    barcode_counts: dict[str, int] = {}

    try:
        bam = pysam.AlignmentFile(bam_path, "rb")

        # Scan mitochondrial reads
        for read in bam.fetch(mito_chr):
            if read.is_unmapped or read.is_duplicate:
                continue

            # Get barcode from tag
            if read.has_tag(barcode_tag):
                barcode = str(read.get_tag(barcode_tag))
                barcode_counts[barcode] = barcode_counts.get(barcode, 0) + 1

        bam.close()

    except Exception as e:
        logger.error("Failed to extract barcodes: %s", e)
        raise

    barcodes = [bc for bc, count in barcode_counts.items() if count >= min_reads]
    barcodes.sort()

    logger.info("  Found %d total barcodes", len(barcode_counts))
    logger.info("  Retained %d barcodes with >= %d reads", len(barcodes), min_reads)

    return barcodes
