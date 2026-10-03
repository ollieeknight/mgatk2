"""Extract barcodes from BAM files"""

import logging
from collections import Counter

import pysam

logger = logging.getLogger(__name__)


def extract_barcodes_from_bam(
    bam_path: str, barcode_tag: str = "CB", mito_chr: str = "chrM", min_reads: int = 10
) -> list[str]:
    """Sorted barcodes carrying at least min_reads non-duplicate reads on mito_chr."""
    logger.info("Extracting '%s' barcodes from %s reads...", barcode_tag, mito_chr)
    counts: Counter[str] = Counter()
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        for read in bam.fetch(mito_chr):
            if not read.is_unmapped and not read.is_duplicate and read.has_tag(barcode_tag):
                counts[str(read.get_tag(barcode_tag))] += 1

    barcodes = sorted(barcode for barcode, n in counts.items() if n >= min_reads)
    logger.info(
        "  Retained %d of %d barcodes with >= %d reads", len(barcodes), len(counts), min_reads
    )
    return barcodes
