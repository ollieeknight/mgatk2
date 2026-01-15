"""BAM reading and deduplication for mgatk2."""

import logging
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam
from tqdm import tqdm

from core.config import PipelineConfig, SimpleRead
from core.exceptions import (
    BAMFormatError,
    BAMReadError,
    NoBarcodeTagsError,
    NoChrMReadsError,
)

logger = logging.getLogger(__name__)


class BAMReader:
    """Reads and filters BAM files for mtDNA analysis."""

    def __init__(self, bam_path: str, config: PipelineConfig, barcodes: set[str]):
        self.bam_path = Path(bam_path)
        self.config = config
        self.barcodes = barcodes

        if not self.bam_path.exists():
            raise BAMReadError(str(bam_path), "File does not exist")

        self._validate_bam_file()

    def _validate_bam_file(self):
        try:
            bam = pysam.AlignmentFile(str(self.bam_path), "rb")
        except Exception as e:
            raise BAMFormatError(str(self.bam_path), f"Cannot open: {e}") from e

        # Check for mitochondrial chromosome
        available = list(bam.references)
        for mito_name in ["chrM", "MT", "M", "chrMT"]:
            if mito_name in available:
                if self.config.mito_chr != mito_name:
                    logger.info(f"Using mitochondrial chromosome: {mito_name}")
                    self.config.mito_chr = mito_name
                break
        else:
            bam.close()
            raise NoChrMReadsError(str(self.bam_path), available)

        # Quick check for barcode tags
        for i, read in enumerate(bam.fetch(self.config.mito_chr)):
            if read.has_tag(self.config.barcode_tag):
                break
            if i >= 1000:
                bam.close()
                raise NoBarcodeTagsError(str(self.bam_path), self.config.barcode_tag, i)

        bam.close()

    def collect_reads_by_barcode(self) -> tuple[dict, dict]:
        """Collect reads grouped by barcode with deduplication.

        If barcodes contains only ['bulk'], performs bulk calling
        without barcode filtering.
        """
        reads_by_barcode = defaultdict(list)

        # Check if we're doing bulk calling
        is_bulk_mode = self.barcodes == {"bulk"}

        # Track seen fragments per barcode for deduplication
        seen_fragments_with_length: dict[str, set] = defaultdict(set)
        seen_fragments_position_only: dict[str, set] = defaultdict(set)

        # Statistics tracking
        total_reads = 0
        filtered_reads = 0
        duplicate_reads_with_length = 0
        duplicate_reads_position_only = 0

        try:
            with pysam.AlignmentFile(str(self.bam_path), "rb") as bam:
                # Stream through chrM reads once
                for read in tqdm(
                    bam.fetch(self.config.mito_chr),
                    desc="Reading chrM",
                    unit=" read",
                    bar_format="{desc}: {n:,} read [{elapsed}, {rate_fmt}]",
                ):
                    total_reads += 1

                    # Skip unmapped/secondary/supplementary
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    if is_bulk_mode:
                        # For bulk calling, assign all reads to 'bulk' barcode
                        barcode = "bulk"
                    else:
                        # Check if read has barcode tag
                        if not read.has_tag(self.config.barcode_tag):
                            continue

                        barcode = read.get_tag(self.config.barcode_tag)

                        # Only keep reads from barcodes in our list
                        if barcode not in self.barcodes:
                            continue

                    # Deduplication: mark duplicates based on alignment
                    # position and strand
                    # - alignment_start: (ref_start, is_reverse) - Picard
                    # - alignment_and_fragment_length: adds template_length
                    # Both in-memory; original mgatk used Picard (ext tool)
                    if not self.config.dedup.skip:
                        ref_start = read.reference_start
                        is_rev = read.is_reverse
                        fragment_length_key = (
                            ref_start,
                            is_rev,
                            abs(read.template_length),
                        )
                        position_only_key = (ref_start, is_rev)

                        # Check both methods for statistics
                        is_fragment_length_dup = (
                            fragment_length_key in seen_fragments_with_length[barcode]
                        )
                        is_position_only_dup = (
                            position_only_key in seen_fragments_position_only[barcode]
                        )

                        # Add to both tracking sets
                        seen_fragments_with_length[barcode].add(fragment_length_key)
                        seen_fragments_position_only[barcode].add(position_only_key)

                        # Count duplicates for both methods
                        if is_fragment_length_dup:
                            duplicate_reads_with_length += 1
                        if is_position_only_dup:
                            duplicate_reads_position_only += 1

                        # Skip read only if duplicate by the SELECTED method
                        if self.config.dedup.use_fragment_length and is_fragment_length_dup:
                            continue  # Skip fragment-length duplicate
                        if not self.config.dedup.use_fragment_length and is_position_only_dup:
                            continue  # Skip position-only duplicate

                    # Convert to lightweight SimpleRead
                    simple_read = SimpleRead(
                        reference_start=read.reference_start,
                        is_reverse=read.is_reverse,
                        mapping_quality=read.mapping_quality,
                        query_sequence=read.query_sequence.encode("ascii"),  # Store as bytes
                        query_qualities=np.array(read.query_qualities, dtype=np.int8),
                        cigar=read.cigartuples if read.cigartuples else [],
                        is_proper_pair=read.is_proper_pair,
                        is_paired=read.is_paired,
                        template_length=read.template_length,
                    )
                    reads_by_barcode[barcode].append(simple_read)
                    filtered_reads += 1

        except Exception as e:
            raise BAMReadError(str(self.bam_path), f"Read error: {e}") from e

        if not self.config.dedup.skip:
            duplicates_removed = (
                duplicate_reads_with_length
                if self.config.dedup.use_fragment_length
                else duplicate_reads_position_only
            )
        logger.info(
            "%s duplicate reads removed (%.1f%%)",
            f"{duplicates_removed:,}",
            duplicates_removed / total_reads * 100,
        )

        logger.info(
            (
                f"Kept {filtered_reads:,} reads from "
                f"{len(reads_by_barcode):,} barcodes at an average of "
                f"{filtered_reads / len(reads_by_barcode):.0f} reads/cell"
            )
            if len(reads_by_barcode) > 1
            else (f"Kept {filtered_reads:,} reads from " f"{len(reads_by_barcode):,} barcodes")
        )

        # Clear deduplication sets to free memory
        del seen_fragments_with_length
        del seen_fragments_position_only

        # Return reads directly in memory (no staging)
        stats = {
            "total_reads": total_reads,
            "filtered_reads": filtered_reads,
            "n_barcodes": len(reads_by_barcode),
            "duplicate_reads_with_length": duplicate_reads_with_length,
            "duplicate_reads_position_only": duplicate_reads_position_only,
        }

        return dict(reads_by_barcode), stats
