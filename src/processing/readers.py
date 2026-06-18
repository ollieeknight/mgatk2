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

        # Track seen fragments per barcode for deduplication.
        # Only build the selected set; second set only when DEBUG logging active.
        use_frag = self.config.dedup.use_fragment_length
        _debug = logger.isEnabledFor(logging.DEBUG)
        seen_primary: dict[str, set] = defaultdict(set)
        seen_secondary: dict[str, set] = defaultdict(set) if _debug else {}

        # Statistics tracking
        total_reads = 0
        filtered_reads = 0
        duplicate_reads_primary = 0
        duplicate_reads_secondary = 0

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

                    # NH/NM filters — matches mgatk filterClipBam.py:36-38
                    if self.config.quality.nh_max > 0 and read.has_tag("NH"):
                        if read.get_tag("NH") > self.config.quality.nh_max:
                            continue
                    if self.config.quality.nm_max > 0:
                        nm = None
                        if read.has_tag("NM"):
                            nm = read.get_tag("NM")
                        elif read.has_tag("nM"):
                            nm = read.get_tag("nM")
                        if nm is not None and nm > self.config.quality.nm_max:
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

                    # Deduplication: single set for the selected strategy.
                    # Second set tracked only in DEBUG mode for comparative stats.
                    if not self.config.dedup.skip:
                        ref_start = read.reference_start
                        is_rev = read.is_reverse
                        tlen = abs(read.template_length) if read.template_length else 0

                        if use_frag:
                            primary_key = (ref_start, is_rev, tlen)
                            is_dup = primary_key in seen_primary[barcode]
                            seen_primary[barcode].add(primary_key)
                            if _debug:
                                alt_key = (ref_start, is_rev)
                                if alt_key in seen_secondary[barcode]:
                                    duplicate_reads_secondary += 1
                                seen_secondary[barcode].add(alt_key)
                        else:
                            primary_key = (ref_start, is_rev)
                            is_dup = primary_key in seen_primary[barcode]
                            seen_primary[barcode].add(primary_key)
                            if _debug:
                                alt_key = (ref_start, is_rev, tlen)
                                if alt_key in seen_secondary[barcode]:
                                    duplicate_reads_secondary += 1
                                seen_secondary[barcode].add(alt_key)

                        if is_dup:
                            duplicate_reads_primary += 1
                            continue

                    # Convert to lightweight SimpleRead
                    if read.query_qualities is not None:
                        qualities = np.array(read.query_qualities, dtype=np.int8)
                    else:
                        # Assign maximum quality (60) if qualities are missing
                        qualities = np.full(len(read.query_sequence), 60, dtype=np.int8)

                    simple_read = SimpleRead(
                        reference_start=read.reference_start,
                        is_reverse=read.is_reverse,
                        mapping_quality=read.mapping_quality,
                        query_sequence=read.query_sequence.encode("ascii"),
                        query_qualities=qualities,
                        cigar=read.cigartuples if read.cigartuples else [],
                        is_proper_pair=read.is_proper_pair,
                        is_paired=read.is_paired,
                        template_length=(
                            read.template_length if read.template_length is not None else 0
                        ),
                    )
                    reads_by_barcode[barcode].append(simple_read)
                    filtered_reads += 1

        except Exception as e:
            raise BAMReadError(str(self.bam_path), f"Read error: {e}") from e

        if not self.config.dedup.skip:
            logger.info(
                "%d duplicate reads removed (%.1f%%)",
                duplicate_reads_primary,
                duplicate_reads_primary / total_reads * 100 if total_reads else 0,
            )
            if _debug and duplicate_reads_secondary:
                logger.debug(
                    "Alt dedup strategy would have removed %d duplicates",
                    duplicate_reads_secondary,
                )

        logger.info(
            f"Kept {filtered_reads:,} reads from {len(reads_by_barcode):,} barcodes at an average of {filtered_reads / len(reads_by_barcode):.0f} reads/cell"
            if len(reads_by_barcode) > 1
            else f"Kept {filtered_reads:,} reads from {len(reads_by_barcode):,} barcodes"
        )

        # Clear deduplication sets to free memory
        del seen_primary
        seen_secondary.clear()

        # Return reads directly in memory (no staging)
        stats = {
            "total_reads": total_reads,
            "filtered_reads": filtered_reads,
            "n_barcodes": len(reads_by_barcode),
            "duplicate_reads_primary": duplicate_reads_primary,
            "duplicate_reads_secondary": duplicate_reads_secondary,
        }

        return dict(reads_by_barcode), stats
