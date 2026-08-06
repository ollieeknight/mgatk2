"""BAM validation and bulk fragment collection for mgatk2.

Single-cell counting does not go through this module: see processing/pileup.py.
"""

import logging
from pathlib import Path

import numpy as np
import pysam

from core.config import PipelineConfig, SimpleRead
from core.exceptions import (
    BAMFormatError,
    BAMReadError,
    InvalidInputError,
    NoBarcodeTagsError,
    NoChrMReadsError,
)
from processing.fragments import (
    deduplicate_fragments,
    group_reads_into_fragments,
)

logger = logging.getLogger(__name__)


class BAMReader:
    """Validates a BAM/CRAM and collects reads for the paired/bulk fragment path."""

    def __init__(
        self,
        bam_path: str,
        config: PipelineConfig,
        barcodes: set[str] | None = None,
        reference_filename: str | None = None,
    ):
        self.bam_path = Path(bam_path)
        self.config = config
        self.barcodes = barcodes or {"bulk"}
        self.reference_filename = reference_filename
        self._is_cram = self.bam_path.suffix.lower() == ".cram"

        if not self.bam_path.exists():
            raise BAMReadError(str(bam_path), "File does not exist")

        self._validate_bam_file()

    def _open(self):
        if self._is_cram:
            if not self.reference_filename:
                raise InvalidInputError("A reference FASTA is required to decode CRAM input")
            return pysam.AlignmentFile(
                str(self.bam_path), "rc", reference_filename=self.reference_filename
            )
        return pysam.AlignmentFile(str(self.bam_path), "rb")

    def _validate_bam_file(self):
        try:
            bam = self._open()
        except Exception as e:
            raise BAMFormatError(str(self.bam_path), f"Cannot open: {e}") from e

        available = list(bam.references)
        # Requested name wins; the aliases are only a fallback.
        for mito_name in [self.config.mito_chr, "chrM", "MT", "M", "chrMT"]:
            if mito_name in available:
                if self.config.mito_chr != mito_name:
                    logger.info(f"Using mitochondrial chromosome: {mito_name}")
                    self.config.mito_chr = mito_name
                break
        else:
            bam.close()
            raise NoChrMReadsError(str(self.bam_path), available)

        if self.barcodes != {"bulk"}:
            reads_checked = 0
            for i, read in enumerate(bam.fetch(self.config.mito_chr)):
                reads_checked = i + 1
                if read.has_tag(self.config.barcode_tag):
                    break
            else:
                bam.close()
                raise NoBarcodeTagsError(str(self.bam_path), self.config.barcode_tag, reads_checked)

        bam.close()

    def collect_bulk_reads(self, deduplication: str) -> tuple[list, dict]:
        """Collect a paired-analysis sample and return fragments plus structured QC."""
        stats = {
            "total_reads": 0,
            "primary_reads": 0,
            "low_mapq_reads": 0,
            "qc_failed_reads": 0,
            "preexisting_duplicate_reads": 0,
            "unmapped_reads": 0,
            "secondary_reads": 0,
            "supplementary_reads": 0,
            "missing_sequence_reads": 0,
            "missing_quality_reads": 0,
            "paired_reads": 0,
            "orphan_reads": 0,
            "improper_pair_reads": 0,
            "clipped_reads": 0,
        }
        reads: list[SimpleRead] = []
        try:
            with self._open() as alignment:
                if not alignment.has_index():
                    raise InvalidInputError(f"Alignment index not found for {self.bam_path}")
                reference_length = alignment.get_reference_length(self.config.mito_chr)
                for read in alignment.fetch(self.config.mito_chr):
                    stats["total_reads"] += 1
                    if read.is_unmapped:
                        stats["unmapped_reads"] += 1
                        continue
                    if read.is_secondary:
                        stats["secondary_reads"] += 1
                        continue
                    if read.is_supplementary:
                        stats["supplementary_reads"] += 1
                        continue
                    stats["primary_reads"] += 1
                    if read.is_qcfail:
                        stats["qc_failed_reads"] += 1
                        continue
                    if read.is_duplicate:
                        stats["preexisting_duplicate_reads"] += 1
                        continue
                    if read.mapping_quality < self.config.quality.min_mapq:
                        stats["low_mapq_reads"] += 1
                        continue
                    if read.query_sequence is None:
                        stats["missing_sequence_reads"] += 1
                        continue
                    if read.query_qualities is None:
                        stats["missing_quality_reads"] += 1
                        continue

                    stats["paired_reads"] += int(read.is_paired)
                    stats["orphan_reads"] += int(not read.is_paired or read.mate_is_unmapped)
                    stats["improper_pair_reads"] += int(read.is_paired and not read.is_proper_pair)
                    cigar = read.cigartuples or []
                    stats["clipped_reads"] += int(any(op in {4, 5} for op, _length in cigar))
                    reads.append(
                        SimpleRead(
                            reference_start=read.reference_start,
                            reference_end=read.reference_end or read.reference_start,
                            is_reverse=read.is_reverse,
                            mapping_quality=read.mapping_quality,
                            query_sequence=read.query_sequence.encode("ascii"),
                            query_qualities=np.array(read.query_qualities, dtype=np.int16),
                            cigar=cigar,
                            is_proper_pair=read.is_proper_pair,
                            is_paired=read.is_paired,
                            template_length=read.template_length or 0,
                            query_name=read.query_name,
                            is_read1=read.is_read1,
                            is_read2=read.is_read2,
                            is_qcfail=read.is_qcfail,
                            is_duplicate=read.is_duplicate,
                            read_group=read.get_tag("RG") if read.has_tag("RG") else None,
                        )
                    )
        except InvalidInputError:
            raise
        except Exception as e:
            raise BAMReadError(str(self.bam_path), f"Read error: {e}") from e

        fragments, collisions = group_reads_into_fragments(reads)
        fragments, duplicate_stats = deduplicate_fragments(fragments, deduplication)
        stats.update(duplicate_stats)
        stats["query_name_collisions"] = collisions
        stats["retained_reads"] = sum(len(fragment.reads) for fragment in fragments)
        stats["retained_fragments"] = len(fragments)
        stats["reference_length"] = reference_length
        if not fragments:
            raise NoChrMReadsError(str(self.bam_path), [self.config.mito_chr])
        return fragments, stats
