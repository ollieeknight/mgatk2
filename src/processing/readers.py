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
    InvalidInputError,
    NoBarcodeTagsError,
    NoChrMReadsError,
)
from processing.fragments import (
    alignment_key,
    deduplicate_fragments,
    group_reads_into_fragments,
)

logger = logging.getLogger(__name__)


class BAMReader:
    """Reads and filters BAM files for mtDNA analysis."""

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
        for mito_name in ["chrM", "MT", "M", "chrMT"]:
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

    def collect_reads_by_barcode(self) -> tuple[dict, dict]:
        """Collect reads grouped by barcode with deduplication.

        If barcodes contains only ['bulk'], performs bulk calling
        without barcode filtering.
        """
        reads_by_barcode = defaultdict(list)

        is_bulk_mode = self.barcodes == {"bulk"}

        # The alternative key is tracked only when debug logging is enabled.
        use_frag = self.config.dedup.use_fragment_length
        primary_mode = "alignment_and_fragment_length" if use_frag else "alignment_start"
        secondary_mode = "alignment_start" if use_frag else "alignment_and_fragment_length"
        _debug = logger.isEnabledFor(logging.DEBUG)
        seen_primary: dict[str, set] = defaultdict(set)
        seen_secondary: dict[str, set] = defaultdict(set) if _debug else {}

        total_reads = 0
        filtered_reads = 0
        duplicate_reads_primary = 0
        duplicate_reads_secondary = 0

        try:
            with self._open() as bam:
                for read in tqdm(
                    bam.fetch(self.config.mito_chr),
                    desc="Reading chrM",
                    unit=" read",
                    bar_format="{desc}: {n:,} read [{elapsed}, {rate_fmt}]",
                ):
                    total_reads += 1

                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue

                    # Keep these filters in step with mgatk's filterClipBam.py.
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
                        barcode = "bulk"
                    else:
                        if not read.has_tag(self.config.barcode_tag):
                            continue

                        barcode = read.get_tag(self.config.barcode_tag)

                        if barcode not in self.barcodes:
                            continue

                    if not self.config.dedup.skip:
                        ref_start = read.reference_start
                        is_rev = read.is_reverse
                        tlen = abs(read.template_length) if read.template_length else 0

                        primary_key = alignment_key(ref_start, is_rev, tlen, primary_mode)
                        is_dup = primary_key in seen_primary[barcode]
                        seen_primary[barcode].add(primary_key)
                        if _debug:
                            alt_key = alignment_key(ref_start, is_rev, tlen, secondary_mode)
                            if alt_key in seen_secondary[barcode]:
                                duplicate_reads_secondary += 1
                            seen_secondary[barcode].add(alt_key)

                        if is_dup:
                            duplicate_reads_primary += 1
                            continue

                    if read.query_qualities is not None:
                        qualities = np.array(read.query_qualities, dtype=np.int8)
                    else:
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

        del seen_primary
        seen_secondary.clear()

        stats = {
            "total_reads": total_reads,
            "filtered_reads": filtered_reads,
            "n_barcodes": len(reads_by_barcode),
            "duplicate_reads_primary": duplicate_reads_primary,
            "duplicate_reads_secondary": duplicate_reads_secondary,
        }

        return dict(reads_by_barcode), stats

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
