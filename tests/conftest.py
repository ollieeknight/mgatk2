"""Synthetic alignment fixtures for paired mitochondrial tests."""

from pathlib import Path

import pysam
import pytest


def _write_alignment(path: Path, reference: Path, reads: list[dict]) -> Path:
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrM", "LN": 40}]}
    mode = "wc" if path.suffix == ".cram" else "wb"
    kwargs = {"reference_filename": str(reference)} if mode == "wc" else {}
    with pysam.AlignmentFile(path, mode, header=header, **kwargs) as output:
        for specification in sorted(reads, key=lambda item: item["start"]):
            read = pysam.AlignedSegment(output.header)
            read.query_name = specification["name"]
            read.query_sequence = specification.get("sequence", "A" * 15)
            read.flag = specification.get("flag", 0)
            read.reference_id = 0
            read.reference_start = specification["start"]
            read.mapping_quality = specification.get("mapq", 60)
            read.cigar = specification.get("cigar", ((0, len(read.query_sequence)),))
            qualities = specification.get("qualities", "D" * len(read.query_sequence))
            read.query_qualities = (
                None if qualities is None else pysam.qualitystring_to_array(qualities)
            )
            for tag, value in specification.get("tags", {}).items():
                read.set_tag(tag, value)
            read.next_reference_id = specification.get("next_reference_id", -1)
            read.next_reference_start = specification.get("next_start", -1)
            read.template_length = specification.get("template_length", 0)
            output.write(read)
    pysam.index(str(path))
    return path


@pytest.fixture
def alignment_factory():
    return _write_alignment


@pytest.fixture
def barcoded_bam(tmp_path, alignment_factory):
    """Two cells on a 40bp chrM. cell-1 has a duplicate pair; cell-2 is reverse."""
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    pysam.faidx(str(reference))

    reads = [
        {"name": "r1", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}},
        # identical alignment to r1: dropped by deduplication
        {"name": "r2", "start": 0, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}},
        {"name": "r3", "start": 4, "sequence": "ACGT" * 5, "tags": {"CB": "cell-1"}},
        {"name": "r4", "start": 0, "sequence": "C" * 20, "flag": 16, "tags": {"CB": "cell-2"}},
        {"name": "r5", "start": 0, "sequence": "G" * 20, "tags": {"CB": "not-a-cell"}},
    ]
    return alignment_factory(tmp_path / "cells.bam", reference, reads)


@pytest.fixture
def paired_files(tmp_path, alignment_factory):
    reference = tmp_path / "reference.fa"
    reference.write_text(">chrM\n" + "A" * 40 + "\n")
    pysam.faidx(str(reference))

    tumor_reads = []
    for index, start in enumerate((3, 4, 5)):
        sequence = list("A" * 15)
        sequence[10 - start] = "C"
        tumor_reads.append({"name": f"tumor{index}", "start": start, "sequence": "".join(sequence)})
    normal_reads = [
        {"name": f"normal{index}", "start": start} for index, start in enumerate((3, 4, 5))
    ]
    tumor_bam = alignment_factory(tmp_path / "tumor.bam", reference, tumor_reads)
    normal_bam = alignment_factory(tmp_path / "normal.bam", reference, normal_reads)
    tumor_cram = alignment_factory(tmp_path / "tumor.cram", reference, tumor_reads)
    normal_cram = alignment_factory(tmp_path / "normal.cram", reference, normal_reads)
    return {
        "reference": reference,
        "tumor_bam": tumor_bam,
        "normal_bam": normal_bam,
        "tumor_cram": tumor_cram,
        "normal_cram": normal_cram,
    }
