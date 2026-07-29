"""Create a tiny deterministic 10x-style BAM for Makefile integration tests."""

import argparse
from pathlib import Path

import pysam


def create_fixture(output_dir: Path) -> Path:
    outs = output_dir / "outs"
    outs.mkdir(parents=True, exist_ok=True)
    bam_path = outs / "possorted_bam.bam"
    header = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chrM", "LN": 16569}]}

    with pysam.AlignmentFile(bam_path, "wb", header=header) as bam:
        for index in range(24):
            read = pysam.AlignedSegment(bam.header)
            read.query_name = f"read-{index:02d}"
            read.query_sequence = "ACGT" * 10
            read.flag = 0
            read.reference_id = 0
            read.reference_start = 100 + index * 10
            read.mapping_quality = 60
            read.cigar = ((0, 40),)
            read.query_qualities = pysam.qualitystring_to_array("I" * 40)
            read.set_tag("CB", f"cell-{index % 2 + 1}")
            bam.write(read)

    pysam.index(str(bam_path))
    (outs / "filtered_peak_bc_matrix").mkdir(exist_ok=True)
    (outs / "filtered_peak_bc_matrix" / "barcodes.tsv").write_text("cell-1\ncell-2\n")
    return outs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", type=Path)
    create_fixture(parser.parse_args().output_dir)
