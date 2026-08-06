"""Build a tiny anonymised 10x Multi (scRNA/CITE-seq) fixture for integration tests.

Derives a small, real-read-shaped fixture from a source CellRanger Multi `outs/`
tree, keeping only a couple of barcodes' reads (matching the scale of the ATAC
fixture in create_integration_fixture.py) and stripping the sample ID / real
cluster paths / username from the BAM header and CSVs.
"""

import argparse
import re
from pathlib import Path

import pysam

SAMPLE_NAME = "test_sample"
KEEP_BARCODES = {"GCATTAGGTCGAGATG-1", "CGGAACCCATCCTTGC-1"}
REDACTED_PREFIX = "/data/test"
CEPHFS_PATH_RE = re.compile(r"/data/cephfs-\d[^\s,]*")


def _redact_path(match: re.Match) -> str:
    segments = match.group(0).split("/")
    return REDACTED_PREFIX + "/" + "/".join(segments[-2:])


def _anonymise(text: str, real_sample_id: str) -> str:
    text = text.replace(real_sample_id, SAMPLE_NAME)
    return CEPHFS_PATH_RE.sub(_redact_path, text)


def build_fixture(source_outs: Path, output_dir: Path, real_sample_id: str) -> Path:
    src_count_dir = source_outs / "per_sample_outs" / real_sample_id / "count"
    dst_outs = output_dir / "outs"
    dst_count_dir = dst_outs / "per_sample_outs" / SAMPLE_NAME / "count"
    dst_count_dir.mkdir(parents=True, exist_ok=True)

    _write_bam(
        src_count_dir / "sample_alignments.bam",
        dst_count_dir / "sample_alignments.bam",
        real_sample_id,
    )
    _write_barcodes(
        src_count_dir / "sample_filtered_barcodes.csv",
        dst_count_dir / "sample_filtered_barcodes.csv",
    )
    _write_text_file(source_outs / "config.csv", dst_outs / "config.csv", real_sample_id)
    _write_text_file(
        source_outs / "per_sample_outs" / real_sample_id / "metrics_summary.csv",
        dst_outs / "per_sample_outs" / SAMPLE_NAME / "metrics_summary.csv",
        real_sample_id,
    )
    return dst_outs


def _write_bam(src_bam: Path, dst_bam: Path, real_sample_id: str) -> None:
    with pysam.AlignmentFile(str(src_bam), "rb") as src:
        header = src.header.to_dict()
        for rg in header.get("RG", []):
            for key in ("ID", "PU", "SM"):
                if key in rg:
                    rg[key] = _anonymise(rg[key], real_sample_id)
        for pg in header.get("PG", []):
            if "CL" in pg:
                pg["CL"] = _anonymise(pg["CL"], real_sample_id)

        with pysam.AlignmentFile(str(dst_bam), "wb", header=header) as dst:
            for read in src.fetch(until_eof=True):
                if not read.has_tag("CB") or read.get_tag("CB") not in KEEP_BARCODES:
                    continue
                if read.has_tag("RG"):
                    read.set_tag("RG", _anonymise(read.get_tag("RG"), real_sample_id))
                dst.write(read)

    pysam.index(str(dst_bam))


def _write_barcodes(src_csv: Path, dst_csv: Path) -> None:
    with open(src_csv) as f_in, open(dst_csv, "w") as f_out:
        for line in f_in:
            barcode = line.strip().rsplit(",", 1)[-1]
            if barcode in KEEP_BARCODES:
                f_out.write(line)


def _write_text_file(src_path: Path, dst_path: Path, real_sample_id: str) -> None:
    dst_path.parent.mkdir(parents=True, exist_ok=True)
    dst_path.write_text(_anonymise(src_path.read_text(), real_sample_id))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("source_outs", type=Path, help="path to real CellRanger Multi outs/")
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("real_sample_id")
    args = parser.parse_args()
    build_fixture(args.source_outs, args.output_dir, args.real_sample_id)
