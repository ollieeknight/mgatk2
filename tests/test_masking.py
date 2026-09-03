"""Coverage for the hardmask-fasta command and its masking helpers."""

import gzip
import importlib

import pytest
from click.testing import CliRunner

from utils.masking import (
    detect_fasta_chr_prefix,
    get_blacklist_path,
    load_blacklist_regions,
    mask_fasta,
    normalise_bed_chromosomes,
    normalise_genome_name,
)


@pytest.mark.parametrize(
    ("supplied", "expected"),
    [
        ("hg38", "hg38"),
        ("GRCh38", "hg38"),
        ("grch37", "hg19"),
        ("GRCm38", "mm10"),
        ("mm9", "mm9"),
    ],
)
def test_genome_aliases_normalise(supplied, expected):
    assert normalise_genome_name(supplied) == expected


def test_unknown_genome_is_rejected():
    with pytest.raises(ValueError, match="not supported"):
        normalise_genome_name("rn6")


@pytest.mark.parametrize("build", ["hg38", "hg19", "mm10", "mm9"])
def test_every_bundled_blacklist_is_present_and_parses(build):
    regions = load_blacklist_regions(get_blacklist_path(build))

    assert regions
    assert all(start < end for _chrom, start, end in regions)


def test_masking_uses_bed_half_open_coordinates(tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">chr1\n" + "ACGT" * 10 + "\n")
    destination = tmp_path / "out.fa"

    stats = mask_fasta(source, destination, [("chr1", 4, 12)])

    sequence = "".join(
        line.strip() for line in destination.read_text().splitlines() if not line.startswith(">")
    )
    assert sequence == "ACGT" + "N" * 8 + "ACGT" * 7
    assert stats == {"chromosomes_processed": 1, "regions_masked": 1, "bases_masked": 8}


def test_masking_leaves_unlisted_contigs_untouched(tmp_path):
    source = tmp_path / "in.fa"
    source.write_text(">chr1\nAAAA\n>chrM\nCCCC\n")
    destination = tmp_path / "out.fa"

    mask_fasta(source, destination, [("chr1", 0, 4)])

    assert destination.read_text() == ">chr1\nNNNN\n>chrM\nCCCC\n"


def test_gzipped_input_and_output_round_trip(tmp_path):
    source = tmp_path / "in.fa.gz"
    with gzip.open(source, "wt") as handle:
        handle.write(">chr1\n" + "A" * 20 + "\n")
    destination = tmp_path / "out.fa.gz"

    mask_fasta(source, destination, [("chr1", 0, 5)])

    with gzip.open(destination, "rt") as handle:
        assert handle.read() == ">chr1\n" + "N" * 5 + "A" * 15 + "\n"


def test_chr_prefix_detection_and_bed_renaming(tmp_path):
    prefixed = tmp_path / "prefixed.fa"
    prefixed.write_text(">chr1\nAAAA\n>chr2\nAAAA\n")
    plain = tmp_path / "plain.fa"
    plain.write_text(">1\nAAAA\n>2\nAAAA\n")

    assert detect_fasta_chr_prefix(prefixed) is True
    assert detect_fasta_chr_prefix(plain) is False
    assert normalise_bed_chromosomes([("1", 0, 4)], True) == [("chr1", 0, 4)]
    assert normalise_bed_chromosomes([("chr1", 0, 4)], False) == [("1", 0, 4)]


def test_hardmask_fasta_command_runs_end_to_end(tmp_path):
    command_module = importlib.import_module("cli.commands.mask")
    source = tmp_path / "genome.fa"
    source.write_text(">chr1\n" + "A" * 100 + "\n>chrM\n" + "C" * 50 + "\n")
    destination = tmp_path / "masked.fa"

    result = CliRunner().invoke(
        command_module.hardmask_fasta,
        ["--input-fasta", str(source), "--output-fasta", str(destination), "--genome", "GRCh38"],
    )

    assert result.exit_code == 0, result.output
    assert destination.exists()
    # Bundled BEDs are nuclear-side only, so chrM must survive untouched.
    assert "C" * 50 in destination.read_text().replace("\n", "")


def test_hardmask_fasta_rejects_an_unknown_genome(tmp_path):
    command_module = importlib.import_module("cli.commands.mask")
    source = tmp_path / "genome.fa"
    source.write_text(">chr1\nAAAA\n")

    result = CliRunner().invoke(
        command_module.hardmask_fasta,
        [
            "--input-fasta",
            str(source),
            "--output-fasta",
            str(tmp_path / "out.fa"),
            "--genome",
            "rn6",
        ],
    )

    assert result.exit_code == 1
