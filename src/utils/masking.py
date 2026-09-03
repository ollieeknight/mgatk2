"""Genome-build naming and hard-masking for the ``hardmask-fasta`` command."""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

from data.blacklists import bundled_bed_path

logger = logging.getLogger(__name__)

# Aliases are matched case-insensitively, so only lowercase keys are stored.
GENOME_ALIASES = {
    "hg38": "hg38",
    "grch38": "hg38",
    "hg19": "hg19",
    "grch37": "hg19",
    "mm10": "mm10",
    "grcm38": "mm10",
    "mm9": "mm9",
    "grcm37": "mm9",
}


def normalise_genome_name(genome: str) -> str:
    """Normalise a genome build name to its canonical bundled-BED form."""
    normalised = GENOME_ALIASES.get(genome.lower())
    if normalised is None:
        raise ValueError(
            f"Genome '{genome}' not supported. Supported genomes: "
            f"{', '.join(sorted(GENOME_ALIASES))}"
        )
    return normalised


def get_blacklist_path(genome: str) -> Path:
    """Path to the bundled NUMT blacklist BED for a genome build."""
    return bundled_bed_path(normalise_genome_name(genome))


def _open_text(path: Path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_blacklist_regions(blacklist_path: Path) -> list[tuple[str, int, int]]:
    """Load ``(chrom, start, end)`` NUMT regions from a BED file."""
    regions = []
    with open(blacklist_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            parts = line.split("\t")
            if len(parts) >= 3:
                regions.append((parts[0], int(parts[1]), int(parts[2])))
    return regions


def detect_fasta_chr_prefix(fasta_path: Path) -> bool:
    """True when the FASTA names contigs ``chr1``/``chrM`` rather than ``1``/``MT``."""
    with_prefix = without_prefix = 0
    with _open_text(fasta_path) as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            contig = line[1:].split()[0]
            if contig.startswith("chr"):
                with_prefix += 1
            elif contig.isdigit() or contig in ("X", "Y", "M", "MT"):
                without_prefix += 1
            if with_prefix + without_prefix >= 5:
                break
    return with_prefix > without_prefix


def normalise_bed_chromosomes(
    regions: list[tuple[str, int, int]], use_chr_prefix: bool
) -> list[tuple[str, int, int]]:
    """Rewrite BED contig names to match the target FASTA's naming convention."""
    normalised = []
    for chrom, start, end in regions:
        if use_chr_prefix and not chrom.startswith("chr"):
            chrom = f"chr{chrom}"
        elif not use_chr_prefix and chrom.startswith("chr"):
            chrom = chrom[3:]
        normalised.append((chrom, start, end))
    return normalised


def mask_fasta(
    input_fasta: Path,
    output_fasta: Path,
    numt_regions: list[tuple[str, int, int]] | None = None,
    line_width: int = 60,
) -> dict:
    """Write ``input_fasta`` with every NUMT region replaced by ``N``."""
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    for chrom, start, end in numt_regions or []:
        by_chrom.setdefault(chrom, []).append((start, end))

    stats = {"chromosomes_processed": 0, "regions_masked": 0, "bases_masked": 0}
    open_output = gzip.open if str(output_fasta).endswith(".gz") else open

    def emit(handle, chrom, sequence):
        if chrom is None:
            return
        # BED is 0-based half-open, which is exactly Python slice semantics.
        masked = list(sequence)
        for start, end in by_chrom.get(chrom, []):
            stop = min(end, len(masked))
            if stop > start:
                masked[start:stop] = "N" * (stop - start)
                stats["regions_masked"] += 1
                stats["bases_masked"] += stop - start
        joined = "".join(masked)
        handle.write(f">{chrom}\n")
        for offset in range(0, len(joined), line_width):
            handle.write(joined[offset : offset + line_width] + "\n")

    with _open_text(input_fasta) as f_in, open_output(output_fasta, "wt") as f_out:
        chrom = None
        sequence: list[str] = []
        for line in f_in:
            if line.startswith(">"):
                emit(f_out, chrom, "".join(sequence))
                chrom = line[1:].split()[0]
                sequence = []
                stats["chromosomes_processed"] += 1
            else:
                sequence.append(line.strip())
        emit(f_out, chrom, "".join(sequence))

    return stats
