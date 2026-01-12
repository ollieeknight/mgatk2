"""
Genome nomenclature utilities for handling different naming conventions.
"""

import sys
from pathlib import Path

if sys.version_info >= (3, 9):
    from importlib.resources import files
else:
    from importlib_resources import files

# Genome nomenclature mapping - bidirectional
GENOME_ALIASES = {
    # Human genomes
    "hg38": "hg38",
    "GRCh38": "hg38",
    "grch38": "hg38",
    "hg19": "hg19",
    "GRCh37": "hg19",
    "grch37": "hg19",
    # Mouse genomes
    "mm10": "mm10",
    "GRCm38": "mm10",
    "grcm38": "mm10",
    "mm9": "mm9",
    "GRCm37": "mm9",
    "grcm37": "mm9",
}

# Canonical genome names
SUPPORTED_GENOMES = ["hg38", "hg19", "mm10", "mm9"]

# Mitochondrial chromosome naming variations
MT_CHROM_NAMES = {
    "hg38": "chrM",
    "hg19": "chrM",
    "mm10": "chrM",
    "mm9": "chrM",
}


def normalize_genome_name(genome: str) -> str:
    """Normalise genome name to canonical form"""
    normalized = GENOME_ALIASES.get(genome.lower())
    if normalized is None:
        supported = ", ".join(sorted(set(GENOME_ALIASES.keys())))
        raise ValueError(f"Genome '{genome}' not supported. Supported genomes: {supported}")
    return normalized


def get_reference_path(genome: str, data_dir: Path | None = None) -> Path:
    """Get path to mitochondrial reference FASTA for genome"""
    normalized = normalize_genome_name(genome)

    if data_dir is None:
        # Use importlib.resources to get package data
        try:
            package_files = files("data.references")
            ref_path = Path(str(package_files / f"{normalized}_MT.fasta"))
        except (ModuleNotFoundError, TypeError):
            # Fallback for development mode
            ref_path = (
                Path(__file__).parent.parent / "data" / "references" / f"{normalized}_MT.fasta"
            )
    else:
        ref_path = data_dir / f"{normalized}_MT.fasta"

    if not ref_path.exists():
        raise FileNotFoundError(
            f"Reference file not found: {ref_path}\nExpected location: {ref_path}"
        )

    return ref_path


def get_blacklist_path(genome: str, data_dir: Path | None = None) -> Path:
    """Get path to NUMT blacklist BED file for genome"""
    normalized = normalize_genome_name(genome)

    if data_dir is None:
        # Use importlib.resources to get package data
        try:
            package_files = files("data.blacklists")
            blacklist_path = Path(str(package_files / f"{normalized}_numts.bed"))
        except (ModuleNotFoundError, TypeError):
            # Fallback for development mode
            blacklist_path = (
                Path(__file__).parent.parent / "data" / "blacklists" / f"{normalized}_numts.bed"
            )
    else:
        blacklist_path = data_dir / f"{normalized}_numts.bed"

    if not blacklist_path.exists():
        raise FileNotFoundError(
            f"Blacklist file not found: {blacklist_path}\nExpected location: {blacklist_path}"
        )

    return blacklist_path


def get_mt_chrom_name(genome: str) -> str:
    """Get standard mitochondrial chromosome name for genome"""
    normalized = normalize_genome_name(genome)
    return MT_CHROM_NAMES.get(normalized, "chrM")


def load_blacklist_regions(blacklist_path: Path) -> list:
    """Load NUMT blacklist regions from BED file"""
    regions = []
    with open(blacklist_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                regions.append((chrom, start, end))
    return regions
