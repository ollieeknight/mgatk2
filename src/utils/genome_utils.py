"""
Genome nomenclature utilities for handling different naming conventions.
"""

import gzip
from importlib.resources import files
from pathlib import Path

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


def normalise_genome_name(genome: str) -> str:
    """Normalise genome name to canonical form"""
    normalised = GENOME_ALIASES.get(genome.lower())
    if normalised is None:
        supported = ", ".join(sorted(set(GENOME_ALIASES.keys())))
        raise ValueError(f"Genome '{genome}' not supported. Supported genomes: {supported}")
    return normalised


def get_reference_path(genome: str, data_dir: Path | None = None) -> Path:
    """Get path to mitochondrial reference FASTA for genome"""
    normalised = normalise_genome_name(genome)

    if data_dir is None:
        # Use importlib.resources to get package data
        try:
            package_files = files("data.references")
            ref_path = Path(str(package_files / f"{normalised}_MT.fasta"))
        except (ModuleNotFoundError, TypeError):
            # Fallback for development mode
            ref_path = (
                Path(__file__).parent.parent / "data" / "references" / f"{normalised}_MT.fasta"
            )
    else:
        ref_path = data_dir / f"{normalised}_MT.fasta"

    if not ref_path.exists():
        raise FileNotFoundError(
            f"Reference file not found: {ref_path}\nExpected location: {ref_path}"
        )

    return ref_path


def get_blacklist_path(genome: str, data_dir: Path | None = None) -> Path:
    """Get path to NUMT blacklist BED file for genome"""
    normalised = normalise_genome_name(genome)

    if data_dir is None:
        # Use importlib.resources to get package data
        try:
            package_files = files("data.blacklists")
            blacklist_path = Path(str(package_files / f"{normalised}_numts.bed"))
        except (ModuleNotFoundError, TypeError):
            # Fallback for development mode
            blacklist_path = (
                Path(__file__).parent.parent / "data" / "blacklists" / f"{normalised}_numts.bed"
            )
    else:
        blacklist_path = data_dir / f"{normalised}_numts.bed"

    if not blacklist_path.exists():
        raise FileNotFoundError(
            f"Blacklist file not found: {blacklist_path}\nExpected location: {blacklist_path}"
        )

    return blacklist_path


def get_mt_chrom_name(genome: str) -> str:
    """Get standard mitochondrial chromosome name for genome"""
    normalised = normalise_genome_name(genome)
    return MT_CHROM_NAMES.get(normalised, "chrM")


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


def detect_fasta_chr_prefix(fasta_path: Path) -> bool:
    """
    Detect whether FASTA uses 'chr' prefix by checking first few chromosome names.

    Returns:
        True if chromosomes use 'chr' prefix (chr1, chr2, chrM, etc.)
        False if chromosomes don't use prefix (1, 2, MT, etc.)
    """
    open_func = gzip.open if str(fasta_path).endswith(".gz") else open
    mode = "rt" if str(fasta_path).endswith(".gz") else "r"

    chr_count = 0
    no_chr_count = 0

    with open_func(fasta_path, mode) as f:
        for line in f:
            if line.startswith(">"):
                chrom = line.strip()[1:].split()[0]
                # Check if it's a standard chromosome (not scaffold/contig)
                if chrom.startswith("chr"):
                    chr_count += 1
                elif chrom.isdigit() or chrom in ("X", "Y", "M", "MT"):
                    no_chr_count += 1

                # Sample first 5 chromosomes
                if chr_count + no_chr_count >= 5:
                    break

    # Return True if majority use chr prefix
    return chr_count > no_chr_count


def normalise_bed_chromosomes(regions: list, use_chr_prefix: bool) -> list:
    """
    Normalise BED file chromosome names to match FASTA naming convention.

    Args:
        regions: List of (chrom, start, end) tuples from BED file
        use_chr_prefix: Whether target FASTA uses 'chr' prefix

    Returns:
        List of regions with normalised chromosome names
    """
    normalised = []
    for chrom, start, end in regions:
        if use_chr_prefix:
            # Add chr prefix if not present
            if not chrom.startswith("chr"):
                chrom = f"chr{chrom}"
        else:
            # Remove chr prefix if present
            if chrom.startswith("chr"):
                chrom = chrom[3:]
        normalised.append((chrom, start, end))
    return normalised
