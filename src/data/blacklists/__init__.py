"""Blacklist files for NUMT regions"""

from pathlib import Path

_BUNDLED: dict[str, Path] = {
    "hg38": Path(__file__).parent / "hg38_numts.bed",
    "hg19": Path(__file__).parent / "hg19_numts.bed",
    "mm10": Path(__file__).parent / "mm10_numts.bed",
    "mm9": Path(__file__).parent / "mm9_numts.bed",
}


def load_blacklist_positions(
    build: str = "hg38",
    custom_bed: str | None = None,
    mito_chr: str = "chrM",
) -> set[int]:
    """Return 1-based positions on mito_chr that fall within the blacklist BED.

    Args:
        build: genome build for bundled BED (hg38/hg19/mm10/mm9) or "none"
        custom_bed: path to a custom BED file (overrides build)
        mito_chr: chromosome name to filter on

    Returns:
        Set of 1-based chrM positions; empty set when build="none" and no custom_bed
    """
    if build.lower() == "none" and custom_bed is None:
        return set()

    bed_path = Path(custom_bed) if custom_bed else _BUNDLED.get(build.lower())
    if bed_path is None:
        raise ValueError(
            f"Unknown blacklist build '{build}'. Choose from: {', '.join(_BUNDLED)} or 'none'"
        )
    if not bed_path.exists():
        raise FileNotFoundError(f"Blacklist BED not found: {bed_path}")

    positions: set[int] = set()
    with open(bed_path) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            if chrom != mito_chr:
                continue
            # BED is 0-based half-open [start, end) → convert to 1-based inclusive
            for pos in range(start + 1, end + 1):
                positions.add(pos)

    return positions
