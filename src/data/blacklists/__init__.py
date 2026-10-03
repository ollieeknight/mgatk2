"""Bundled nuclear NUMT BEDs and a BED position reader."""

from pathlib import Path

BLACKLIST_DIR = Path(__file__).parent


def load_bed_positions(bed_path: str, mito_chr: str) -> set[int]:
    """1-based positions on mito_chr covered by a BED file."""
    positions: set[int] = set()
    with open(bed_path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != mito_chr:
                continue
            # BED is 0-based half-open; positions here are 1-based inclusive.
            positions.update(range(int(parts[1]) + 1, int(parts[2]) + 1))
    return positions
