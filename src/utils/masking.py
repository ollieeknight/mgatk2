"""
FASTA masking utilities for creating hard-masked reference genomes.
"""

import gzip
import logging
from pathlib import Path

logger = logging.getLogger(__name__)


def mask_fasta(
    input_fasta: Path,
    output_fasta: Path,
    numt_regions: list | None = None,
    mt_chrom: str = "chrM",
    mt_blacklist_positions: set[int] | None = None,
) -> dict:
    """Hard-mask reference genome FASTA with N's at specified regions"""
    stats = {
        "chromosomes_processed": 0,
        "numt_regions_masked": 0,
        "numt_bases_masked": 0,
        "mt_positions_masked": 0,
        "total_output_size": 0,
    }

    # Create mask dictionaries for faster lookup
    numt_mask: dict[str, list[tuple[int, int]]] = {}
    if numt_regions:
        for chrom, start, end in numt_regions:
            if chrom not in numt_mask:
                numt_mask[chrom] = []
            numt_mask[chrom].append((start, end))

    # Handle both gzipped and uncompressed FASTA files
    open_func = gzip.open if str(input_fasta).endswith(".gz") else open
    mode = "rt" if str(input_fasta).endswith(".gz") else "r"

    with open_func(input_fasta, mode) as f_in, open(output_fasta, "w") as f_out:
        current_chrom = None
        sequence: list[str] = []

        for line in f_in:
            if line.startswith(">"):
                # Process previous chromosome before starting new one
                if current_chrom is not None and sequence:
                    masked_seq = _process_chromosome(
                        current_chrom,
                        "".join(sequence),
                        numt_mask.get(current_chrom, []),
                        mt_chrom,
                        mt_blacklist_positions,
                        stats,
                    )
                    _write_fasta_sequence(f_out, current_chrom, masked_seq)
                    stats["total_output_size"] += len(masked_seq)

                # Start new chromosome
                current_chrom = line.strip()[1:].split()[0]
                sequence = []
                stats["chromosomes_processed"] += 1
                if stats["chromosomes_processed"] <= 25:
                    logger.debug("Processing chromosome: %s", current_chrom)
            else:
                sequence.append(line.strip())

        # Process last chromosome
        if current_chrom is not None and sequence:
            masked_seq = _process_chromosome(
                current_chrom,
                "".join(sequence),
                numt_mask.get(current_chrom, []),
                mt_chrom,
                mt_blacklist_positions,
                stats,
            )
            _write_fasta_sequence(f_out, current_chrom, masked_seq)
            stats["total_output_size"] += len(masked_seq)

    return stats


def _process_chromosome(
    chrom: str,
    sequence: str,
    numt_regions: list,
    mt_chrom: str,
    mt_blacklist_positions: set[int] | None,
    stats: dict,
) -> str:
    """Process a single chromosome sequence, applying masks"""
    seq_list = list(sequence)

    # Mask NUMT regions on nuclear chromosomes
    if numt_regions:
        for start, end in numt_regions:
            # BED format is 0-indexed, half-open [start, end)
            for i in range(start, min(end, len(seq_list))):
                seq_list[i] = "N"
            stats["numt_regions_masked"] += 1
            stats["numt_bases_masked"] += min(end, len(seq_list)) - start

    # Mask blacklist positions on mitochondrial chromosome
    if chrom == mt_chrom and mt_blacklist_positions:
        for pos in mt_blacklist_positions:
            # Positions are 1-indexed
            idx = pos - 1
            if 0 <= idx < len(seq_list):
                seq_list[idx] = "N"
                stats["mt_positions_masked"] += 1

    return "".join(seq_list)


def _write_fasta_sequence(file_handle, chrom: str, sequence: str, line_width: int = 60):
    """Write FASTA sequence with proper line wrapping"""
    file_handle.write(f">{chrom}\n")
    for i in range(0, len(sequence), line_width):
        file_handle.write(sequence[i : i + line_width] + "\n")


def extract_mt_sequence(fasta_path: Path, mt_chrom: str = "MT") -> str:
    """Extract mitochondrial sequence from FASTA file"""
    current_chrom = None
    sequence: list[str] = []

    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                # Found new chromosome
                if current_chrom == mt_chrom and sequence:
                    # Already collected MT sequence
                    return "".join(sequence)

                # Start new chromosome
                current_chrom = line.strip()[1:].split()[0]
                sequence = []
            elif current_chrom == mt_chrom:
                sequence.append(line.strip())

    # Check if we found it at the end
    if current_chrom == mt_chrom and sequence:
        return "".join(sequence)

    raise ValueError(f"Mitochondrial chromosome '{mt_chrom}' not found in {fasta_path}")
