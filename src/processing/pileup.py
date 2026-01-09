"""Pileup generation and filtering for mtDNA analysis."""

from __future__ import annotations

import numpy as np

from core.config import PipelineConfig, SimpleRead


class PileupGenerator:
    """Generate pileup data from aligned reads with quality filtering."""

    def __init__(self, config: PipelineConfig):
        self.config = config
        self.bases = ["A", "C", "G", "T"]
        self.base_to_idx = {"A": 0, "C": 1, "G": 2, "T": 3}

    def generate_pileup(self, reads: list[SimpleRead]) -> dict[int, dict[str, int]]:
        """Count bases at each position, stratified by strand."""
        if not reads:
            return {}

        mito_length = self.config.mito_length
        base_counts = np.zeros((mito_length, 4, 2), dtype=np.uint32)
        tn5_cuts = np.zeros((mito_length, 2), dtype=np.uint32)

        min_mapq = self.config.quality.min_mapq
        min_baseq = self.config.quality.min_baseq
        min_dist_from_end = self.config.quality.min_distance_from_end
        base_to_idx = self.base_to_idx

        for read in reads:
            if read.mapping_quality < min_mapq:
                continue

            is_reverse = read.is_reverse
            strand_idx = 1 if is_reverse else 0
            qualities = read.query_qualities
            sequence = read.query_sequence
            sequence_view = memoryview(sequence)
            read_length = len(sequence)

            if is_reverse:
                start_pos = read.reference_start + read_length - 1
                if 0 <= start_pos < mito_length:
                    tn5_cuts[start_pos, 1] += 1
            else:
                start_pos = read.reference_start
                if 0 <= start_pos < mito_length:
                    tn5_cuts[start_pos, 0] += 1

            ref_pos = read.reference_start
            query_pos = 0

            for op, length in read.cigar:
                if op in [0, 7, 8]:
                    start_refpos = max(0, ref_pos)
                    end_refpos = min(ref_pos + length, mito_length)

                    if start_refpos >= end_refpos:
                        query_pos += length
                        ref_pos += length
                        continue

                    offset = start_refpos - ref_pos

                    if min_dist_from_end > 0:
                        valid_q_start = min_dist_from_end
                        valid_q_end = read_length - min_dist_from_end
                    else:
                        valid_q_start = 0
                        valid_q_end = read_length

                    for i in range(end_refpos - start_refpos):
                        current_qpos = query_pos + offset + i

                        if not (valid_q_start <= current_qpos < valid_q_end):
                            continue

                        if qualities[current_qpos] < min_baseq:
                            continue

                        base = chr(sequence_view[current_qpos]).upper()
                        base_idx = base_to_idx.get(base)
                        if base_idx is None:
                            continue

                        base_counts[start_refpos + i, base_idx, strand_idx] += 1

                    query_pos += length
                    ref_pos += length
                elif op in [2, 3]:
                    ref_pos += length
                elif op == 4:
                    query_pos += length

        pileup = {}
        bases = self.bases

        for pos in range(mito_length):
            total_fwd = base_counts[pos, :, 0].sum()
            total_rev = base_counts[pos, :, 1].sum()
            depth = total_fwd + total_rev

            if depth == 0 and tn5_cuts[pos].sum() == 0:
                continue

            pos_counts = {
                "depth": depth.item(),
                "tn5_cuts_fwd": tn5_cuts[pos, 0].item(),
                "tn5_cuts_rev": tn5_cuts[pos, 1].item(),
            }

            for base_idx in range(4):
                base = bases[base_idx]
                fwd_count = base_counts[pos, base_idx, 0].item()
                rev_count = base_counts[pos, base_idx, 1].item()
                total_count = fwd_count + rev_count

                pos_counts[base] = total_count
                pos_counts[f"{base}_fwd"] = fwd_count
                pos_counts[f"{base}_rev"] = rev_count

            pileup[pos] = pos_counts

        return pileup

    def filter_strand_bias(self, pileup: dict[int, dict[str, int]]) -> dict[int, dict[str, int]]:
        """Remove positions where most reads come from a single strand."""
        filtered = {}
        max_bias = self.config.quality.max_strand_bias
        bases = self.bases

        for pos, counts in pileup.items():
            filtered_counts = counts.copy()

            for base in bases:
                fwd = counts[f"{base}_fwd"]
                rev = counts[f"{base}_rev"]
                total = fwd + rev

                if total > 0:
                    bias = max(fwd, rev) / total
                    if bias > max_bias:
                        filtered_counts[base] = 0
                        filtered_counts[f"{base}_fwd"] = 0
                        filtered_counts[f"{base}_rev"] = 0

            filtered_counts["depth"] = sum(filtered_counts[b] for b in bases)

            if filtered_counts["depth"] > 0:
                filtered[pos] = filtered_counts

        return filtered
