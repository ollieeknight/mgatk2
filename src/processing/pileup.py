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
        if getattr(self.config, "pileup_mode", "classic") == "fast":
            return self._generate_pileup_fast(reads)
        return self._generate_pileup_classic(reads)

    def _generate_pileup_classic(self, reads: list[SimpleRead]) -> dict[int, dict[str, int]]:
        """Original per-base Python loop. Reference implementation — do not modify."""
        mito_length = self.config.mito_length
        base_counts = np.zeros((mito_length, 4, 2), dtype=np.uint32)
        compute_tn5 = getattr(self.config, "compute_tn5", True)
        tn5_cuts = np.zeros((mito_length, 2), dtype=np.uint32) if compute_tn5 else None

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

            if tn5_cuts is not None:
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

        return self._build_pileup_dict(base_counts, tn5_cuts)

    def _generate_pileup_fast(self, reads: list[SimpleRead]) -> dict[int, dict[str, int]]:
        """Vectorized numpy pileup. Replaces inner per-base loop with boolean masking."""
        mito_length = self.config.mito_length
        base_counts = np.zeros((mito_length, 4, 2), dtype=np.uint32)
        compute_tn5 = getattr(self.config, "compute_tn5", True)
        tn5_cuts = np.zeros((mito_length, 2), dtype=np.uint32) if compute_tn5 else None

        min_mapq = self.config.quality.min_mapq
        min_baseq = self.config.quality.min_baseq
        min_dist = self.config.quality.min_distance_from_end

        # ASCII values for A C G T
        base_ascii = ((65, 0), (67, 1), (71, 2), (84, 3))

        for read in reads:
            if read.mapping_quality < min_mapq:
                continue

            strand_idx = 1 if read.is_reverse else 0
            sequence = read.query_sequence
            qualities = read.query_qualities
            read_length = len(sequence)

            if tn5_cuts is not None:
                if read.is_reverse:
                    tp = read.reference_start + read_length - 1
                    if 0 <= tp < mito_length:
                        tn5_cuts[tp, 1] += 1
                else:
                    tp = read.reference_start
                    if 0 <= tp < mito_length:
                        tn5_cuts[tp, 0] += 1

            ref_pos = read.reference_start
            query_pos = 0

            for op, length in read.cigar:
                if op in [0, 7, 8]:
                    r_start = max(0, ref_pos)
                    r_end = min(ref_pos + length, mito_length)

                    if r_start >= r_end:
                        query_pos += length
                        ref_pos += length
                        continue

                    # query slice bounds accounting for ref clamping
                    q_off = r_start - ref_pos
                    run_len = r_end - r_start
                    q_start = query_pos + q_off
                    q_end = q_start + run_len

                    # zero-copy numpy views over the bytes/array slices
                    q_bytes = np.frombuffer(sequence[q_start:q_end], dtype=np.uint8)
                    q_quals = qualities[q_start:q_end]

                    # quality mask
                    valid = q_quals >= min_baseq

                    # distance-from-end mask (absolute query positions)
                    if min_dist > 0:
                        abs_qpos = np.arange(q_start, q_end, dtype=np.int32)
                        dist_mask = (abs_qpos >= min_dist) & (abs_qpos < read_length - min_dist)
                        valid = valid & dist_mask

                    # per-base scatter: 4 boolean ops instead of run_len Python iterations
                    for b_ascii, b_idx in base_ascii:
                        hits = valid & (q_bytes == b_ascii)
                        if hits.any():
                            base_counts[r_start:r_end][hits, b_idx, strand_idx] += 1

                    query_pos += length
                    ref_pos += length
                elif op in [2, 3]:
                    ref_pos += length
                elif op == 4:
                    query_pos += length

        return self._build_pileup_dict(base_counts, tn5_cuts)

    def _build_pileup_dict(
        self,
        base_counts: np.ndarray,
        tn5_cuts: np.ndarray | None,
    ) -> dict[int, dict[str, int]]:
        """Convert base_counts array to pileup dict. Shared by both engines."""
        mito_length = self.config.mito_length
        bases = self.bases
        pileup = {}

        for pos in range(mito_length):
            total_fwd = base_counts[pos, :, 0].sum()
            total_rev = base_counts[pos, :, 1].sum()
            depth = total_fwd + total_rev

            tn5_fwd = tn5_cuts[pos, 0].item() if tn5_cuts is not None else 0
            tn5_rev = tn5_cuts[pos, 1].item() if tn5_cuts is not None else 0

            if depth == 0 and tn5_fwd == 0 and tn5_rev == 0:
                continue

            pos_counts = {
                "depth": depth.item(),
                "tn5_cuts_fwd": tn5_fwd,
                "tn5_cuts_rev": tn5_rev,
            }

            for base_idx in range(4):
                base = bases[base_idx]
                fwd_count = base_counts[pos, base_idx, 0].item()
                rev_count = base_counts[pos, base_idx, 1].item()
                pos_counts[base] = fwd_count + rev_count
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
