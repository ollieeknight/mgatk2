"""Shared fragment grouping, deduplication, and mate-overlap resolution."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

import numpy as np

from core.config import SimpleRead
from processing.pileup import BASE_LUT


@dataclass
class Fragment:
    """Primary alignments belonging to one sequenced fragment."""

    query_name: str
    reads: list[SimpleRead] = field(default_factory=list)


def fragment_key(fragment: Fragment, mode: str) -> tuple:
    """Coordinate key for a whole fragment, so a mate pair is one deduplication unit."""
    keys = []
    for read in fragment.reads:
        key = (read.reference_start, read.is_reverse)
        if mode == "alignment_and_fragment_length":
            key += (abs(read.template_length),)
        keys.append(key)
    return tuple(sorted(keys))


def representative_score(fragment: Fragment) -> tuple[int, int, int, int]:
    """Score duplicate representatives without depending on input order."""
    return (
        sum(sum(read.query_qualities) for read in fragment.reads),
        sum(read.mapping_quality for read in fragment.reads),
        sum(read.is_proper_pair for read in fragment.reads),
        len(fragment.reads),
    )


def group_reads_into_fragments(reads: list[SimpleRead]) -> tuple[list[Fragment], int]:
    """Group mates by query name and report groups containing extra alignments."""
    grouped: dict[str, list[SimpleRead]] = defaultdict(list)
    for index, read in enumerate(reads):
        grouped[read.query_name or f"__unnamed_{index}"].append(read)
    collisions = sum(max(0, len(group) - 2) for group in grouped.values())
    return [Fragment(name, group) for name, group in grouped.items()], collisions


def deduplicate_fragments(
    fragments: list[Fragment], mode: str
) -> tuple[list[Fragment], dict[str, int]]:
    """Select one deterministic representative for each compatibility key."""
    if mode == "none":
        return fragments, {"duplicate_groups": 0, "duplicate_reads": 0}

    groups: dict[tuple, list[Fragment]] = defaultdict(list)
    for fragment in fragments:
        groups[fragment_key(fragment, mode)].append(fragment)

    retained: list[Fragment] = []
    duplicate_groups = duplicate_reads = 0
    for group in groups.values():
        group.sort(key=lambda f: tuple(-v for v in representative_score(f)) + (f.query_name,))
        retained.append(group[0])
        if len(group) > 1:
            duplicate_groups += 1
            duplicate_reads += sum(len(fragment.reads) for fragment in group[1:])
    retained.sort(key=lambda f: f.query_name)
    return retained, {"duplicate_groups": duplicate_groups, "duplicate_reads": duplicate_reads}


def _orientation(fragment: Fragment) -> int:
    """0 for F1R2, 1 for F2R1, -1 unless both mates are present on opposite strands."""
    read1 = next((read for read in fragment.reads if read.is_read1), None)
    read2 = next((read for read in fragment.reads if read.is_read2), None)
    if read1 is None or read2 is None or read1.is_reverse == read2.is_reverse:
        return -1
    return 1 if read1.is_reverse else 0


def _reduce(ufunc, values: np.ndarray, starts: np.ndarray) -> np.ndarray:
    """Per-group reduction; reduceat cannot take an empty array."""
    return ufunc.reduceat(values, starts) if len(starts) else values[:0]


def fragment_observations(
    fragments: list[Fragment], min_baseq: int, min_distance_from_end: int
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    """SNV observations with mate overlaps collapsed to one per fragment and position.

    Returns parallel arrays (position, allele, strand, base_quality,
    mapping_quality, distance, orientation) plus overlap counts. Where mates
    overlap and agree, the forward, higher-MAPQ observation of the best base
    quality is kept. Where they disagree, the higher-quality base wins, and a
    quality tie masks the position. Indels are deliberately not represented.
    """
    reads = [read for fragment in fragments for read in fragment.reads]
    fragment_of = np.repeat(np.arange(len(fragments)), [len(f.reads) for f in fragments])
    orientation = np.array([_orientation(f) for f in fragments], dtype=np.int64)[fragment_of]

    # Every aligned block as (read, query start, reference start, length).
    blocks = []
    for index, read in enumerate(reads):
        query_pos, ref_pos = 0, read.reference_start
        for op, length in read.cigar:
            if op in (0, 7, 8):
                blocks.append((index, query_pos, ref_pos, length))
                query_pos += length
                ref_pos += length
            elif op in (1, 4):
                query_pos += length
            elif op in (2, 3):
                ref_pos += length
    block_read, block_query, block_ref, block_length = (
        np.array(blocks, dtype=np.int64).reshape(-1, 4).T
    )

    # Expand blocks to one row per aligned base.
    read_of = np.repeat(block_read, block_length)
    offset = np.arange(len(read_of)) - np.repeat(
        np.cumsum(block_length) - block_length, block_length
    )
    query_pos = np.repeat(block_query, block_length) + offset
    ref_pos = np.repeat(block_ref, block_length) + offset

    read_length = np.array([len(read.query_sequence) for read in reads], dtype=np.int64)
    read_start = np.cumsum(read_length) - read_length
    sequence = np.frombuffer(b"".join(read.query_sequence for read in reads), dtype=np.uint8)
    qualities = np.frombuffer(b"".join(read.query_qualities for read in reads), dtype=np.uint8)

    length_of = read_length[read_of]
    usable = (query_pos >= min_distance_from_end) & (query_pos < length_of - min_distance_from_end)
    index = read_start[read_of] + np.where(usable, query_pos, 0)
    allele = BASE_LUT[sequence[index]].astype(np.int64)
    base_quality = qualities[index].astype(np.int64)
    usable &= (base_quality >= min_baseq) & (allele < 4)

    read_of, query_pos, ref_pos = read_of[usable], query_pos[usable], ref_pos[usable]
    allele, base_quality, length_of = allele[usable], base_quality[usable], length_of[usable]
    is_reverse = np.array([read.is_reverse for read in reads], dtype=np.int64)[read_of]
    mapping_quality = np.array([read.mapping_quality for read in reads], dtype=np.int64)[read_of]
    distance = np.minimum(query_pos, length_of - 1 - query_pos)

    # One group per (fragment, position); a stable sort keeps read order within it.
    key = fragment_of[read_of] * (int(ref_pos.max(initial=0)) + 1) + ref_pos
    order = np.argsort(key, kind="stable")
    key = key[order]
    starts = np.flatnonzero(np.r_[True, key[1:] != key[:-1]])
    size = np.diff(np.r_[starts, len(key)])
    group = np.repeat(np.arange(len(starts)), size)

    if not len(key):
        starts = size = group = np.zeros(0, dtype=np.int64)
    sorted_allele, sorted_quality = allele[order], base_quality[order]
    best = sorted_quality == _reduce(np.maximum, sorted_quality, starts)[group]
    agree = _reduce(np.minimum, sorted_allele, starts) == _reduce(np.maximum, sorted_allele, starts)
    best_low = _reduce(np.minimum, np.where(best, sorted_allele, 4), starts)
    best_high = _reduce(np.maximum, np.where(best, sorted_allele, -1), starts)
    keep = agree | (best_low == best_high)

    # Within each group: best quality first, then forward, higher MAPQ, and (on
    # disagreement only) further from the read end; ties fall to read order.
    tiebreak = np.where(agree[group], 0, -distance[order])
    ranked = np.lexsort(
        (
            np.arange(len(order)),
            tiebreak,
            -mapping_quality[order],
            is_reverse[order],
            ~best,
            group,
        )
    )
    chosen = order[ranked[starts]][keep]

    overlapping = size > 1
    stats = {
        "overlap_positions": int(overlapping.sum()),
        "overlap_agreements": int((overlapping & agree).sum()),
        "overlap_disagreements": int((overlapping & ~agree).sum()),
    }
    arrays = {
        "position": ref_pos[chosen],
        "allele": allele[chosen],
        "strand": is_reverse[chosen],
        "base_quality": base_quality[chosen],
        "mapping_quality": mapping_quality[chosen],
        "distance": distance[chosen],
        "orientation": orientation[read_of[chosen]],
    }
    return arrays, stats
