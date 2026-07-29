"""Shared fragment grouping, deduplication, and mate-overlap resolution."""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field

from core.config import SimpleRead


@dataclass
class Fragment:
    """Primary alignments belonging to one sequenced fragment."""

    query_name: str
    reads: list[SimpleRead] = field(default_factory=list)


@dataclass(frozen=True)
class Observation:
    """One quality-qualified base observation."""

    base: str
    base_quality: int
    mapping_quality: int
    read_position: int
    is_reverse: bool
    clipped: bool
    orientation: str | None


def alignment_key(start: int, is_reverse: bool, template_length: int, mode: str) -> tuple:
    """Return the historical mgatk2 coordinate deduplication key."""
    key = (start, is_reverse)
    if mode == "alignment_and_fragment_length":
        return (*key, abs(template_length))
    if mode == "alignment_start":
        return key
    raise ValueError(f"Deduplication mode has no compatibility key: {mode}")


def compatibility_key(read: SimpleRead, mode: str) -> tuple:
    """Return a compatibility key for a lightweight read."""
    return alignment_key(read.reference_start, read.is_reverse, read.template_length, mode)


def canonical_fragment_key(fragment: Fragment, mode: str) -> tuple:
    """Combine mate keys deterministically so a pair is one deduplication unit."""
    return tuple(sorted(compatibility_key(read, mode) for read in fragment.reads))


def representative_score(fragment: Fragment) -> tuple[int, int, int, int]:
    """Score duplicate representatives without depending on input order."""
    return (
        sum(int(q) for read in fragment.reads for q in read.query_qualities),
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
        groups[canonical_fragment_key(fragment, mode)].append(fragment)

    retained: list[Fragment] = []
    duplicate_groups = duplicate_reads = 0
    for group in groups.values():
        group.sort(key=lambda f: tuple(-v for v in representative_score(f)) + (f.query_name,))
        retained.append(group[0])
        if len(group) > 1:
            duplicate_groups += 1
            duplicate_reads += sum(len(fragment.reads) for fragment in group[1:])
    retained.sort(key=lambda f: f.query_name)
    return retained, {
        "duplicate_groups": duplicate_groups,
        "duplicate_reads": duplicate_reads,
    }


def fragment_orientation(fragment: Fragment) -> str | None:
    """Return VCF-compatible pair orientation when both mates are present."""
    read1 = next((read for read in fragment.reads if read.is_read1), None)
    read2 = next((read for read in fragment.reads if read.is_read2), None)
    if read1 is None or read2 is None:
        return None
    if not read1.is_reverse and read2.is_reverse:
        return "F1R2"
    if read1.is_reverse and not read2.is_reverse:
        return "F2R1"
    return None


def read_observations(
    read: SimpleRead, min_baseq: int, min_distance_from_end: int, orientation: str | None
) -> dict[int, Observation]:
    """Extract SNV observations; indels are deliberately not represented."""
    observations = {}
    clipped = any(op in {4, 5} for op, _length in read.cigar)
    read_length = len(read.query_sequence)
    for query_pos, ref_pos in read.get_aligned_pairs():
        if query_pos < min_distance_from_end or query_pos >= read_length - min_distance_from_end:
            continue
        quality = int(read.query_qualities[query_pos])
        if quality < min_baseq:
            continue
        base = chr(read.query_sequence[query_pos]).upper()
        if base not in "ACGT":
            continue
        observations[ref_pos] = Observation(
            base=base,
            base_quality=quality,
            mapping_quality=read.mapping_quality,
            read_position=query_pos,
            is_reverse=read.is_reverse,
            clipped=clipped,
            orientation=orientation,
        )
    return observations


def resolve_fragment_observations(
    fragment: Fragment, min_baseq: int, min_distance_from_end: int
) -> tuple[dict[int, Observation], dict[str, object]]:
    """Collapse mate overlaps to at most one observation per reference position."""
    orientation = fragment_orientation(fragment)
    by_position: dict[int, list[Observation]] = defaultdict(list)
    for read in fragment.reads:
        for position, observation in read_observations(
            read, min_baseq, min_distance_from_end, orientation
        ).items():
            by_position[position].append(observation)

    resolved = {}
    stats: dict[str, object] = {
        "overlap_positions": 0,
        "overlap_agreements": 0,
        "overlap_disagreements": 0,
        "disagreement_positions": set(),
    }
    for position, observations in by_position.items():
        if len(observations) == 1:
            resolved[position] = observations[0]
            continue
        stats["overlap_positions"] += 1
        best_quality = max(observation.base_quality for observation in observations)
        best = [
            observation for observation in observations if observation.base_quality == best_quality
        ]
        bases = {observation.base for observation in observations}
        if len(bases) == 1:
            stats["overlap_agreements"] += 1
            resolved[position] = sorted(
                best, key=lambda observation: (observation.is_reverse, -observation.mapping_quality)
            )[0]
        elif len({observation.base for observation in best}) == 1:
            stats["overlap_disagreements"] += 1
            stats["disagreement_positions"].add(position)
            resolved[position] = sorted(
                best,
                key=lambda observation: (
                    observation.is_reverse,
                    -observation.mapping_quality,
                    observation.read_position,
                ),
            )[0]
        else:
            stats["overlap_disagreements"] += 1
            stats["disagreement_positions"].add(position)
    return resolved, stats
