import numpy as np

from core.config import SimpleRead
from processing.fragments import (
    Fragment,
    canonical_fragment_key,
    deduplicate_fragments,
    group_reads_into_fragments,
    resolve_fragment_observations,
)


def _read(
    name,
    sequence="AAAA",
    quality=30,
    start=10,
    reverse=False,
    template_length=100,
    **kwargs,
):
    qualities = [quality] * len(sequence) if isinstance(quality, int) else quality
    return SimpleRead(
        reference_start=start,
        reference_end=start + len(sequence),
        is_reverse=reverse,
        mapping_quality=kwargs.pop("mapping_quality", 60),
        query_sequence=sequence.encode(),
        query_qualities=np.array(qualities),
        cigar=[(0, len(sequence))],
        template_length=template_length,
        query_name=name,
        **kwargs,
    )


def test_deduplication_modes_choose_the_best_read():
    low = Fragment("z", [_read("z", quality=10)])
    high = Fragment("a", [_read("a", quality=40)])
    different_length = Fragment("length", [_read("length", template_length=120)])

    retained, stats = deduplicate_fragments(
        [low, high, different_length], "alignment_and_fragment_length"
    )
    assert [fragment.query_name for fragment in retained] == ["a", "length"]
    assert stats == {"duplicate_groups": 1, "duplicate_reads": 1}

    retained, _ = deduplicate_fragments([low, high, different_length], "alignment_start")
    assert [fragment.query_name for fragment in retained] == ["a"]

    retained, _ = deduplicate_fragments([low, high], "none")
    assert retained == [low, high]


def test_deduplication_is_stable_for_paired_reads():
    first = Fragment(
        "a",
        [
            _read("a", is_read1=True, is_paired=True),
            _read("a", start=20, reverse=True, is_read2=True, is_paired=True),
        ],
    )
    second = Fragment(
        "b",
        [
            _read("b", is_read1=True, is_paired=True),
            _read("b", start=20, reverse=True, is_read2=True, is_paired=True),
        ],
    )

    retained, _ = deduplicate_fragments([second, first], "alignment_and_fragment_length")

    assert [fragment.query_name for fragment in retained] == ["a"]
    assert len(canonical_fragment_key(first, "alignment_and_fragment_length")) == 2


def test_grouping_reports_extra_alignments():
    fragments, collisions = group_reads_into_fragments(
        [_read("pair"), _read("pair", start=20), _read("pair", start=30), _read("orphan")]
    )

    assert sorted(len(fragment.reads) for fragment in fragments) == [1, 3]
    assert collisions == 1


def test_agreeing_mates_count_once():
    fragment = Fragment(
        "pair",
        [
            _read("pair", start=0, is_paired=True),
            _read("pair", start=0, reverse=True, is_paired=True),
        ],
    )

    observations, stats = resolve_fragment_observations(fragment, 20, 0)

    assert len(observations) == 4
    assert stats["overlap_agreements"] == 4


def test_overlap_disagreements_use_quality_or_are_masked():
    high_quality = Fragment(
        "pair",
        [
            _read("pair", "A", 30, start=0, is_paired=True),
            _read("pair", "C", 20, start=0, reverse=True, is_paired=True),
        ],
    )
    observations, stats = resolve_fragment_observations(high_quality, 0, 0)
    assert observations[0].base == "A"
    assert stats["overlap_disagreements"] == 1

    tied = Fragment(
        "pair",
        [
            _read("pair", "A", 30, start=0, is_paired=True),
            _read("pair", "C", 30, start=0, reverse=True, is_paired=True),
        ],
    )
    observations, stats = resolve_fragment_observations(tied, 0, 0)
    assert observations == {}
    assert stats["disagreement_positions"] == {0}


def test_separate_mates_are_unchanged():
    fragment = Fragment(
        "pair",
        [
            _read("pair", "AA", 30, start=0, is_paired=True),
            _read("pair", "CC", 30, start=3, reverse=True, is_paired=True),
        ],
    )

    observations, stats = resolve_fragment_observations(fragment, 0, 0)

    assert set(observations) == {0, 1, 3, 4}
    assert stats["overlap_positions"] == 0
