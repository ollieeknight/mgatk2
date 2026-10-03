from core.config import SimpleRead
from processing.fragments import (
    Fragment,
    deduplicate_fragments,
    fragment_observations,
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
        is_reverse=reverse,
        mapping_quality=kwargs.pop("mapping_quality", 60),
        query_sequence=sequence.encode(),
        query_qualities=bytes(qualities),
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


A, C = 0, 1


def test_agreeing_mates_count_once():
    fragment = Fragment(
        "pair",
        [
            _read("pair", start=0, is_paired=True),
            _read("pair", start=0, reverse=True, is_paired=True),
        ],
    )

    observations, stats = fragment_observations([fragment], 20, 0)

    assert observations["position"].tolist() == [0, 1, 2, 3]
    assert observations["strand"].tolist() == [0, 0, 0, 0]  # forward mate kept
    assert stats["overlap_agreements"] == 4


def test_overlap_disagreements_use_quality_or_are_masked():
    def disagreeing(quality_a, quality_c):
        fragment = Fragment(
            "pair",
            [
                _read("pair", "A", quality_a, start=0, is_paired=True),
                _read("pair", "C", quality_c, start=0, reverse=True, is_paired=True),
            ],
        )
        return fragment_observations([fragment], 0, 0)

    observations, stats = disagreeing(30, 20)
    assert observations["allele"].tolist() == [A]
    assert stats["overlap_disagreements"] == 1

    observations, _ = disagreeing(20, 30)
    assert observations["allele"].tolist() == [C]

    observations, _ = disagreeing(30, 30)
    assert observations["position"].size == 0


def test_end_trimming_and_insertions_keep_register():
    read = _read("solo", "ACGTTTACGT", start=100)
    read.cigar = [(4, 2), (0, 3), (1, 2), (0, 3)]  # soft clip, match, insertion, match
    observations, _ = fragment_observations([Fragment("solo", [read])], 0, 1)

    # Query 0-1 clipped, 9 trimmed as a read end; the insertion (5-6) advances
    # the query without consuming the reference.
    assert observations["position"].tolist() == [100, 101, 102, 103, 104]
    assert observations["distance"].tolist() == [2, 3, 4, 2, 1]
