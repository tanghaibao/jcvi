import pytest

from jcvi.utils.range import Range, range_closest


@pytest.mark.parametrize(
    "input,expected",
    [("chr1:1000-1", Range(seqid="chr1", start=1, end=1000, score=0, id=0))],
)
def test_range_parse(input, expected):
    from jcvi.utils.range import range_parse

    assert range_parse(input) == expected


@pytest.mark.parametrize(
    "a,b,expected", [((30, 45), (55, 65), None), ((48, 65), (45, 55), [48, 55])]
)
def test_range_intersect(a, b, expected):
    from jcvi.utils.range import range_intersect

    assert range_intersect(a, b) == expected


@pytest.mark.parametrize(
    "a,b,ratio,expected",
    [
        (("1", 30, 45), ("1", 41, 55), False, 5),
        (("1", 21, 45), ("1", 41, 75), True, 0.2),
        (("1", 30, 45), ("1", 15, 55), False, 16),
        (("1", 30, 45), ("1", 15, 55), True, 1.0),
        (("1", 30, 45), ("1", 57, 68), False, 0),
        (("1", 30, 45), ("2", 42, 55), False, 0),
        (("1", 30, 45), ("2", 42, 55), True, 0.0),
    ],
)
def test_range_overlap(a, b, ratio, expected):
    from jcvi.utils.range import range_overlap

    assert range_overlap(a, b, ratio) == expected


@pytest.mark.parametrize(
    "a,b,distmode,expected",
    [
        (("1", 30, 45, "+"), ("1", 45, 55, "+"), "ss", (26, "++")),
        (("1", 30, 45, "-"), ("1", 57, 68, "-"), "ss", (39, "--")),
        (("1", 30, 42, "-"), ("1", 45, 55, "+"), "ss", (26, "-+")),
        (("1", 30, 42, "+"), ("1", 45, 55, "-"), "ee", (2, "+-")),
    ],
)
def test_range_distance(a, b, distmode, expected):
    from jcvi.utils.range import range_distance

    assert range_distance(a, b, distmode) == expected


@pytest.mark.parametrize(
    "ranges,expected", [([(30, 45), (40, 50), (10, 100)], (10, 100))]
)
def test_range_minmax(ranges, expected):
    from jcvi.utils.range import range_minmax

    assert range_minmax(ranges) == expected


@pytest.mark.parametrize(
    "ranges,target,left,expected",
    [
        (
            [("1", 30, 40), ("1", 33, 35), ("1", 10, 20)],
            ("1", 22, 25),
            True,
            ("1", 10, 20),
        ),
        (
            [("1", 30, 40), ("1", 33, 35), ("1", 10, 20)],
            ("1", 22, 25),
            False,
            ("1", 33, 35),
        ),
        ([("1", 30, 40), ("1", 33, 35), ("1", 10, 20)], ("1", 2, 5), True, None),
    ],
)
def test_range_closest(ranges, target, left, expected):
    from jcvi.utils.range import range_closest

    assert range_closest(ranges, target, left) == expected


@pytest.mark.parametrize(
    "ranges,sizes,expected",
    [
        ([("1", 30, 40), ("1", 45, 50), ("1", 10, 30)], {}, [("1", 41, 44)]),
        ([("1", 30, 40), ("1", 42, 50)], {}, [("1", 41, 41)]),
        (
            [("1", 30, 40), ("1", 42, 50)],
            {"1": 70},
            [("1", 1, 29), ("1", 41, 41), ("1", 51, 70)],
        ),
    ],
)
def test_range_interleave(ranges, sizes, expected):
    from jcvi.utils.range import range_interleave

    assert range_interleave(ranges, sizes) == expected


@pytest.mark.parametrize(
    "ranges,dist,expected",
    [
        ([("1", 30, 45), ("1", 40, 50), ("1", 10, 50)], 0, [("1", 10, 50)]),
        ([("1", 30, 40), ("1", 45, 50)], 0, [("1", 30, 40), ("1", 45, 50)]),
        ([("1", 30, 40), ("1", 45, 50)], 5, [("1", 30, 50)]),
    ],
)
def test_range_merge(ranges, dist, expected):
    from jcvi.utils.range import range_merge

    assert range_merge(ranges, dist) == expected


@pytest.mark.parametrize(
    "ranges,expected",
    [
        ([("1", 30, 45), ("1", 40, 50), ("1", 10, 50)], 41),
        ([("1", 30, 45), ("2", 40, 50)], 27),
        ([("1", 30, 45), ("1", 45, 50)], 21),
        ([], 0),
    ],
)
def test_range_union(ranges, expected):
    from jcvi.utils.range import range_union

    assert range_union(ranges) == expected


@pytest.mark.parametrize(
    "ranges,expected",
    [
        ([("1", 30, 45), ("1", 40, 50), ("1", 10, 50)], 41),
        ([("1", 30, 45), ("2", 40, 50)], 27),
        ([("1", 30, 45), ("1", 45, 50)], 21),
        ([], 0),
    ],
)
def test_range_span(ranges, expected):
    from jcvi.utils.range import range_span

    assert range_span(ranges) == expected


@pytest.mark.parametrize(
    "ranges,expected",
    [
        (
            [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)],
            [[0, 1], [2]],
        ),
    ],
)
def test_range_piles(ranges, expected):
    from jcvi.utils.range import range_piles

    assert list(range_piles(ranges)) == expected


@pytest.mark.parametrize(
    "ranges,expected",
    [
        (
            [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)],
            [(0, 1)],
        ),
    ],
)
def test_range_conflict(ranges, expected):
    from jcvi.utils.range import range_conflict

    assert list(range_conflict(ranges)) == expected


@pytest.mark.parametrize(
    "ranges,expected",
    [
        (
            [
                Range("1", 0, 9, 22, 0),
                Range("1", 3, 18, 24, 1),
                Range("1", 10, 28, 20, 2),
            ],
            (
                [
                    Range(seqid="1", start=0, end=9, score=22, id=0),
                    Range(seqid="1", start=10, end=28, score=20, id=2),
                ],
                42,
            ),
        ),
        (
            [Range("2", 0, 1, 3, 0), Range("2", 1, 4, 3, 1), Range("3", 5, 7, 3, 2)],
            (
                [
                    Range(seqid="2", start=0, end=1, score=3, id=0),
                    Range(seqid="3", start=5, end=7, score=3, id=2),
                ],
                6,
            ),
        ),
    ],
)
def test_range_chain(ranges, expected):
    from jcvi.utils.range import range_chain

    assert range_chain(ranges) == expected
