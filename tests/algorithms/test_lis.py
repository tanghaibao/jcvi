#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "input_array,expected",
    [((4, 5, 1, 2, 3), (3, 1)), ((1, 2, 3, 5, 4), (4, 2)), ((1, 2, 1), (2, 0))],
)
def test_longest_monotonic_subseq_length(input_array, expected):
    from jcvi.algorithms.lis import longest_monotonic_subseq_length

    assert longest_monotonic_subseq_length(input_array) == expected


@pytest.mark.parametrize("input_array,expected", [(range(3), 3), ([3, 1, 2, 0], 2)])
def test_longest_increasing_subseq_length(input_array, expected):
    from jcvi.algorithms.lis import longest_increasing_subseq_length

    assert longest_increasing_subseq_length(input_array) == expected


@pytest.mark.parametrize(
    "input_array,expected", [(range(3), [0, 1, 2]), ([3, 1, 2, 0], [1, 2])]
)
def test_longest_increasing_subsequence(input_array, expected):
    from jcvi.algorithms.lis import longest_increasing_subsequence

    assert longest_increasing_subsequence(input_array) == expected


@pytest.mark.parametrize(
    "input_array,expected", [([23, 19, 97, 16, 37, 44, 88, 77, 26], [97, 88, 77, 26])]
)
def test_longest_decreasing_subsequence(input_array, expected):
    from jcvi.algorithms.lis import longest_decreasing_subsequence

    assert longest_decreasing_subsequence(input_array) == expected


@pytest.mark.parametrize(
    "input_array,expected", [([(3, 3), (2, 2), (1, 1), (0, 5)], ([(0, 5)], 5))]
)
def test_heaviest_increasing_subsequence(input_array, expected):
    from jcvi.algorithms.lis import heaviest_increasing_subsequence

    assert heaviest_increasing_subsequence(input_array) == expected
