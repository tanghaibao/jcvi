#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest

from jcvi.algorithms.formula import calc_ldscore


@pytest.mark.parametrize(
    "a,b,expected",
    [
        ("AAA", "AAA", 0.0),
        ("AAB", "ABB", 0.25),
        ("AAB", "BBB", 0.0),
        ("AABB", "BBAA", 1.0),
    ],
)
def test_calc_ldscore(a: str, b: str, expected: float):
    assert calc_ldscore(a, b) == expected
