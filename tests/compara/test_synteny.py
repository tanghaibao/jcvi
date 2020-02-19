#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "ia,ib,expected",
    [
        ((), (), "+"),
        ((1,), (2,), "+"),
        ((1, 2, 3), (4, 5, 6), "+"),
        ((3, 2), (1, 2), "-"),
        (list(range(10)), list(reversed(range(10))), "-"),
    ],
)
def test_get_orientation(ia, ib, expected):
    from jcvi.compara.synteny import get_orientation

    assert get_orientation(ia, ib) == expected
