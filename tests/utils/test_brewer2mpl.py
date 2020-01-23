#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "map_type,number,expected",
    [("Qualitative", 9, 4), ("Sequential", 9, 18), ("Diverging", 9, 9)],
)
def test_get_maps(map_type, number, expected):
    from jcvi.utils.brewer2mpl import get_maps

    got = len(get_maps(map_type, number))
    assert got == expected
