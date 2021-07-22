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


@pytest.mark.parametrize(
    "hintfile,bed_filenames",
    [
        ("grape.peach.last", ("grape.bed", "peach.bed")),
        ("./grape.peach.last", ("./grape.bed", "./peach.bed")),
        (
            "../../test/opt/grape.peach.last",
            ("../../test/opt/grape.bed", "../../test/opt/peach.bed"),
        ),
    ],
)
def test_get_bed_filenames(hintfile, bed_filenames):
    from types import SimpleNamespace
    from jcvi.compara.synteny import get_bed_filenames

    opts = SimpleNamespace()
    opts.qbed, opts.sbed = None, None
    assert get_bed_filenames(hintfile, None, opts) == bed_filenames
