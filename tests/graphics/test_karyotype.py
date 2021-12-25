#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "sid,rev,expected",
    [
        ("chr123abc456", {}, "123"),
        ("scaffold_55", {"scaffold_55"}, "55-"),
        ("PMOO", {}, "P"),
        ("contig_xvv", {"contig_xvv"}, "X-"),
        ("chr01", {}, "1"),
        ("Ttru_Chr1", {}, "1"),
        ("pseudochromosome_1", {}, "1"),
    ],
)
def test_make_circle_name(sid, rev, expected):
    from jcvi.graphics.karyotype import make_circle_name

    assert make_circle_name(sid, rev) == expected, "Expect {}".format(expected)
