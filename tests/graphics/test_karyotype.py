#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import pytest

from jcvi.apps.base import cleanup
from jcvi.graphics.karyotype import main as karyotype_main, make_circle_name


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
        ("chrZ", {}, "Z"),
        ("ChrX", {"ChrX"}, "X-"),
    ],
)
def test_make_circle_name(sid, rev, expected):
    assert make_circle_name(sid, rev) == expected, "Expect {}".format(expected)


def test_main():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    cleanup("karyotype.pdf")
    image_name = karyotype_main(["seqids", "layout"])
    assert op.exists(image_name)
    cleanup("karyotype_with_comments.pdf")
    image_name = karyotype_main(
        ["seqids_with_comments", "layout", "-o", "karyotype_with_comments.pdf"]
    )
    assert op.exists(image_name)
    cleanup("karyotype_with_empty_lines.pdf")
    image_name = karyotype_main(
        ["seqids_with_empty_lines", "layout", "-o", "karyotype_with_empty_lines.pdf"]
    )
    os.chdir(cwd)
