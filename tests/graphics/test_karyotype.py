#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import pytest

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


def test_main(copy_data):
    copy_data(op.join(op.dirname(__file__), "data"))

    image_name = karyotype_main(["seqids", "layout"])
    assert op.exists(image_name)
    image_name = karyotype_main(
        ["seqids_with_comments", "layout", "-o", "karyotype_with_comments.pdf"]
    )
    assert op.exists(image_name)
    image_name = karyotype_main(
        ["seqids_with_empty_lines", "layout", "-o", "karyotype_with_empty_lines.pdf"]
    )
    assert op.exists(image_name)


TRACKS = """\
 .6, .1, .8, 0, , Grape, top, grape.bed
 .4, .1, .8, 0, , Peach, top, peach.bed
"""


@pytest.mark.parametrize(
    "seqids,layout,error",
    [
        # https://github.com/tanghaibao/jcvi/issues/771
        ("chr1\n", TRACKS + "e, 0, 1, grape.peach.anchors.simple\n", "has 1 rows"),
        # https://github.com/tanghaibao/jcvi/issues/557
        ("chr1\nPp01\nPp02\n", TRACKS, "has 3 rows"),
        # https://github.com/tanghaibao/jcvi/issues/769
        ("chr1\nPp01\n", TRACKS + "e, 1, 2, grape.peach.anchors.simple\n", "0-1"),
        # https://github.com/tanghaibao/jcvi/issues/788
        ("chr1\nPp01\n", " .8, .1, .8, 0, , Fake, top\n" + TRACKS, "8 columns"),
        (
            "chrF\nPp01\n",
            TRACKS.replace("grape.bed", "fake.bed")
            + "e, 0, 1, grape.peach.anchors.simple\n",
            "not in",
        ),
        # https://github.com/tanghaibao/jcvi/issues/768
        ("chr1\nPp01\n", TRACKS + "e, 0, 1, three.anchors\n", "6 columns"),
    ],
)
def test_main_invalid_input(copy_data, seqids, layout, error):
    copy_data(op.join(op.dirname(__file__), "data"))
    with open("fake.bed", "w") as fw:
        print("chrF\t0\t100\tfake", file=fw)
    with open("three.anchors", "w") as fw:
        print("GSVIVT01012255001\tPrupe.1G281700.1\t110", file=fw)
    with open("seqids_bad", "w") as fw:
        fw.write(seqids)
    with open("layout_bad", "w") as fw:
        fw.write(layout)

    with pytest.raises(ValueError, match=error):
        karyotype_main(["seqids_bad", "layout_bad", "--notex"])


def test_main_reversed_track_label(copy_data):
    from jcvi.graphics.base import plt

    copy_data(op.join(op.dirname(__file__), "data"))
    # https://github.com/tanghaibao/jcvi/issues/774, peach runs right to left
    with open("layout_rev", "w") as fw:
        fw.write(" .6, .1, .8, 0, , Grape, top, grape.bed\n")
        fw.write(" .4, .8, .1, 0, , Peach, top, peach.bed\n")
        fw.write("e, 0, 1, grape.peach.anchors.simple\n")
    karyotype_main(["seqids", "layout_rev", "--notex"])
    texts = {t.get_text(): t for t in plt.gcf().axes[0].texts}
    plt.close("all")
    # Labels sit left of the track, aligned away from it
    for label in ("Grape", "Peach"):
        assert texts[label].get_position()[0] < 0.1
        assert texts[label].get_ha() == "right"


@pytest.mark.parametrize("xstart,xend", [(0.1, 0.8), (0.8, 0.1)])
def test_track_gaps(xstart, xend):
    from types import SimpleNamespace

    from jcvi.graphics.base import plt
    from jcvi.graphics.karyotype import Track

    seqids = ["chr1", "chr2", "chr3"]
    t = SimpleNamespace(
        empty=False, y=0.5, xstart=xstart, xend=xend, rotation=0, va="top",
        label="t", label_va="center", color="k", seqids=seqids, rev=set(),
        sizes=dict(zip(seqids, (100, 200, 300))), bed=None, order=None,
        order_in_chr=None,
    )  # fmt: skip
    ax = plt.figure().add_axes((0, 0, 1, 1))
    tr = Track(ax, t, gap=0.01, draw=False)
    plt.close("all")

    sign = 1 if xend > xstart else -1
    ends = [tr.offsets[x] + tr.ratio * tr.sizes[x] for x in seqids]
    # Chromosomes are separated by the gap, also on reversed tracks (#774)
    for end, nxt in zip(ends, seqids[1:]):
        assert sign * (tr.offsets[nxt] - end) == pytest.approx(0.01)
    assert tr.offsets["chr1"] == xstart
    assert ends[-1] == pytest.approx(xend)
