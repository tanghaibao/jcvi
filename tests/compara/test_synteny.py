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


def test_process_blast_for_liftover(tmp_path):
    from jcvi.compara.synteny import process_blast_for_liftover
    from jcvi.formats.bed import Bed

    bed = tmp_path / "a.bed"
    bed.write_text(
        "".join(
            f"chr{c}\t{i * 100}\t{i * 100 + 50}\tg{c}{i}\n"
            for c in (1, 2)
            for i in range(5)
        )
    )
    order = Bed(str(bed)).order
    hit = "{}\t{}\t90\t100\t1\t0\t1\t100\t1\t100\t1e-50\t{}\n"
    blast = tmp_path / "a.b.last"
    blast.write_text(
        hit.format("g10", "g11", 50)  # anchor pair
        + hit.format("g10", "g11", 80)  # second HSP, ignored
        + hit.format("g12", "g13", 70)  # same chr pair, near anchor
        + hit.format("g12", "g13", 90)  # second HSP, ignored
        + hit.format("g20", "g10", 60)  # chr pair without anchors
    )

    accepted, all_hits, nblast = process_blast_for_liftover(
        str(blast), order, order, {("g10", "g11")}, {("chr1", "chr1")}, ostrip=False
    )
    assert nblast == 5
    assert accepted == {("g10", "g11"): "50"}
    assert list(all_hits) == [("chr1", "chr1")]
    assert all_hits[("chr1", "chr1")].tolist() == [[0, 1, 50], [2, 3, 70]]
