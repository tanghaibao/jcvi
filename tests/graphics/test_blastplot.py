#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import pytest


@pytest.mark.parametrize(
    "seqid,expected",
    [
        ("chr123abc456", "chr123abc456"),
        ("contig100", "c100"),
        ("scaffold_55", "c55"),
        ("PMOO", "PMOO"),
        ("contig_xvv", "xvv"),
        ("chr01", "chr01"),
    ],
)
def test_rename_seqid(seqid, expected):
    from jcvi.graphics.blastplot import rename_seqid

    assert rename_seqid(seqid) == expected, "Expect {}".format(expected)
