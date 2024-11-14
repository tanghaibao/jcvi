#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op


def test_anchorfile():
    from jcvi.compara.base import AnchorFile

    test_file = op.join(op.dirname(__file__), "test.anchors")
    anchor_file = AnchorFile(test_file)
    assert not anchor_file.is_empty
    assert len(list(anchor_file.iter_blocks())) == 3
    assert len(list(anchor_file.iter_pairs())) == 14
