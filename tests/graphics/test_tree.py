#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op


def remove_ifexists(filename):
    if op.exists(filename):
        os.remove(filename)


def test_tree_main():
    from jcvi.graphics.tree import main as tree_main

    demopdf = "demo.pdf"
    remove_ifexists(demopdf)

    tree_main(["demo"])
    assert op.exists(demopdf)
    remove_ifexists(demopdf)
