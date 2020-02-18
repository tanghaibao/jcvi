#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op


def remove_ifexists(filename):
    if op.exists(filename):
        os.remove(filename)


def disabled_test_tree_main():
    from jcvi.graphics.tree import main as tree_main

    demo = "demo.png"
    remove_ifexists(demo)

    tree_main(["demo", "--format", "png"])
    assert op.exists(demo)
    remove_ifexists(demo)
