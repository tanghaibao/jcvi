#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op


def disabled_test_tree_main():
    from jcvi.apps.base import cleanup
    from jcvi.graphics.tree import main as tree_main

    demo = "demo.png"
    cleanup(demo)

    tree_main(["demo", "--format", "png"])
    assert op.exists(demo)
    cleanup(demo)
