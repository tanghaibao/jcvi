#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os
import os.path as op
import tempfile

import matplotlib

matplotlib.use("Agg")

import pytest


def test_tree_main(tmp_path):
    """Render the built-in demo tree to a PNG and check the file is created."""
    from jcvi.apps.base import cleanup
    from jcvi.graphics.tree import main as tree_main

    demo = str(tmp_path / "demo.png")
    cleanup(demo)

    os.chdir(tmp_path)
    tree_main(["demo", "--format", "png"])
    assert op.exists(str(tmp_path / "demo.png"))


def test_tree_parse():
    """Test that the ete4 Tree can be created from a newick string."""
    from ete4 import Tree

    newick = "(((A:0.1,B:0.2)90:0.3,(C:0.15,D:0.25)80:0.2)100:0.05,E:0.4);"
    t = Tree(newick)
    assert t is not None
    leaves = list(t.leaf_names())
    assert sorted(leaves) == ["A", "B", "C", "D", "E"]


def test_tree_parse_tree_function():
    """Test parse_tree handles a plain newick file (no HPD)."""
    from jcvi.graphics.tree import parse_tree

    newick = "(((A:0.1,B:0.2)90:0.3,(C:0.15,D:0.25)80:0.2)100:0.05,E:0.4);"
    with tempfile.NamedTemporaryFile(mode="w", suffix=".nwk", delete=False) as f:
        f.write(newick)
        fname = f.name
    try:
        t, hpd = parse_tree(fname)
        assert hpd is None
        leaves = list(t.leaf_names())
        assert sorted(leaves) == ["A", "B", "C", "D", "E"]
    finally:
        os.unlink(fname)


def test_tree_truncate_name():
    """Test the truncate_name helper function."""
    from jcvi.graphics.tree import truncate_name

    assert truncate_name("ABCDEF", rule=None) == "ABCDEF"
    assert truncate_name("ABCDEF", rule="head3") == "DEF"
    assert truncate_name("ABCDEF", rule="ohead3") == "ABC"
    assert truncate_name("ABCDEF", rule="tail2") == "ABCD"
    assert truncate_name("ABCDEF", rule="otail2") == "EF"


def test_tree_draw(tmp_path):
    """Test that draw_tree runs without error on a simple tree."""
    import matplotlib.pyplot as plt
    from ete4 import Tree

    from jcvi.graphics.tree import draw_tree

    newick = "(((A:0.1,B:0.2)0.9:0.3,(C:0.15,D:0.25)0.8:0.2)1.0:0.05,E:0.4);"
    t = Tree(newick)

    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_axes([0, 0, 1, 1])
    draw_tree(ax, t)
    out = str(tmp_path / "tree_draw.png")
    fig.savefig(out)
    plt.close(fig)
    assert op.exists(out)
