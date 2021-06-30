#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Script to plot diagrams of assembly graph in polyploids.
"""

from collections import defaultdict
from random import choice, sample

from brewer2mpl import get_map
from graphviz import Graph
from matplotlib.colors import to_hex
from more_itertools import pairwise


def make_sequence(seq, name="S"):
    """
    Make unique nodes for sequence graph.
    """
    return ["{}_{}_{}".format(name, i, x) for i, x in enumerate(seq)]


def sequence_to_graph(G, seq, color="black"):
    """
    Automatically construct graph given a sequence of characters.
    """
    for x in seq:
        if x.endswith("_1"):  # Mutation
            G.node(x, color=color, width="0.1", shape="circle", label="")
        else:
            G.node(x, color=color)
    for a, b in pairwise(seq):
        G.edge(a, b, color=color)


def zip_sequences(G, allseqs):
    """
    Fuse certain nodes together, if they contain same data except for the
    sequence name.
    """
    for s in zip(*allseqs):
        groups = defaultdict(list)
        for x in s:
            part = x.split("_", 1)[1]
            groups[part].append(x)
        for part, g in groups.items():
            with G.subgraph(name="cluster_" + part) as c:
                for x in g:
                    c.node(x)
                c.attr(style="invis")


def main():
    SIZE = 20
    PLOIDY = 2
    MUTATIONS = 2

    indices = range(SIZE)
    # Build fake data
    seqA = list("0" * SIZE)
    allseqs = [seqA[:] for x in range(PLOIDY)]  # Hexaploid
    for s in allseqs:
        for i in [choice(indices) for x in range(MUTATIONS)]:
            s[i] = "1"

    allseqs = [
        make_sequence(s, name=name)
        for (s, name) in zip(allseqs, [str(x) for x in range(PLOIDY)])
    ]

    # Build graph structure
    G = Graph("Assembly graph", filename="graph")
    G.attr(rankdir="LR", fontname="Helvetica", splines="true")
    G.attr(ranksep=".2", nodesep="0.02")
    G.attr("node", shape="point")
    G.attr("edge", dir="none", penwidth="4")

    colorset = get_map("Set2", "qualitative", 8).mpl_colors
    colorset = [to_hex(x) for x in colorset]
    colors = sample(colorset, PLOIDY)
    for s, color in zip(allseqs, colors):
        sequence_to_graph(G, s, color=color)
    zip_sequences(G, allseqs)

    # Output graph
    G.view()


if __name__ == "__main__":
    main()
