#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Script to plot diagrams of assembly graph in polyploids.
"""

from itertools import combinations
from random import choice
from graphviz import Digraph

from jcvi.utils.iter import pairwise


def make_sequence(seq, name="S"):
    """
    Make unique nodes for sequence graph.
    """
    return ["{}_{}_{}".format(name, i, x) for i, x in enumerate(seq)]


def sequence_to_graph(G, seq, color='black'):
    """
    Automatically construct graph given a sequence of characters.
    """
    with G.subgraph(name=color) as c:
        c.attr('node', color=color)
        c.attr('edge', color=color)
        for x in seq:
            c.node(x)
        for a, b in pairwise(seq):
            c.edge(a, b)


def zip_sequence(G, seqA, seqB, color="lightgray"):
    """
    Fuse certain nodes together, if they contain same data except for the
    sequence name.
    """
    for a, b in zip(seqA, seqB):
        if a.split('_', 1)[1] != b.split('_', 1)[1]:
            continue
        G.edge(a, b, color=color)


def main():
    size = 30
    indices = range(size)
    # Build fake data
    seqA = list("0" * size)
    seqB = seqA[:]
    seqC = seqA[:]
    allseqs = [seqA, seqB, seqC]
    for s in allseqs:
        for i in [choice(indices) for x in range(5)]:
            s[i] = "1"

    allseqs = [make_sequence(s, name=name) for (s, name) in \
                zip(allseqs, "ABC")]

    # Build graph structure
    G = Digraph("Assembly graph", filename="graph")
    G.attr(rankdir='LR')
    G.attr('node', shape='point')
    G.attr('edge', dir='none', penwidth='4')

    colors = ['green', 'magenta', 'lightslategray']
    for s, color in zip(allseqs, colors):
        sequence_to_graph(G, s, color=color)
    for sA, sB in combinations(allseqs, 2):
        zip_sequence(G, sA, sB)

    # Output graph
    G.view()


if __name__ == '__main__':
    main()
