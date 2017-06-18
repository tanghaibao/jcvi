#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Script to plot diagrams of assembly graph in polyploids.
"""

import networkx as nx

from jcvi.utils.iter import pairwise
from jcvi.utils.webcolors import name_to_rgb


def make_sequence(seq, name="S"):
    """
    Make unique nodes for sequence graph.
    """
    return ["{}_{}_{}".format(name, i, x) for i, x in enumerate(seq)]


def make_rgb_color(color):
    """
    Convert web color name into a dictionary of color compatible with GEXF.
    """
    r, g, b = name_to_rgb(color)
    color = {'color': {'r': r, 'g': g, 'b': b, 'a': 0}}
    return color


def sequence_to_graph(G, seq, color='black'):
    """
    Automatically construct graph given a sequence of characters.
    """
    color = make_rgb_color(color)
    for x in seq:
        G.add_node(x)
        G.node[x]['viz'] = color
    for a, b in pairwise(seq):
        G.add_edge(a, b)
        G.edge[a][b]['viz'] = color


def zip_sequence(G, seqA, seqB, color="white"):
    """
    Fuse certain nodes together, if they contain same data except for the
    sequence name.
    """
    color = make_rgb_color(color)
    for a, b in zip(seqA, seqB):
        if a.split('_', 1)[1] == b.split('_', 1)[1]:
            G.add_edge(a, b)
            G.edge[a][b]['viz'] = color


def main():
    # Build fake data
    seqA = list("0" * 20)
    seqB = list("1" * 20)
    seqB[1] = "0"
    seqB[10] = "0"
    seqA = make_sequence(seqA, name="A")
    seqB = make_sequence(seqB, name="B")

    # Build graph structure
    G = nx.Graph()
    sequence_to_graph(G, seqA, color='green')
    sequence_to_graph(G, seqB, color='magenta')
    zip_sequence(G, seqA, seqB)

    # Output graph
    gexf = "graph.gexf"
    nx.write_gexf(G, gexf)


if __name__ == '__main__':
    main()
