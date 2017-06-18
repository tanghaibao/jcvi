#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Script to plot diagrams of assembly graph in polyploids.
"""

import networkx as nx
from itertools import combinations
from random import choice

from jcvi.utils.iter import pairwise
from jcvi.utils.webcolors import name_to_rgb


def make_rgb_color(color):
    """
    Convert web color name into a dictionary of color compatible with GEXF.
    """
    r, g, b = name_to_rgb(color)
    return {'r': r, 'g': g, 'b': b, 'a': 0}


def make_sequence(seq, name="S"):
    """
    Make unique nodes for sequence graph.
    """
    return ["{}_{}_{}".format(name, i, x) for i, x in enumerate(seq)]


def make_node_viz(color):
    """
    Make GEXF-compatible node viz.
    """
    return {'color': make_rgb_color(color)}


def make_edge_viz(color, thickness=4):
    """
    Make GEXF-compatible edge viz.
    """
    return {'color': make_rgb_color(color),
            'thickness': thickness}


def sequence_to_graph(G, seq, color='black'):
    """
    Automatically construct graph given a sequence of characters.
    """
    node_viz = make_node_viz(color)
    edge_viz = make_edge_viz(color)
    for x in seq:
        G.add_node(x)
        G.node[x]['viz'] = node_viz
    for a, b in pairwise(seq):
        G.add_edge(a, b)
        G.edge[a][b]['viz'] = edge_viz


def zip_sequence(G, seqA, seqB, color="white"):
    """
    Fuse certain nodes together, if they contain same data except for the
    sequence name.
    """
    edge_viz = make_edge_viz(color)
    for a, b in zip(seqA, seqB):
        if a.split('_', 1)[1] != b.split('_', 1)[1]:
            continue
        G.add_edge(a, b, weight=2)
        G.edge[a][b]['viz'] = edge_viz


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
    G = nx.Graph()
    colors = ['green', 'magenta', 'lightslategray']
    for s, color in zip(allseqs, colors):
        sequence_to_graph(G, s, color=color)
    for sA, sB in combinations(allseqs, 2):
        zip_sequence(G, sA, sB)

    # Output graph
    gexf = "graph.gexf"
    nx.write_gexf(G, gexf)


if __name__ == '__main__':
    main()
