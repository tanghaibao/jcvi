#!/usr/bin/env python
# -*- coding: UTF-8 -*-


"""
Script to plot diagrams of assembly graph in polyploids.
"""

import networkx as nx

from jcvi.utils.iter import pairwise
from jcvi.utils.webcolors import name_to_rgb


def sequence_to_graph(G, seq, color='black'):
    """
    Automatically construct graph given a sequence of characters.
    """
    seq = ["{}|{}".format(i, x) for i, x in enumerate(seq)]
    r, g, b = name_to_rgb(color)
    color = {'color': {'r': r, 'g': g, 'b': b, 'a': 0}}
    for x in seq:
        G.add_node(x)
        G.node[x]['viz'] = color
    for a, b in pairwise(seq):
        G.add_edge(a, b)
        G.edge[a][b]['viz'] = color


def main():
    seqA = "0" * 20
    seqB = "1" * 20
    G = nx.Graph()
    sequence_to_graph(G, seqA, color='green')
    sequence_to_graph(G, seqB, color='magenta')
    gexf = "graph.gexf"
    nx.write_gexf(G, gexf)


if __name__ == '__main__':
    main()
