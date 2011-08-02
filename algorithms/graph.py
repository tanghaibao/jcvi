#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for the common graph algorithms. Common usages are:

>>> edges = [(1, 2), (2, 3), (4, 5)]
>>> g = nx.DiGraph(edges)
>>> c = weakly_connected_components(g)
>>> print c
[[1, 2, 3], [4, 5]]
>>> sub = g.subgraph(c[0])
>>> topological_sort(sub)
[1, 2, 3]
"""

import networkx as nx

from networkx.algorithms.dag import topological_sort
from networkx.algorithms.components.weakly_connected import \
        weakly_connected_components
from networkx.algorithms.components.connected import connected_components
from networkx.algorithms.shortest_paths.generic import shortest_path


if __name__ == '__main__':
    import doctest
    doctest.testmod()
