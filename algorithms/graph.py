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
from collections import deque

from networkx.algorithms.dag import topological_sort
from networkx.algorithms.components.weakly_connected import \
        weakly_connected_components
from networkx.algorithms.components.connected import connected_components
from networkx.algorithms.shortest_paths.generic import shortest_path

"""
Bidirectional graph.
"""
dirs = (">", "<")


class BiNode (object):

    def __init__(self, v):
        self.v = v
        self.ins = []
        self.outs = []

    def get_next(self, tag="<"):
        """
        This function is tricky and took me a while to figure out.

        The tag specifies the direction where the current edge came from.

         tag     ntag
        ---> V >----> U
             cur      next

        This means the next vertex should follow the outs since this tag is
        inward '<'. Check if there are multiple branches if len(L) == 1, and
        also check if the next it finds has multiple incoming edges though if
        len(B) == 1.
        """
        next, ntag = None, None

        L = self.outs if tag == "<" else self.ins

        if len(L) == 1:
            e, = L
            if e.v1.v == self.v:
                next, ntag = e.v2, e.o2
                ntag = "<" if ntag == ">" else ">"  # Flip tag if on other end
            else:
                next, ntag = e.v1, e.o1

        if next:  # Validate the next vertex
            B = next.ins if ntag == "<" else next.outs
            if len(B) > 1:
                return None, None

        return next, ntag

    def __str__(self):
        return str(self.v)

    __repr__ = __str__


class BiEdge (object):

    def __init__(self, v1, v2, o1, o2):

        self.v1 = v1
        self.v2 = v2

        assert o1 in dirs and o2 in dirs
        self.o1 = o1
        self.o2 = o2

        if v1 > v2:
            self.flip()

    def __str__(self):
        return "".join(str(x) for x in \
                 (self.v1, self.o1, "--", self.o2, self.v2))

    def flip(self):
        self.v2, self.v1 = self.v1, self.v2
        o1, o2 = self.o1, self.o2
        self.o1 = ">" if o2 == "<" else "<"
        self.o2 = ">" if o1 == "<" else "<"


class BiGraph (object):

    def __init__(self):
        self.nodes = {}
        self.edges = {}

    def add_edge(self, e):
        v1, v2 = e.v1, e.v2

        assert isinstance(e, BiEdge)
        for v in (v1, v2):
            if v not in self.nodes:
                self.nodes[v] = BiNode(v)
        n1 = self.nodes.get(v1)
        n2 = self.nodes.get(v2)
        l = n1.outs if e.o1 == ">" else n1.ins
        r = n2.ins if e.o2 == ">" else n2.outs
        l.append(e)
        r.append(e)
        e.v1, e.v2 = n1, n2
        self.edges[(v1, v2)] = e

    def iter_paths(self):

        discovered = set()
        for v, vv in self.nodes.items():
            if v in discovered:
                continue

            path = deque([vv])
            #print "cur", vv

            prev, ptag = vv.get_next(tag=">")
            while prev:
                #print "prev", prev, ptag
                if prev.v in discovered:
                    break
                path.appendleft(prev)
                prev, ptag = prev.get_next(tag=ptag)

            next, ntag = vv.get_next(tag="<")
            while next:
                #print "next", next, ntag
                if next.v in discovered:
                    break
                path.append(next)
                next, ntag = next.get_next(tag=ntag)

            discovered |= set(x.v for x in path)
            yield path

    def path(self, path):
        from jcvi.utils.iter import pairwise

        oo = []
        isCurrentPlusOrientation = True
        oo.append((path[0], isCurrentPlusOrientation))
        if len(path) == 1:
            m = "Singleton {0}".format(path[0])
            return m, oo

        edges = []
        for a, b in pairwise(path):
            av, bv = a.v, b.v
            flip = False
            if av > bv:
                av, bv = bv, av
                flip = True
            e = self.edges[(av, bv)]
            if e.o1 != e.o2:
                isCurrentPlusOrientation = not isCurrentPlusOrientation
            oo.append((b.v, isCurrentPlusOrientation))

            if flip:
                e.flip()
                se = str(e)
                e.flip()
            else:
                se = str(e)
            edges.append(se)

        return "|".join(edges), oo

    def __str__(self):
        return "BiGraph with {0} nodes and {1} edges".\
                format(len(self.nodes), len(self.edges))


def bigraph_test():
    g = BiGraph()
    g.add_edge(BiEdge(1, 2, ">", "<"))
    g.add_edge(BiEdge(2, 3, "<", "<"))
    g.add_edge(BiEdge(5, 3, ">", ">"))
    g.add_edge(BiEdge(4, 3, "<", ">"))
    g.add_edge(BiEdge(4, 6, ">", ">"))
    g.add_edge(BiEdge(7, 1, ">", ">"))
    g.add_edge(BiEdge(7, 5, "<", ">"))
    g.add_edge(BiEdge(8, 6, ">", "<"))
    print g
    for e in g.edges.values():
        print e
    for path in g.iter_paths():
        p, oo = g.path(path)
        print p
        print oo


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    bigraph_test()
