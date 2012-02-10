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

import sys
import logging

import networkx as nx
from collections import deque

from networkx.algorithms.dag import topological_sort
from networkx.algorithms.components.weakly_connected import \
        weakly_connected_components
from networkx.algorithms.components.connected import connected_components
from networkx.algorithms.shortest_paths.generic import shortest_path

from jcvi.formats.base import must_open
from jcvi.apps.base import debug
debug()

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

    def __init__(self, v1, v2, o1, o2, color="black"):

        self.v1 = v1
        self.v2 = v2

        assert o1 in dirs and o2 in dirs
        self.o1 = o1
        self.o2 = o2

        if v1 > v2:
            self.flip()

        self.color = color

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

            #print "cur", v
            discovered.add(v)
            prev, ptag = vv.get_next(tag=">")
            while prev:
                #print prev, ptag
                if prev.v in discovered:
                    break
                path.appendleft(prev)
                discovered.add(prev.v)
                prev, ptag = prev.get_next(tag=ptag)

            next, ntag = vv.get_next(tag="<")
            while next:
                #print next, ntag
                if next.v in discovered:
                    break
                path.append(next)
                discovered.add(next.v)
                next, ntag = next.get_next(tag=ntag)

            #discovered |= set(x.v for x in path)
            yield path

    def path(self, path):
        from jcvi.utils.iter import pairwise

        oo = []
        if len(path) == 1:
            m = "Singleton {0}".format(path[0])
            oo.append((path[0].v, True))
            return m, oo

        edges = []
        for a, b in pairwise(path):
            av, bv = a.v, b.v
            flip = False
            if av > bv:
                av, bv = bv, av
                flip = True
            e = self.edges[(av, bv)]
            if flip:
                e.flip()

            if not oo:  # First edge imports two nodes
                oo.append((e.v1.v, e.o1 == ">"))
            last = oo[-1]
            assert last == (e.v1.v, e.o1 == ">")
            oo.append((e.v2.v, e.o2 == ">"))

            if flip:
                se = str(e)
                e.flip()
            else:
                se = str(e)
            edges.append(se)

        return "|".join(edges), oo

    def __str__(self):
        return "BiGraph with {0} nodes and {1} edges".\
                format(len(self.nodes), len(self.edges))

    def read(self, filename, color="black"):
        fp = open(filename)
        nedges = 0
        for row in fp:
            a, b = row.strip().split("--")
            oa = a[-1]
            ob = b[0]
            a, b = a.strip("<>"), b.strip("<>")
            self.add_edge(BiEdge(a, b, oa, ob, color=color))
            nedges += 1
        logging.debug("A total of {0} edges imported from `{1}` (color={2}).".
                      format(nedges, filename, color))

    def write(self, filename="stdout"):

        fw = must_open(filename, "w")
        for e in self.edges.values():
            print >> fw, e
        logging.debug("Graph written to `{0}`.".format(filename))

    def draw(self, pngfile, dpi=96, verbose=False, namestart=0,
                nodehighlight=None, prog="circo"):
        import pygraphviz as pgv

        G = pgv.AGraph()
        for e in self.edges.values():
            arrowhead = (e.o1 == ">")
            arrowtail = (e.o2 == "<")
            if e.o1 != e.o2:  # Not sure why this is necessary
                arrowhead = not arrowhead
                arrowtail = not arrowtail
            arrowhead = "normal" if arrowhead else "inv"
            arrowtail = "normal" if arrowtail else "inv"
            v1, v2 = e.v1, e.v2
            v1, v2 = str(v1)[namestart:], str(v2)[namestart:]
            G.add_edge(v1, v2, color=e.color,
                       arrowhead=arrowhead, arrowtail=arrowtail)

        if nodehighlight:
            for n in nodehighlight:
                n = n[namestart:]
                n = G.get_node(n)
                n.attr["shape"] = "box"

        G.graph_attr.update(dpi=str(dpi))
        if verbose:
            G.write(sys.stderr)
        G.draw(pngfile, prog=prog)
        logging.debug("Graph written to `{0}`.".format(pngfile))


def bigraph_test():
    g = BiGraph()
    g.add_edge(BiEdge(1, 2, ">", "<"))
    g.add_edge(BiEdge(2, 3, "<", "<", color="red"))
    g.add_edge(BiEdge(2, 3, ">", ">", color="blue"))
    g.add_edge(BiEdge(5, 3, ">", ">"))
    g.add_edge(BiEdge(4, 3, "<", ">"))
    g.add_edge(BiEdge(4, 6, ">", ">"))
    g.add_edge(BiEdge(7, 1, ">", ">"))
    g.add_edge(BiEdge(7, 5, "<", ">"))
    g.add_edge(BiEdge(8, 6, ">", "<"))
    print g
    g.write()
    for path in g.iter_paths():
        p, oo = g.path(path)
        print p
        print oo

    g.draw("demo.png", verbose=True)


if __name__ == '__main__':
    import doctest
    doctest.testmod()

    bigraph_test()
