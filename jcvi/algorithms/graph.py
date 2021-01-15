#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for the common graph algorithms.
"""
import sys
import logging

import networkx as nx
from collections import deque
from more_itertools import pairwise

from jcvi.formats.base import must_open


"""
Bidirectional graph.
"""
dirs = (">", "<")
trans = str.maketrans("+?-", ">><")


class BiNode(object):
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
            (e,) = L
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


class BiEdge(object):
    def __init__(self, v1, v2, o1, o2, color="black", length=None):
        o1 = o1.translate(trans)
        o2 = o2.translate(trans)
        assert o1 in dirs and o2 in dirs
        self.o1 = o1
        self.o2 = o2

        self.color = color
        self.length = length

    def __str__(self):
        return "".join(str(x) for x in (self.v1, self.o1, "--", self.o2, self.v2))

    def flip(self):
        self.v2, self.v1 = self.v1, self.v2
        o1, o2 = self.o1, self.o2
        self.o1 = ">" if o2 == "<" else "<"
        self.o2 = ">" if o1 == "<" else "<"


class BiGraph(object):
    def __init__(self):
        self.nodes = {}
        self.edges = {}

    def __str__(self):
        return "BiGraph with {0} nodes and {1} edges".format(
            len(self.nodes), len(self.edges)
        )

    def add_node(self, v):
        if v not in self.nodes:
            self.nodes[v] = BiNode(v)

    def add_edge(self, v1, v2, o1, o2, color="black", length=None):
        for v in (v1, v2):
            self.add_node(v)
        n1 = self.nodes.get(v1)
        n2 = self.nodes.get(v2)

        if (v1, v2) in self.edges or (v2, v1) in self.edges:
            return

        e = BiEdge(v1, v2, o1, o2, color=color, length=length)
        l = n1.outs if e.o1 == ">" else n1.ins
        r = n2.ins if e.o2 == ">" else n2.outs
        l.append(e)
        r.append(e)
        e.v1, e.v2 = n1, n2
        if v1 > v2:
            v1, v2 = v2, v1
            e.flip()
        self.edges[(v1, v2)] = e

    def get_node(self, v):
        return self.nodes[v]

    def get_edge(self, av, bv):
        flip = False
        if av > bv:
            av, bv = bv, av
            flip = True
        e = self.edges[(av, bv)]
        if flip:
            e.flip()
        return e

    def iter_paths(self):

        discovered = set()
        for v, vv in self.nodes.items():
            if v in discovered:
                continue

            path = deque([vv])

            # print "cur", v
            discovered.add(v)
            prev, ptag = vv.get_next(tag=">")
            while prev:
                # print prev, ptag
                if prev.v in discovered:
                    break
                path.appendleft(prev)
                discovered.add(prev.v)
                prev, ptag = prev.get_next(tag=ptag)

            next, ntag = vv.get_next(tag="<")
            while next:
                # print next, ntag
                if next.v in discovered:
                    break
                path.append(next)
                discovered.add(next.v)
                next, ntag = next.get_next(tag=ntag)

            # discovered |= set(x.v for x in path)
            yield path

    def path(self, path, flip=False):
        oo = []
        if len(path) == 1:
            m = "Singleton {0}".format(path[0])
            oo.append((path[0].v, True))
            return m, oo

        edges = []
        for a, b in pairwise(path):
            av, bv = a.v, b.v
            e = self.get_edge(av, bv)

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

    def read(self, filename, color="black"):
        fp = open(filename)
        nedges = 0
        for row in fp:
            a, b = row.strip().split("--")
            oa = a[-1]
            ob = b[0]
            a, b = a.strip("<>"), b.strip("<>")
            self.add_edge(a, b, oa, ob, color=color)
            nedges += 1
        logging.debug(
            "A total of {0} edges imported from `{1}` (color={2}).".format(
                nedges, filename, color
            )
        )

    def write(self, filename="stdout"):

        fw = must_open(filename, "w")
        for e in self.edges.values():
            print(e, file=fw)
        logging.debug("Graph written to `{0}`.".format(filename))

    def draw(
        self,
        pngfile,
        dpi=96,
        verbose=False,
        namestart=0,
        nodehighlight=None,
        prog="circo",
    ):
        import pygraphviz as pgv

        G = pgv.AGraph()
        for e in self.edges.values():
            arrowhead = e.o1 == ">"
            arrowtail = e.o2 == "<"
            if e.o1 != e.o2:  # Not sure why this is necessary
                arrowhead = not arrowhead
                arrowtail = not arrowtail
            arrowhead = "normal" if arrowhead else "inv"
            arrowtail = "normal" if arrowtail else "inv"
            v1, v2 = e.v1, e.v2
            v1, v2 = str(v1)[namestart:], str(v2)[namestart:]
            G.add_edge(v1, v2, color=e.color, arrowhead=arrowhead, arrowtail=arrowtail)

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

    def get_next(self, node, tag="<"):
        return self.get_node(node).get_next(tag)

    def get_path(self, n1, n2, tag="<"):
        # return all intermediate nodes on path n1 -> n2
        path = deque()
        next, ntag = self.get_next(n1, tag=tag)
        while next:
            if next.v == n2:
                return path
            path.append((next, ntag))
            next, ntag = next.get_next(tag=ntag)
        return path if n2 is None else None


def graph_stats(G, diameter=False):
    logging.debug("Graph stats: |V|={0}, |E|={1}".format(len(G), G.size()))
    if diameter:
        d = max(nx.diameter(H) for H in nx.connected_component_subgraphs(G))
        logging.debug("Graph diameter: {0}".format(d))


def graph_local_neighborhood(G, query, maxdegree=10000, maxsize=10000):
    c = [k for k, d in G.degree().iteritems() if d > maxdegree]
    if c:
        logging.debug("Remove {0} nodes with deg > {1}".format(len(c), maxdegree))
    G.remove_nodes_from(c)

    logging.debug("BFS search from {0}".format(query))

    queue = set(query)
    # BFS search of max depth
    seen = set(query)
    coresize = len(query)
    depth = 0
    while True:
        neighbors = set()
        for q in queue:
            if q not in G:
                continue
            neighbors |= set(G.neighbors(q))
        queue = neighbors - seen
        if not queue:
            break

        if len(seen | queue) > maxsize + coresize:
            break

        seen |= queue
        print(
            "iter: {0}, graph size={1} ({2} excluding core)".format(
                depth, len(seen), len(seen) - coresize
            ),
            file=sys.stderr,
        )
        depth += 1

    return G.subgraph(seen)


def graph_simplify(G):
    """
    Simplify big graphs: remove spurs and contract unique paths.
    """
    spurs = []
    path_nodes = []
    for k, d in G.degree().iteritems():
        if d == 1:
            spurs.append(k)
        elif d == 2:
            path_nodes.append(k)

    logging.debug("Remove {0} spurs.".format(len(spurs)))
    G.remove_nodes_from(spurs)

    SG = G.subgraph(path_nodes)
    cc = nx.connected_components(SG)
    for c in cc:
        if len(c) == 1:
            continue
        c = set(c)
        neighbors = set()
        for x in c:
            neighbors |= set(G.neighbors(x))
        neighbors -= c
        newtag = list(c)[0] + "*"
        for n in neighbors:
            G.add_edge(newtag, n)
        G.remove_nodes_from(c)
    logging.debug(
        "Contract {0} path nodes into {1} nodes.".format(len(path_nodes), len(cc))
    )


def bigraph_test():
    g = BiGraph()
    g.add_edge(1, 2, ">", "<")
    g.add_edge(2, 3, "<", "<", color="red")
    g.add_edge(2, 3, ">", ">", color="blue")
    g.add_edge(5, 3, ">", ">")
    g.add_edge(4, 3, "<", ">")
    g.add_edge(4, 6, ">", ">")
    g.add_edge(7, 1, ">", ">")
    g.add_edge(7, 5, "<", ">")
    g.add_edge(8, 6, ">", "<")
    print(g)
    g.write()
    for path in g.iter_paths():
        p, oo = g.path(path)
        print(p)
        print(oo)

    # g.draw("demo.png", verbose=True)


def update_weight(G, a, b, w):
    if G.has_edge(a, b):  # Parallel edges found!
        G[a][b]["weight"] += w
    else:
        G.add_edge(a, b, weight=w)


def make_paths(paths, weights=None):
    """
    Zip together paths. Called by merge_paths().
    """
    npaths = len(paths)
    weights = weights or [1] * npaths
    assert len(paths) == len(weights)

    G = nx.DiGraph()
    for path, w in zip(paths, weights):
        for a, b in pairwise(path):
            update_weight(G, a, b, w)
    return G


def reduce_paths(G):
    """
    Make graph into a directed acyclic graph (DAG).
    """
    from jcvi.algorithms.lpsolve import min_feedback_arc_set

    while not nx.is_directed_acyclic_graph(G):
        edges = []
        for a, b, w in G.edges_iter(data=True):
            w = w["weight"]
            edges.append((a, b, w))
        mf, mf_score = min_feedback_arc_set(edges)
        for a, b, w in mf:
            G.remove_edge(a, b)

    assert nx.is_directed_acyclic_graph(G)
    G = transitive_reduction(G)
    return G


def draw_graph(G, pngfile, prog="dot"):
    G = nx.to_agraph(G)
    G.draw(pngfile, prog=prog)
    logging.debug("Graph written to `{0}`.".format(pngfile))


def transitive_reduction(G):
    """
    Returns a transitive reduction of a graph.  The original graph
    is not modified.

    A transitive reduction H of G has a path from x to y if and
    only if there was a path from x to y in G.  Deleting any edge
    of H destroys this property.  A transitive reduction is not
    unique in general.  A transitive reduction has the same
    transitive closure as the original graph.

    A transitive reduction of a complete graph is a tree.  A
    transitive reduction of a tree is itself.

    >>> G = nx.DiGraph([(1, 2), (1, 3), (2, 3), (2, 4), (3, 4)])
    >>> H = transitive_reduction(G)
    >>> H.edges()
    [(1, 2), (2, 3), (3, 4)]
    """
    H = G.copy()
    for a, b, w in G.edges_iter(data=True):
        # Try deleting the edge, see if we still have a path
        # between the vertices
        H.remove_edge(a, b)
        if not nx.has_path(H, a, b):  # we shouldn't have deleted it
            H.add_edge(a, b, w)
    return H


def merge_paths(paths, weights=None):
    """
    Zip together sorted lists.

    >>> paths = [[1, 2, 3], [1, 3, 4], [2, 4, 5]]
    >>> G = merge_paths(paths)
    >>> nx.topological_sort(G)
    [1, 2, 3, 4, 5]
    >>> paths = [[1, 2, 3, 4], [1, 2, 3, 2, 4]]
    >>> G = merge_paths(paths, weights=(1, 2))
    >>> nx.topological_sort(G)
    [1, 2, 3, 4]
    """
    G = make_paths(paths, weights=weights)
    G = reduce_paths(G)
    return G


def longest_path_weighted_nodes(G, source, target, weights=None):
    """
    The longest path problem is the problem of finding a simple path of maximum
    length in a given graph. While for general graph, this problem is NP-hard,
    but if G is a directed acyclic graph (DAG), longest paths in G can be found
    in linear time with dynamic programming.

    >>> G = nx.DiGraph([(1, 2), (1, 3), (2, "M"), (3, "M")])
    >>> longest_path_weighted_nodes(G, 1, "M", weights={1: 1, 2: 1, 3: 2, "M": 1})
    ([1, 3, 'M'], 4)
    """
    assert nx.is_directed_acyclic_graph(G)

    tree = nx.topological_sort(G)
    node_to_index = dict((t, i) for i, t in enumerate(tree))

    nnodes = len(tree)
    weights = [weights.get(x, 1) for x in tree] if weights else [1] * nnodes
    score, fromc = weights[:], [-1] * nnodes
    si = node_to_index[source]
    ti = node_to_index[target]
    for a in tree[si:ti]:
        ai = node_to_index[a]
        for b, w in G[a].items():
            bi = node_to_index[b]
            w = w.get("weight", 1)
            d = score[ai] + weights[bi] * w  # Favor heavier edges
            if d <= score[bi]:
                continue
            score[bi] = d  # Update longest distance so far
            fromc[bi] = ai

    # Backtracking
    path = []
    while ti != -1:
        path.append(ti)
        ti = fromc[ti]

    path = [tree[x] for x in path[::-1]]
    return path, score[ti]


if __name__ == "__main__":
    import doctest

    doctest.testmod()
    # bigraph_test()
