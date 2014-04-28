#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Multiple TSPs problem are cases where several salesmen travel together. As they
travel between cities, they each incur different cost. In addition, they visit
different sets cities, but are constrained by a common visiting order. For
illustration:

            A B C D E F G H
Salesman 1  x   x   x x   x
Salesman 2    x x x     x x

The goal is to minimize total combined cost for all salesmen. This is sometimes
also called "synchronized TSP". In the above example, final cost is:

Salesman 1 cost: AC + CE + EF + FH
Salesman 2 cost: BC + CD + DG + GH
"""

import numpy as np
from random import sample
from itertools import combinations

from jcvi.algorithms.lpsolve import hamiltonian
from jcvi.algorithms.lpsolve import node_to_edge, \
            edges_to_path, MINIMIZE, BINARY, GENERNAL, LPInstance, summation


def make_mtsp_data(POINTS, SALESMEN):
    """
    Make test data that has known solutions. All points lay on single line of
    length 1, and we simulate two salesmen. The score of the optimal solution is
    therefore 2.
    """
    scaffolds = ["P{0}".format(x) for x in xrange(POINTS)]
    positions = zip(scaffolds, np.random.rand(POINTS))
    positions += [("START", 0), ("END", 1)]
    positions = dict(positions)
    print sorted((v, k) for k, v in positions.items())
    salesmen = []
    for i in xrange(SALESMEN):
        subscaffolds = sample(scaffolds, POINTS * 3 / 5)
        edges = []
        for a, b in combinations(subscaffolds, 2):
            d = abs(positions[a] - positions[b])
            edges.append((a, b, d))
            edges.append((b, a, d))
        for p in subscaffolds:
            d = positions[p]
            edges.append(("START", p, d))
            edges.append((p, "END", 1 - d))
        print hamiltonian(edges, directed=True)
        salesmen.append(edges)
    return salesmen


def mhamiltonian(edges):
    incident, nodes = node_to_edge(edges, directed=False)
    DUMMY = "DUMMY"
    dummy_edges = edges + [(DUMMY, x, 0) for x in nodes]
    # Make graph symmetric
    all_edges = dummy_edges[:]
    for e in dummy_edges:  # flip source and link
        new_edge = tuple([e[1], e[0]] + list(e[2:]))
        all_edges.append(new_edge)

    results, obj_val = mtsp(all_edges)
    if results:
        results = [x for x in results if DUMMY not in x]
        results = edges_to_path(results)
    return results, obj_val


def mtsp(edges):
    """
    Calculates shortest cycle that traverses each node exactly once. Also known
    as the Traveling Salesman Problem (TSP).
    """
    incoming, outgoing, nodes = node_to_edge(edges)

    nedges, nnodes = len(edges), len(nodes)
    L = LPInstance()

    L.print_objective(edges, objective=MINIMIZE)
    constraints = []
    # For each node, select exactly 1 incoming and 1 outgoing edge
    for v in nodes:
        incoming_edges = incoming[v]
        outgoing_edges = outgoing[v]
        icc = summation(incoming_edges)
        occ = summation(outgoing_edges)
        constraints.append("{0} = 1".format(icc))
        constraints.append("{0} = 1".format(occ))

    # Subtour elimination - Miller-Tucker-Zemlin (MTZ) formulation
    # <http://en.wikipedia.org/wiki/Travelling_salesman_problem>
    # Desrochers and laporte, 1991 (DFJ) has a stronger constraint
    # See also:
    # G. Laporte / The traveling salesman problem: Overview of algorithms
    start_step = nedges + 1
    u0 = nodes[0]
    nodes_to_steps = dict((n, start_step + i) for i, n in enumerate(nodes[1:]))
    edge_store = dict((e[:2], i) for i, e in enumerate(edges))
    for i, e in enumerate(edges):
        a, b = e[:2]
        if u0 in (a, b):
            continue
        na, nb = nodes_to_steps[a], nodes_to_steps[b]
        con_ab = " x{0} - x{1} + {2}x{3}".format(na, nb, nnodes - 1, i + 1)
        if (b, a) in edge_store:  # This extra term is the stronger DFJ formulation
            j = edge_store[(b, a)]
            con_ab += " + {0}x{1}".format(nnodes - 3, j + 1)
        con_ab += " <= {0}".format(nnodes - 2)
        constraints.append(con_ab)

    L.print_constraints(constraints)

    # Step variables u_i bound between 1 and n, as additional variables
    bounds = []
    for i in xrange(start_step, nedges + nnodes):
        bounds.append("1 <= x{0} <= {1}".format(i, nnodes - 1))
    L.print_bounds(bounds)

    L.print_vars(nedges, vars=BINARY)
    L.print_vars(nnodes - 1, offset=start_step, vars=GENERNAL)
    L.print_end()

    selected, obj_val = L.lpsolve()
    results = sorted(x for i, x in enumerate(edges) if i in selected) \
                    if selected else None

    return results, obj_val


def main():
    make_mtsp_data(10, 2)


if __name__ == '__main__':
    main()
