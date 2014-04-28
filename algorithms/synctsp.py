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

from jcvi.algorithms.lpsolve import node_to_edge, edges_to_path, \
            MINIMIZE, BINARY, GENERNAL, LPInstance, summation
from jcvi.utils.iter import flatten


def make_data(POINTS, SALESMEN):
    """
    Make test data that has known solutions. All points lay on single line of
    length 1, and we simulate two salesmen. The score of the optimal solution is
    therefore 2.
    """
    scaffolds = ["P{0}".format(x) for x in xrange(POINTS)]
    positions = zip(scaffolds, np.random.rand(POINTS))
    positions += [("START", 0), ("END", 1)]
    positions = dict(positions)
    answer = sorted((v, k) for k, v in positions.items())
    dd, answer = zip(*answer)
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
        salesmen.append(edges)
    return salesmen, answer


def sync_hamiltonian(salesmen):
    """
    Salesmen is simply a list of edges. All edges are assumed to be directed, so
    caller is responsible for setting up distances in both directions.
    """
    DUMMY = "DUMMY"
    dummy_salesmen = []
    for edges in salesmen:
        incident, nodes = node_to_edge(edges, directed=False)
        dummy_edges = edges + [(DUMMY, x, 0) for x in nodes] + \
                              [(x, DUMMY, 0) for x in nodes]
        dummy_salesmen.append(dummy_edges)

    results, obj_val = sync_tsp(dummy_salesmen)
    if results:
        results = [x for x in results if DUMMY not in x]
        results = edges_to_path(results)
    return results, obj_val


def sync_tsp(salesmen):
    """
    Calculates shortest cycle that traverses each node exactly once. Also known
    as the Traveling Salesman Problem (TSP).
    """
    L = LPInstance()
    all_edges = list(flatten(salesmen))
    incident, all_nodes = node_to_edge(all_edges, directed=False)
    nedges, nnodes = len(all_edges), len(all_nodes)

    L.print_objective(all_edges, objective=MINIMIZE)
    constraints = []
    current_nedges = 0
    for edges in salesmen:
        incoming, outgoing, nodes = node_to_edge(edges)
        # For each node, select exactly 1 incoming and 1 outgoing edge
        for v in nodes:
            incoming_edges = [x + current_nedges for x in incoming[v]]
            outgoing_edges = [x + current_nedges for x in outgoing[v]]
            icc = summation(incoming_edges)
            occ = summation(outgoing_edges)
            constraints.append("{0} = 1".format(icc))
            constraints.append("{0} = 1".format(occ))
        current_nedges += len(edges)

    assert current_nedges == nedges

    # Subtour elimination - Miller-Tucker-Zemlin (MTZ) formulation
    start_step = nedges + 1
    u0 = all_nodes[0]
    nodes_to_steps = dict((n, start_step + i) for i, n in enumerate(all_nodes[1:]))
    current_nedges = 0
    for edges in salesmen:
        edge_store = dict((e[:2], i + current_nedges) for i, e in enumerate(edges))
        for a, b, w in edges:
            if u0 in (a, b):
                continue
            na, nb = nodes_to_steps[a], nodes_to_steps[b]
            i = edge_store[a, b]
            con_ab = " x{0} - x{1} + {2}x{3}".format(na, nb, nnodes - 1, i + 1)
            con_ab += " <= {0}".format(nnodes - 2)
            constraints.append(con_ab)
        current_nedges += len(edges)

    L.print_constraints(constraints)

    # Step variables u_i bound between 1 and n, as additional variables
    bounds = []
    for i in xrange(start_step, nedges + nnodes):
        bounds.append("1 <= x{0} <= {1}".format(i, nnodes - 1))
    L.print_bounds(bounds)

    L.print_vars(nedges, vars=BINARY)
    L.print_vars(nnodes - 1, offset=start_step, vars=GENERNAL)
    L.print_end()

    #print L.handle.getvalue()
    selected, obj_val = L.lpsolve(clean=False)
    results = sorted(x for i, x in enumerate(all_edges) if i in selected) \
                    if selected else None

    return results, obj_val


def main():
    salesmen, answer = make_data(100, 2)
    tour, val = sync_hamiltonian(salesmen)
    print "Solution found:", tour, val
    print "Truth:", answer


if __name__ == '__main__':
    main()
