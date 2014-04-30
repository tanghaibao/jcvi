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
            MINIMIZE, LPInstance, summation
from jcvi.utils.iter import flatten, pairwise


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
        subscaffolds = sorted((positions[k], k) for k in subscaffolds)
        dd, subscaffolds = zip(*subscaffolds)
        print i, subscaffolds
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

    #results, obj_val = sync_tsp(dummy_salesmen)
    results, obj_val = sync_tsp_gurobi(dummy_salesmen)
    if results:
        results = [x for x in results if DUMMY not in x]
        results = edges_to_path(results)
    return results, obj_val


def sync_tsp_gurobi(salesmen):
    from gurobipy import Model, GRB, quicksum

    all_edges = list(flatten(salesmen))
    incident, all_nodes = node_to_edge(all_edges, directed=False)
    idx = dict((n, i) for i, n in enumerate(all_nodes))
    nedges, n = len(all_edges), len(all_nodes)

    m = Model()

    step = lambda x: "u_{0}".format(x)
    # Create variables
    vars = {}
    for i, (a, b, w) in enumerate(all_edges):
        vars[i] = m.addVar(obj=w, vtype=GRB.BINARY, name=str(i))
    for u in all_nodes[1:]:
        u = step(u)
        vars[u] = m.addVar(obj=0, vtype=GRB.INTEGER, name=u)
    m.update()

    current_nedges = 0
    # Add degree constraint
    for edges in salesmen:
        incoming, outgoing, nodes = node_to_edge(edges)
        # For each node, select exactly 1 incoming and 1 outgoing edge
        for v in nodes:
            incoming_edges = [x + current_nedges for x in incoming[v]]
            outgoing_edges = [x + current_nedges for x in outgoing[v]]
            m.addConstr(quicksum(vars[x] for x in incoming_edges) == 1)
            m.addConstr(quicksum(vars[x] for x in outgoing_edges) == 1)
        current_nedges += len(edges)

    assert current_nedges == nedges

    # Bounds for step variables
    for u in nodes[1:]:
        u = step(u)
        vars[u].lb = 1
        vars[u].ub = n - 1

    # Subtour elimination - Miller-Tucker-Zemlin (MTZ) formulation
    u0 = nodes[0]
    current_nedges = 0
    for edges in salesmen:
        for i, (a, b, w) in enumerate(edges):
            if u0 in (a, b):
                continue
            a, b = step(a), step(b)
            na, nb, ne = vars[a], vars[b], vars[i + current_nedges]
            m.addConstr(na - nb + (n - 1) * ne <= n - 2)
        current_nedges += len(edges)
    m.update()

    # Given a list of edges, finds the shortest subtour
    def subtour(s_edges):
        visited = [False] * n
        cycles = []
        lengths = []
        selected = [[] for i in range(n)]
        for x, y in s_edges:
            selected[x].append(y)
        while True:
            current = visited.index(False)
            thiscycle = [current]
            while True:
                visited[current] = True
                neighbors = [x for x in selected[current] if not visited[x]]
                if len(neighbors) == 0:
                    break
                current = neighbors[0]
                thiscycle.append(current)
            cycles.append(thiscycle)
            lengths.append(len(thiscycle))
            if sum(lengths) == n:
                break
        return cycles[lengths.index(min(lengths))]

    def subtourelim(model, where):
        if where != GRB.callback.MIPSOL:
            return
        selected = []
        # make a list of edges selected in the solution
        sol = model.cbGetSolution([model._vars[i] for i in range(nedges)])
        selected = [all_edges[i] for i, x in enumerate(sol) if x > .5]
        # find the shortest cycle in the selected edge list
        current_nedges = 0
        for salesman in salesmen:
            incoming, outgoing, nodes = node_to_edge(edges)
            edge_store = dict(((a, b), i + current_nedges) \
                               for i, (a, b, w) in salesman)
            salesmen_selected = [(idx[a], idx[b]) for a, b, w in selected \
                                 if a in nodes and b in nodes]
            tour = subtour(salesmen_selected)
            if len(tour) == len(nodes):
                return
            # add a subtour elimination constraint
            c = tour
            incident = [edge_store[a, b] for a, b in pairwise(c + [c[0]])]
            model.cbLazy(quicksum(model._vars[x] for x in incident) <= len(tour) - 1)
            current_nedges += len(salesman)

    m._vars = vars
    #m.params.LazyConstraints = 1
    m.optimize()

    selected = [v.varName for v in m.getVars() if v.x > .5]
    selected = [int(x) for x in selected if x[:2] != "u_"]
    results = sorted(x for i, x in enumerate(all_edges) if i in selected) \
                    if selected else None
    return results, m.objVal


def sync_tsp(salesmen):
    """
    Calculates shortest cycle that traverses each node exactly once. Also known
    as the Traveling Salesman Problem (TSP).
    """
    L = LPInstance()
    all_edges = list(flatten(salesmen))
    incident, all_nodes = node_to_edge(all_edges, directed=False)
    nedges, nnodes = len(all_edges), len(all_nodes)

    L.add_objective(all_edges, objective=MINIMIZE)
    balance = []
    current_nedges = 0
    for edges in salesmen:
        incoming, outgoing, nodes = node_to_edge(edges)
        # For each node, select exactly 1 incoming and 1 outgoing edge
        for v in nodes:
            incoming_edges = [x + current_nedges for x in incoming[v]]
            outgoing_edges = [x + current_nedges for x in outgoing[v]]
            icc = summation(incoming_edges)
            occ = summation(outgoing_edges)
            balance.append("{0} = 1".format(icc))
            balance.append("{0} = 1".format(occ))
        current_nedges += len(edges)

    assert current_nedges == nedges

    # Subtour elimination - Miller-Tucker-Zemlin (MTZ) formulation
    mtz = []
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
            mtz.append(con_ab)
        current_nedges += len(edges)

    L.constraints = balance + mtz

    # Step variables u_i bound between 1 and n, as additional variables
    bounds = []
    for i in xrange(start_step, nedges + nnodes):
        bounds.append(" 1 <= x{0} <= {1}".format(i, nnodes - 1))
    L.bounds = bounds

    L.add_vars(nedges)
    L.add_vars(nnodes - 1, offset=start_step, binary=False)

    selected, obj_val = L.lpsolve(clean=False)
    results = sorted(x for i, x in enumerate(all_edges) if i in selected) \
                    if selected else None

    return results, obj_val


def main():
    salesmen, answer = make_data(200, 2)
    tour, val = sync_hamiltonian(salesmen)
    print "Solution found:", tour, val
    print "Truth:", answer


if __name__ == '__main__':
    main()
