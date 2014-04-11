#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
TSP solver using Concorde. This is much faster than the LP-formulation in
algorithms.lpsolve.tsp().
"""

import os.path as op
import os
import logging

import numpy as np

from jcvi.formats.base import must_open
from jcvi.algorithms.lpsolve import populate_edge_weights, node_to_edge
from jcvi.apps.base import mkdir, debug, which, sh
debug()

NEG_INF = -1000
Work_dir = "tsp_work"


class Concorde (object):

    def __init__(self, edges, work_dir=Work_dir, clean=False, verbose=False,
                       precision=0):

        self.work_dir = work_dir
        self.clean = clean
        self.verbose = verbose

        mkdir(work_dir)
        tspfile = op.join(work_dir, "data.tsp")
        self.print_to_tsplib(edges, tspfile, precision=precision)
        retcode, outfile = self.run_concorde(tspfile)
        self.tour = self.parse_output(outfile)

    def print_to_tsplib(self, edges, tspfile, precision=0):
        """
        See TSPlib format:
        <https://www.iwr.uni-heidelberg.de/groups/comopt/software/TSPLIB95/>

        NAME: bayg29
        TYPE: TSP
        COMMENT: 29 Cities in Bavaria, geographical distances
        DIMENSION: 29
        EDGE_WEIGHT_TYPE: EXPLICIT
        EDGE_WEIGHT_FORMAT: UPPER_ROW
        DISPLAY_DATA_TYPE: TWOD_DISPLAY
        EDGE_WEIGHT_SECTION
        (... numbers ...)
        """
        fw = must_open(tspfile, "w")
        incident, nodes = node_to_edge(edges, directed=False)
        self.nodes = nodes
        nodes_indices = dict((n, i) for i, n in enumerate(nodes))
        self.nnodes = nnodes = len(nodes)

        # TSPLIB requires explicit weights to be integral, and non-negative
        weights = [x[-1] for x in edges]
        max_x, min_x = max(weights), min(weights)
        inf = 2 * max(abs(max_x), abs(min_x))
        factor = 10 ** precision

        print >> fw, "NAME: data"
        print >> fw, "TYPE: TSP"
        print >> fw, "DIMENSION: {0}".format(nnodes)

        D = np.ones((nnodes, nnodes), dtype=float) * inf
        for a, b, w in edges:
            ia, ib = nodes_indices[a], nodes_indices[b]
            D[ia, ib] = D[ib, ia] = w
        D = (D - min_x) * factor
        D = D.astype(int)

        print >> fw, "EDGE_WEIGHT_TYPE: EXPLICIT"
        print >> fw, "EDGE_WEIGHT_FORMAT: FULL_MATRIX"
        print >> fw, "EDGE_WEIGHT_SECTION"
        for row in D:  # Dump the full matrix
            print >> fw, " " + " ".join(str(x) for x in row)

        print >> fw, "EOF"
        fw.close()
        logging.debug("Write TSP instance to `{0}`".format(tspfile))

    def run_concorde(self, tspfile):
        outfile = op.join(self.work_dir, "data.sol")
        if op.exists(outfile):
            os.remove(outfile)

        cc = "concorde"
        assert which(cc), "You must install `concorde` on your PATH" + \
                          " [http://www.math.uwaterloo.ca/tsp/concorde.html]"
        cmd = "{0} -x -o {1} {2}".format(cc, outfile, tspfile)

        outf = None if self.verbose else "/dev/null"
        retcode = sh(cmd, outfile=outf, errfile=outf)
        return retcode, outfile

    def parse_output(self, outfile):
        fp = open(outfile)
        dimension = int(fp.next().strip())  # header
        assert dimension == self.nnodes
        tour = []
        for row in fp:
            tour += [int(x) for x in row.split()]
        tour = [self.nodes[x] for x in tour]
        return tour


def hamiltonian(edges, symmetric=True, precision=0):
    """
    Calculates shortest path that traverses each node exactly once. Convert
    Hamiltonian path problem to TSP by adding one dummy point that has a distance
    of zero to all your other points. Solve the TSP and get rid of the dummy
    point - what remains is the Hamiltonian Path.

    >>> g = [(1,2), (2,3), (3,4), (4,2), (3,5)]
    >>> hamiltonian(g)
    [1, 2, 4, 3, 5]
    >>> hamiltonian([(1, 2), (2, 3)], symmetric=False)
    [1, 2, 3]
    """
    edges = populate_edge_weights(edges)
    incident, nodes = node_to_edge(edges, directed=False)
    DUMMY = "DUMMY"
    dummy_edges = edges + [(DUMMY, x, 0) for x in nodes]
    if not symmetric:
        dummy_edges += [(x, DUMMY, 0) for x in nodes]
        dummy_edges = reformulate_atsp_as_tsp(dummy_edges)

    tour = tsp(dummy_edges, precision=precision)
    dummy_index = tour.index(DUMMY)
    tour = tour[dummy_index:] + tour[:dummy_index]
    if symmetric:
        path = tour[1:]
    else:
        dummy_star_index = tour.index((DUMMY, '*'))
        assert dummy_star_index in (1, len(tour) - 1), tour
        if dummy_star_index == len(tour) - 1:  # need to flip
            tour = tour[1:] + tour[:1]
            tour = tour[::-1]
        path = tour[1:]
        path = [x for x in path if not isinstance(x, tuple)]

    return path


def tsp(edges, precision=0):
    c = Concorde(edges, precision=precision)
    return c.tour


def reformulate_atsp_as_tsp(edges):
    """
    To reformulate the ATSP as a TSP, for each city a dummy city (e.g, for New
    York, a dummy city New York* is added. Between each city and its
    corresponding dummy city a negative or very small distance with value cheap
    is used. This makes sure that each cities always occurs in the solution
    together with its dummy city. The original distances are used between the
    cities and the dummy cities, where each city is responsible for the distance
    going to the city and the dummy city is responsible for the distance coming
    from the city. The distances between all cities and the distances between
    all dummy cities are set to infeasible.
    """
    incident, nodes = node_to_edge(edges, directed=False)
    new_edges = []
    for a, b, w in edges:
        new_edges.append(((a, '*'), b, w))
    for n in nodes:
        new_edges.append((n, (n, '*'), NEG_INF))  # A negative weight
    return new_edges


def demo():
    from itertools import combinations
    from jcvi.utils.iter import pairwise
    from jcvi.graphics.base import plt, savefig

    POINTS = 100
    x = np.random.randn(POINTS)
    y = np.random.randn(POINTS)
    edges = []
    xy = zip(x, y)
    for ia, ib in combinations(range(POINTS), 2):
        ax, ay = xy[ia]
        bx, by = xy[ib]
        dist = ((ax - bx) ** 2 + (ay - by) ** 2) ** .5
        edges.append((ia, ib, dist))

    tour = hamiltonian(edges)
    plt.plot(x, y, "ro")
    for ia, ib in pairwise(tour):
        plt.plot((x[ia], x[ib]), (y[ia], y[ib]), "r-")

    savefig("demo.pdf")


if __name__ == '__main__':

    import doctest
    doctest.testmod()
    #demo()
