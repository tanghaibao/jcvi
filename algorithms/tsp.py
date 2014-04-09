#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
TSP solver using Concorde. This is much faster than the LP-formulation in
algorithms.lpsolve.tsp().
"""

import os.path as op
import os
import logging
import math

import numpy as np

from jcvi.formats.base import must_open
from jcvi.utils.iter import pairwise
from jcvi.algorithms.lpsolve import populate_edge_weights, node_to_edge
from jcvi.apps.base import mkdir, debug, which, sh
debug()

Work_dir = "tsp_work"


class Concorde (object):

    def __init__(self, edges, work_dir=Work_dir, clean=False, verbose=False):

        self.work_dir = work_dir
        self.clean = clean
        self.verbose = verbose

        mkdir(work_dir)
        tspfile = op.join(work_dir, "data.tsp")
        self.print_to_tsplib(edges, tspfile)
        retcode, outfile = self.run_concorde(tspfile)
        self.tour = self.parse_output(outfile)

    def print_to_tsplib(self, edges, tspfile):
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

        # TSPLIB requires explicit weights to be integral
        weights = [x[-1] for x in edges if x[-1]]
        min_weight = min(weights) if weights else 0
        factor = 1. / min_weight if min_weight < 1 else 1
        factor = min(factor, 1000)
        edges = [(a, b, int(w * factor)) for a, b, w in edges]
        Maxdist = 10 ** len(str(sum(x[-1] for x in edges))) - 1

        print >> fw, "NAME: data"
        print >> fw, "TYPE: TSP"
        print >> fw, "DIMENSION: {0}".format(nnodes)

        D = np.ones((nnodes, nnodes), dtype=int) * Maxdist
        for a, b, w in edges:
            ia, ib = nodes_indices[a], nodes_indices[b]
            D[ia, ib] = D[ib, ia] = w

        print >> fw, "EDGE_WEIGHT_TYPE: EXPLICIT"
        print >> fw, "EDGE_WEIGHT_FORMAT: FULL_MATRIX"
        print >> fw, "EDGE_WEIGHT_SECTION"
        for row in D:  # Dump the full matrix
            print >> fw, " " + " ".join(str(int(x)) for x in row)

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
        retcode = sh(cmd, outfile=outf)
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


def hamiltonian(edges):
    """
    Calculates shortest path that traverses each node exactly once. Convert
    Hamiltonian path problem to TSP by adding one dummy point that has a distance
    of zero to all your other points. Solve the TSP and get rid of the dummy
    point - what remains is the Hamiltonian Path.

    >>> g = [(1,2), (2,3), (3,4), (4,2), (3,5)]
    >>> hamiltonian(g)
    [1, 2, 4, 3, 5]
    >>> g = [(1,2), (2,3), (1,4), (2,5), (3,6), (5,6)]
    >>> hamiltonian(g)
    [4, 1, 2, 3, 6, 5]
    """
    edges = populate_edge_weights(edges)
    incident, nodes = node_to_edge(edges, directed=False)
    DUMMY = "DUMMY"
    dummy_edges = edges + [(DUMMY, x, 0) for x in nodes]

    tour = tsp(dummy_edges)
    dummy_index = tour.index(DUMMY)
    path = tour[dummy_index + 1:] + tour[:dummy_index]
    path = min(path, path[::-1])

    return path


def tsp(edges):
    c = Concorde(edges)
    return c.tour


if __name__ == '__main__':

    import doctest
    doctest.testmod()
