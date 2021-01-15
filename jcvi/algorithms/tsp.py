#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
TSP solver using Concorde or OR-tools. This is much faster than the LP-formulation in
algorithms.lpsolve.tsp(). See also:
https://developers.google.com/optimization/routing/tsp
"""
import os.path as op
import os
import logging
import shutil
import numpy as np

from collections import defaultdict
from dataclasses import dataclass
from itertools import combinations
from more_itertools import pairwise

from jcvi.formats.base import FileShredder, must_open
from jcvi.apps.base import mkdir, which, sh


INF = 10000
NEG_INF = -INF
Work_dir = "tsp_work"


@dataclass
class TSPDataModel:
    edges: list  # List of tuple (source, target, weight)

    def distance_matrix(self, precision=0) -> tuple:
        """Compute the distance matrix

        Returns:
            np.array: Numpy square matrix with integer entries as distance
        """
        _, _, nodes = node_to_edge(self.edges, directed=False)
        nodes_indices = dict((n, i) for i, n in enumerate(nodes))
        nnodes = len(nodes)

        # TSPLIB requires explicit weights to be integral, and non-negative
        weights = [x[-1] for x in self.edges]
        max_x, min_x = max(weights), min(weights)
        inf = 2 * max(abs(max_x), abs(min_x))
        factor = 10 ** precision
        logging.debug(
            "TSP rescale: max_x=%d, min_x=%d, inf=%d, factor=%d",
            max_x,
            min_x,
            inf,
            factor,
        )

        D = np.ones((nnodes, nnodes), dtype=float) * inf
        for a, b, w in self.edges:
            ia, ib = nodes_indices[a], nodes_indices[b]
            D[ia, ib] = D[ib, ia] = w
        D = (D - min_x) * factor
        D = D.astype(int)
        return D, nodes

    def solve(self, time_limit=5, concorde=False, precision=0) -> list:
        """Solve the TSP instance.

        Args:
            time_limit (int, optional): Time limit to run. Default to 5 seconds.
            concorde (bool, optional): Shall we run concorde? Defaults to False.
            precision (int, optional): Float precision of distance. Defaults to 0.

        Returns:
            list: Ordered list of node indices to visit
        """
        if concorde:
            return Concorde(self, precision=precision).tour

        # Use OR-tools
        from ortools.constraint_solver import routing_enums_pb2
        from ortools.constraint_solver import pywrapcp

        D, nodes = self.distance_matrix(precision)
        nnodes = len(nodes)

        # Create the routing index manager
        manager = pywrapcp.RoutingIndexManager(nnodes, 1, 0)

        # Create routing model
        routing = pywrapcp.RoutingModel(manager)

        def distance_callback(from_index, to_index):
            """Returns the distance between the two nodes."""
            from_node = manager.IndexToNode(from_index)
            to_node = manager.IndexToNode(to_index)
            return D[from_node, to_node]

        transit_callback_index = routing.RegisterTransitCallback(distance_callback)

        # Define cost of each arc
        routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

        # Search strategy
        search_parameters = pywrapcp.DefaultRoutingSearchParameters()
        search_parameters.local_search_metaheuristic = (
            routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
        )
        search_parameters.time_limit.seconds = time_limit

        # Solve the problem
        solution = routing.SolveWithParameters(search_parameters)

        tour = []
        logging.info("Objective: %d", solution.ObjectiveValue())
        index = routing.Start(0)
        route_distance = 0
        while not routing.IsEnd(index):
            tour.append(manager.IndexToNode(index))
            previous_index = index
            index = solution.Value(routing.NextVar(index))
            route_distance = routing.GetArcCostForVehicle(previous_index, index, 0)
        logging.info("Route distance: %d", route_distance)

        return [nodes[x] for x in tour]


class Concorde(object):
    def __init__(
        self,
        data: TSPDataModel,
        work_dir=Work_dir,
        clean=True,
        verbose=False,
        precision=0,
        seed=666,
    ):
        """Run concorde on TSP instance

        Args:
            datamodel (TSPDataModel): TSP instance with edge weights
            work_dir ([type], optional): Path to the work dir. Defaults to Work_dir.
            clean (bool, optional): Clean up intermediate results. Defaults to True.
            verbose (bool, optional): Show verbose messages. Defaults to False.
            precision (int, optional): Float precision of distance. Defaults to 0.
            seed (int, optional): Random seed. Defaults to 666.
        """
        self.data = data
        self.work_dir = work_dir
        self.clean = clean
        self.verbose = verbose

        mkdir(work_dir)
        tspfile = op.join(work_dir, "data.tsp")
        self.print_to_tsplib(tspfile, precision=precision)
        _, outfile = self.run_concorde(tspfile, seed=seed)
        self.tour = self.parse_output(outfile)

        if clean:
            shutil.rmtree(work_dir)
            residual_output = ["data.sol", "data.res", "Odata.res"]
            FileShredder(residual_output, verbose=False)

    def print_to_tsplib(self, tspfile, precision=0):
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
        D, nodes = self.data.distance_matrix(precision)
        self.nodes = nodes
        self.nnodes = len(nodes)

        print("NAME: data", file=fw)
        print("TYPE: TSP", file=fw)
        print("DIMENSION: {}".format(self.nnodes), file=fw)
        print("EDGE_WEIGHT_TYPE: EXPLICIT", file=fw)
        print("EDGE_WEIGHT_FORMAT: FULL_MATRIX", file=fw)
        print("EDGE_WEIGHT_SECTION", file=fw)

        for row in D:  # Dump the full matrix
            print(" " + " ".join(str(x) for x in row), file=fw)

        print("EOF", file=fw)
        fw.close()
        logging.debug("Write TSP instance to `%s`", tspfile)

    def run_concorde(self, tspfile, seed=666):
        outfile = op.join(self.work_dir, "data.sol")
        if op.exists(outfile):
            os.remove(outfile)

        cc = "concorde"
        assert which(cc), (
            "You must install `concorde` on your PATH"
            + " [http://www.math.uwaterloo.ca/tsp/concorde.html]"
        )
        cmd = "{0} -s {1} -x -o {2} {3}".format(cc, seed, outfile, tspfile)

        outf = None if self.verbose else "/dev/null"
        retcode = sh(cmd, outfile=outf, errfile=outf)
        return retcode, outfile

    def parse_output(self, outfile):
        fp = open(outfile)
        dimension = int(next(fp).strip())  # header
        assert dimension == self.nnodes
        tour = []
        for row in fp:
            tour += [int(x) for x in row.split()]
        tour = [self.nodes[x] for x in tour]
        return tour


def node_to_edge(edges, directed=True):
    """
    From list of edges, record per node, incoming and outgoing edges
    """
    outgoing = defaultdict(set)
    incoming = defaultdict(set) if directed else outgoing
    nodes = set()
    for i, edge in enumerate(edges):
        (
            a,
            b,
        ) = edge[:2]
        outgoing[a].add(i)
        incoming[b].add(i)
        nodes.add(a)
        nodes.add(b)
    nodes = list(nodes)
    return outgoing, incoming, nodes


def populate_edge_weights(edges):
    # assume weight is 1 if not specified
    new_edges = []
    for e in edges:
        assert len(e) in (2, 3)
        if len(e) == 2:
            a, b = e
            w = 1
        else:
            a, b, w = e
        new_edges.append((a, b, w))
    return new_edges


def hamiltonian(edges, directed=False, time_limit=5, concorde=False, precision=0):
    """
    Calculates shortest path that traverses each node exactly once. Convert
    Hamiltonian path problem to TSP by adding one dummy point that has a distance
    of zero to all your other points. Solve the TSP and get rid of the dummy
    point - what remains is the Hamiltonian Path.

    >>> g = [(1,2), (2,3), (3,4), (4,2), (3,5)]
    >>> hamiltonian(g)
    [1, 2, 4, 3, 5]
    >>> hamiltonian([(1, 2), (2, 3)], directed=True)
    [1, 2, 3]
    """
    edges = populate_edge_weights(edges)
    _, _, nodes = node_to_edge(edges, directed=False)
    DUMMY = "DUMMY"
    dummy_edges = edges + [(DUMMY, x, 0) for x in nodes]
    if directed:
        dummy_edges += [(x, DUMMY, 0) for x in nodes]
        dummy_edges = reformulate_atsp_as_tsp(dummy_edges)

    tour = tsp(
        dummy_edges, time_limit=time_limit, concorde=concorde, precision=precision
    )

    dummy_index = tour.index(DUMMY)
    tour = tour[dummy_index:] + tour[:dummy_index]
    if directed:
        dummy_star_index = tour.index((DUMMY, "*"))
        assert dummy_star_index in (1, len(tour) - 1), tour
        if dummy_star_index == len(tour) - 1:  # need to flip
            tour = tour[1:] + tour[:1]
            tour = tour[::-1]
        path = tour[1:]
        path = [x for x in path if not isinstance(x, tuple)]
    else:
        path = tour[1:]

    return path


def tsp(edges, time_limit=5, concorde=False, precision=0) -> list:
    """Compute TSP solution

    Args:
        edges (list): List of tuple (source, target, weight)
        time_limit (int, optional): Time limit to run. Default to 5 seconds.
        concorde (bool, optional): Shall we run concorde? Defaults to False.
        precision (int, optional): Float precision of distance. Defaults to 0.

    Returns:
        list: List of nodes to visit
    """
    data = TSPDataModel(edges)
    return data.solve(time_limit=time_limit, concorde=concorde, precision=precision)


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
    _, _, nodes = node_to_edge(edges, directed=False)
    new_edges = []
    for a, b, w in edges:
        new_edges.append(((a, "*"), b, w))
    for n in nodes:
        new_edges.append((n, (n, "*"), NEG_INF))  # A negative weight
    return new_edges


def make_data(N, directed=False):
    x = np.random.randn(N)
    y = np.random.randn(N)
    xy = list(zip(x, y))
    M = np.zeros((N, N), dtype=float)
    for ia, ib in combinations(range(N), 2):
        ax, ay = xy[ia]
        bx, by = xy[ib]
        d = ((ax - bx) ** 2 + (ay - by) ** 2) ** 0.5
        M[ia, ib] = M[ib, ia] = d

    edges = []
    for ia, ib in combinations(range(N), 2):
        edges.append((ia, ib, M[ia, ib]))
        if directed:
            edges.append((ib, ia, M[ib, ia]))

    return x, y, M, edges


def evaluate(tour, M):
    score = 0
    for ia, ib in pairwise(tour):
        score += M[ia, ib]
    return score


def plot_data(x, y, tour, M):
    from jcvi.graphics.base import plt, savefig

    plt.plot(x, y, "ro")
    for ia, ib in pairwise(tour):
        plt.plot((x[ia], x[ib]), (y[ia], y[ib]), "r-")

    score = evaluate(tour, M)
    plt.title("Score={0:.2f}".format(score))

    savefig("demo.pdf")


def concorde_demo(POINTS=100):
    x, y, M, edges = make_data(POINTS)
    ctour = hamiltonian(edges, precision=3)
    plot_data(x, y, ctour, M)


def compare_lpsolve_to_concorde(POINTS=80, directed=False):
    from jcvi.algorithms.lpsolve import hamiltonian as lhamiltonian

    _, _, M, edges = make_data(POINTS, directed=directed)
    ltour = lhamiltonian(edges, directed=directed)
    print(ltour, evaluate(ltour, M))

    ctour = hamiltonian(edges, directed=directed, precision=3)
    print(ctour, evaluate(ctour, M))


if __name__ == "__main__":
    import doctest

    doctest.testmod()
