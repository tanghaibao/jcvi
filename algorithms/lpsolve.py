#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Implement a few MIP solvers, based on benchmark found on <http://scip.zib.de/>
SCIP solver is ~16x faster than GLPK solver.  However, I found in rare cases
it will segfault. Therefore the default is SCIP, the program will switch to
GLPK solver for crashed cases.

The input lp_data is assumed in .lp format, see below

>>> lp_data = '''
... Maximize
...  5 x1 + 3 x2 + 2 x3
... Subject to
...  x2 + x3 <= 1
... Binary
...  x1
...  x2
...  x3
... End'''
>>> print SCIPSolver(lp_data).results
[0, 1]
>>> print GLPKSolver(lp_data).results
[0, 1]
"""

import os
import os.path as op
import sys
import shutil
import logging
import cStringIO

from collections import defaultdict

from jcvi.utils.cbook import fill
from jcvi.apps.base import sh, debug
debug()


Work_dir = "lpsolve_work"

# CPLEX LP format
# <http://lpsolve.sourceforge.net/5.0/CPLEX-format.htm>
MAXIMIZE = "Maximize"
MINIMIZE = "Minimize"
SUBJECTTO = "Subject To"
BOUNDS = "Bounds"
BINARY = "Binary"
END = "End"


class AbstractMIPSolver(object):
    """
    Base class for LP solvers
    """
    def __init__(self, lp_data, work_dir=Work_dir, clean=True, verbose=False):

        self.work_dir = work_dir
        self.clean = clean
        self.verbose = verbose

        if not os.path.isdir(work_dir):
            os.mkdir(work_dir)

        lpfile = op.join(work_dir, "data.lp")  # problem instance
        logging.debug("write MIP instance to `{0}`".format(lpfile))

        fw = open(lpfile, "w")
        fw.write(lp_data)
        fw.close()

        retcode, outfile = self.run(lpfile)
        if retcode < 0:
            self.results = []
        else:
            self.results = self.parse_output(outfile)

        if self.results:
            logging.debug("optimized objective value ({0})".\
                    format(self.obj_val))

    def run(self, lp_data, work_dir):
        raise NotImplementedError

    def parse_output(self):
        raise NotImplementedError

    def cleanup(self):
        shutil.rmtree(self.work_dir)


class GLPKSolver(AbstractMIPSolver):
    """
    GNU Linear Programming Kit (GLPK) solver, wrapper for calling GLPSOL
    """
    def run(self, lpfile):

        outfile = op.join(self.work_dir, "data.lp.out")  # verbose output
        listfile = op.join(self.work_dir, "data.lp.list")  # simple output
        # cleanup in case something wrong happens
        for f in (outfile, listfile):
            if op.exists(f):
                os.remove(f)

        cmd = "glpsol --cuts --fpump --lp {0} -o {1} -w {2}".format(lpfile,
                outfile, listfile)

        outf = None if self.verbose else "/dev/null"
        retcode = sh(cmd, outfile=outf)

        if retcode == 127:
            logging.error("You need to install program `glpsol` " + \
                          "[http://www.gnu.org/software/glpk/]")
            return -1, None

        return retcode, listfile

    def parse_output(self, listfile, clean=False):

        filtered_list = []

        fp = open(listfile)
        header = fp.readline()
        columns, rows = header.split()
        rows = int(rows)
        data = fp.readlines()
        self.obj_val = int(data[0].split()[-1])
        # the info are contained in the last several lines
        results = [int(x) for x in data[-rows:]]
        results = [i for i, x in enumerate(results) if x == 1]

        fp.close()

        if self.clean:
            self.cleanup()

        return results


class SCIPSolver(AbstractMIPSolver):
    """
    SCIP solver, wrapper for calling SCIP executable
    """
    def run(self, lpfile):

        outfile = self.work_dir + "/data.lp.out"  # verbose output
        if op.exists(outfile):
            os.remove(outfile)

        cmd = "scip -f {0} -l {1}".format(lpfile, outfile)

        outf = None if self.verbose else "/dev/null"
        retcode = sh(cmd, outfile=outf)

        if retcode == 127:
            logging.error("You need to install program `scip` " +\
                          "[http://scip.zib.de/]")
            return -1, None

        return retcode, outfile

    def parse_output(self, outfile):

        fp = open(outfile)
        for row in fp:
            if row.startswith("objective value"):
                obj_row = row
                break

        results = []
        for row in fp:
            """
            objective value:               8
            x1                             1   (obj:5)
            x2                             1   (obj:3)
            """
            if row.strip() == "":  # blank line ends the section
                break
            x = row.split()[0]
            results.append(int(x[1:]) - 1)  # 0-based indexing

        if results:
            self.obj_val = int(obj_row.split(":")[1])

        fp.close()

        if self.clean:
            self.cleanup()

        return results


def node_to_edge(edges):
    """
    From list of edges, record per node, incoming and outgoing edges
    """
    outgoing = defaultdict(set)
    incoming = defaultdict(set)
    nodes = set()
    for i, edge in enumerate(edges):
        a, b, w = edge
        outgoing[a].add(i)
        incoming[b].add(i)
        nodes.add(a)
        nodes.add(b)
    return outgoing, incoming, nodes


def print_objective(lp_handle, edges, objective=MAXIMIZE):
    """
    CPLEX LP format commonly contains three blocks:
    objective, constraints, vars
    spec <http://lpsolve.sourceforge.net/5.0/CPLEX-format.htm>
    """
    print >> lp_handle, objective
    items = [" + {0}x{1}".format(w, i + 1) \
            for i, (a, b, w) in enumerate(edges)]
    sums = fill(items, width=10)
    print >> lp_handle, sums


def print_constraints(lp_handle, constraints):
    print >> lp_handle, SUBJECTTO
    for cc in constraints:
        print >> lp_handle, cc


def print_vars(lp_handle, nedges, vars=BINARY):
    print >> lp_handle, vars
    for i in xrange(nedges):
        print >> lp_handle, " x{0}".format(i + 1)

    print >> lp_handle, END


def lpsolve(lp_handle, solver="scip", clean=True):

    solver = SCIPSolver if solver == "scip" else GLPKSolver
    lp_data = lp_handle.getvalue()
    lp_handle.close()

    g = solver(lp_data, clean=clean)
    selected = set(g.results)
    return selected, g.obj_val


def path(edges, source, sink, flavor="longest"):
    """
    Calculates shortest/longest path from list of edges in a graph

    >>> g = [(1,2,1),(2,3,9),(2,4,3),(2,5,2),(3,6,8),(4,6,10),(4,7,4)]
    >>> g += [(6,8,7),(7,9,5),(8,9,6),(9,10,11)]
    >>> print path(g, 1, 8, flavor="shortest")
    ([(1, 2, 1), (2, 4, 3), (4, 6, 10), (6, 8, 7)], 21)
    >>> print path(g, 1, 8, flavor="longest")
    ([(1, 2, 1), (2, 3, 9), (3, 6, 8), (6, 8, 7)], 25)
    """
    outgoing, incoming, nodes = node_to_edge(edges)

    nedges = len(edges)
    lp_handle = cStringIO.StringIO()

    assert flavor in ("longest", "shortest")

    objective = MAXIMIZE if flavor == "longest" else MINIMIZE
    print_objective(lp_handle, edges, objective=objective)

    # Balancing constraint, incoming edges equal to outgoing edges except
    # source and sink

    constraints = []
    for v in nodes:
        incoming_edges = incoming[v]
        outgoing_edges = outgoing[v]
        icc = "".join(" + x{0}".format(i + 1) for i in incoming_edges)
        occ = "".join(" + x{0}".format(i + 1) for i in outgoing_edges)

        if v == source:
            if not outgoing_edges:
                return None
            constraints.append("{0} = 1".format(occ))
        elif v == sink:
            if not incoming_edges:
                return None
            constraints.append("{0} = 1".format(icc))
        else:
            # Balancing
            constraints.append("{0}{1} = 0".format(icc, occ.replace('+', '-')))
            # Simple path
            if incoming_edges:
                constraints.append("{0} <= 1".format(icc))
            if outgoing_edges:
                constraints.append("{0} <= 1".format(occ))

    print_constraints(lp_handle, constraints)
    print_vars(lp_handle, nedges, vars=BINARY)

    selected, obj_val = lpsolve(lp_handle)
    results = sorted(x for i, x in enumerate(edges) if i in selected)
    if not results:
        results = None

    return results, obj_val


if __name__ == '__main__':

    import doctest
    doctest.testmod()
