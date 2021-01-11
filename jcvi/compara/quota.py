#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Quota synteny alignment (QUOTA-ALIGN)

%prog [options] anchorsfile --qbed=qbedfile --sbed=sbedfile

This python program does the following:
1. merge 2D-overlapping blocks (now skipped, but existed in original version)
2. build constraints that represent 1D-overlap among blocks
3. feed the data into the linear programming solver

The algorithm is described in Tang et al. BMC Bioinformatics 2011.
"Screening synteny blocks in pairwise genome comparisons through integer
programming."
"""

import logging
import os.path as op
import sys

from dataclasses import dataclass
from io import StringIO

from jcvi.algorithms.lpsolve import GLPKSolver, SCIPSolver
from jcvi.compara.synteny import AnchorFile, _score, check_beds
from jcvi.formats.base import must_open
from jcvi.utils.console import printf
from jcvi.apps.base import OptionParser


@dataclass
class MIPDataModel:
    """Data model for use with OR-tools. Modeled after the tutorial."""

    constraint_coeffs: list
    bounds: list
    obj_coeffs: list
    num_vars: int
    num_constraints: int


def get_1D_overlap(eclusters, depth=1):
    """
    Find blocks that are 1D overlapping,
    returns cliques of block ids that are in conflict
    """
    overlap_set = set()
    active = set()

    ends = []
    for i, (chr, left, right) in enumerate(eclusters):
        ends.append((chr, left, 0, i))  # 0/1 for left/right-ness
        ends.append((chr, right, 1, i))
    ends.sort()

    chr_last = ""
    for chr, _, left_right, i in ends:
        if chr != chr_last:
            active.clear()
        if left_right == 0:
            active.add(i)
        else:
            active.remove(i)

        if len(active) > depth:
            overlap_set.add(tuple(sorted(active)))

        chr_last = chr

    return overlap_set


def make_range(clusters, extend=0):
    """
    Convert to interval ends from a list of anchors
    extend modifies the xmax, ymax boundary of the box,
    which can be positive or negative
    very useful when we want to make the range as fuzzy as we specify
    """
    eclusters = []
    for cluster in clusters:
        xlist, ylist, _ = zip(*cluster)
        score = _score(cluster)

        xchr, xmin = min(xlist)
        xchr, xmax = max(xlist)
        ychr, ymin = min(ylist)
        ychr, ymax = max(ylist)

        # allow fuzziness to the boundary
        xmax += extend
        ymax += extend
        # because extend can be negative values, we don't want it to be less than min
        if xmax < xmin:
            xmin, xmax = xmax, xmin
        if ymax < ymin:
            ymin, ymax = ymax, ymin

        eclusters.append(((xchr, xmin, xmax), (ychr, ymin, ymax), score))

    return eclusters


def get_constraints(clusters, quota=(1, 1), Nmax=0):
    """
    Check pairwise cluster comparison, if they overlap then mark edge as conflict
    """
    qa, qb = quota
    eclusters = make_range(clusters, extend=-Nmax)
    # (1-based index, cluster score)
    nodes = [(i + 1, c[-1]) for i, c in enumerate(eclusters)]

    eclusters_x, eclusters_y, _ = zip(*eclusters)

    # represents the contraints over x-axis and y-axis
    constraints_x = get_1D_overlap(eclusters_x, qa)
    constraints_y = get_1D_overlap(eclusters_y, qb)

    return nodes, constraints_x, constraints_y


def create_data_model(nodes, constraints_x, qa, constraints_y, qb):
    """
    Maximize
     4 x1 + 2 x2 + 3 x3 + x4
    Subject To
     x1 + x2 <= 1
    End
    """
    num_vars = len(nodes)
    obj_coeffs = [0] * num_vars
    for i, score in nodes:
        obj_coeffs[i - 1] = score

    num_constraints = 0
    constraint_coeffs = []
    bounds = []
    for c in constraints_x:
        constraint = [0] * num_vars
        for x in c:
            constraint[x] = 1
        constraint_coeffs.append(constraint)
        bounds.append(qa)
    num_constraints = len(constraints_x)

    # non-self
    if not (constraints_x is constraints_y):
        for c in constraints_y:
            constraint = [0] * num_vars
            for x in c:
                constraint[x] = 1
            constraint_coeffs.append(constraint)
            bounds.append(qb)
        num_constraints += len(constraints_y)

    return MIPDataModel(
        constraint_coeffs, bounds, obj_coeffs, num_vars, num_constraints
    )


def format_lp(data: MIPDataModel):
    """Format data dictionary into MIP formatted string.

    Args:
        data (DataModel): Data model that contains vars and constraints

    Returns:
        str: MIP formatted string
    """
    lp_handle = StringIO()

    lp_handle.write("Maximize\n ")
    records = 0
    for i, score in enumerate(data.obj_coeffs):
        lp_handle.write("+ %d x%d " % (score, i))
        # SCIP does not like really long string per row
        records += 1
        if records % 10 == 0:
            lp_handle.write("\n")
    lp_handle.write("\n")

    lp_handle.write("Subject To\n")
    for constraint, bound in zip(data.constraint_coeffs, data.bounds):
        additions = " + ".join("x{}".format(i) for (i, x) in enumerate(constraint) if x)
        lp_handle.write(" %s <= %d\n" % (additions, bound))

    printf(
        "number of variables ({}), number of constraints ({})".format(
            data.num_vars, data.num_constraints
        ),
        file=sys.stderr,
    )

    lp_handle.write("Binary\n")
    for i in range(data.num_vars):
        lp_handle.write(" x{}\n".format(i))

    lp_handle.write("End\n")

    lp_data = lp_handle.getvalue()
    lp_handle.close()

    return lp_data


def create_solver(data: MIPDataModel, backend: str = "SCIP"):
    """
    Create OR-tools solver instance. See also:
    https://developers.google.com/optimization/mip/mip_var_array

    Args:
        data (MIPDataModel): Data model that contains vars and constraints
        backend (str, optional): Backend for the MIP solver. Defaults to "SCIP".

    Returns:
        OR-tools solver instance
    """
    from ortools.linear_solver import pywraplp

    solver = pywraplp.Solver.CreateSolver(backend)
    x = {}
    for j in range(data.num_vars):
        x[j] = solver.IntVar(0, 1, "x[%i]" % j)
    printf("Number of variables =", solver.NumVariables())

    for bound, constraint_coeff in zip(data.bounds, data.constraint_coeffs):
        constraint = solver.RowConstraint(0, bound, "")
        for j, coeff in enumerate(constraint_coeff):
            if coeff:
                constraint.SetCoefficient(x[j], coeff)
    printf("Number of constraints =", solver.NumConstraints())

    objective = solver.Objective()
    for j, score in enumerate(data.obj_coeffs):
        objective.SetCoefficient(x[j], score)
    objective.SetMaximization()

    return solver, x


def has_ortools() -> bool:
    """Do we have an installation of OR-tools?

    Returns:
        bool: True if installed
    """
    try:
        from ortools.linear_solver import pywraplp

        return True
    except ImportError:
        return False


def solve_lp(
    clusters,
    quota,
    work_dir="work",
    Nmax=0,
    self_match=False,
    verbose=False,
):
    """
    Solve the formatted LP instance
    """
    qb, qa = quota  # flip it
    nodes, constraints_x, constraints_y = get_constraints(clusters, (qa, qb), Nmax=Nmax)

    if self_match:
        constraints_x = constraints_y = constraints_x | constraints_y

    data = create_data_model(nodes, constraints_x, qa, constraints_y, qb)
    filtered_list = []

    if has_ortools():
        # Use OR-tools
        from ortools.linear_solver import pywraplp

        solver, x = create_solver(data)
        status = solver.Solve()
        if status == pywraplp.Solver.OPTIMAL:
            printf("Objective value =", solver.Objective().Value())
            filtered_list = [
                j for j in range(data.num_vars) if x[j].solution_value() == 1
            ]
            printf("Problem solved in {} milliseconds".format(solver.wall_time()))
            printf("Problem solved in {} iterations".format(solver.iterations()))
            printf("Problem solved in {} branch-and-bound nodes".format(solver.nodes()))
    else:
        # Use custom formatter
        lp_data = format_lp(data)
        filtered_list = SCIPSolver(lp_data, work_dir, verbose=verbose).results
        if not filtered_list:
            printf("SCIP fails... trying GLPK", file=sys.stderr)
            filtered_list = GLPKSolver(lp_data, work_dir, verbose=verbose).results

    return filtered_list


def read_clusters(qa_file, qorder, sorder):
    """Read in the clusters from anchors file

    Args:
        qa_file (str): Path to input file
        qorder (dict): Dictionary to find position of feature in query
        sorder (dict): Dictionary to find position of feature in subject

    Returns:
        List: List of matches and scores
    """
    af = AnchorFile(qa_file)
    blocks = af.blocks
    clusters = []
    for block in blocks:
        cluster = []
        for a, b, score in block:
            ia, oa = qorder[a]
            ib, ob = sorder[b]
            ca, cb = oa.seqid, ob.seqid
            cluster.append(((ca, ia), (cb, ib), score))
        clusters.append(cluster)

    return clusters


def main(args):
    p = OptionParser(__doc__)

    p.set_beds()
    p.add_option(
        "--quota",
        default="1:1",
        help="`quota mapping` procedure -- screen blocks to constrain mapping"
        " (useful for orthology), "
        "put in the format like (#subgenomes expected for genome X):"
        "(#subgenomes expected for genome Y)",
    )
    p.add_option(
        "--Nm",
        dest="Nmax",
        type="int",
        default=10,
        help="distance cutoff to tolerate two blocks that are "
        "slightly overlapping (cutoff for `quota mapping`) "
        "[default: %default units (gene or bp dist)]",
    )

    p.add_option(
        "--self",
        dest="self_match",
        action="store_true",
        default=False,
        help="you might turn this on when screening paralogous blocks, "
        "esp. if you have reduced mirrored blocks into non-redundant set",
    )
    p.set_verbose(help="Show verbose solver output")

    p.add_option(
        "--screen",
        default=False,
        action="store_true",
        help="generate new anchors file",
    )

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (qa_file,) = args
    _, _, qorder, sorder, _ = check_beds(qa_file, p, opts)

    # sanity check for the quota
    if opts.quota:
        try:
            qa, qb = opts.quota.split(":")
            qa, qb = int(qa), int(qb)
        except ValueError:
            printf(
                "quota string should be the form x:x (2:4, 1:3, etc.)", file=sys.stderr
            )
            sys.exit(1)

        if opts.self_match and qa != qb:
            raise Exception(
                "when comparing genome to itself, "
                "quota must be the same number "
                "(like 1:1, 2:2) you have %s" % opts.quota
            )
        quota = (qa, qb)

    self_match = opts.self_match

    clusters = read_clusters(qa_file, qorder, sorder)
    for cluster in clusters:
        assert len(cluster) > 0

    # below runs `quota mapping`
    work_dir = op.join(op.dirname(op.abspath(qa_file)), "work")

    selected_ids = solve_lp(
        clusters,
        quota,
        work_dir=work_dir,
        Nmax=opts.Nmax,
        self_match=self_match,
        verbose=opts.verbose,
    )

    logging.debug("Selected %d blocks", len(selected_ids))
    prefix = qa_file.rsplit(".", 1)[0]
    suffix = "{}x{}".format(qa, qb)
    outfile = ".".join((prefix, suffix))
    fw = must_open(outfile, "w")
    print(",".join(str(x) for x in selected_ids), file=fw)
    fw.close()
    logging.debug("Screened blocks ids written to `%s`", outfile)

    if opts.screen:
        from jcvi.compara.synteny import screen

        new_qa_file = ".".join((prefix, suffix, "anchors"))
        largs = [qa_file, new_qa_file, "--ids", outfile]
        if opts.qbed and opts.sbed:
            largs += ["--qbed={0}".format(opts.qbed)]
            largs += ["--sbed={0}".format(opts.sbed)]
        screen(largs)


if __name__ == "__main__":
    main(sys.argv[1:])
