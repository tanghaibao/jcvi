#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Syntenic path assembly.
"""

import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.agp import order_to_agp
from jcvi.formats.blast import BlastSlow
from jcvi.formats.sizes import Sizes
from jcvi.utils.iter import pairwise
from jcvi.algorithms.graph import BiGraph, BiEdge
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('fromblast', 'Generate path from BLAST file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fromblast(args):
    """
    %prog fromblast blastfile subject.fasta

    Generate path from BLAST file. If multiple subjects map to the same query,
    an edge is constructed between them (with the link provided by the query).

    The BLAST file MUST be filtered, chained, supermapped, and sorted by --query.
    """
    p = OptionParser(fromblast.__doc__)
    p.add_option("--verbose", default=False, action="store_true",
                 help="Print verbose reports to stdout [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blastfile, subjectfasta = args
    blast = BlastSlow(blastfile, sorted=True)
    g = BiGraph()
    for query, blines in groupby(blast, key=lambda x: x.query):
        blines = list(blines)
        for a, b in pairwise(blines):
            asub, bsub = a.subject, b.subject
            if asub == bsub:
                continue

            atag = ">" if a.orientation == "+" else "<"
            btag = ">" if b.orientation == "+" else "<"
            g.add_edge(BiEdge(asub, bsub, atag, btag))

    logging.debug(str(g))
    paths = []
    for path in g.iter_paths():
        m, oo = g.path(path)
        if len(oo) == 1:  # Singleton path
            continue
        paths.append(oo)
        if opts.verbose:
            print m
            print oo

    npaths = len(paths)
    ntigs = sum(len(x) for x in paths)
    logging.debug("Graph decomposed to {0} paths with {1} components.".\
                  format(npaths, ntigs))

    agpfile = blastfile + ".agp"
    sizes = Sizes(subjectfasta)
    fwagp = open(agpfile, "w")
    scaffolded = set()
    for i, oo in enumerate(paths):
        ctgorder = [(str(ctg), ("+" if strand else "-")) \
                     for ctg, strand in oo]
        scaffolded |= set(ctg for ctg, strand in ctgorder)
        object = "pmol_{0:04d}".format(i)
        order_to_agp(object, ctgorder, sizes.mapping, fwagp)

    # Get the singletons as well
    nsingletons = 0
    for ctg, size in sizes.iter_sizes():
        if ctg in scaffolded:
            continue

        ctgorder = [(ctg, "+")]
        object = ctg
        order_to_agp(object, ctgorder, sizes.mapping, fwagp)
        nsingletons += 1
    logging.debug("Written {0} unscaffolded singletons.".format(nsingletons))

    fwagp.close()
    logging.debug("AGP file written to `{0}`.".format(agpfile))


if __name__ == '__main__':
    main()
