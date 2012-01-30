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
        ('happy', 'Make graph from happy mapping data'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def happy(args):
    """
    %prog happy happy.txt

    Make bi-directed graph from HAPPY mapping data. JCVI encodes uncertainties
    in the order of the contigs / scaffolds.

    : separates scaffolds
    + means telomere (though the telomere repeats may not show because the
    telomere-adjacent sequence is missing)
    - means that the scaffold is in reverse orientation to that shown in the 2003
    TIGR scaffolds.

    Ambiguities are represented as follows, using Paul Dear.s description:
    [ ] means undetermined orientation. error quite possible (70% confidence?)
    ( ) means uncertain orientation. small chance of error (90% confidence?)
    { } means uncertain order.

    Example:
    +-8254707:8254647:-8254690:{[8254694]:[8254713]:[8254531]:[8254797]}:8254802:8254788+
    """
    from string import maketrans
    from jcvi.utils.iter import pairwise

    p = OptionParser(happy.__doc__)
    p.add_option("--prefix", help="Add prefix to the name [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    happyfile, = args
    prefix = opts.prefix

    certain = "certain.graph"
    uncertain = "uncertain.graph"
    fw1 = open(certain, "w")
    fw2 = open(uncertain, "w")
    n_certain = n_uncertain = 0

    fp = open(happyfile)
    trans = maketrans("[](){}", "      ")
    for row in fp:
        row = row.strip().strip("+")
        row = row.translate(trans)
        scfs = [x.strip("+") for x in row.split(":")]
        for a, b in pairwise(scfs):
            oa = '<' if a.strip()[0] == '-' else '>'
            ob = '<' if b.strip()[0] == '-' else '>'

            is_uncertain = a[-1] == ' ' or b[0] == ' '
            if is_uncertain:
                n_uncertain += 1
            else:
                n_certain += 1

            a = a.strip().strip('-')
            b = b.strip().strip('-')

            if prefix:
                a = prefix + a
                b = prefix + b

            e = BiEdge(a, b, oa, ob)
            fw = fw2 if is_uncertain else fw1
            print >> fw, e

    logging.debug("Certain edges: {0}, Uncertain edges: {1} written to `{2}`".\
                  format(n_certain, n_uncertain, ",".join((certain, uncertain))))


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

    g.write("graph.txt")
    g.draw("graph.pdf")

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
