#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SSPACE scaffolding-related operations.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.fasta import gaps
from jcvi.formats.sizes import Sizes
from jcvi.formats.base import BaseFile, read_block
from jcvi.formats.agp import AGP
from jcvi.utils.iter import pairwise
from jcvi.algorithms.graph import BiGraph, BiEdge
from jcvi.apps.base import ActionDispatcher, debug
debug()


class EvidenceLine (object):

    def __init__(self, row, sizes):
        # f_tig3222|size7922|links348|gaps-109|merged16
        args = row.strip().split("|")
        nargs = len(args)

        tig = args[0]
        o, mtig = tig.split("_")
        tig = int(mtig.replace("tig", ""))
        assert o in ('f', 'r')
        self.o = ">" if o == 'f' else '<'

        name, size = sizes[tig]
        self.tig = name
        self.size = int(args[1].replace("size", ""))
        assert self.size == size, "{0} and {1} size mismatch".\
                format(mtig, name)

        if nargs > 2:
            self.links = int(args[2].replace("links", ""))
            self.gaps = int(args[3].replace("gaps", ""))
        if nargs > 4:
            self.merged = int(args[4].replace("merged", ""))


class EvidenceFile (BaseFile):

    def __init__(self, filename, fastafile):
        super(EvidenceFile, self).__init__(filename)
        sz = Sizes(fastafile)
        sizes = [None]  # tig-list starts at 1
        for name, size in sz.iter_sizes():
            sizes.append((name, size))
        self.sizes = sizes

    @property
    def graph(self):
        filename = self.filename
        sizes = self.sizes
        g = BiGraph()
        fp = open(filename)
        for header, lines in read_block(fp, ">"):
            lines = [EvidenceLine(x, sizes) for x in lines if x.strip()]

            for a, b in pairwise(lines):
                e = BiEdge(a.tig, b.tig, a.o, b.o, color=a.gaps)
                g.add_edge(e)

            if len(lines) == 1:  # Singleton scaffold
                a = lines[0]
                g.add_node(a.tig)

        return g


def main():

    actions = (
        ('anchor', 'anchor contigs to upgrade existing structure'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_target(p, name):
    before, before_tag = p.get_next(name, ">")
    if not before:  # Start of a scaffold
        return (None, ">")
    next, next_tag = p.get_next(name)
    if not next:  # End of a scaffold
        return (None, "<")
    # Internal to a scaffold
    return (next.v, "<")


def anchor(args):
    """
    %prog anchor evidencefile scaffolds.fasta contigs.fasta

    Use SSPACE evidencefile to scaffold contigs into existing scaffold
    structure, as in `scaffolds.fasta`. Contigs.fasta were used by SSPACE
    directly to scaffold.

    Rules:
    1. Only update existing structure by anchoring contigs in
    2. Promote singleton contigs only if they are >= 10Kb.
    """
    p = OptionParser(anchor.__doc__)
    p.add_option("--mingap", default=10, type="int",
                 help="Option -minGap used with gapSplit [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    evidencefile, scaffolds, contigs = args
    splitfasta, oagp, cagp = gaps([scaffolds, "--split"])

    agp = AGP(cagp)
    p = agp.graph

    ef = EvidenceFile(evidencefile, contigs)
    q = ef.graph

    logging.debug("Reference graph: {0}".format(p))
    logging.debug("Patch graph: {0}".format(q))

    newagp = list(agp)
    for a in agp:
        if a.is_gap:
            continue

        name = a.component_id
        target_name, tag = get_target(p, name)
        path = q.get_path(name, target_name, tag=tag)
        if tag == ">":
            path.reverse()

        if path and len(path) > 3:  # Heuristic, the patch must not be too long
            path = None

        print name, target_name, path


if __name__ == '__main__':
    main()
