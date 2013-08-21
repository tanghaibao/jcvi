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
        self.tig = tig
        assert o in ('f', 'r')
        self.o = ">" if o == 'f' else '<'

        name, size = sizes[tig]
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

        g = BiGraph()
        fp = open(filename)
        scf = {}  # Composition of SSPACE scfs
        for header, lines in read_block(fp, ">"):
            z = scf[header] = []
            for r in lines:
                if r.strip() == "":
                    continue
                r = EvidenceLine(r, sizes)
                z.append(r)

            for a, b in pairwise(z):
                e = BiEdge(a.tig, b.tig, a.o, b.o, color=a.gaps)
                g.add_edge(e)

        logging.debug(g)
        self.graph = g

    def get_next(self, tig, o="<"):
        return self.graph.nodes[tig].get_next(o)


def main():

    actions = (
        ('anchor', 'anchor contigs to upgrade existing structure'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def anchor(args):
    """
    %prog anchor evidencefile scaffolds.fasta contigs.fasta

    Use SSPACE evidencefile to scaffold contigs into existing scaffold
    structure, as in `scaffolds.fasta`. Contigs.fasta were used by SSPACE
    directly to scaffold.
    """
    p = OptionParser(anchor.__doc__)
    p.add_option("--mingap", default=10, type="int",
                 help="Option -minGap used with gapSplit [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    evidencefile, scaffolds, contigs = args
    splitfasta, oagp, cagp = gaps([scaffolds, "--split"])

    ef = EvidenceFile(evidencefile, contigs)


if __name__ == '__main__':
    main()
