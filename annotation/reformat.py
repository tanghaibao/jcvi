#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert common output files from gene prediction softwares into gff3 format.

Similar to the utilities in DAWGPAWS.
<http://dawgpaws.sourceforge.net/man.html>
"""

import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.gff import GffLine, Gff
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('rename', 'rename genes for annotation release'),
        ('augustus', 'convert augustus output into gff3'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rename(args):
    """
    %prog rename genes.bed gaps.bed

    Rename genes for annotation release.

    For genes on chromosomes (e.g. the 12th gene on C1):
    Bo1g00120

    For genes on scaffolds (e.g. the 12th gene on unplaced Scaffold00285):
    Bo00285s120

    The genes identifiers will increment by 10. So assuming no gap, these are
    the consecutive genes:
    Bo1g00120, Bo1g00130, Bo1g00140...
    Bo00285s120, Bo00285s130, Bo00285s140...

    When we encounter gaps, we would like the increment to be larger. For example,
    Bo1g00120, <gap>, Bo1g01120...
    """
    import string
    from jcvi.formats.bed import Bed, BedLine

    p = OptionParser(rename.__doc__)
    p.add_option("-a", dest="gene_increment", default=10, type="int",
                 help="Increment for continuous genes [default: %default]")
    p.add_option("-b", dest="gap_increment", default=1000, type="int",
                 help="Increment for gaps [default: %default]")
    p.add_option("--pad0", default=6, type="int",
                 help="Pad gene identifiers with 0 [default: %default]")
    p.add_option("--prefix", default="Bo",
                 help="Genome prefix [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    genebed, gapbed = args
    prefix = opts.prefix
    gene_increment = opts.gene_increment
    gap_increment = opts.gap_increment

    genes = Bed(genebed)
    fp = open(gapbed)
    for row in fp:
        genes.append(BedLine(row))

    genes.sort(key=genes.key)
    idsfile = prefix + ".ids"
    newbedfile = prefix + ".bed"
    fw = open(idsfile, "w")

    for chr, lines in groupby(genes, key=lambda x: x.seqid):
        lines = list(lines)
        pad0 = opts.pad0 if len(lines) > 1000 else 3
        isChr = chr[0].upper() == 'C'
        digits = "".join(x for x in chr if x in string.digits)
        gs = "g" if isChr else "s"
        pp = prefix + digits + gs
        idx = 0
        for r in lines:
            isGap = r.strand not in ("+", "-")
            if isGap:
                idx += gap_increment - gene_increment
                continue
            else:
                idx += gene_increment
            accn = pp + "{0:0{1}d}".format(idx, pad0)
            oldaccn = r.accn
            print >> fw, "\t".join((oldaccn, accn))
            r.accn = accn

    genes.print_to_file(newbedfile)
    logging.debug("Converted IDs written to `{0}`.".format(idsfile))
    logging.debug("Converted bed written to `{0}`.".format(newbedfile))


def augustus(args):
    """
    %prog augustus augustus.gff3 > reformatted.gff3

    AUGUSTUS does generate a gff3 (--gff3=on) but need some refinement.
    """
    p = OptionParser(augustus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    ingff3, = args
    gff = Gff(ingff3)
    for g in gff:
        if g.type not in ("gene", "transcript", "CDS"):
            continue

        if g.type == "transcript":
            g.type = "mRNA"

        print g


if __name__ == '__main__':
    main()
