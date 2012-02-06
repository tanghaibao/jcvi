#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert common output files from gene prediction softwares into gff3 format.

Similar to the utilities in DAWGPAWS.
<http://dawgpaws.sourceforge.net/man.html>
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.gff import GffLine
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('augustus', 'convert augustus output into gff3'),
        ('rename', 'change the IDs within the gff3'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rename(args):
    """
    %prog reindex in.gff3 > reindexed.gff3

    Change the IDs within the gff3.
    """
    p = OptionParser(rename.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    ingff3, = args
    fp = open(ingff3)
    for row in fp:
        if row[0] == '#':
            continue

        g = GffLine(row)
        id, = g.attributes["ID"]
        prefix = "{0}.{1}.".format(g.source, g.seqid)
        g.attributes["ID"] = [prefix + id]

        if "Parent" in g.attributes:
            parents = g.attributes["Parent"]
            g.attributes["Parent"] = [prefix + x for x in parents]

        g.attributes_text = g._attributes_text()
        print g


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
    fp = open(ingff3)
    for row in fp:
        if row[0] == '#':
            continue

        g = GffLine(row)
        if g.type not in ("gene", "transcript", "CDS"):
            continue

        if g.type == "transcript":
            g.type = "mRNA"

        print g


if __name__ == '__main__':
    main()
