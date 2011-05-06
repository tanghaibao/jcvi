#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
btab format, used by BLAST and NUCMER, found spec here:

<http://www.agcol.arizona.edu/~matthew/formats.html>
"""

import sys
from optparse import OptionParser

from jcvi.formats.blast import BlastLine
from jcvi.apps.base import ActionDispatcher, debug
debug()


class BtabLine (object):

    def __init__(self, row):
        args = row.strip().split("\t")
        self.query = args[0]
        #self.date = args[1]
        self.querylen = int(args[2])
        #self.algorithm = args[3]
        #self.subjectfilename = args[4]
        self.subject = args[5]
        self.qstart = int(args[6])
        self.qstop = int(args[7])
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        self.pctid = float(args[10])
        self.score = float(args[13])
        #self.description = args[14]
        #self.strand = "-" if args[17]=="Minus" else "Plus"
        self.subjectlen = int(args[18])
        self.evalue = float(args[19])

    @property
    def blastline(self):
        # some fields are not represented so ignore
        return "\t".join((self.query, self.subject,
                "%.2f" % self.pctid,
                "0", "0", "0",
                "%d" % self.qstart, "%d" % self.qstop,
                "%d" % self.sstart, "%d" % self.sstop,
                "%.1g" % self.evalue, "%.1f" % self.score))


def main():

    actions = (
        ('blast', 'convert back to BLAST -m8 format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def blast(args):
    """
    %prog blast btabfile

    convert back to BLAST -m8 format
    """
    p = OptionParser(blast.__doc__)

    opts, args = p.parse_args(args)

    if len(args) == 1:
        btabfile = args[0]
    else:
        sys.exit(p.print_help())

    fp = open(btabfile)
    for row in fp:
        b = BtabLine(row)
        print b.blastline


if __name__ == '__main__':
    main()
