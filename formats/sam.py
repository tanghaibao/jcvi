"""
SAM alignment format. There are other tools that handles better SAM and BAM.
This script simply parses the lines in SAM into human readable fields.

http://samtools.sourceforge.net/SAM1.pdf
"""

import sys

from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.apps.base import ActionDispatcher


class SamLine (object):

    def __init__(self, row):

        args = row.strip().split("\t")
        self.qname = args[0]
        self.flag = args[1]
        self.rname = args[2]
        self.pos = args[3]
        self.mapq = args[4]
        self.cigar = args[5]
        self.mrnm = args[6]
        self.mpos = args[7]
        self.isize = args[8]
        self.seq = args[9]
        self.qual = args[10]

    @property
    def pairline(self):
        qpos = self.cigar.split('H', 1)[0]
        return "%s:%s\t%s:%s" % (self.qname, qpos, self.rname, self.pos)


class Sam (LineFile):

    def __init__(self, filename, callback=None):

        fp = open(filename)
        for row in fp:
            if row[0]=='@': continue
            s = SamLine(row)
            if callback: callback(s)


def main():

    actions = (
        ('pair', 'parse sam file and get pairs'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pair(args):
    '''
    %prog pair sam_file [--options]

    Parses the sam file and retrieve in pairs format,
    query:pos ref:pos
    '''
    p = OptionParser(pair.__doc__)

    opts, args = p.parse_args(args)
    if len(args)!=1:
        sys.exit(p.print_help())

    def callback(s):
        print s.pairline
    sam = Sam (args[0], callback=callback)


if __name__ == '__main__':
    main()
