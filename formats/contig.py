#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
TIGR contig format, see spec:

<http://www.cbcb.umd.edu/research/contig_representation.shtml#contig>
"""

import sys

from optparse import OptionParser

from jcvi.formats.base import BaseFile, read_block
from jcvi.apps.base import ActionDispatcher, debug
debug()


class ReadLine (object):

    def __init__(self, row):
        # '#16(0) [RC] 3046 bases, 00000000 checksum. {3046 1} <1 3046>'
        assert row[0] == '#'
        self.id = row.strip("#").split('(')[0]

    def __str__(self):
        return self.id

    __repr__ = __str__


class ContigLine (object):

    def __init__(self, row):
        # '##1 6 8914 bases, 00000000 checksum.'
        assert row[:2] == '##'
        self.id = row.strip("#").split()[0]
        self.reads = []

    def __str__(self):
        return ":".join((self.id, str(self.reads)))

    __repr__ = __str__


class ContigFile (BaseFile):

    def __init__(self, filename):
        super(ContigFile, self).__init__(filename)
        self.fp = open(filename)

    def iter_records(self):
        c = None
        for a, b in read_block(self.fp, "#"):
            if a[:2] == '##':
                if c:
                    yield c
                c = ContigLine(a)
            else:
                c.reads.append(ReadLine(a))
        if c:  # last one
            yield c


def main():
    """
    %prog contigfile

    Prints out the contigs and their associated reads.
    """
    p = OptionParser(main.__doc__)
    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(not p.print_help())

    contigfile, = args
    c = ContigFile(contigfile)
    for rec in c.iter_records():
        print rec


if __name__ == '__main__':
    main()
