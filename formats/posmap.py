#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
POSMAP (POSitional MAPping) files are part of the Celera Assembler output. 

Specs:
http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=POSMAP
"""

import sys
import csv

from collections import namedtuple
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.blast import report_pairs
from jcvi.apps.base import ActionDispatcher, debug
debug()

class Library (object):
    
    def __init__(self, library, minsize, maxsize):
        self.library = library
        self.minsize = minsize
        self.maxsize = maxsize

    def __str__(self):
        return "\t".join(str(x) for x in \
                (self.library, self.minsize, self.maxsize))


class Mate (object):

    def __init__(self, read1, read2, library):
        self.read1 = read1
        self.read2 = read2
        self.library = library

    def __str__(self):
        return "\t".join((self.read1, self.read2, self.library))


class Frags (object):
    pass


MatesLine = namedtuple("MatesLine", 
        "firstReadID secondReadID mateStatus")

class Mates (object):

    def __init__(self, filename):
        fp = csv.reader(open(filename), delimiter='\t')
        for row in fp:
            b = MatesLine._make(row)


class FrgScfLine (object):

    def __init__(self, row):
        atoms = row.split()
        self.fragmentID = atoms[0]
        self.scaffoldID = atoms[1]
        self.begin = int(atoms[2])
        self.end = int(atoms[3])
        self.orientation = '+' if atoms[4]=='f' else '-'


class FrgScf (object):

    def __init__(self, filename):
        fp = csv.reader(open(filename), delimiter='\t')
        for row in fp:
            b = FrgScfLine(row)


class Posmap (LineFile):

    # dispatch based on filename
    mapping = {
            "frags": Frags,
            "mates": Mates,
            "frgscf": FrgScf,
            }

    def __init__(self, filename):
        super(Posmap, self).__init__(filename)


    def parse(self):

        filename = self.filename
        suffix = filename.rsplit(".", 1)[-1]
        assert suffix in self.mapping, \
                "`{0}` unknown format".format(filename)

        # dispatch to the proper handler 
        klass = self.mapping[suffix]
        return klass(filename)


def main():

    actions = (
        ('parse', 'parse posmap file'),
        ('pairs', 'report insert statistics for read pairs')
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def parse(args):
    """
    %prog parse posmap

    """
    p = OptionParser(parse.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    posmapfile = args[0]
    p = Posmap(posmapfile)

    p.parse()


def pairs(args):
    """
    %prog pairs frgscffile 
    
    report summary of the frgscf, how many paired ends mapped, avg
    distance between paired ends, etc. Reads have to be in the form of
    `READNAME{/1,/2}`
    """
    p = OptionParser(pairs.__doc__)
    p.add_option("--cutoff", dest="cutoff", default=1e9, type="int",
            help="distance to call valid links between PE [default: %default]")
    p.add_option("--pairs", dest="pairsfile", 
            default=True, action="store_true",
            help="write valid pairs to pairsfile")
    p.add_option("--inserts", dest="insertsfile", default=True, 
            help="write insert sizes to insertsfile and plot distribution " + \
            "to insertsfile.pdf")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    cutoff = opts.cutoff
    if cutoff < 0: cutoff = 1e9
    frgscffile, = args

    basename = frgscffile.split(".")[0]
    pairsfile = ".".join((basename, "pairs")) if opts.pairsfile else None
    insertsfile = ".".join((basename, "inserts")) if opts.insertsfile else None

    fp = open(frgscffile)
    data = [FrgScfLine(row) for row in fp]
    data.sort(key=lambda x: x.fragmentID)

    report_pairs(data, cutoff, dialect="frgscf", pairsfile=pairsfile,
           insertsfile=insertsfile)


if __name__ == '__main__':
    main()
