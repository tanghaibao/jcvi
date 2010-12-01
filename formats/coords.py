#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
parses JCVI software NUCMER (http://mummer.sourceforge.net/manual/)
output - mostly as *.coords file.

when run as commandline,
$ python -m jcvi.formats.coords test.coords

will produce a list of BACs, mapped position and orientation (needs to 
be >95% query coverage) in bed format
"""

import sys
import itertools
import logging

from base import LineFile


class CoordsLine (object):

    __slots__ = ("start1", "end1", "start2", "end2", "len1", "len2", "identity",
            "reflen", "querylen", "refcov", "querycov", "ref", "query",
            "quality", "orientation")

    """
    151317 151396 | 79 1 | 80 79 | 98.75 | Mt3.5.1Chr1 contig_173759
    """
    def __init__(self, row):
        
        row = row.replace("|", "")
        atoms = row.split()
        
        self.start1 = int(atoms[0])
        self.end1 = int(atoms[1])
        
        self.start2 = int(atoms[2])
        self.end2 = int(atoms[3])
        self.orientation = '+' if self.end2 >= self.start2 else '-'
        
        self.len1 = int(atoms[4])
        self.len2 = int(atoms[5])

        self.identity = float(atoms[6])
        
        self.reflen = int(atoms[7])
        self.querylen = int(atoms[8])
        
        self.refcov = float(atoms[9]) / 100.
        self.querycov = float(atoms[10]) / 100.
        # this is taken from CoGeBlast:
        # the coverage of the hit muliplied by percent seq identity
        # range from 0-100
        self.quality = self.identity * self.querycov
        
        self.ref = atoms[11]
        self.query = atoms[12]

    def __str__(self):
        # bed formatted line
        return '\t'.join((self.ref, str(self.start1), str(self.end1), 
                self.query, self.orientation))


class Coords (LineFile):
    
    """
    when parsing the .coords file, first skip first 5 lines
    [S1] [E1] | [S2] [E2] | [LEN 1] [LEN 2] | [% IDY] | [TAGS]

    then each row would be composed as this
    """
    def __init__(self, filename):
        super(Coords, self).__init__(filename)

        fp = open(filename)
        self.cmd = fp.next()
        
        for x in xrange(4): fp.next()
        for row in fp:
            self.append(CoordsLine(row))

    @property
    def hits(self):
        """
        returns a dict with query => blastline
        """
        # sort descending with score=identity * coverage
        self.sort(key=lambda x: (x.query, -x.quality))

        hits = dict((query, list(blines)) for (query, blines) in \
                itertools.groupby(self, lambda x: x.query))

        # sort back to order by reference positions
        self.sort(key=lambda x: (x.ref, x.start1))

        return hits


    @property
    def best_hits(self):
        """
        returns a dict with query => best mapped position 
        """
        # sort descending with score=identity * coverage
        self.sort(key=lambda x: (x.query, -x.quality))

        best_hits = dict((query, blines.next()) for (query, blines) in \
                itertools.groupby(self, lambda x: x.query))
        
        # sort back to order by reference positions
        self.sort(key=lambda x: (x.ref, x.start1))

        return best_hits


if __name__ == '__main__':
    
    from optparse import OptionParser

    p = OptionParser(__doc__)

    opts, args = p.parse_args()

    try:
        coordsfile, = args
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    coords = Coords(coordsfile)

    for c in coords:
        if c.quality < 80: continue
        print c
