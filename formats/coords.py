#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
parses JCVI software NUCMER (http://mummer.sourceforge.net/manual/)
output - mostly as *.coords file.

when run as commandline,
$ python -m jcvi.formats.coords test.coords quality_cutoff

will produce a list of BACs, mapped position and orientation (needs to 
be beyond quality cutoff, say 50) in bed format
"""

import sys
import itertools
import logging

from jcvi.formats.base import LineFile


class CoordsLine (object):

    """
    The coords line looks like (in one line):
        2953     4450  |      525     2023  |     1498     1499  |    98.07  |
        8046     2023  |    18.62    74.10  | gi|270341414|gb|AC182814.30|_2    contig_100476
    """
    def __init__(self, row):
        
        row = row.replace(" | ", "")
        atoms = row.split()
        assert len(atoms)==13, "expecting 13 columns"
        
        self.start1 = int(atoms[0])
        self.end1 = int(atoms[1])
        
        self.start2 = int(atoms[2])
        self.end2 = int(atoms[3])
        
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
        score = int(self.quality * 10)
        return '\t'.join((self.ref, str(self.start1-1), str(self.end1), 
                self.query, str(score), self.orientation))

    @property
    def orientation(self):
        """
        the orientation of the alignment between query and ref (ref is always
        plus strand)
        """
        return '+' if self.end2 >= self.start2 else '-'


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
        
        #for x in xrange(4): fp.next()
        for row in fp:
            try:
                self.append(CoordsLine(row))
            except AssertionError, e:
                pass

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
        coordsfile = args[0]
        if len(args) > 1:
            quality_cutoff = int(args[1])
        else:
            quality_cutoff = 50
    except Exception, e:
        logging.error(str(e))
        sys.exit(p.print_help())

    coords = Coords(coordsfile)

    for c in coords:
        if c.quality < quality_cutoff: continue
        print c
