#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Python implementation of BlastLine, an alternative Cython implementation is
available in .cblast.BlastLine, which may be up to 2x faster
"""

class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si', 'orientation')

    def __init__(self, sline):
        args = sline.split("\t")
        self.query = args[0]
        self.subject = args[1]
        self.pctid = float(args[2])
        self.hitlen = int(args[3])
        self.nmismatch = int(args[4])
        self.ngaps = int(args[5])
        self.qstart = int(args[6])
        self.qstop = int(args[7])
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        if len(args) > 10:
            self.evalue = float(args[10])
            self.score = float(args[11])

        self.orientation = '+'
        if self.qstart > self.qstop:
            self.qstart, self.qstop = self.qstop, self.qstart
            self.orientation = '-'
        if self.sstart > self.sstop:
            self.sstart, self.sstop = self.sstop, self.sstart
            self.orientation = '-'

    @property
    def has_score(self):
        return hasattr(self, "score")

    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        if self.has_score:
            args = [getattr(self, attr) for attr in BlastLine.__slots__[:12]]
        else:
            args = [getattr(self, attr) for attr in BlastLine.__slots__[:10]]
        if self.orientation == '-':
            args[8], args[9] = args[9], args[8]
        return "\t".join(str(x) for x in args)

    @property
    def swapped(self):
        """
        Swap query and subject.
        """
        args = [getattr(self, attr) for attr in BlastLine.__slots__[:12]]
        args[0:2] = [self.subject, self.query]
        args[6:10] = [self.sstart, self.sstop, self.qstart, self.qstop]
        if self.orientation == '-':
            args[8], args[9] = args[9], args[8]
        b = "\t".join(str(x) for x in args)
        return BlastLine(b)

    @property
    def bedline(self):
        return "\t".join(str(x) for x in \
                (self.subject, self.sstart - 1, self.sstop, self.query,
                 self.score, self.orientation))


