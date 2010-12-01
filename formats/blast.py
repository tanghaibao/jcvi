
"""
parses tabular BLAST -m8 (-format 6 in BLAST+) format
"""

import itertools
from base import LineFile


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si')
 
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
        self.evalue = float(args[10])
        self.score = float(args[11])
 
    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(map(str, [getattr(self, attr) \
                for attr in BlastLine.__slots__[:-4]]))


class Blast (LineFile):
    """
    Collection of BlastLine
    """
    def __init__(self, filename):
        super(Blast, self).__init__(filename)

        fp = open(filename)
        for row in fp:
            self.append(BlastLine(row))

        self.sort(key=lambda x: (x.query, -x.score))

    @property
    def hits(self):
        """
        returns a dict with query => blastline
        """
        return dict((query, list(blines)) for (query, blines) in \
                itertools.groupby(self, lambda x: x.query))

    @property
    def best_hits(self):
        """
        returns a dict with query => best blasthit
        """
        return dict((query, blines.next()) for (query, blines) in \
                itertools.groupby(self, lambda x: x.query))

