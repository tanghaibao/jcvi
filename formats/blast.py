"""
parses tabular BLAST -m8 (-format 6 in BLAST+) format
"""

import sys
import logging
import itertools

from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.coords import print_stats
from jcvi.apps.base import ActionDispatcher, debug
debug()


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


def get_stats(blastfile):
    
    from jcvi.utils.range import range_union

    logging.debug("report stats on `%s`" % blastfile)
    fp = open(blastfile)
    ref_ivs = []
    qry_ivs = []
    identicals = 0
    alignlen = 0

    for row in fp:
        c = BlastLine(row)
        qstart, qstop = c.qstart, c.qstop
        if qstart > qstop: qstart, qstop = qstop, qstart
        qry_ivs.append((c.query, qstart, qstop))

        sstart, sstop = c.sstart, c.sstop
        if sstart > sstop: sstart, sstop = sstop, sstart
        ref_ivs.append((c.subject, sstart, sstop))
        
        alen = sstop - sstart
        alignlen += alen
        identicals += c.pctid / 100. * alen 

    qrycovered = range_union(qry_ivs)
    refcovered = range_union(ref_ivs)
    id_pct = identicals * 100. / alignlen 

    return qrycovered, refcovered, id_pct 


def main():
    
    actions = (
        ('summary', 'provide summary on id% and cov%'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary blastfile 
    
    provide summary on id% and cov%, for both query and reference
    """
    p = OptionParser(summary.__doc__)

    opts, args = p.parse_args(args)

    if len(args)==1:
        blastfile = args[0]
    else:
        sys.exit(p.print_help())

    qrycovered, refcovered, id_pct = get_stats(blastfile)
    print_stats(qrycovered, refcovered, id_pct)


if __name__ == '__main__':
    main()
