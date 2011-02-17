"""
parses tabular BLAST -m8 (-format 6 in BLAST+) format
"""

import sys
import logging

from itertools import groupby
from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.formats.coords import print_stats
from jcvi.utils.range import range_distance
from jcvi.apps.base import ActionDispatcher, debug
debug()


class BlastLine(object):
    __slots__ = ('query', 'subject', 'pctid', 'hitlen', 'nmismatch', 'ngaps', \
                 'qstart', 'qstop', 'sstart', 'sstop', 'evalue', 'score', \
                 'qseqid', 'sseqid', 'qi', 'si', 'qstrand', 'sstrand')
 
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

        self.qstrand = self.sstrand = '+'
        if self.sstart > self.sstop:
            self.sstrand = '-'
            self.sstart, self.sstop = self.sstop, self.sstart
 
    def __repr__(self):
        return "BlastLine('%s' to '%s', eval=%.3f, score=%.1f)" % \
                (self.query, self.subject, self.evalue, self.score)

    def __str__(self):
        return "\t".join(str(x) for x in \
                [getattr(self, attr) for attr in BlastLine.__slots__[:-4]])


class Blast (LineFile):
    """
    We can have a Blast class that loads entire file into memory, this is
    not very efficient for big files; when the BLAST file is generated
    by BLAST or BLAT, the file is already sorted, with same query grouped
    """
    def __init__(self, filename):
        super(Blast, self).__init__(filename)
        self.fp = open(filename)

    def iter_hits(self):
        self.fp.seek(0)
        for query, blines in groupby(self.fp, key=lambda x: BlastLine(x).query):
            blines = [BlastLine(x) for x in blines]
            blines.sort(key=lambda x: -x.score) # descending score
            yield query, blines 

    def iter_best_hit(self, N=1):
        self.fp.seek(0)
        for query, blines in groupby(self.fp, key=lambda x: BlastLine(x).query):
            blines = [BlastLine(x) for x in blines]
            blines.sort(key=lambda x: -x.score)
            for x in blines[:N]:
                yield query, x 
    
    @property
    def hits(self):
        """
        returns a dict with query => blastline
        """
        return dict(self.iter_hits())

    @property
    def best_hits(self):
        """
        returns a dict with query => best blasthit
        """
        return dict(self.iter_best_hit())


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


def filter(args):
    """
    %prog filter test.blast > new.blast

    produce a new blast file and filter based on score
    """
    p = OptionParser(filter.__doc__)
    p.add_option("--score", dest="score", default=0., type="float",
            help="score cutoff [default: %default]")
    p.add_option("--pctid", dest="pctid", default=0., type="float",
            help="pctid cutoff [default: %default]")
    p.add_option("--hitlen", dest="hitlen", default=0., type="float",
            help="pctid cutoff [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    fp = open(args[0])
    for row in fp:
        c = BlastLine(row)

        if c.score < opts.score: continue
        if c.pctid < opts.pctid: continue
        if c.hitlen < opts.hitlen: continue

        print row.rstrip()


def main():
    
    actions = (
        ('summary', 'provide summary on id% and cov%'),
        ('filter', 'filter BLAST file (based on e.g. score)'),
        ('best', 'get best BLAST hit per query'),
        ('pairs', 'print paired-end reads of BLAST tabular output'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def report_pairs(data, cutoff=300000, dialect="blast", print_pairs=False,
        print_inserts=False):
    """
    This subroutine is used by the pairs function in blast.py and cas.py.
    Reports number of fragments and pairs as well as linked pairs
    """
    dialect = dialect[0]
    num_fragments, num_pairs = 0, 0
    linked_dist = []
    # +- (forward-backward) is `innie`, -+ (backward-forward) is `outie`
    orientations = defaultdict(int)

    if dialect=="b":
        key = lambda x: x.query.split("/")[0] 
    else:
        key = lambda x: x.readname.split("/")[0] 

    for pe, lines in groupby(data, key=key):   
        lines = list(lines)
        if len(lines)==2: 
            num_pairs += 1
            a, b = lines

            if dialect=="b":
                asubject, astart, astop = a.subject, a.sstart, a.sstop
                bsubject, bstart, bstop = b.subject, b.sstart, b.sstop
     
                astrand = a.sstrand
                bstrand = b.sstrand

            else:
                asubject, astart, astop = a.refnum, a.refstart, a.refstop
                bsubject, bstart, bstop = b.refnum, b.refstart, b.refstop
                if -1 in (astart, bstart): continue

                astrand = a.strand
                bstrand = b.strand

            dist, orientation = range_distance(\
                    (asubject, astart, astop, astrand), 
                    (bsubject, bstart, bstop, bstrand))

            if 0 <= dist <= cutoff:
                linked_dist.append(dist)
                if print_pairs:
                    for b in lines: print b
                orientations[orientation] += 1
        else:
            num_fragments += 1
    
    import numpy as np
    print >>sys.stderr, "%d fragments, %d pairs" % (num_fragments, num_pairs)
    num_links = len(linked_dist)
    print >>sys.stderr, "%d pairs (%.1f%%) are linked (cutoff=%d)" % \
            (num_links, num_links*100./num_pairs, cutoff)
    print >>sys.stderr, "median distance between PE: %d" % np.median(linked_dist) 
    print >>sys.stderr, "\nOrientations:"
    for orientation, count in sorted(orientations.items()):
        print >>sys.stderr, "{0}: {1}".format(orientation, count)

    if print_inserts:
        print "\n".join(str(x) for x in linked_dist)


def pairs(args):
    """
    %prog pairs blastfile 
    
    report summary of the cas tabular results, how many paired ends mapped, avg
    distance between paired ends, etc. Reads have to be in the form of
    `READNAME{/1,/2}`
    """
    p = OptionParser(pairs.__doc__)
    p.add_option("--cutoff", dest="cutoff", default=0, type="int",
            help="distance to call valid links between PE [default: %default]")
    p.add_option("--pairs", dest="pairs", default=False, action="store_true",
            help="write valid pairs to stdout [default: %default]")
    p.add_option("--inserts", dest="inserts", default=False, action="store_true",
            help="write insert sizes to stdout [default: %default]")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    cutoff = opts.cutoff
    if cutoff <= 0: cutoff = 1e10
    print_pairs = opts.pairs
    print_inserts = opts.inserts

    blastfile = args[0]
    fp = open(blastfile)
    data = [BlastLine(row) for row in fp]
    data.sort(key=lambda x: x.query)

    report_pairs(data, cutoff, dialect="blast", print_pairs=print_pairs,
            print_inserts=print_inserts)

    
def best(args):
    """
    %prog best blastfile

    print the best hit for each query in the blastfile
    """
    p = OptionParser(best.__doc__)

    p.add_option("-N", dest="N", default=1, type="int",
            help="get best N hits [default: %default]")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    blastfile = args[0]
    b = Blast(blastfile)
    for q, bline in b.iter_best_hit(N=opts.N):
        print bline


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
