#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
POSMAP (POSitional MAPping) files are part of the Celera Assembler output.

Specs:
http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=POSMAP
"""

import os.path as op
import sys
import csv
import logging

from collections import namedtuple, defaultdict
from optparse import OptionParser
from itertools import groupby

from jcvi.formats.base import LineFile
from jcvi.formats.blast import set_options_pairs, report_pairs
from jcvi.apps.base import ActionDispatcher, sh, debug
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
        self.begin = int(atoms[2]) + 1  # convert to 1-based
        self.end = int(atoms[3])
        self.orientation = '+' if atoms[4] == 'f' else '-'

    @property
    def bedline(self):
        s = '\t'.join(str(x) for x in (self.scaffoldID, self.begin - 1, self.end,
            self.fragmentID, "na", self.orientation))
        return s


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
        ('bed', 'convert to bed format'),
        ('index', 'index the frgscf.posmap.sorted file'),
        ('query', 'query region from frgscf index'),
        ('dup', 'estimate level of redundancy based on position collision'),
        ('reads', 'report read counts per scaffold (based on frgscf)'),
        ('pairs', 'report insert statistics for read pairs')
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def index(args):
    """
    %prog index frgscf.sorted

    Compress frgscffile.sorted and index it using `tabix`.
    """
    p = OptionParser(index.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    frgscffile, = args
    gzfile = frgscffile + ".gz"
    cmd = "bgzip -c {0}".format(frgscffile)

    if not op.exists(gzfile):
        sh(cmd, outfile=gzfile)

    tbifile = gzfile + ".tbi"
    # Sequence, begin, end in 2, 3, 4-th column, respectively
    cmd = "tabix -s 2 -b 3 -e 4 {0}".format(gzfile)

    if not op.exists(tbifile):
        sh(cmd)


def query(args):
    """
    %prog query frgscf.sorted scfID:start-end

    Query certain region to get frg placement, using random access. Build index
    if not present.
    """
    p = OptionParser(query.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    frgscffile, region = args
    gzfile = frgscffile + ".gz"
    tbifile = gzfile + ".tbi"
    outfile = region + ".posmap"

    if not (op.exists(gzfile) and op.exists(tbifile)):
        index(frgscffile)

    assert op.exists(gzfile) and op.exists(tbifile)

    cmd = "tabix {0} {1}".format(gzfile, region)
    sh(cmd, outfile=outfile)


def reads(args):
    """
    %prog reads frgscffile

    Report read counts per scaffold (based on frgscf).
    """
    p = OptionParser(reads.__doc__)
    p.add_option("-p", dest="prefix_length", default=4, type="int",
            help="group the reads based on the first N chars [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    frgscffile, = args
    prefix_length = opts.prefix_length

    fp = open(frgscffile)
    keyfn = lambda: defaultdict(int)
    counts = defaultdict(keyfn)
    for row in fp:
        f = FrgScfLine(row)
        fi = f.fragmentID[:prefix_length]
        counts[f.scaffoldID][fi] += 1

    for scf, count in sorted(counts.items()):
        print "{0}\t{1}".format(scf,
                ", ".join("{0}:{1}".format(*x) for x in sorted(count.items())))


def bed(args):
    """
    %prog bed frgscffile

    Convert the frgscf posmap file to bed format.
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    frgscffile, = args
    bedfile = frgscffile.rsplit(".", 1)[0] + ".bed"
    fw = open(bedfile, "w")

    fp = open(frgscffile)
    for row in fp:
        f = FrgScfLine(row)
        print >> fw, f.bedline

    logging.debug("File written to {0}.".format(bedfile))


def dup(args):
    """
    %prog dup frgscffile

    Use the frgscf posmap file as an indication of the coverage of the library.
    Large insert libraries are frequently victims of high levels of redundancy.
    """
    p = OptionParser(dup.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    frgscffile, = args

    fp = open(frgscffile)
    data = [FrgScfLine(row) for row in fp]
    # we need to separate forward and reverse reads, because the position
    # collisions are handled differently
    forward_data = [x for x in data if x.orientation == '+']
    reverse_data = [x for x in data if x.orientation == '-']

    counts = defaultdict(int)
    key = lambda x: (x.scaffoldID, x.begin)
    forward_data.sort(key=key)
    for k, data in groupby(forward_data, key=key):
        data = list(data)
        count = len(data)
        counts[count] += 1

    key = lambda x: (x.scaffoldID, x.end)
    reverse_data.sort(key=key)
    for k, data in groupby(forward_data, key=key):
        data = list(data)
        count = len(data)
        counts[count] += 1

    prefix = frgscffile.split(".")[0]
    print >> sys.stderr, "Duplication level in `{0}`".format(prefix)
    print >> sys.stderr, "=" * 40
    for c, v in sorted(counts.items()):
        if c > 10:
            break
        label = "unique" if c == 1 else "{0} copies".format(c)
        print >> sys.stderr, "{0}: {1}".format(label, v)


def pairs(args):
    """
    See __doc__ for set_options_pairs().
    """
    p = set_options_pairs()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    frgscffile, = args

    basename = frgscffile.split(".")[0]
    pairsfile = ".".join((basename, "pairs")) if opts.pairsfile else None
    insertsfile = ".".join((basename, "inserts")) if opts.insertsfile else None

    fp = open(frgscffile)
    data = [FrgScfLine(row) for i, row in enumerate(fp) if i < opts.nrows]

    ascii = not opts.pdf
    return report_pairs(data, opts.cutoff, opts.mateorientation,
           dialect="frgscf", pairsfile=pairsfile, insertsfile=insertsfile,
           rclip=opts.rclip, ascii=ascii, bins=opts.bins)


if __name__ == '__main__':
    main()
