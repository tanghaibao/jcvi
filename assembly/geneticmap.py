#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use genetic map to break chimeric scaffolds, order and orient scaffolds onto
chromosomes.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.base import BaseFile, LineFile, must_open, read_block
from jcvi.apps.base import ActionDispatcher, debug
debug()


class BinMap (BaseFile, dict):

    def __init__(self, filename):
        super(BinMap, self).__init__(filename)

        fp = open(filename)
        for header, seq in read_block(fp, "group "):
            lg = header.split()[-1]
            self[lg] = []
            for s in seq:
                if s.strip() == '' or s[0] == ';':
                    continue
                marker, pos = s.split()
                pos = int(float(pos) * 1000)
                self[lg].append((marker, pos))

    def print_to_bed(self, filename="stdout"):
        fw = must_open(filename, "w")
        for lg, markers in sorted(self.items()):
            for marker, pos in markers:
                print >> fw, "\t".join(str(x) for x in \
                        (lg, pos, pos + 1, marker))
        fw.close()


class MSTMapLine (object):

    def __init__(self, row):
        args = row.split()
        self.id = args[0]
        self.seqid, pos = self.id.split(".")
        self.pos = int(pos)
        self.genotype = "".join(args[1:])

    def __len__(self):
        return len(self.genotype)

    def __str__(self):
        return "{0}: {1}".format(self.id, self.genotype)

    @property
    def bedline(self):
        return "\t".join(str(x) for x in \
                (self.seqid, self.pos - 1, self.pos, self.id))


class MSTMap (LineFile):

    def __init__(self, filename):
        super(MSTMap, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row.startswith("locus_name"):
                break

        for row in fp:
            self.append(MSTMapLine(row))


def hamming_distance(a, b, ignore=None):
    dist = 0
    for x, y in zip(a, b):
        if ignore and ignore in (x, y):
            continue
        if x != y:
            dist += 1
    return dist


def main():

    actions = (
        ('breakpoint', 'find scaffold breakpoints using genetic map'),
        ('fasta', 'extract markers based on map'),
        ('placeone', 'attempt to place one scaffold'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fasta(args):
    """
    %prog fasta map.out scaffolds.fasta markers.fasta

    Extract marker sequences based on map.
    """
    from jcvi.formats.bed import Bed, fastaFromBed
    from jcvi.formats.sizes import Sizes

    p = OptionParser(fasta.__doc__)
    p.add_option("--extend", default=1000, type="int",
                 help="Extend seq flanking the gaps [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mapout, sfasta = args
    Flank = opts.extend
    pf = mapout.split(".")[0]
    mapbed = pf + ".bed"
    bm = BinMap(mapout)
    bm.print_to_bed(mapbed)

    bed = Bed(mapbed, sorted=False)
    markersbed = pf + ".markers.bed"
    fw = open(markersbed, "w")
    sizes = Sizes(sfasta).mapping
    for b in bed:
        accn = b.accn
        scf, pos = accn.split(".")
        pos = int(pos)
        start = max(0, pos - Flank)
        end = min(pos + Flank, sizes[scf])
        print >> fw, "\t".join(str(x) for x in \
                    (scf, start, end, accn))

    fw.close()

    fastaFromBed(markersbed, sfasta, name=True)


def place(genotype, map, exclude):
    d = [(hamming_distance(genotype, x.genotype, ignore="-"), x) \
            for x in map if x.seqid != exclude]
    return min(d)


def placeone(args):
    """
    %prog placeone map scaffold_id

    Attempt to place one scaffold.
    """
    p = OptionParser(placeone.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mstmap, scf = args
    data = MSTMap(mstmap)
    for d in data:
        if d.seqid != scf:
            continue

        dist, besthit = place(d.genotype, data, scf)
        width = len(d) + len(d.id) + 10
        print >> sys.stderr, "distance = {0}".format(dist)
        print >> sys.stderr, str(d).rjust(width)
        print >> sys.stderr, str(besthit).rjust(width)
        print >> sys.stderr, "-" * width


OK, BREAK, END = range(3)

def check_markers(a, b, maxdiff):

    if a.seqid != b.seqid:
        return END, None
    diff = hamming_distance(a.genotype, b.genotype, ignore="-")
    max_allowed = len(a) * maxdiff
    if diff <= max_allowed:
        return OK, None

    return BREAK, (a.seqid, a.pos, b.pos)


def breakpoint(args):
    """
    %prog breakpoint mstmap.input > breakpoints.bed

    Find scaffold breakpoints using genetic map. Use formats.vcf.mstmap() to
    generate the input for this routine.
    """
    from jcvi.utils.iter import pairwise

    p = OptionParser(breakpoint.__doc__)
    p.add_option("--diff", default=.1, type="float",
                 help="Maximum ratio of differences allowed [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    mstmap, = args
    diff = opts.diff
    data = MSTMap(mstmap)

    # Remove singleton markers (avoid double cross-over)
    good = []
    nsingletons = 0
    for i in xrange(1, len(data) - 1):
        a = data[i]
        left_label, left_rr = check_markers(data[i - 1], a, diff)
        right_label, right_rr = check_markers(a, data[i + 1], diff)

        if left_label == BREAK and right_label == BREAK:
            nsingletons += 1
            continue

        good.append(a)

    logging.debug("A total of {0} singleton markers removed.".format(nsingletons))

    for a, b in pairwise(good):
        label, rr = check_markers(a, b, diff)
        if label == BREAK:
            print "\t".join(str(x) for x in rr)


if __name__ == '__main__':
    main()
