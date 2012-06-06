#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Patch the sequences of one assembly using sequences from another assembly. This
is tested on merging the medicago WGS assembly with the clone-by-clone assembly.
"""

import os.path as op
import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.bed import Bed, mergeBed, fastaFromBed
from jcvi.formats.fasta import Fasta
from jcvi.formats.sizes import Sizes
from jcvi.utils.range import range_parse
from jcvi.formats.base import must_open, FileMerger
from jcvi.apps.base import ActionDispatcher, debug, sh, mkdir
debug()


def main():

    actions = (
        ('refine', 'find gaps within or near breakpoint regions'),
        ('fill', 'perform gap filling using one assembly vs the other'),
        ('patcher', 'given om alignment, prepare the patchers'),
        ('install', 'install patches into backbone'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fill(args):
    """
    %prog fill gaps.bed bad.fasta

    Perform gap filling of one assembly (bad) using sequences from another.
    """
    p = OptionParser(fill.__doc__)
    p.add_option("--extend", default=2000, type="int",
                 help="Extend seq flanking the gaps [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gapsbed, badfasta = args
    Ext = opts.extend

    gapsbed = mergeBed(gapsbed, d=Ext, nms=True)

    bed = Bed(gapsbed)
    sizes = Sizes(badfasta).mapping
    pf = gapsbed.rsplit(".", 1)[0]
    extbed = pf + ".ext.bed"
    fw = open(extbed, "w")
    for b in bed:
        gapname = b.accn
        start, end = max(0, b.start - Ext - 1), b.start - 1
        print >> fw, "\t".join(str(x) for x in \
                             (b.seqid, start, end, gapname + "L"))
        start, end = b.end, min(sizes[b.seqid], b.end + Ext)
        print >> fw, "\t".join(str(x) for x in \
                             (b.seqid, start, end, gapname + "R"))
    fw.close()

    fastaFromBed(extbed, badfasta, name=True)


def install(args):
    """
    %prog install patchers.bed patchers.fasta backbone.fasta alt.fasta

    Install patches into backbone, using sequences from alternative assembly.
    The patches sequences are generated via jcvi.assembly.patch.fill().
    """
    from jcvi.apps.base import blast
    from jcvi.formats.blast import BlastSlow

    p = OptionParser(install.__doc__)
    p.add_option("--rclip", default=1, type="int",
            help="Pair ID is derived from rstrip N chars [default: %default]")
    p.add_option("--maxsize", default=250000, type="int",
            help="Maximum size of patchers to be replaced [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    pbed, pfasta, bbfasta, altfasta = args
    Max = opts.maxsize  # Max DNA size to replace gap
    rclip = opts.rclip

    blastfile = blast([altfasta, pfasta,"--wordsize=100", "--pctid=99"])
    order = Bed(pbed).order

    key1 = lambda x: x.query
    key2 = lambda x: x.query[:-rclip] if rclip else key1
    data = BlastSlow(blastfile)

    for pe, lines in groupby(data, key=key2):
        lines = list(lines)
        if len(lines) != 2:
            continue

        a, b = lines

        aquery, bquery = a.query, b.query
        asubject, bsubject = a.subject, b.subject
        if asubject != bsubject:
            continue

        astrand, bstrand = a.orientation, b.orientation
        assert aquery[-1] == 'L' and bquery[-1] == 'R', str((aquery, bquery))

        ai, ax = order[aquery]
        bi, bx = order[bquery]
        qstart, qstop = ax.start + a.qstart - 1, bx.start + b.qstop - 1

        if astrand == '+' and bstrand == '+':
            sstart, sstop = a.sstart, b.sstop

        elif astrand == '-' and bstrand == '-':
            sstart, sstop = b.sstart, a.sstop

        else:
            continue

        if sstart > sstop:
            continue

        if sstop > sstart + Max:
            continue

        name = aquery[:-1] + "LR"
        print "\t".join(str(x) for x in \
                (name, qstart, qstop,
                 asubject, sstart, sstop, astrand))


def refine(args):
    """
    %prog refine breakpoints.bed gaps.bed

    Find gaps within or near breakpoint region.
    """
    p = OptionParser(refine.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    breakpointsbed, gapsbed = args
    cmd = "intersectBed -wao -a {0} -b {1}".format(breakpointsbed, gapsbed)

    pf = "{0}.{1}".format(breakpointsbed.split(".")[0], gapsbed.split(".")[0])
    ingapsbed = pf + ".bed"
    sh(cmd, outfile=ingapsbed)

    fp = open(ingapsbed)
    data = [x.split() for x in fp]
    nogapsbed = pf + ".nogaps.bed"
    largestgapsbed = pf + ".largestgaps.bed"
    nogapsfw = open(nogapsbed, "w")
    largestgapsfw = open(largestgapsbed, "w")
    for b, gaps in groupby(data, key=lambda x: x[:3]):
        gaps = list(gaps)
        if len(gaps) == 1 and gaps[0][3] == ".":
            gap = gaps[0]
            assert gap[4] == "-1"
            print >> nogapsfw, "\t".join(b)
            continue

        gaps = [(int(x[-1]), x) for x in gaps]
        maxgap = max(gaps)[1]
        print >> largestgapsfw, "\t".join(maxgap[3:])

    nogapsfw.close()
    largestgapsfw.close()

    closestgapsbed = pf + ".closestgaps.bed"
    closestgapsfw = open(closestgapsbed, "w")
    cmd = "closestBed -a {0} -b {1}".format(nogapsbed, gapsbed)
    cmd += " | cut -f4-7"
    sh(cmd, outfile=closestgapsbed)

    refinedbed = pf + ".refined.bed"
    FileMerger([largestgapsbed, closestgapsbed], outfile=refinedbed).merge()

    # Clean-up
    cmd = "rm -f " + " ".join([nogapsbed, largestgapsbed, closestgapsbed])
    sh(cmd)


def merge_ranges(beds):

    m = [x.accn for x in beds]

    mr = [range_parse(x) for x in m]
    mc = set(x.seqid for x in mr)
    if len(mc) != 1:
        logging.error("Multiple seqid found in pocket. Aborted.")
        return

    mc = list(mc)[0]
    ms = min(x.start for x in mr)
    me = max(x.end for x in mr)

    neg_strands = sum(1 for x in beds if x.strand == '-')
    pos_strands = len(beds) - neg_strands
    strand = '-' if neg_strands > pos_strands else '+'

    return mc, ms, me, strand


def patcher(args):
    """
    %prog patcher backbone.bed other.bed

    Given optical map alignment, prepare the patchers. Use --backbone to suggest
    which assembly is the major one, and the patchers will be extracted from
    another assembly.
    """
    from jcvi.formats.bed import uniq

    p = OptionParser(patcher.__doc__)
    p.add_option("--backbone", default="OM",
                 help="Prefix of the backbone assembly [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    backbonebed, otherbed = args
    backbonebed = uniq([backbonebed])
    otherbed = uniq([otherbed])

    bb = opts.backbone
    pf = backbonebed.split(".")[0]
    key = lambda x: (x.seqid, x.start, x.end)
    is_bb = lambda x: x.startswith(bb)

    # Make a uniq bed keeping backbone at redundant intervals
    cmd = "intersectBed -v -wa"
    cmd += " -a {0} -b {1}".format(otherbed, backbonebed)
    outfile = otherbed.rsplit(".", 1)[0] + ".not." + backbonebed
    sh(cmd, outfile=outfile)

    uniqbed = Bed()
    uniqbedfile = pf + ".merged.bed"
    uniqbed.extend(Bed(backbonebed))
    uniqbed.extend(Bed(outfile))
    uniqbed.print_to_file(uniqbedfile, sorted=True)

    # Condense adjacent intervals, allow some chaining
    bed = uniqbed
    key = lambda x: range_parse(x.accn).seqid

    bed_fn = pf + ".patchers.bed"
    bed_fw = open(bed_fn, "w")

    for k, sb in groupby(bed, key=key):
        sb = list(sb)
        chr, start, end, strand = merge_ranges(sb)

        id = "{0}:{1}-{2}".format(chr, start, end)
        print >> bed_fw, "\t".join(str(x) for x in \
                (chr, start, end, 1000, strand))

    bed_fw.close()


if __name__ == '__main__':
    main()
