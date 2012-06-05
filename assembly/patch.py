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
    %prog install patch.blast.bed backbone.fasta alt.fasta

    Install patches into backbone, using sequences from alternative assembly.
    """
    p = OptionParser(install.__doc__)
    p.add_option("--rclip", default=1, type="int",
            help="pair ID is derived from rstrip N chars [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, bbfasta, altfasta = args
    rclip = opts.rclip
    data = Bed(bedfile)

    key1 = lambda x: x.accn
    key2 = lambda x: x.accn[:-rclip] if rclip else key1
    data.sort(key=key1)
    Max = 200000  # Max DNA size to replace gap

    for pe, lines in groupby(data, key=key2):
        lines = list(lines)
        if len(lines) != 2:
            continue

        a, b = lines

        asubject, astart, astop = a.seqid, a.start, a.end
        bsubject, bstart, bstop = b.seqid, b.start, b.end

        if asubject != bsubject:
            continue

        aquery, bquery = a.accn, b.accn
        astrand, bstrand = a.strand, b.strand
        assert aquery[-1] == 'L' and bquery[-1] == 'R', str((aquery, bquery))

        """
        Case 1: L matches +, R matches +, L_stop < R_start
        LLLLL.......RRRRR
        >>>>>>>>>>>>>>>>>

        Case 2: L matches -, R matches -, L_stop < R_start
        LLLLL.......RRRRR
        <<<<<<<<<<<<<<<<<
        """
        if astrand == '+' and bstrand == '+':
            start, end = astop + 1, bstart - 1

        elif astrand == '-' and bstrand == '-':
            start, end = bstop + 1, astart - 1

        else:
            continue

        if start > end:
            continue

        if end > start + Max:
            continue

        name = aquery[:-1] + "LR"
        print "\t".join(str(x) for x in (asubject, start - 1, end, name, astrand))



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
    %prog patcher backbone.bed other.bed other.fasta

    Given optical map alignment, prepare the patchers. Use --backbone to suggest
    which assembly is the major one, and the patchers will be extracted from
    another assembly.
    """
    from jcvi.formats.bed import uniq

    p = OptionParser(patcher.__doc__)
    p.add_option("--backbone", default="OM",
                 help="Prefix of the backbone assembly [default: %default]")
    p.add_option("--flank", default=50000, type="int",
                 help="Extend flanks for patchers [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    backbonebed, otherbed, fastafile = args
    backbonebed = uniq([backbonebed])
    otherbed = uniq([otherbed])

    bb = opts.backbone
    pf = backbonebed.split(".")[0]
    key = lambda x: (x.seqid, x.start, x.end)
    is_bb = lambda x: x.startswith(bb)

    # Make a uniq bed keeping backbone at redundant intervals
    cmd = "intersectBed -v"
    cmd += " -a {0} -b {1}".format(otherbed, backbonebed)
    outfile = otherbed.rsplit(".", 1)[0] + ".not." + backbonebed
    sh(cmd, outfile=outfile)

    uniqbed = Bed()
    uniqbedfile = pf + ".merged.bed"
    uniqbed.extend(Bed(backbonebed))
    uniqbed.extend(Bed(outfile))
    uniqbed.sort(key=uniqbed.nullkey)
    uniqbed.print_to_file(uniqbedfile)

    # Condense adjacent intervals, allow some chaining
    bed = uniqbed
    key = lambda x: range_parse(x.accn).seqid
    Flank = opts.flank
    sizes = Sizes(fastafile).mapping

    bed_fn = pf + ".patchers.bed"
    bed_fw = open(bed_fn, "w")

    Flank = 0
    for k, sb in groupby(bed, key=key):
        sb = list(sb)
        chr, start, end, strand = merge_ranges(sb)
        size = sizes[chr]
        start = max(start - Flank, 0)
        end = min(end + Flank, size)
        if is_bb(chr):
            continue

        id = "{0}:{1}-{2}".format(chr, start, end)
        print >> bed_fw, "\t".join(str(x) for x in (chr, start, end))

    bed_fw.close()

    #fastaFromBed(bed_fn, fastafile)


if __name__ == '__main__':
    main()
