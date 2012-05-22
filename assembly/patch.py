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

from jcvi.formats.bed import Bed
from jcvi.formats.fasta import Fasta
from jcvi.utils.range import range_parse
from jcvi.formats.base import must_open, FileMerger
from jcvi.apps.base import ActionDispatcher, debug, sh, mkdir
debug()


def main():

    actions = (
        ('refine', 'find gaps within or near breakpoint regions'),
        ('prepare', 'given om alignment, prepare the patchers'),
        ('install', 'install patches into backbone'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def install(args):
    """
    %prog install patchers.fasta backbone.fasta

    Install one patch into backbone
    """
    from jcvi.apps.command import run_megablast

    p = OptionParser(install.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    patchfasta, bbfasta = args

    blastfile = patchfasta + ".blast"
    run_megablast(infile=patchfasta, outfile=blastfile, db=bbfasta,
                  pctid=99, hitlen=1000)
    #cmd = "blastn -query {0} -subject {1}".format(patchfasta, bbfasta)
    #cmd += " -evalue 0.01 -outfmt 6 -perc_identity 99"


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


def prepare(args):
    """
    %prog prepare om_alignment.bed seq.fasta

    Given optical map alignment, prepare the patchers. Use --backbone to suggest
    which assembly is the major one, and the patchers will be extracted from
    another assembly.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(prepare.__doc__)
    p.add_option("--backbone", default="Scaffold",
                 help="Prefix of the backbone assembly [default: %default]")
    p.add_option("--flank", default=50000, type="int",
                 help="Extend flanks for patchers [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    bb = opts.backbone
    pf = bedfile.split(".")[0]

    # Make a uniq bed keeping backbone at redundant intervals
    bed = Bed(bedfile)
    uniqbed = Bed()
    key = lambda x: (x.seqid, x.start, x.end)
    is_bb = lambda x: x.startswith(bb)
    for k, sb in groupby(bed, key=key):
        sb = list(sb)
        backbone = [x for x in sb if is_bb(x.accn)]
        others = [x for x in sb if not is_bb(x.accn)]
        if backbone and others:
            uniqbed.extend(backbone)
        else:
            uniqbed.extend(sb)

    if not bedfile.endswith(".uniq.bed"):
        uniqbedfile = bedfile.rsplit(".", 1)[0] + ".uniq.bed"
        uniqbed.print_to_file(uniqbedfile)

    # Condense adjacent intervals, allow some chaining
    bed = uniqbed
    key = lambda x: range_parse(x.accn).seqid
    Flank = opts.flank
    sizes = Sizes(fastafile).mapping

    bed_fn = pf + ".patchers.bed"
    bed_fw = open(bed_fn, "w")

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

    fastafn = pf + ".patchers.fasta"
    cmd = "fastaFromBed -fi {0} -bed {1} -fo {2}".\
            format(fastafile, bed_fn, fastafn)
    sh(cmd)


if __name__ == '__main__':
    main()
