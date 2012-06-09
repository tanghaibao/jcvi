#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Patch the sequences of one assembly using sequences from another assembly. This
is tested on merging the medicago WGS assembly with the clone-by-clone assembly.

There are a few techniques, used in curating medicago assembly.

1. Split chimeric scaffolds based on genetic map and then refine breakpoints
2. Create patchers by mix-and-max guided by optical map
3. Find gaps and fill N's using alternative assembly
4. Add telomeric sequences
5. Find gaps in optical map
6. Insert unplaced scaffolds using mates
"""

import os.path as op
import sys
import logging

from itertools import groupby
from optparse import OptionParser
from collections import defaultdict

from jcvi.formats.bed import Bed, BedLine, complementBed, mergeBed, fastaFromBed
from jcvi.formats.fasta import Fasta
from jcvi.formats.sizes import Sizes
from jcvi.utils.range import range_parse, range_distance, ranges_depth, \
            range_minmax, range_overlap
from jcvi.formats.base import must_open, FileMerger, FileShredder
from jcvi.apps.base import ActionDispatcher, debug, sh, mkdir
debug()


def main():

    actions = (
        ('refine', 'find gaps within or near breakpoint regions'),
        ('patcher', 'given om alignment, prepare the patchers'),

        ('fill', 'perform gap filling using one assembly vs the other'),
        ('install', 'install patches into backbone'),

        ('tips', 'append telomeric sequences based on patchers and complements'),
        ('gaps', 'create patches around OM gaps'),

        ('bambus', 'find candidate scaffolds to insert based on mates'),
        ('insert', 'insert scaffolds into assembly'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def insert(args):
    """
    %prog insert candidates.bed gaps.bed chrs.fasta unplaced.fasta

    Insert scaffolds into assembly.
    """
    from jcvi.formats.agp import mask, bed
    from jcvi.formats.sizes import agp

    p = OptionParser(insert.__doc__)
    p.add_option("--prefix", default="unplaced",
                 help="Prefix of the unplaced scaffolds [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    candidates, gapsbed, chrfasta, unplacedfasta = args
    refinedbed = refine([candidates, gapsbed])
    sizes = Sizes(unplacedfasta).mapping
    cbed = Bed(candidates)
    corder = cbed.order
    gbed = Bed(gapsbed)
    gorder = gbed.order

    gpbed = Bed()
    gappositions = {}  # (chr, start, end) => gapid

    fp = open(refinedbed)
    gap_to_scf = defaultdict(list)
    seen = set()
    for row in fp:
        atoms = row.split()
        unplaced = atoms[3]
        strand = atoms[5]
        gapid = atoms[9]
        if gapid not in seen:
            seen.add(gapid)
            gi, gb = gorder[gapid]
            gpbed.append(gb)
            gappositions[(gb.seqid, gb.start, gb.end)] = gapid
        gap_to_scf[gapid].append((unplaced, strand))

    gpbedfile = "candidate.gaps.bed"
    gpbed.print_to_file(gpbedfile, sorted=True)

    agpfile = agp([chrfasta])
    maskedagpfile = mask([agpfile, gpbedfile])
    maskedbedfile = maskedagpfile.rsplit(".", 1)[0] + ".bed"
    bed([maskedagpfile, "--outfile={0}".format(maskedbedfile)])

    mbed = Bed(maskedbedfile)
    beds = []
    for b in mbed:
        sid = b.seqid
        key = (sid, b.start, b.end)
        if key not in gappositions:
            beds.append(b)
            continue

        gapid = gappositions[key]
        scfs = gap_to_scf[gapid]

        # For scaffolds placed in the same gap, sort according to positions
        scfs.sort(key=lambda x: corder[x[0]][1].start + corder[x[0]][1].end)
        for scf, strand in scfs:
            size = sizes[scf]
            beds.append(BedLine("\t".join(str(x) for x in \
                    (scf, 0, size, sid, 1000, strand))))

    finalbed = Bed()
    finalbed.extend(beds)
    finalbedfile = "final.bed"
    finalbed.print_to_file(finalbedfile)

    # Clean-up
    toclean = [gpbedfile, agpfile, maskedagpfile, maskedbedfile]
    FileShredder(toclean)


def bambus(args):
    """
    %prog bambus bambus.bed bambus.mates total.fasta

    Insert unplaced scaffolds based on mates.
    """
    from jcvi.utils.iter import pairwise
    from jcvi.formats.posmap import MatesFile

    p = OptionParser(bambus.__doc__)
    p.add_option("--prefix", default="unplaced",
                 help="Prefix of the unplaced scaffolds [default: %default]")
    p.add_option("--minlinks", default=3, type="int",
                 help="Minimum number of links to place [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, matesfile, fastafile = args
    pf = matesfile.rsplit(".", 1)[0]
    logfile = pf + ".log"
    log = open(logfile, "w")

    mf = MatesFile(matesfile)
    maxdist = max(x.max for x in mf.libraries.values())
    logging.debug("Max separation: {0}".format(maxdist))

    prefix = opts.prefix
    minlinks = opts.minlinks

    is_unplaced = lambda x: x.startswith(prefix)
    bed = Bed(bedfile, sorted=False)
    beds = []
    unplaced = defaultdict(list)

    for a, b in pairwise(bed):
        aname, bname = a.accn, b.accn
        aseqid, bseqid = a.seqid, b.seqid

        pa, la = mf[aname]
        if pa != bname:
            continue

        ia = is_unplaced(aseqid)
        ib = is_unplaced(bseqid)
        if ia == ib:
            continue

        if ia:
            a, b = b, a

        unplaced[b.seqid].append((a, b))
        beds.extend([a, b])

    sizes = Sizes(fastafile)
    candidatebed = Bed()
    cbeds = []
    # For each unplaced scaffold, find most likely placement and orientation
    for scf, beds in sorted(unplaced.items()):
        print >> log
        ranges = []
        for a, b in beds:
            aname, astrand = a.accn, a.strand
            bname, bstrand = b.accn, b.strand
            aseqid, bseqid = a.seqid, b.seqid
            pa, lib = mf[aname]

            print >> log, a
            print >> log, b

            flip_b = (astrand == bstrand)
            fbstrand = '-' if flip_b else '+'
            if flip_b:
                b.reverse_complement(sizes)

            lmin, lmax = lib.min, lib.max

            L = sizes.get_size(scf)
            assert astrand in ('+', '-')
            if astrand == '+':
                offset = a.start - b.end
                sstart, sstop = offset + lmin, offset + lmax
            else:
                offset = a.end - b.start + L
                sstart, sstop = offset - lmax, offset - lmin

            # Prevent out of range error
            size = sizes.get_size(aseqid)
            sstart = max(0, sstart)
            sstop = max(0, sstop)
            sstart = min(size - 1, sstart)
            sstop = min(size - 1, sstop)

            start_range = (aseqid, sstart, sstop, scf, 1, fbstrand)
            print >> log, "*" + "\t".join(str(x) for x in start_range)
            ranges.append(start_range)

        mranges = [x[:3] for x in ranges]
        # Determine placement by finding the interval with the most support
        rd = ranges_depth(mranges, sizes.mapping, verbose=False)
        alldepths = []
        for depth in rd:
            alldepths.extend(depth)
        print >> log, alldepths

        maxdepth = max(alldepths, key=lambda x: x[-1])[-1]
        if maxdepth < minlinks:
            print >> log, "Insufficient links ({0} < {1})".format(maxdepth, minlinks)
            continue

        candidates = [x for x in alldepths if x[-1] == maxdepth]
        nseqids = len(set(x[0] for x in candidates))
        msg = "Multiple conflicting candidates found"
        if nseqids != 1:
            print >> log, msg
            continue

        seqid, mmin, mmax, depth = candidates[0]
        mmin, mmax = range_minmax([x[1:3] for x in candidates])
        if (mmax - mmin) > maxdist:
            print >> log, msg
            continue

        # Determine orientation by voting
        nplus, nminus = 0, 0
        arange = (seqid, mmin, mmax)
        for sid, start, end, sf, sc, fbstrand in ranges:
            brange = (sid, start, end)
            if range_overlap(arange, brange):
                if fbstrand == '+':
                    nplus += 1
                else:
                    nminus += 1

        fbstrand = '+' if nplus >= nminus else '-'

        candidate = (seqid, mmin, mmax, scf, depth, fbstrand)
        bedline = BedLine("\t".join((str(x) for x in candidate)))
        cbeds.append(bedline)
        print >> log, "Plus: {0}, Minus: {1}".format(nplus, nminus)
        print >> log, candidate

    candidatebed.extend(cbeds)
    logging.debug("A total of {0} scaffolds can be placed.".\
                    format(len(candidatebed)))
    log.close()

    candidatebedfile = pf + ".candidate.bed"
    candidatebed.print_to_file(candidatebedfile, sorted=True)


def gaps(args):
    """
    %prog gaps OM.bed fastafile

    Create patches around OM gaps.
    """
    from jcvi.formats.bed import uniq
    from jcvi.utils.iter import pairwise

    p = OptionParser(gaps.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ombed, fastafile = args
    ombed = uniq([ombed])
    bed = Bed(ombed)

    for a, b in pairwise(bed):
        om_a = (a.seqid, a.start, a.end, "+")
        om_b = (b.seqid, b.start, b.end, "+")
        ch_a = range_parse(a.accn)
        ch_b = range_parse(b.accn)
        ch_a = (ch_a.seqid, ch_a.start, ch_a.end, "+")
        ch_b = (ch_b.seqid, ch_b.start, ch_b.end, "+")

        om_dist, x = range_distance(om_a, om_b, distmode="ee")
        ch_dist, x = range_distance(ch_a, ch_b, distmode="ee")

        if om_dist <= 0 and ch_dist <= 0:
            continue

        print a
        print b
        print om_dist, ch_dist


def tips(args):
    """
    %prog tips patchers.bed complements.bed original.fasta backbone.fasta

    Append telomeric sequences based on patchers and complements.
    """
    p = OptionParser(tips.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    pbedfile, cbedfile, sizesfile, bbfasta = args

    pbed = Bed(pbedfile, sorted=False)
    cbed = Bed(cbedfile, sorted=False)

    complements = dict()
    for object, beds in groupby(cbed, key=lambda x: x.seqid):
        beds = list(beds)
        complements[object] = beds

    sizes = Sizes(sizesfile).mapping
    bbsizes = Sizes(bbfasta).mapping
    tbeds = []

    for object, beds in groupby(pbed, key=lambda x: x.accn):
        beds = list(beds)
        startbed, endbed = beds[0], beds[-1]
        start_id, end_id = startbed.seqid, endbed.seqid
        if startbed.start == 1:
            start_id = None
        if endbed.end == sizes[end_id]:
            end_id = None
        print >> sys.stderr, object, start_id, end_id
        if start_id:
            b = complements[start_id][0]
            b.accn = object
            tbeds.append(b)
        tbeds.append(BedLine("\t".join(str(x) for x in \
                        (object, 0, bbsizes[object], object, 1000, "+"))))
        if end_id:
            b = complements[end_id][-1]
            b.accn = object
            tbeds.append(b)

    tbed = Bed()
    tbed.extend(tbeds)

    tbedfile = "tips.bed"
    tbed.print_to_file(tbedfile)


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

    gapsbed = mergeBed(gapsbed, d=2 * Ext, nms=True)

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

    The output is a bedfile that can be converted to AGP using
    jcvi.formats.agp.frombed().
    """
    from jcvi.apps.base import blast
    from jcvi.formats.blast import BlastSlow
    from jcvi.formats.fasta import SeqIO

    p = OptionParser(install.__doc__)
    p.add_option("--rclip", default=1, type="int",
            help="Pair ID is derived from rstrip N chars [default: %default]")
    p.add_option("--maxsize", default=1000000, type="int",
            help="Maximum size of patchers to be replaced [default: %default]")
    p.add_option("--prefix", help="Prefix of the new object [default: %default]")
    p.add_option("--strict", default=False, action="store_true",
            help="Only update if replacement has no gaps [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    pbed, pfasta, bbfasta, altfasta = args
    Max = opts.maxsize  # Max DNA size to replace gap
    rclip = opts.rclip
    prefix = opts.prefix

    blastfile = blast([altfasta, pfasta,"--wordsize=100", "--pctid=99"])
    order = Bed(pbed).order

    beforebed, afterbed = "before.bed", "after.bed"
    fwa = open(beforebed, "w")
    fwb = open(afterbed, "w")

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
        print >> fwa, "\t".join(str(x) for x in \
                    (ax.seqid, qstart - 1, qstop, name, 1000, "+"))
        print >> fwb, "\t".join(str(x) for x in \
                    (asubject, sstart - 1, sstop, name, 1000, astrand))

    fwa.close()
    fwb.close()

    beforefasta = fastaFromBed(beforebed, bbfasta, name=True, stranded=True)
    afterfasta = fastaFromBed(afterbed, altfasta, name=True, stranded=True)

    # Exclude the replacements that contain more Ns than before
    ah = SeqIO.parse(beforefasta, "fasta")
    bh = SeqIO.parse(afterfasta, "fasta")
    count_Ns = lambda x: x.seq.count('n') + x.seq.count('N')
    exclude = set()
    for arec, brec in zip(ah, bh):
        an = count_Ns(arec)
        bn = count_Ns(brec)
        if opts.strict:
            if bn == 0:
                continue

        elif bn < an:
            continue

        id = arec.id
        exclude.add(id)

    logging.debug("Ignore {0} updates because of decreasing quality."\
                    .format(len(exclude)))

    abed = Bed(beforebed, sorted=False)
    bbed = Bed(afterbed, sorted=False)
    abed = [x for x in abed if x.accn not in exclude]
    bbed = [x for x in bbed if x.accn not in exclude]

    abedfile = "before.filtered.bed"
    bbedfile = "after.filtered.bed"
    afbed = Bed()
    afbed.extend(abed)
    bfbed = Bed()
    bfbed.extend(bbed)

    afbed.print_to_file(abedfile)
    bfbed.print_to_file(bbedfile)

    # Shuffle the two bedfiles together
    sizesfile = Sizes(bbfasta).filename
    shuffled = "shuffled.bed"

    abedfile = complementBed(abedfile, sizesfile)
    afbed = Bed(abedfile)
    all = []
    j = 0  # index in bfbed
    cj = 0  # index of new molecule

    for chr, beds in groupby(afbed, key=lambda x: x.seqid):
        beds = list(beds)
        cj += 1
        if prefix:
            accn = prefix + str(cj)

        b = beds[0]
        if prefix:
            b.accn = accn
        all.append(b)
        for b in beds[1:]:
            nb = bfbed[j]
            if prefix:
                b.accn = nb.accn = accn
            all.append(nb)
            j += 1
            all.append(b)

    shuffledbed = Bed()
    shuffledbed.extend(all)
    shuffledbed.print_to_file(shuffled)


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
    ncols = len(open(breakpointsbed).next().split())
    logging.debug("File {0} contains {1} columns.".format(breakpointsbed, ncols))
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
    for b, gaps in groupby(data, key=lambda x: x[:ncols]):
        gaps = list(gaps)
        gap = gaps[0]
        if len(gaps) == 1 and gap[-1] == "0":
            assert gap[-2] == "."
            print >> nogapsfw, "\t".join(b)
            continue

        gaps = [(int(x[-1]), x) for x in gaps]
        maxgap = max(gaps)[1]
        print >> largestgapsfw, "\t".join(maxgap)

    nogapsfw.close()
    largestgapsfw.close()

    closestgapsbed = pf + ".closestgaps.bed"
    closestgapsfw = open(closestgapsbed, "w")
    cmd = "closestBed -a {0} -b {1} -d".format(nogapsbed, gapsbed)
    sh(cmd, outfile=closestgapsbed)

    refinedbed = pf + ".refined.bed"
    FileMerger([largestgapsbed, closestgapsbed], outfile=refinedbed).merge()

    # Clean-up
    toclean = [nogapsbed, largestgapsbed, closestgapsbed]
    FileShredder(toclean)

    return refinedbed


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
    p.add_option("--object", default="object",
                 help="New object name [default: %default]")
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
                (chr, start, end, opts.object, 1000, strand))

    bed_fw.close()


if __name__ == '__main__':
    main()
