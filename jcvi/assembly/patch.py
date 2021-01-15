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
import sys
import math
import logging

from collections import defaultdict
from itertools import groupby
from more_itertools import pairwise, roundrobin

from jcvi.formats.bed import Bed, complementBed, mergeBed, fastaFromBed, summary
from jcvi.formats.blast import BlastSlow
from jcvi.formats.sizes import Sizes
from jcvi.utils.range import (
    range_parse,
    range_distance,
    ranges_depth,
    range_minmax,
    range_overlap,
    range_merge,
    range_closest,
    range_interleave,
)
from jcvi.formats.base import FileMerger, FileShredder
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (
        # OM guided approach
        ("refine", "find gaps within or near breakpoint regions"),
        ("patcher", "given om alignment, prepare the patchers"),
        # Gap filling through sequence matching
        ("fill", "perform gap filling using one assembly vs the other"),
        ("install", "install patches into backbone"),
        # Placement through mates and manual insertions and deletions
        ("bambus", "find candidate scaffolds to insert based on mates"),
        ("insert", "insert scaffolds into assembly"),
        ("eject", "eject scaffolds from assembly"),
        ("closest", "find the nearest gaps flanking suggested regions"),
        # Misc
        ("tips", "append telomeric sequences based on patchers and complements"),
        ("gaps", "create patches around OM gaps"),
        # Touch-up
        ("pasteprepare", "prepare sequences for paste"),
        ("paste", "paste in good sequences in the final assembly"),
        ("pastegenes", "paste in zero or low coverage genes"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pastegenes(args):
    """
    %prog pastegenes coverage.list old.genes.bed new.genes.bed old.assembly

    Paste in zero or low coverage genes.  For a set of neighboring genes
    missing, add the whole cassette as unplaced scaffolds. For singletons the
    program will try to make a patch.
    """
    from jcvi.formats.base import DictFile
    from jcvi.utils.cbook import gene_name

    p = OptionParser(pastegenes.__doc__)
    p.add_option(
        "--cutoff",
        default=90,
        type="int",
        help="Coverage cutoff to call gene missing",
    )
    p.add_option(
        "--flank",
        default=2000,
        type="int",
        help="Get the seq of size on two ends",
    )
    p.add_option(
        "--maxsize",
        default=50000,
        type="int",
        help="Maximum size of patchers to be replaced",
    )
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    coveragefile, oldbed, newbed, oldassembly = args
    cutoff = opts.cutoff
    flank = opts.flank
    maxsize = opts.maxsize

    coverage = DictFile(coveragefile, valuepos=2, cast=float)

    obed = Bed(oldbed)
    order = obed.order
    bed = [x for x in obed if x.accn in coverage]
    key = lambda x: coverage[x.accn] >= cutoff

    extrabed = "extra.bed"
    extendbed = "extend.bed"
    pastebed = "paste.bed"

    fw = open(extrabed, "w")
    fwe = open(extendbed, "w")
    fwp = open(pastebed, "w")
    fw_ids = open(extendbed + ".ids", "w")

    singletons, large, large_genes = 0, 0, 0
    for chr, chrbed in groupby(bed, key=lambda x: x.seqid):
        chrbed = list(chrbed)
        for good, beds in groupby(chrbed, key=key):
            if good:
                continue

            beds = list(beds)
            blocksize = len(set([gene_name(x.accn) for x in beds]))
            if blocksize == 1:
                singletons += 1
                accn = beds[0].accn
                gi, gb = order[accn]
                leftb = obed[gi - 1]
                rightb = obed[gi + 1]
                leftr = leftb.range
                rightr = rightb.range
                cur = gb.range
                distance_to_left, oo = range_distance(leftr, cur)
                distance_to_right, oo = range_distance(cur, rightr)
                span, oo = range_distance(leftr, rightr)

                if distance_to_left <= distance_to_right and distance_to_left > 0:
                    label = "LEFT"
                else:
                    label = "RIGHT"

                if 0 < span <= maxsize:
                    print(
                        "\t".join(
                            str(x) for x in (chr, leftb.start, rightb.end, gb.accn)
                        ),
                        file=fwp,
                    )

                print(leftb, file=fwe)
                print(gb, file=fwe)
                print(rightb, file=fwe)
                print(
                    "L:{0} R:{1} [{2}]".format(
                        distance_to_left, distance_to_right, label
                    ),
                    file=fwe,
                )
                print(gb.accn, file=fw_ids)
                continue

            large += 1
            large_genes += blocksize

            ranges = [(x.start, x.end) for x in beds]
            rmin, rmax = range_minmax(ranges)
            rmin -= flank
            rmax += flank

            name = "-".join((beds[0].accn, beds[-1].accn))
            print("\t".join(str(x) for x in (chr, rmin - 1, rmax, name)), file=fw)

    fw.close()
    fwe.close()

    extrabed = mergeBed(extrabed, d=flank, nms=True)
    fastaFromBed(extrabed, oldassembly, name=True)
    summary([extrabed])

    logging.debug("Singleton blocks : {0}".format(singletons))
    logging.debug("Large blocks : {0} ({1} genes)".format(large, large_genes))


def pasteprepare(args):
    """
    %prog pasteprepare bacs.fasta

    Prepare sequences for paste.
    """
    p = OptionParser(pasteprepare.__doc__)
    p.add_option(
        "--flank",
        default=5000,
        type="int",
        help="Get the seq of size on two ends",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (goodfasta,) = args
    flank = opts.flank
    pf = goodfasta.rsplit(".", 1)[0]
    extbed = pf + ".ext.bed"

    sizes = Sizes(goodfasta)
    fw = open(extbed, "w")
    for bac, size in sizes.iter_sizes():
        print("\t".join(str(x) for x in (bac, 0, min(flank, size), bac + "L")), file=fw)
        print(
            "\t".join(str(x) for x in (bac, max(size - flank, 0), size, bac + "R")),
            file=fw,
        )
    fw.close()

    fastaFromBed(extbed, goodfasta, name=True)


def paste(args):
    """
    %prog paste flanks.bed flanks_vs_assembly.blast backbone.fasta

    Paste in good sequences in the final assembly.
    """
    from jcvi.formats.bed import uniq

    p = OptionParser(paste.__doc__)
    p.add_option(
        "--maxsize",
        default=300000,
        type="int",
        help="Maximum size of patchers to be replaced",
    )
    p.add_option("--prefix", help="Prefix of the new object")
    p.set_rclip(rclip=1)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    pbed, blastfile, bbfasta = args
    maxsize = opts.maxsize  # Max DNA size to replace gap
    order = Bed(pbed).order

    beforebed, afterbed = blast_to_twobeds(
        blastfile, order, log=True, rclip=opts.rclip, maxsize=maxsize, flipbeds=True
    )
    beforebed = uniq([beforebed])

    afbed = Bed(beforebed)
    bfbed = Bed(afterbed)

    shuffle_twobeds(afbed, bfbed, bbfasta, prefix=opts.prefix)


def eject(args):
    """
    %prog eject candidates.bed chr.fasta

    Eject scaffolds from assembly, using the range identified by closest().
    """
    p = OptionParser(eject.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    candidates, chrfasta = args
    sizesfile = Sizes(chrfasta).filename
    cbedfile = complementBed(candidates, sizesfile)

    cbed = Bed(cbedfile)
    for b in cbed:
        b.accn = b.seqid
        b.score = 1000
        b.strand = "+"

    cbed.print_to_file()


def closest(args):
    """
    %prog closest candidates.bed gaps.bed fastafile

    Identify the nearest gaps flanking suggested regions.
    """
    p = OptionParser(closest.__doc__)
    p.add_option(
        "--om",
        default=False,
        action="store_true",
        help="The bedfile is OM blocks",
    )
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    candidates, gapsbed, fastafile = args
    sizes = Sizes(fastafile).mapping
    bed = Bed(candidates)
    ranges = []
    for b in bed:
        r = range_parse(b.accn) if opts.om else b
        ranges.append([r.seqid, r.start, r.end])

    gapsbed = Bed(gapsbed)
    granges = [(x.seqid, x.start, x.end) for x in gapsbed]

    ranges = range_merge(ranges)
    for r in ranges:
        a = range_closest(granges, r)
        b = range_closest(granges, r, left=False)
        seqid = r[0]

        if a is not None and a[0] != seqid:
            a = None
        if b is not None and b[0] != seqid:
            b = None

        mmin = 1 if a is None else a[1]
        mmax = sizes[seqid] if b is None else b[2]

        print("\t".join(str(x) for x in (seqid, mmin - 1, mmax)))


def insert(args):
    """
    %prog insert candidates.bed gaps.bed chrs.fasta unplaced.fasta

    Insert scaffolds into assembly.
    """
    from jcvi.formats.agp import mask, bed
    from jcvi.formats.sizes import agp

    p = OptionParser(insert.__doc__)
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
        if len(atoms) <= 6:
            continue
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
    finalbed = Bed()
    for b in mbed:
        sid = b.seqid
        key = (sid, b.start, b.end)
        if key not in gappositions:
            finalbed.add("{0}\n".format(b))
            continue

        gapid = gappositions[key]
        scfs = gap_to_scf[gapid]

        # For scaffolds placed in the same gap, sort according to positions
        scfs.sort(key=lambda x: corder[x[0]][1].start + corder[x[0]][1].end)
        for scf, strand in scfs:
            size = sizes[scf]
            finalbed.add("\t".join(str(x) for x in (scf, 0, size, sid, 1000, strand)))

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
    from jcvi.formats.bed import BedLine
    from jcvi.formats.posmap import MatesFile

    p = OptionParser(bambus.__doc__)
    p.add_option(
        "--prefix",
        default="scaffold",
        help="Prefix of the unplaced scaffolds",
    )
    p.add_option(
        "--minlinks",
        default=3,
        type="int",
        help="Minimum number of links to place",
    )
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

        if aname not in mf:
            continue

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
        print(file=log)
        ranges = []
        for a, b in beds:
            aname, astrand = a.accn, a.strand
            bname, bstrand = b.accn, b.strand
            aseqid, bseqid = a.seqid, b.seqid
            pa, lib = mf[aname]

            print(a, file=log)
            print(b, file=log)

            flip_b = astrand == bstrand
            fbstrand = "-" if flip_b else "+"
            if flip_b:
                b.reverse_complement(sizes)

            lmin, lmax = lib.min, lib.max

            L = sizes.get_size(scf)
            assert astrand in ("+", "-")
            if astrand == "+":
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
            print("*" + "\t".join(str(x) for x in start_range), file=log)
            ranges.append(start_range)

        mranges = [x[:3] for x in ranges]
        # Determine placement by finding the interval with the most support
        rd = ranges_depth(mranges, sizes.mapping, verbose=False)
        alldepths = []
        for depth in rd:
            alldepths.extend(depth)
        print(alldepths, file=log)

        maxdepth = max(alldepths, key=lambda x: x[-1])[-1]
        if maxdepth < minlinks:
            print("Insufficient links ({0} < {1})".format(maxdepth, minlinks), file=log)
            continue

        candidates = [x for x in alldepths if x[-1] == maxdepth]
        nseqids = len(set(x[0] for x in candidates))
        if nseqids != 1:
            msg = "Multiple conflicting candidates found"
            print(msg, file=log)
            continue

        seqid, mmin, mmax, depth = candidates[0]
        mmin, mmax = range_minmax([x[1:3] for x in candidates])

        if mmin >= mmax:
            msg = "Invalid (min, max) range"
            print("Invalid (min, max) range", file=log)
            continue

        if (mmax - mmin) > maxdist:
            msg = "(min, max) distance greater than library maxdist"
            print(msg, file=log)
            continue

        # Determine orientation by voting
        nplus, nminus = 0, 0
        arange = (seqid, mmin, mmax)
        for sid, start, end, sf, sc, fbstrand in ranges:
            brange = (sid, start, end)
            if range_overlap(arange, brange):
                if fbstrand == "+":
                    nplus += 1
                else:
                    nminus += 1

        fbstrand = "+" if nplus >= nminus else "-"

        candidate = (seqid, mmin, mmax, scf, depth, fbstrand)
        bedline = BedLine("\t".join((str(x) for x in candidate)))
        cbeds.append(bedline)
        print("Plus: {0}, Minus: {1}".format(nplus, nminus), file=log)
        print(candidate, file=log)

    candidatebed.extend(cbeds)
    logging.debug("A total of {0} scaffolds can be placed.".format(len(candidatebed)))
    log.close()

    candidatebedfile = pf + ".candidate.bed"
    candidatebed.print_to_file(candidatebedfile, sorted=True)


def gaps(args):
    """
    %prog gaps OM.bed fastafile

    Create patches around OM gaps.
    """
    from jcvi.formats.bed import uniq

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

        print(a)
        print(b)
        print(om_dist, ch_dist)


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
        print(object, start_id, end_id, file=sys.stderr)
        if start_id:
            b = complements[start_id][0]
            b.accn = object
            tbeds.append(b)
        tbeds.append(
            BedLine(
                "\t".join(
                    str(x) for x in (object, 0, bbsizes[object], object, 1000, "+")
                )
            )
        )
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
    p.add_option(
        "--extend",
        default=2000,
        type="int",
        help="Extend seq flanking the gaps",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gapsbed, badfasta = args
    Ext = opts.extend

    gapdist = 2 * Ext + 1  # This is to prevent to replacement ranges intersect
    gapsbed = mergeBed(gapsbed, d=gapdist, nms=True)

    bed = Bed(gapsbed)
    sizes = Sizes(badfasta).mapping
    pf = gapsbed.rsplit(".", 1)[0]
    extbed = pf + ".ext.bed"
    fw = open(extbed, "w")
    for b in bed:
        gapname = b.accn
        start, end = max(0, b.start - Ext - 1), b.start - 1
        print("\t".join(str(x) for x in (b.seqid, start, end, gapname + "L")), file=fw)
        start, end = b.end, min(sizes[b.seqid], b.end + Ext)
        print("\t".join(str(x) for x in (b.seqid, start, end, gapname + "R")), file=fw)
    fw.close()

    fastaFromBed(extbed, badfasta, name=True)


def blast_to_twobeds(
    blastfile, order, log=False, rclip=1, maxsize=300000, flipbeds=False
):

    abed, bbed = "before.bed", "after.bed"
    beforebed, afterbed = abed, bbed
    if flipbeds:
        beforebed, afterbed = afterbed, beforebed

    fwa = open(beforebed, "w")
    fwb = open(afterbed, "w")
    if log:
        logfile = "problems.log"
        log = open(logfile, "w")

    key1 = lambda x: x.query
    key2 = lambda x: x.query[:-rclip] if rclip else key1
    data = BlastSlow(blastfile)
    OK = "OK"

    seen = set()
    for pe, lines in groupby(data, key=key2):
        label = OK
        lines = list(lines)
        if len(lines) != 2:
            label = "Singleton"

        else:
            a, b = lines

            aquery, bquery = a.query, b.query
            asubject, bsubject = a.subject, b.subject
            if asubject != bsubject:
                label = "Different chr {0}|{1}".format(asubject, bsubject)

            else:
                astrand, bstrand = a.orientation, b.orientation
                assert aquery[-1] == "L" and bquery[-1] == "R", str((aquery, bquery))

                ai, ax = order[aquery]
                bi, bx = order[bquery]
                qstart, qstop = ax.start + a.qstart - 1, bx.start + b.qstop - 1

                if astrand == "+" and bstrand == "+":
                    sstart, sstop = a.sstart, b.sstop

                elif astrand == "-" and bstrand == "-":
                    sstart, sstop = b.sstart, a.sstop

                else:
                    label = "Strand {0}|{1}".format(astrand, bstrand)

                if sstart > sstop:
                    label = "Start beyond stop"

                if sstop > sstart + maxsize:
                    label = "Stop beyond start plus {0}".format(maxsize)

        aquery = lines[0].query
        bac_name = aquery[:-1]
        seen.add(bac_name)
        name = bac_name + "LR"

        if label != OK:
            if log:
                print("\t".join((name, label)), file=log)
            continue

        print(
            "\t".join(str(x) for x in (ax.seqid, qstart - 1, qstop, name, 1000, "+")),
            file=fwa,
        )
        print(
            "\t".join(
                str(x) for x in (asubject, sstart - 1, sstop, name, 1000, astrand)
            ),
            file=fwb,
        )

    # Missing
    if log:
        label = "Missing"
        for k in order.keys():
            k = k[:-1]
            if k not in seen:
                seen.add(k)
                k += "LR"
                print("\t".join((k, label)), file=log)
        log.close()

    fwa.close()
    fwb.close()

    return abed, bbed


def shuffle_twobeds(afbed, bfbed, bbfasta, prefix=None):
    # Shuffle the two bedfiles together
    sz = Sizes(bbfasta)
    sizes = sz.mapping
    shuffled = "shuffled.bed"
    border = bfbed.order

    all = []
    afbed.sort(key=afbed.nullkey)
    totalids = len(sizes)
    pad = int(math.log10(totalids)) + 1
    cj = 0
    seen = set()
    accn = lambda x: "{0}{1:0{2}d}".format(prefix, x, pad)

    for seqid, aa in afbed.sub_beds():
        cj += 1
        abeds, bbeds, beds = [], [], []
        size = sizes[seqid]
        ranges = [(x.seqid, x.start, x.end) for x in aa]
        cranges = range_interleave(ranges, sizes={seqid: size}, empty=True)
        for crange in cranges:
            if crange:
                seqid, start, end = crange
                bedline = "\t".join(str(x) for x in (seqid, start - 1, end))
                abeds.append(BedLine(bedline))
            else:
                abeds.append(None)

        for a in aa:
            gapid = a.accn
            bi, b = border[gapid]
            if a.strand == "-":
                b.extra[1] = b.strand = "-" if b.strand == "+" else "+"

            bbeds.append(b)

        n_abeds = len(abeds)
        n_bbeds = len(bbeds)
        assert n_abeds - n_bbeds == 1, "abeds: {0}, bbeds: {1}".format(n_abeds, n_bbeds)

        beds = [x for x in roundrobin(abeds, bbeds) if x]
        if prefix:
            for b in beds:
                b.accn = accn(cj)

        all.extend(beds)
        seen.add(seqid)

    # Singletons
    for seqid, size in sz.iter_sizes():
        if seqid in seen:
            continue

        bedline = "\t".join(str(x) for x in (seqid, 0, size, accn(cj)))
        b = BedLine(bedline)

        cj += 1
        if prefix:
            b.accn = accn(cj)

        all.append(b)

    shuffledbed = Bed()
    shuffledbed.extend(all)
    shuffledbed.print_to_file(shuffled)

    return shuffledbed


def install(args):
    """
    %prog install patchers.bed patchers.fasta backbone.fasta alt.fasta

    Install patches into backbone, using sequences from alternative assembly.
    The patches sequences are generated via jcvi.assembly.patch.fill().

    The output is a bedfile that can be converted to AGP using
    jcvi.formats.agp.frombed().
    """
    from jcvi.apps.align import blast
    from jcvi.formats.fasta import SeqIO

    p = OptionParser(install.__doc__)
    p.set_rclip(rclip=1)
    p.add_option(
        "--maxsize",
        default=300000,
        type="int",
        help="Maximum size of patchers to be replaced",
    )
    p.add_option("--prefix", help="Prefix of the new object")
    p.add_option(
        "--strict",
        default=False,
        action="store_true",
        help="Only update if replacement has no gaps",
    )
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    pbed, pfasta, bbfasta, altfasta = args
    maxsize = opts.maxsize  # Max DNA size to replace gap
    rclip = opts.rclip

    blastfile = blast([altfasta, pfasta, "--wordsize=100", "--pctid=99"])
    order = Bed(pbed).order
    beforebed, afterbed = blast_to_twobeds(
        blastfile, order, rclip=rclip, maxsize=maxsize
    )

    beforefasta = fastaFromBed(beforebed, bbfasta, name=True, stranded=True)
    afterfasta = fastaFromBed(afterbed, altfasta, name=True, stranded=True)

    # Exclude the replacements that contain more Ns than before
    ah = SeqIO.parse(beforefasta, "fasta")
    bh = SeqIO.parse(afterfasta, "fasta")
    count_Ns = lambda x: x.seq.count("n") + x.seq.count("N")
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

    logging.debug(
        "Ignore {0} updates because of decreasing quality.".format(len(exclude))
    )

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

    shuffle_twobeds(afbed, bfbed, bbfasta, prefix=opts.prefix)


def refine(args):
    """
    %prog refine breakpoints.bed gaps.bed

    Find gaps within or near breakpoint region.

    For breakpoint regions with no gaps, there are two options:
    - Break in the middle of the region
    - Break at the closest gap (--closest)
    """
    p = OptionParser(refine.__doc__)
    p.add_option(
        "--closest",
        default=False,
        action="store_true",
        help="In case of no gaps, use closest",
    )
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
            assert gap[-3] == "."
            print("\t".join(b), file=nogapsfw)
            continue

        gaps = [(int(x[-1]), x) for x in gaps]
        maxgap = max(gaps)[1]
        print("\t".join(maxgap), file=largestgapsfw)

    nogapsfw.close()
    largestgapsfw.close()
    beds = [largestgapsbed]
    toclean = [nogapsbed, largestgapsbed]

    if opts.closest:
        closestgapsbed = pf + ".closestgaps.bed"
        cmd = "closestBed -a {0} -b {1} -d".format(nogapsbed, gapsbed)
        sh(cmd, outfile=closestgapsbed)
        beds += [closestgapsbed]
        toclean += [closestgapsbed]
    else:
        pointbed = pf + ".point.bed"
        pbed = Bed()
        bed = Bed(nogapsbed)
        for b in bed:
            pos = (b.start + b.end) / 2
            b.start, b.end = pos, pos
            pbed.append(b)
        pbed.print_to_file(pointbed)
        beds += [pointbed]
        toclean += [pointbed]

    refinedbed = pf + ".refined.bed"
    FileMerger(beds, outfile=refinedbed).merge()

    # Clean-up
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

    neg_strands = sum(1 for x in beds if x.strand == "-")
    pos_strands = len(beds) - neg_strands
    strand = "-" if neg_strands > pos_strands else "+"

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
    p.add_option(
        "--backbone",
        default="OM",
        help="Prefix of the backbone assembly",
    )
    p.add_option("--object", default="object", help="New object name")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    backbonebed, otherbed = args
    backbonebed = uniq([backbonebed])
    otherbed = uniq([otherbed])

    pf = backbonebed.split(".")[0]
    key = lambda x: (x.seqid, x.start, x.end)

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

        print(
            "\t".join(str(x) for x in (chr, start, end, opts.object, 1000, strand)),
            file=bed_fw,
        )

    bed_fw.close()


if __name__ == "__main__":
    main()
