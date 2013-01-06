#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert common output files from gene prediction softwares into gff3 format.

Similar to the utilities in DAWGPAWS.
<http://dawgpaws.sourceforge.net/man.html>
"""

import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.gff import GffLine, Gff
from jcvi.utils.cbook import number
from jcvi.apps.base import ActionDispatcher, debug, need_update
debug()


FRAME, RETAIN, OVERLAP, NEW = "FRAME", "RETAIN", "OVERLAP", "NEW"
PRIORITY = (FRAME, RETAIN, OVERLAP, NEW)


def main():

    actions = (
        ('rename', 'rename genes for annotation release'),
        # Medicago gene renumbering
        ('renumber', 'renumber genes for annotation updates'),
        ('instantiate', 'instantiate NEW genes tagged by renumber'),
        # External gene prediction programs
        ('augustus', 'convert augustus output into gff3'),
        ('tRNAscan', 'convert tRNAscan-SE output into gff3'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def instantiate(args):
    """
    %prog instantiate tagged.bed blacklist.ids

    instantiate NEW genes tagged by renumber.
    """
    from jcvi.formats.base import SetFile

    p = OptionParser(instantiate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    taggedbed, blacklist = args
    black = SetFile(blacklist)
    black = set(atg_name(x) for x in black)

    # Run through the bed, identify stretch of NEW ids to instantiate,
    # identify the flanking FRAMEs, interpolate!
    bed = Bed(taggedbed)
    errorslog = open("errors.log", "w")
    tagkey = lambda x: x.rsplit("|", 1)[-1]
    for chr, sbed in bed.sub_beds():
        if "chr" not in chr:
            continue

        current_chr = number(chr)

        ranks = []
        for i, s in enumerate(sbed):
            nametag = s.extra[2]
            tag = tagkey(nametag)

            if tag in (NEW, FRAME):
                ranks.append((i, nametag))

        blocks = []
        for tag, names in groupby(ranks, key=lambda x: tagkey(x[-1])):
            names = list(names)
            if tag == NEW:
                blocks.append((tag, [x[0] for x in names]))
            else:
                start, end = names[0][-1], names[-1][-1]
                start, end = atg_name(start)[-1], atg_name(end)[-1]
                blocks.append((tag, [start, end]))

        for i, (tag, info) in enumerate(blocks):
            if tag != NEW:
                continue

            needed = len(info)
            start_id = 0 if i == 0 else blocks[i - 1][1][-1]
            end_id = start_id + 10000 if i == len(blocks) -1 \
                        else blocks[i + 1][1][0]

            assert end_id > start_id, \
                "end ({0}) > start ({1})".format(end_id, start_id)

            spots = end_id - start_id - 1
            available = len([x for x in xrange(start_id + 1, end_id) if
                                (current_chr, x) not in black])
            message = "chr{0} need {1} ids, has {2} spots ({3} available)".\
                    format(current_chr, needed, spots, available)

            start_gene = gene_name(current_chr, start_id)
            end_gene = gene_name(current_chr, end_id)
            message += " between {0} - {1}".format(start_gene, end_gene)

            if available < needed:
                print >> errorslog, message
            else:
                print message

    errorslog.close()


def atg_name(name):

    name = name.upper().rsplit(".", 1)[0]
    if "G" in name:
        first, second = name.rsplit("G", 1)
    elif "TE" in name:
        first, second = name.rsplit("TE", 1)
    else:
        return None, None

    chr = number(first)
    rank = number(second)

    return chr, rank


def gene_name(current_chr, x, prefix="Medtr", sep="g", pad0=6):
    return "{0}{1}{2}{3:0{4}}".format(prefix, current_chr, sep, x, pad0)


def prepare(bedfile):
    """
    Remove prepended tags in gene names.
    """
    pf = bedfile.rsplit(".", 1)[0]
    abedfile = pf + ".a.bed"
    bbedfile = pf + ".b.bed"
    fwa = open(abedfile, "w")
    fwb = open(bbedfile, "w")

    bed = Bed(bedfile)
    seen = set()
    for b in bed:
        accns = b.accn.split(";")
        new_accns = []
        for accn in accns:
            if ":" in accn:
                method, a = accn.split(":", 1)
                if method in ("liftOver", "GMAP", ""):
                    accn = a
            if accn in seen:
                logging.error("Duplicate id {0} found. Ignored.".format(accn))
                continue

            new_accns.append(accn)
            b.accn = accn
            print >> fwa, b
            seen.add(accn)

        b.accn = ";".join(new_accns)
        print >> fwb, b
    fwa.close()
    fwb.close()


def renumber(args):
    """
    %prog renumber Mt35.consolidated.bed > tagged.bed

    Renumber genes for annotation updates.
    """
    from jcvi.algorithms.lis import longest_increasing_subsequence
    from jcvi.formats.bed import mergeBed
    from jcvi.utils.grouper import Grouper

    p = OptionParser(renumber.__doc__)
    p.add_option("--pad0", default=6, type="int",
                 help="Pad gene identifiers with 0 [default: %default]")
    p.add_option("--prefix", default="Medtr",
                 help="Genome prefix [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    pad0 = opts.pad0
    prefix = opts.prefix

    pf = bedfile.rsplit(".", 1)[0]
    abedfile = pf + ".a.bed"
    bbedfile = pf + ".b.bed"
    if need_update(bedfile, (abedfile, bbedfile)):
        prepare(bedfile)

    mbed = Bed(bbedfile)
    g = Grouper()
    for s in mbed:
        accn = s.accn
        g.join(*accn.split(";"))

    bed = Bed(abedfile)
    for chr, sbed in bed.sub_beds():
        if "chr" not in chr:
            continue

        current_chr = number(chr)
        ranks = []

        gg = set()
        for s in sbed:
            accn = s.accn
            achr, arank = atg_name(accn)
            if achr != current_chr:
                continue
            ranks.append(arank)
            gg.add(accn)

        lranks = longest_increasing_subsequence(ranks)
        print >> sys.stderr, current_chr, len(sbed), "==>", len(ranks), \
                    "==>", len(lranks)

        granks = set(gene_name(current_chr, x) for x in lranks) | \
                 set(gene_name(current_chr, x, sep="te") for x in lranks)

        tagstore = {}
        for s in sbed:
            achr, arank = atg_name(s.accn)
            accn = s.accn
            if accn in granks:
                tag = (accn, FRAME)
            elif accn in gg:
                tag = (accn, RETAIN)
            else:
                tag = (".", NEW)

            tagstore[accn] = tag

        # Find cases where genes overlap
        for s in sbed:
            accn = s.accn
            gaccn = g[accn]
            tags = [((tagstore[x][-1] if x in tagstore else NEW), x) for x in gaccn]
            group = [(PRIORITY.index(tag), x) for tag, x in tags]
            best = min(group)[-1]

            if accn != best:
                tag = (best, OVERLAP)
            else:
                tag = tagstore[accn]

            print "\t".join((str(s), "|".join(tag)))


def rename(args):
    """
    %prog rename genes.bed gaps.bed

    Rename genes for annotation release.

    For genes on chromosomes (e.g. the 12th gene on C1):
    Bo1g00120

    For genes on scaffolds (e.g. the 12th gene on unplaced Scaffold00285):
    Bo00285s120

    The genes identifiers will increment by 10. So assuming no gap, these are
    the consecutive genes:
    Bo1g00120, Bo1g00130, Bo1g00140...
    Bo00285s120, Bo00285s130, Bo00285s140...

    When we encounter gaps, we would like the increment to be larger. For example,
    Bo1g00120, <gap>, Bo1g01120...
    """
    import string

    p = OptionParser(rename.__doc__)
    p.add_option("-a", dest="gene_increment", default=10, type="int",
                 help="Increment for continuous genes [default: %default]")
    p.add_option("-b", dest="gap_increment", default=1000, type="int",
                 help="Increment for gaps [default: %default]")
    p.add_option("--pad0", default=6, type="int",
                 help="Pad gene identifiers with 0 [default: %default]")
    p.add_option("--prefix", default="Bo",
                 help="Genome prefix [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    genebed, gapbed = args
    prefix = opts.prefix
    gene_increment = opts.gene_increment
    gap_increment = opts.gap_increment

    genes = Bed(genebed)
    fp = open(gapbed)
    for row in fp:
        genes.append(BedLine(row))

    genes.sort(key=genes.key)
    idsfile = prefix + ".ids"
    newbedfile = prefix + ".bed"
    gap_increment -= gene_increment
    assert gap_increment >= 0

    fw = open(idsfile, "w")
    for chr, lines in groupby(genes, key=lambda x: x.seqid):
        lines = list(lines)
        pad0 = opts.pad0 if len(lines) > 1000 else 3
        isChr = chr[0].upper() == 'C'
        digits = "".join(x for x in chr if x in string.digits)
        gs = "g" if isChr else "s"
        pp = prefix + digits + gs
        idx = 0
        if isChr:
            idx += gap_increment

        for r in lines:
            isGap = r.strand not in ("+", "-")
            if isGap:
                idx += gap_increment
                continue
            else:
                idx += gene_increment
            accn = pp + "{0:0{1}d}".format(idx, pad0)
            oldaccn = r.accn
            print >> fw, "\t".join((oldaccn, accn))
            r.accn = accn

    genes.print_to_file(newbedfile)
    logging.debug("Converted IDs written to `{0}`.".format(idsfile))
    logging.debug("Converted bed written to `{0}`.".format(newbedfile))


def augustus(args):
    """
    %prog augustus augustus.gff3 > reformatted.gff3

    AUGUSTUS does generate a gff3 (--gff3=on) but need some refinement.
    """
    p = OptionParser(augustus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    ingff3, = args
    gff = Gff(ingff3)
    for g in gff:
        if g.type not in ("gene", "transcript", "CDS"):
            continue

        if g.type == "transcript":
            g.type = "mRNA"

        print g


def tRNAscan(args):
    """
    %prog tRNAscan all.trna > all.trna.gff3

    Convert tRNAscan-SE output into gff3 format.

    Sequence                tRNA            Bounds          tRNA    Anti    Intron Bounds   Cove
    Name            tRNA #  Begin           End             Type    Codon   Begin   End     Score
    --------        ------  ----            ------          ----    -----   -----   ----    ------
    23231           1       335355          335440          Tyr     GTA     335392  335404  69.21
    23231           2       1076190         1076270         Leu     AAG     0       0       66.33

    Conversion based on PERL one-liner in:
    <https://github.com/sujaikumar/assemblage/blob/master/README-annotation.md>
    """
    from jcvi.formats.gff import sort

    p = OptionParser(tRNAscan.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    trnaout, = args
    gffout = trnaout + ".gff3"
    fp = open(trnaout)
    fw = open(gffout, "w")

    fp.next()
    fp.next()
    row = fp.next()
    assert row.startswith("--------")

    for row in fp:
        atoms = [x.strip() for x in row.split("\t")]
        contig, trnanum, start, end, aa, codon, \
                intron_start, intron_end, score = atoms

        start, end = int(start), int(end)
        orientation = '+'
        if start > end:
            start, end = end, start
            orientation = '-'

        source = "tRNAscan"
        type = "tRNA"
        if codon == "???":
            codon = "XXX"

        comment = "ID={0}.tRNA.{1};Name=tRNA-{2} (anticodon: {3})".\
                format(contig, trnanum, aa, codon)

        print >> fw, "\t".join(str(x) for x in (contig, source, type, start,\
                            end, score, orientation, ".", comment))

    fw.close()
    sort([gffout, "-i"])


if __name__ == '__main__':
    main()
