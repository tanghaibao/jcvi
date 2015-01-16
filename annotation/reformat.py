#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert common output files from gene prediction softwares into gff3 format.

Similar to the utilities in DAWGPAWS.
<http://dawgpaws.sourceforge.net/man.html>
"""

import os
import sys
import logging
import re

from itertools import groupby, product

from jcvi.formats.bed import Bed, BedLine, sort
from jcvi.formats.base import SetFile, must_open, get_number, flexible_cast
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, \
            need_update, popen, sh


FRAME, RETAIN, OVERLAP, NEW = "FRAME", "RETAIN", "OVERLAP", "NEW"
PRIORITY = (FRAME, RETAIN, OVERLAP, NEW)
new_id_pat = re.compile(r"^\d+\.[cemtx]+\S+")


class Stride (object):
    """
    Allows four basic strides:

        0 10
       0 5 10
      0 3 7 10
     0 2 5 8 10

    We have main parameters, # we need, # available go through all possible
    numbers excluding everything in black.
    """
    def __init__(self, needed, available):
        configurations = ("0", "05", "037", "0258")
        nneeded = len(needed)
        self.conf = None
        self.available = None
        for c in configurations:
            a = [x for x in available if str(x)[-1] in c]
            if len(a) >= nneeded:
                self.conf = c
                self.available = a
                break


class NameRegister (object):

    def __init__(self, prefix="Medtr", pad0="6", uc=False):
        self.black = set()
        self.gaps = []
        self.prefix = prefix
        self.pad0 = pad0
        self.uc = uc

    def get_blacklist(self, filename):
        black = SetFile(filename)
        black = set(atg_name(x) for x in black)
        self.black.update(black)

    def get_gaps(self, filename):
        self.gapfile = filename

    def allocate(self, info, chr, start_id, end_id, id_table):

        start_bp = info[0].start
        end_bp = info[-1].end

        current_chr = chr_number(chr)
        needed = info
        assert end_id > start_id, \
            "end ({0}) > start ({1})".format(end_id, start_id)

        spots = end_id - start_id - 1
        available = [x for x in xrange(start_id + 1, end_id) if
                            (current_chr, x) not in self.black]

        message = "{0} need {1} ids, has {2} spots ({3} available)".\
                format(chr, len(needed), spots, len(available))

        start_gene = gene_name(current_chr, start_id, prefix=self.prefix, \
                pad0=self.pad0, uc=self.uc)
        end_gene = gene_name(current_chr, end_id, prefix=self.prefix,
                pad0=self.pad0, uc=self.uc)
        message += " between {0} - {1}\n".format(start_gene, end_gene)

        assert end_bp > start_bp

        b = "\t".join(str(x) for x in (chr, start_bp - 1, end_bp))
        cmd = "echo '{0}' |".format(b)
        cmd += " intersectBed -a {0} -b stdin".format(self.gapfile)
        gaps = list(BedLine(x) for x in popen(cmd, debug=False))
        ngaps = len(gaps)

        gapsexpanded = []
        GeneDensity = 10000.  # assume 10Kb per gene
        for gap in gaps:
            gap_bp = int(gap.score)
            gap_ids = int(round(gap_bp / GeneDensity))
            gapsexpanded += [gap] * gap_ids

        lines = sorted(info + gapsexpanded, key=lambda x: x.start)

        message += "between bp: {0} - {1}, there are {2} gaps (total {3} ids)".\
                format(start_bp, end_bp, ngaps, len(lines))

        needed = lines
        stride = Stride(needed, available)
        conf = stride.conf
        message += " stride: {0}".format(conf)
        print >> sys.stderr, message

        nneeded = len(needed)
        if conf is None: # prefix rule - prepend version number for spills
            magic = 400000  # version 4
            firstdigit = 100000
            step = 10  # stride for the prefixed ids
            rank = start_id + magic
            if rank > magic + firstdigit:
                rank -= firstdigit
            available = []
            while len(available) != nneeded:
                rank += step
                if (current_chr, rank) in self.black:  # avoid blacklisted ids
                    continue
                available.append(rank)

        else: # follow the best stride
            available = stride.available
            if start_id == 0:  # follow right flank at start of chr
                available = available[- nneeded:]
            else:  # follow left flank otherwise
                available = available[:nneeded]

        # Finally assign the ids
        assert len(needed) == len(available)
        for b, rank in zip(needed, available):
            name = gene_name(current_chr, rank, prefix=self.prefix, \
                    pad0=self.pad0, uc=self.uc)
            print >> sys.stderr, "\t".join((str(b), name))
            id_table[b.accn] = name
            self.black.add((current_chr, rank))
        print >> sys.stderr


def main():

    actions = (
        ('rename', 'rename genes for annotation release'),
        # perform following actions on list files
        ('reindex', 'reindex isoforms per gene locus'),
        ('publocus', 'create pub_locus identifiers according to GenBank specs'),
        # Medicago gene renumbering
        ('annotate', 'annotation new bed file with features from old'),
        ('renumber', 'renumber genes for annotation updates'),
        ('instantiate', 'instantiate NEW genes tagged by renumber'),
        ('plot', 'plot gene identifiers along certain chromosome'),
        # External gene prediction programs
        ('augustus', 'convert augustus output into gff3'),
        ('tRNAscan', 'convert tRNAscan-SE output into gff3'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def plot(args):
    """
    %prog plot tagged.new.bed chr1

    Plot gene identifiers along a particular chromosome, often to illustrate the
    gene id assignment procedure.
    """
    from jcvi.graphics.base import plt, savefig
    from jcvi.graphics.chromosome import ChromosomeMap

    p = OptionParser(plot.__doc__)
    p.add_option("--firstn", type="int", help="Only plot the first N genes")
    p.add_option("--ymax", type="int", help="Y-axis max value")
    p.add_option("--log", action="store_true",
                help="Write plotting data [default: %default]")
    opts, args, iopts = p.set_image_options(args, figsize="6x4")

    if len(args) != 2:
        sys.exit(not p.print_help())

    taggedbed, chr = args
    bed = Bed(taggedbed)
    beds = list(bed.sub_bed(chr))
    old, new = [], []
    i = 0
    for b in beds:
        accn = b.extra[0]
        if "te" in accn:
            continue

        accn, tag = accn.split("|")
        if tag == "OVERLAP":
            continue

        c, r = atg_name(accn)
        if tag == "NEW":
            new.append((i, r))
        else:
            old.append((i, r))
        i += 1

    ngenes = i
    assert ngenes == len(new) + len(old)

    logging.debug("Imported {0} ranks on {1}.".format(ngenes, chr))
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    xstart, xend = .2, .8
    ystart, yend = .2, .8
    pad = .02

    ngenes = opts.firstn or ngenes
    ymax = opts.ymax or 500000

    title = "Assignment of Medtr identifiers"
    if opts.ymax:
        subtitle = "{0}, first {1} genes".format(chr, ngenes)
    else:
        subtitle = "{0}, {1} genes ({2} new)".format(chr, ngenes, len(new))

    chr_map = ChromosomeMap(fig, root, xstart, xend, ystart, yend, pad, 0,
                        ymax, 5, title, subtitle)

    ax = chr_map.axes

    if opts.log:
        from jcvi.utils.table import write_csv
        header = ["x", "y"]
        write_csv(header, new, filename=chr + ".new")
        write_csv(header, old, filename=chr + ".old")

    x, y = zip(*new)
    ax.plot(x, y, "b,")
    x, y = zip(*old)
    ax.plot(x, y, "r,")

    # Legends
    ymid = (ystart + yend) / 2
    y = ymid + pad
    root.plot([.2], [y], "r.", lw=2)
    root.text(.2 + pad, y, "Existing Medtr ids", va="center", size=10)
    y = ymid - pad
    root.plot([.2], [y], "b.", lw=2)
    root.text(.2 + pad, y, "Newly instantiated ids", va="center", size=10)

    ax.set_xlim(0, ngenes)
    ax.set_ylim(0, ymax)
    ax.set_axis_off()

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + ".identifiers." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def instantiate(args):
    """
    %prog instantiate tagged.bed blacklist.ids big_gaps.bed

    instantiate NEW genes tagged by renumber.
    """
    p = OptionParser(instantiate.__doc__)
    p.set_annot_reformat_opts()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    taggedbed, blacklist, gapsbed = args
    r = NameRegister(prefix=opts.prefix, pad0=opts.pad0, uc=opts.uc)
    r.get_blacklist(blacklist)
    r.get_gaps(gapsbed)

    # Run through the bed, identify stretch of NEW ids to instantiate,
    # identify the flanking FRAMEs, interpolate!
    bed = Bed(taggedbed)
    outputbed = taggedbed.rsplit(".", 1)[0] + ".new.bed"
    fw = open(outputbed, "w")

    tagkey = lambda x: x.rsplit("|", 1)[-1]
    for chr, sbed in bed.sub_beds():
        current_chr = chr_number(chr)
        if not current_chr:
            continue

        sbed = list(sbed)

        ranks = []
        for i, s in enumerate(sbed):
            nametag = s.extra[0]
            tag = tagkey(nametag)

            if tag in (NEW, FRAME):
                ranks.append((i, nametag))

        blocks = []
        for tag, names in groupby(ranks, key=lambda x: tagkey(x[-1])):
            names = list(names)
            if tag == NEW:
                blocks.append((tag, [sbed[x[0]] for x in names]))
            else:
                start, end = names[0][-1], names[-1][-1]
                start, end = atg_name(start, retval="rank"), atg_name(end, retval="rank")
                blocks.append((tag, [start, end]))

        id_table = {}  # old to new name conversion
        for i, (tag, info) in enumerate(blocks):
            if tag != NEW:
                continue

            start_id = 0 if i == 0 else blocks[i - 1][1][-1]
            end_id = start_id + 10000 if i == len(blocks) -1 \
                        else blocks[i + 1][1][0]

            r.allocate(info, chr, start_id, end_id, id_table)

        # Output new names
        for i, s in enumerate(sbed):
            nametag = s.extra[0]
            name, tag = nametag.split("|")

            if tag == NEW:
                assert name == '.'
                name = id_table[s.accn]
            elif tag == OVERLAP:
                if name in id_table:
                    name = id_table[name]

            s.extra[0] = "|".join((name, tag))
            print >> fw, s

    fw.close()


def atg_name(name, retval="chr,rank", trimpad0=True):
    atg_name_pat = re.compile(r"""
            ^(?P<locus>
                (?:(?P<prefix>\D+[\D\d\D])\.?)(?P<chr>[\d|C|M]+)(?P<sep>[A-z]+)(?P<rank>\d+)
            )
            \.?(?P<iso>\d+)?
            """, re.VERBOSE)

    seps = ["g", "te", "trna", "s", "u"]
    pad0s = ["rank"]

    if name is not None:
        m = re.match(atg_name_pat, name)
        if m is not None and m.group('sep').lower() in seps:
            retvals = []
            for grp in retval.split(","):
                if grp == 'chr':
                    val = chr_number(m.group(grp))
                else:
                    val = get_number(m.group(grp)) \
                            if trimpad0 and grp in pad0s \
                            else m.group(grp)
                retvals.append(val)

            return (x for x in retvals) if len(retvals) > 1 \
                    else retvals[0]

    return (None for x in retval.split(","))


def gene_name(current_chr, x, prefix="Medtr", sep="g", pad0=6, uc=False):
    identifier = "{0}{1}{2}{3:0{4}}".format(prefix, current_chr, sep, x, pad0)
    if uc: identifier = identifier.upper()
    return identifier


def chr_number(chr):
    chr_pat = re.compile(r"(?P<prefix>\D*)(?P<chr>[\d|C|M]+)$", re.VERBOSE | re.IGNORECASE)

    if chr is not None:
        m = re.match(chr_pat, chr)
        if m is not None:
            return flexible_cast(m.group('chr'))

    return None


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
    from jcvi.utils.grouper import Grouper

    p = OptionParser(renumber.__doc__)
    p.set_annot_reformat_opts()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args

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
        current_chr = chr_number(chr)
        if not current_chr:
            continue

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

        granks = set(gene_name(current_chr, x, prefix=opts.prefix, \
                     pad0=opts.pad0, uc=opts.uc) for x in lranks) | \
                 set(gene_name(current_chr, x, prefix=opts.prefix, \
                     pad0=opts.pad0, sep="te", uc=opts.uc) for x in lranks)

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


def annotate(args):
    """
    %prog annotate new.bed old.bed 2> log

    Annotate the `new.bed` with features from `old.bed` for the purpose of
    gene numbering.

    Ambiguity in ID assignment can be resolved by either of the following 2 methods:
    - `alignment`: make use of global sequence alignment score (calculated by `needle`)
    - `overlap`: make use of overlap length (calculated by `intersectBed`)

    Transfer over as many identifiers as possible while following guidelines:
    http://www.arabidopsis.org/portals/nomenclature/guidelines.jsp#editing

    Note: Following RegExp pattern describes the structure of the identifier
    assigned to features in the `new.bed` file.

    new_id_pat = re.compile(r"^\d+\.[cemtx]+\S+")

    Examples: 23231.m312389, 23231.t004898, 23231.tRNA.144
    Adjust the value of `new_id_pat` manually as per your ID naming conventions.
    """
    from jcvi.utils.grouper import Grouper

    valid_resolve_choices = ["alignment", "overlap"]

    p = OptionParser(annotate.__doc__)
    p.add_option("--resolve", default="alignment", choices=valid_resolve_choices,
                 help="Resolve ID assignment based on a certain metric" \
                        + " [default: %default]")
    p.add_option("--atg_name", default=False, action="store_true",
                help="Specify is locus IDs in `new.bed` file follow ATG nomenclature" \
                        + " [default: %default]")

    g1 = OptionGroup(p, "Optional parameters (alignment):\n" \
            + "Use if resolving ambiguities based on sequence `alignment`")
    g1.add_option("--pid", dest="pid", default=35., type="float",
            help="Percent identity cutoff [default: %default]")
    g1.add_option("--score", dest="score", default=250., type="float",
            help="Alignment score cutoff [default: %default]")
    p.add_option_group(g1)

    g2 = OptionGroup(p, "Optional parameters (overlap):\n" \
            + "Use if resolving ambiguities based on `overlap` length\n" \
            + "Parameters equivalent to `intersectBed`")
    g2.add_option("-f", dest="f", default=0.5, type="float",
            help="Minimum overlap fraction (0.0 - 1.0) [default: %default]")
    g2.add_option("-r", dest="r", default=False, action="store_true",
            help="Require fraction overlap to be reciprocal [default: %default]")
    g2.add_option("-s", dest="s", default=True, action="store_true",
            help="Require same strandedness [default: %default]")
    p.add_option_group(g2)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    nbedfile, obedfile = args
    npf, opf = nbedfile.rsplit(".", 1)[0], obedfile.rsplit(".", 1)[0]

    # Make consolidated.bed
    cbedfile = "consolidated.bed"
    if not os.path.isfile(cbedfile):
        consolidate(nbedfile, obedfile, cbedfile)
    else:
        logging.warning("`{0}` already exists. Skipping step".format(cbedfile))

    logging.warning("Resolving ID assignment ambiguity based on `{0}`".\
            format(opts.resolve))

    if opts.resolve == "alignment":
        # Get pairs and prompt to run needle
        pairsfile = "nw.pairs"
        scoresfile = "nw.scores"
        if not os.path.isfile(pairsfile):
            get_pairs(cbedfile, pairsfile)
        else:
            logging.warning("`{0}` already exists. Checking for needle output".\
                    format(pairsfile))

        # If needle scores do not exist, prompt user to run needle
        if not os.path.isfile(scoresfile):
            logging.error("`{0}` does not exist. Please process {1} using `needle`".\
                    format(scoresfile, pairsfile))
            sys.exit()
    else:
        scoresfile = "ovl.scores"
        # Calculate overlap length using intersectBed
        calculate_ovl(nbedfile, obedfile, opts, scoresfile)

    logging.warning("`{0}' exists. Storing scores in memory".\
            format(scoresfile))
    scores = read_scores(scoresfile, opts)

    # Iterate through consolidated bed and
    # filter piles based on score
    abedline = {}

    cbed = Bed(cbedfile)
    g = Grouper()
    for c in cbed:
        accn = c.accn
        g.join(*accn.split(";"))

    nbedline = {}
    nbed = Bed(nbedfile)
    for line in nbed: nbedline[line.accn] = line

    splits = set()
    for chr, chrbed in nbed.sub_beds():
        abedline, splits = annotate_chr(chr, chrbed, g, scores, nbedline, abedline, opts, splits)

    if splits is not None:
        abedline = process_splits(splits, scores, nbedline, abedline)

    abedfile = npf + ".annotated.bed"
    afh = open(abedfile, "w")
    for accn in abedline:
        print >> afh, abedline[accn]
    afh.close()

    sort([abedfile, "-i"])


def calculate_ovl(nbedfile, obedfile, opts, scoresfile):
    from pybedtools import BedTool
    nbedtool = BedTool(nbedfile)
    obedtool = BedTool(obedfile)

    ab = nbedtool.intersect(obedtool, wao=True, f=opts.f, r=opts.r, s=opts.s)
    cmd = """cut -f4,5,10,13 | awk -F $'\t' 'BEGIN { OFS = FS } ($3 != "."){ print $1,$3,$2,$4; }'"""
    sh(cmd, infile=ab.fn, outfile=scoresfile)


def read_scores(scoresfile, opts):
    scores = {}
    fp = must_open(scoresfile)
    for row in fp:
        (new, old, identity, score) = row.strip().split("\t")
        old = re.sub('\.\d+$', '', old)
        if opts.resolve == "alignment":
            match = re.search("\d+\/\d+\s+\(\s*(\d+\.\d+)%\)", identity)
            pid = match.group(1)
            if float(pid) < opts.pid or float(score) < opts.score:
                continue
        else:
            pid = identity

        if new not in scores:
            scores[new] = []

        scores[new].append((new, old, pid, score))

    return scores


def annotate_chr(chr, chrbed, g, scores, nbedline, abedline, opts, splits):
    current_chr = chr_number(chr)

    for line in chrbed:
        accn = line.accn
        if accn not in g or (opts.atg_name and not current_chr):
            abedline[accn] = line
            continue

        gaccns = g[accn]
        new = [a for a in gaccns if re.search(new_id_pat, a)]
        newgrp = ";".join(sorted(new))

        if accn in scores:
            scores[accn] = sorted(scores[accn], key=lambda x: x[1])
            scores[accn] = sorted(scores[accn], key=lambda x: float(x[3]), reverse=True)

            accns = []
            print >> sys.stderr, accn
            for elem in scores[accn]:
                print >> sys.stderr, "\t" + ", ".join([str(x)\
                        for x in elem[1:]])
                if opts.atg_name:
                    achr, arank = atg_name(elem[1])
                    if not achr or achr != current_chr:
                        continue

                accns.append(elem[1])
                if len(new) > 1:
                    if newgrp not in scores: scores[newgrp] = []
                    scores[newgrp].append(elem)
                else:
                    accns[0:0] = [accn]
                    line.accn = ";".join([str(x) for x in accns])
                if len(scores[accn]) > 1: break

        if len(new) > 1:
            splits.add(newgrp)
        else:
            abedline[line.accn] = line

    return abedline, splits


def process_splits(splits, scores, nbedline, abedline):
    for newgrp in splits:
        new = newgrp.split(";")
        print >> sys.stderr, new
        if newgrp in scores:
            best = {}
            scores[newgrp] = sorted(scores[newgrp], key=lambda x: (x[0], x[1]))
            scores[newgrp] = sorted(scores[newgrp], key=lambda x: float(x[3]), reverse=True)

            for elem in scores[newgrp]:
                if elem[1] not in best:
                    best[elem[1]] = elem[0]

            for n in new:
                line = nbedline[n]
                if n in scores:
                    accns = set()
                    scores[n] = sorted(scores[n], key=lambda x: x[1])
                    scores[n] = sorted(scores[n], key=lambda x: float(x[3]), reverse=True)
                    accns.add(n)
                    print >> sys.stderr, "\t" + n
                    for elem in scores[n]:
                        if not elem[0] == n:
                            continue
                        print >> sys.stderr, "\t\t" + ", ".join([str(x)\
                                for x in elem[1:]])
                        if elem[1] in best and n == best[elem[1]]:
                            accns.add(elem[1])
                            accns = sorted(accns)
                            line.accn = ";".join([str(x) for x in accns])
                            break
                abedline[line.accn] = line
        else:
            for n in new:
                abedline[n] = nbedline[n]

    return abedline


def get_pairs(cbedfile, pairsfile):
    fp = open(pairsfile, "w")
    bed = Bed(cbedfile)
    for b in bed:
        if ";" in b.accn:
            genes = b.accn.split(";")
            new = [x for x in genes if re.search(new_id_pat, x)]
            old = [x for x in genes if not re.search(new_id_pat, x)]
            for a, b in product(new, old):
                print >> fp, "\t".join((a, b))

    fp.close()


def consolidate(nbedfile, obedfile, cbedfile):
    from pybedtools import BedTool
    nbedtool = BedTool(nbedfile)
    obedtool = BedTool(obedfile)

    ab = nbedtool.intersect(obedtool, s=True, u=True)
    ba = obedtool.intersect(nbedtool, s=True, u=True)

    cmd = "cat {0} {1} | sort -k1,1 -k2,2n".format(ab.fn, ba.fn)
    fp = popen(cmd)
    ovl = BedTool(fp.readlines())

    abmerge = ovl.merge(s=True, nms=True, scores="mean").sort()
    cmd = "cat {0}".format(abmerge.fn)
    fp = popen(cmd, debug=False)
    ovl = BedTool(fp.readlines())

    notovl = nbedtool.intersect(ovl.sort(), s=True, v=True)

    infile = "{0} {1}".format(notovl.fn, ovl.fn)
    tmpfile = "/tmp/reformat.{0}.bed".format(os.getpid())
    cmd = "sort -k1,1 -k2,2n"
    sh(cmd, infile=infile, outfile=tmpfile)

    fp = open(cbedfile, "w")
    bed = Bed(tmpfile)
    for b in bed:
        if ";" in b.accn:
            accns = set()
            for accn in b.accn.split(";"):
                accns.add(accn)
            b.accn = ";".join(accns)
        print >> fp, b
    fp.close()
    os.remove(tmpfile)

    sort([cbedfile, "-i"])


def rename(args):
    """
    %prog rename genes.bed [gaps.bed]

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

    Gaps bed file is optional.
    """
    import string

    p = OptionParser(rename.__doc__)
    p.add_option("-a", dest="gene_increment", default=10, type="int",
                 help="Increment for continuous genes [default: %default]")
    p.add_option("-b", dest="gap_increment", default=1000, type="int",
                 help="Increment for gaps [default: %default]")
    p.add_option("--pad0", default=6, type="int",
                 help="Pad gene identifiers with 0 [default: %default]")
    p.add_option("--spad0", default=4, type="int",
                 help="Pad gene identifiers on small scaffolds [default: %default]")
    p.add_option("--prefix", default="Bo",
                 help="Genome prefix [default: %default]")
    p.add_option("--jgi", default=False, action="store_true",
                 help="Create JGI style identifier PREFIX.NN[G|TE]NNNNN.1" + \
                      " [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    genebed = args[0]
    gapbed = args[1] if len(args) == 2 else None
    prefix = opts.prefix
    gene_increment = opts.gene_increment
    gap_increment = opts.gap_increment

    genes = Bed(genebed)
    if gapbed:
        fp = open(gapbed)
        for row in fp:
            genes.append(BedLine(row))

    genes.sort(key=genes.key)
    idsfile = prefix + ".ids"
    newbedfile = prefix + ".bed"
    gap_increment -= gene_increment
    assert gap_increment >= 0

    if opts.jgi:
        prefix += "."
    fw = open(idsfile, "w")
    for chr, lines in groupby(genes, key=lambda x: x.seqid):
        lines = list(lines)
        pad0 = opts.pad0 if len(lines) > 1000 else opts.spad0
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


def reindex(args):
    """
    %prog reindex idsfile > idsfiles.reindex

    Given a list of gene model identifier(following ATG naming conventions), reindex
    all the model IDs with new isoform indices

    Example output:
	Medtr2g100130.1         Medtr2g100130.1		SAME
	Medtr2g100130.3         Medtr2g100130.2		CHANGED
	Medtr2g100130.4         Medtr2g100130.3		CHANGED
	Medtr2g100130.5         Medtr2g100130.4		CHANGED
	Medtr2g100130.9         Medtr2g100130.5		CHANGED
	Medtr2g100130.11        Medtr2g100130.6		CHANGED
	Medtr2g100130.12        Medtr2g100130.7		CHANGED
	Medtr2g100130.13        Medtr2g100130.8		CHANGED

    Last column has 3-flags: "SAME" (Isoform index remains the same),
                             "CHANGED" (Isoform index has changed),
                         and "CHECK" (Isoform index for primary isoform is missing).
    """
    p = OptionParser(reindex.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    index = {}
    idsfile, = args
    fp = must_open(idsfile)
    for row in fp:
        locus, iso = atg_name(row, retval="locus,iso")
        if None in (locus, iso):
            logging.warning("{0} is not a valid gene model identifier".format(row))
            continue
        if locus not in index:
            index[locus] = set()

        index[locus].add(int(iso))

    for locus in index:
        index[locus] = sorted(index[locus])
        l = len(index[locus])
        new = range(1,l+1)
        for idx, (i, ni) in enumerate(zip(index[locus], new)):
            flag = "SAME" if i == ni else "CHANGED"
            if idx == 0 and i != 1:
                flag = "CHECK"
            print "\t".join(x for x in ("{0}.{1}".format(locus, i), \
                                        "{0}.{1}".format(locus, ni), \
                                        flag))


def publocus(args):
    """
    %prog publocus idsfile > idsfiles.publocus

    Given a list of model identifiers, convert each into a GenBank approved
    pub_locus.

    Example output:
    Medtr1g007020.1		MTR_1g007020
    Medtr1g007030.1		MTR_1g007030
    Medtr1g007060.1		MTR_1g007060A
    Medtr1g007060.2		MTR_1g007060B
    """
    from jcvi.utils.cbook import AutoVivification

    p = OptionParser(publocus.__doc__)
    p.add_option("--locus_tag", default="MTR_",
                 help="GenBank locus tag [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    locus_tag = opts.locus_tag

    index = AutoVivification()
    idsfile, = args
    fp = must_open(idsfile)
    for row in fp:
        locus, chrom, sep, rank, iso = atg_name(row, retval="locus,chr,sep,rank,iso")
        if None in (locus, chrom, sep, rank, iso):
            logging.warning("{0} is not a valid gene model identifier".format(row))
            continue
        if locus not in index.keys():
            pub_locus = gene_name(chrom, rank, prefix=locus_tag, sep=sep)
            index[locus]['pub_locus'] = pub_locus
            index[locus]['isos'] = set()

        index[locus]['isos'].add(int(iso))

    for locus in index:
        pub_locus = index[locus]['pub_locus']
        Index[locus]['isos'] = sorted(index[locus]['isos'])
        if len(index[locus]['isos']) > 1:
            new = [chr(n+64) for n in index[locus]['isos'] if n < 27]
            for i, ni in zip(index[locus]['isos'], new):
                print "\t".join(x for x in ("{0}.{1}".format(locus, i), \
                                            "{0}{1}".format(pub_locus, ni)))
        else:
            print "\t".join(x for x in ("{0}.{1}".format(locus, index[locus]['isos'][0]), \
                                        pub_locus))


def augustus(args):
    """
    %prog augustus augustus.gff3 > reformatted.gff3

    AUGUSTUS does generate a gff3 (--gff3=on) but need some refinement.
    """
    from jcvi.formats.gff import Gff

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
