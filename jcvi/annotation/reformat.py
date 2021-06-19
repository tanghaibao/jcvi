#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert common output files from gene prediction software into gff3 format.

Similar to the utilities in DAWGPAWS.
<http://dawgpaws.sourceforge.net/man.html>
"""
import os
import sys
import logging
import re

from collections import defaultdict
from itertools import groupby, product

from jcvi.utils.cbook import AutoVivification
from jcvi.formats.bed import Bed, BedLine, sort
from jcvi.formats.base import SetFile, must_open, get_number, flexible_cast
from jcvi.apps.base import (
    OptionParser,
    OptionGroup,
    ActionDispatcher,
    need_update,
    popen,
    sh,
)


FRAME, RETAIN, OVERLAP, NEW = "FRAME", "RETAIN", "OVERLAP", "NEW"
PRIORITY = (FRAME, RETAIN, OVERLAP, NEW)

new_id_pat = re.compile(r"^\d+\.[cemtx]+\S+")
atg_name_pat = re.compile(
    r"""
        ^(?P<locus>
            (?:(?P<prefix>\w+[\D\d\D])\.?)(?P<chr>[\d|C|M]+)(?P<sep>[A-z]+)(?P<rank>\d+)
        )
        \.?(?P<iso>\d+)?
        """,
    re.VERBOSE,
)


class Stride(object):
    """
      Allows four basic strides and three extended strides:
                      __.
          0 10          |
         0 5 10         | basic set of strides
        0 3 7 10        |
       0 2 5 8 10     __|
      0 2 4 6 8 10      |
     0 1 3 5 7 9 10     | extended set of strides
    01 23 45 67 89 10 __|

      We have main parameters, # we need, # available go through all possible
      numbers excluding everything in black.
    """

    def __init__(self, needed, available, extended=False):
        configurations = ("0", "05", "037", "0258")
        if extended:
            configurations += ("02468", "013579", "0123456789")
        nneeded = len(needed)
        self.conf = None
        self.available = None
        for c in configurations:
            a = [x for x in available if str(x)[-1] in c]
            if len(a) >= nneeded:
                self.conf = c
                self.available = a
                break


class NameRegister(object):
    def __init__(self, prefix="Medtr", pad0=6, uc=False):
        self.black = set()
        self.gaps = []
        self.prefix = prefix
        self.pad0 = pad0
        self.uc = uc

    def get_blacklist(self, filename):
        black = SetFile(filename)
        for x in black:
            chr, rank = atg_name(x)
            self.black.add((chr, rank))

    def get_gaps(self, filename):
        self.gapfile = filename

    def allocate(self, info, chr, start_id, end_id, id_table, extended_stride=False):

        start_bp = info[0].start
        end_bp = info[-1].end

        current_chr = chr_number(chr)
        needed = info
        assert end_id > start_id, "end ({0}) > start ({1})".format(end_id, start_id)

        spots = end_id - start_id - 1
        available = [
            x for x in range(start_id + 1, end_id) if (current_chr, x) not in self.black
        ]

        message = "{0} need {1} ids, has {2} spots ({3} available)".format(
            chr, len(needed), spots, len(available)
        )

        start_gene = gene_name(
            current_chr, start_id, prefix=self.prefix, pad0=self.pad0, uc=self.uc
        )
        end_gene = gene_name(
            current_chr, end_id, prefix=self.prefix, pad0=self.pad0, uc=self.uc
        )
        message += " between {0} - {1}\n".format(start_gene, end_gene)

        assert end_bp > start_bp

        b = "\t".join(str(x) for x in (chr, start_bp - 1, end_bp))
        cmd = "echo '{0}' |".format(b)
        cmd += " intersectBed -a {0} -b stdin".format(self.gapfile)
        gaps = list(BedLine(x) for x in popen(cmd, debug=False))
        ngaps = len(gaps)

        gapsexpanded = []
        GeneDensity = 10000.0  # assume 10Kb per gene
        for gap in gaps:
            gap_bp = int(gap.score)
            gap_ids = int(round(gap_bp / GeneDensity))
            gapsexpanded += [gap] * gap_ids

        lines = sorted(info + gapsexpanded, key=lambda x: x.start)

        message += "between bp: {0} - {1}, there are {2} gaps (total {3} ids)".format(
            start_bp, end_bp, ngaps, len(lines)
        )

        needed = lines
        stride = Stride(needed, available, extended=extended_stride)
        conf = stride.conf
        message += " stride: {0}".format(conf)
        print(message, file=sys.stderr)

        nneeded = len(needed)
        if conf is None:  # prefix rule - prepend version number for spills
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

        else:  # follow the best stride
            available = stride.available
            if start_id == 0:  # follow right flank at start of chr
                available = available[-nneeded:]
            else:  # follow left flank otherwise
                available = available[:nneeded]

        # Finally assign the ids
        assert len(needed) == len(available)
        for b, rank in zip(needed, available):
            name = gene_name(
                current_chr, rank, prefix=self.prefix, pad0=self.pad0, uc=self.uc
            )
            print("\t".join((str(b), name)), file=sys.stderr)
            id_table[b.accn] = name
            self.black.add((current_chr, rank))
        print(file=sys.stderr)


def main():

    actions = (
        ("rename", "rename genes for annotation release"),
        # perform following actions on list files
        ("reindex", "reindex isoforms per gene locus"),
        ("publocus", "create pub_locus identifiers according to GenBank specs"),
        # Medicago gene renumbering
        ("annotate", "annotation new bed file with features from old"),
        ("renumber", "renumber genes for annotation updates"),
        ("instantiate", "instantiate NEW genes tagged by renumber"),
        ("plot", "plot gene identifiers along certain chromosome"),
        # External gene prediction programs
        ("augustus", "convert augustus output into gff3"),
        ("tRNAscan", "convert tRNAscan-SE output into gff3"),
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
    p.add_option("--log", action="store_true", help="Write plotting data")
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

    xstart, xend = 0.2, 0.8
    ystart, yend = 0.2, 0.8
    pad = 0.02

    ngenes = opts.firstn or ngenes
    ymax = opts.ymax or 500000

    title = "Assignment of Medtr identifiers"
    if opts.ymax:
        subtitle = "{0}, first {1} genes".format(chr, ngenes)
    else:
        subtitle = "{0}, {1} genes ({2} new)".format(chr, ngenes, len(new))

    chr_map = ChromosomeMap(
        fig, root, xstart, xend, ystart, yend, pad, 0, ymax, 5, title, subtitle
    )

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
    root.plot([0.2], [y], "r.", lw=2)
    root.text(0.2 + pad, y, "Existing Medtr ids", va="center", size=10)
    y = ymid - pad
    root.plot([0.2], [y], "b.", lw=2)
    root.text(0.2 + pad, y, "Newly instantiated ids", va="center", size=10)

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
    p.add_option(
        "--extended_stride",
        default=False,
        action="store_true",
        help="Toggle extended strides for gene numbering",
    )
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
                start, end = (
                    atg_name(start, retval="rank"),
                    atg_name(end, retval="rank"),
                )
                blocks.append((tag, [start, end]))

        id_table = {}  # old to new name conversion
        for i, (tag, info) in enumerate(blocks):
            if tag != NEW:
                continue

            start_id = 0 if i == 0 else blocks[i - 1][1][-1]
            end_id = start_id + 10000 if i == len(blocks) - 1 else blocks[i + 1][1][0]

            r.allocate(
                info,
                chr,
                start_id,
                end_id,
                id_table,
                extended_stride=opts.extended_stride,
            )

        # Output new names
        for i, s in enumerate(sbed):
            nametag = s.extra[0]
            name, tag = nametag.split("|")

            if tag == NEW:
                assert name == "."
                name = id_table[s.accn]
            elif tag == OVERLAP:
                if name in id_table:
                    name = id_table[name]

            s.extra[0] = "|".join((name, tag))
            print(s, file=fw)

    fw.close()


def atg_name(name, retval="chr,rank", trimpad0=True):
    seps = ["g", "te", "trna", "s", "u", "nc"]
    pad0s = ["rank"]

    if name is not None:
        m = re.match(atg_name_pat, name)
        if m is not None and m.group("sep").lower() in seps:
            retvals = []
            for grp in retval.split(","):
                if grp == "chr":
                    val = chr_number(m.group(grp))
                else:
                    val = (
                        get_number(m.group(grp))
                        if trimpad0 and grp in pad0s
                        else m.group(grp)
                    )
                retvals.append(val)

            return (x for x in retvals) if len(retvals) > 1 else retvals[0]

    return (None for _ in retval.split(","))


def gene_name(current_chr, x, prefix="Medtr", sep="g", pad0=6, uc=False):
    identifier = "{0}{1}{2}{3:0{4}}".format(prefix, current_chr, sep, x, pad0)
    if uc:
        identifier = identifier.upper()
    return identifier


def chr_number(chr):
    chr_pat = re.compile(
        r"(?P<prefix>\D*)(?P<chr>[\d|C|M]+)$", re.VERBOSE | re.IGNORECASE
    )

    if chr is not None:
        m = re.match(chr_pat, chr)
        if m is not None:
            return flexible_cast(m.group("chr"))

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
            print(b, file=fwa)
            seen.add(accn)

        b.accn = ";".join(new_accns)
        print(b, file=fwb)
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

    (bedfile,) = args

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
        print(
            current_chr,
            len(sbed),
            "==>",
            len(ranks),
            "==>",
            len(lranks),
            file=sys.stderr,
        )

        granks = set(
            gene_name(current_chr, x, prefix=opts.prefix, pad0=opts.pad0, uc=opts.uc)
            for x in lranks
        ) | set(
            gene_name(
                current_chr, x, prefix=opts.prefix, pad0=opts.pad0, sep="te", uc=opts.uc
            )
            for x in lranks
        )

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

            print("\t".join((str(s), "|".join(tag))))


def annotate(args):
    r"""
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
    p.add_option(
        "--resolve",
        default="alignment",
        choices=valid_resolve_choices,
        help="Resolve ID assignment based on a certain metric",
    )
    p.add_option(
        "--atg_name",
        default=False,
        action="store_true",
        help="Specify is locus IDs in `new.bed` file follow ATG nomenclature",
    )

    g1 = OptionGroup(
        p,
        "Optional parameters (alignment):\n"
        + "Use if resolving ambiguities based on sequence `alignment`",
    )
    g1.add_option(
        "--pid",
        dest="pid",
        default=35.0,
        type="float",
        help="Percent identity cutoff",
    )
    g1.add_option(
        "--score",
        dest="score",
        default=250.0,
        type="float",
        help="Alignment score cutoff",
    )
    p.add_option_group(g1)

    g2 = OptionGroup(
        p,
        "Optional parameters (overlap):\n"
        + "Use if resolving ambiguities based on `overlap` length\n"
        + "Parameters equivalent to `intersectBed`",
    )
    g2.add_option(
        "-f",
        dest="f",
        default=0.5,
        type="float",
        help="Minimum overlap fraction (0.0 - 1.0)",
    )
    g2.add_option(
        "-r",
        dest="r",
        default=False,
        action="store_true",
        help="Require fraction overlap to be reciprocal",
    )
    g2.add_option(
        "-s",
        dest="s",
        default=True,
        action="store_true",
        help="Require same strandedness",
    )
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

    logging.warning(
        "Resolving ID assignment ambiguity based on `{0}`".format(opts.resolve)
    )

    if opts.resolve == "alignment":
        # Get pairs and prompt to run needle
        pairsfile = "nw.pairs"
        scoresfile = "nw.scores"
        if not os.path.isfile(pairsfile):
            get_pairs(cbedfile, pairsfile)
        else:
            logging.warning(
                "`{0}` already exists. Checking for needle output".format(pairsfile)
            )

        # If needle scores do not exist, prompt user to run needle
        if not os.path.isfile(scoresfile):
            logging.error(
                "`{0}` does not exist. Please process {1} using `needle`".format(
                    scoresfile, pairsfile
                )
            )
            sys.exit()
    else:
        scoresfile = "ovl.scores"
        # Calculate overlap length using intersectBed
        calculate_ovl(nbedfile, obedfile, opts, scoresfile)

    logging.warning("`{0}' exists. Storing scores in memory".format(scoresfile))
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
    for line in nbed:
        nbedline[line.accn] = line

    splits = set()
    for chr, chrbed in nbed.sub_beds():
        abedline, splits = annotate_chr(chr, chrbed, g, scores, abedline, opts, splits)

    if splits is not None:
        abedline = process_splits(splits, scores, nbedline, abedline)

    abedfile = npf + ".annotated.bed"
    afh = open(abedfile, "w")
    for accn in abedline:
        print(abedline[accn], file=afh)
    afh.close()

    sort([abedfile, "-i"])


def calculate_ovl(nbedfile, obedfile, opts, scoresfile):
    from pybedtools import BedTool

    nbedtool = BedTool(nbedfile)
    obedtool = BedTool(obedfile)

    ab = nbedtool.intersect(obedtool, wao=True, f=opts.f, r=opts.r, s=opts.s)
    cmd = """cut -f4,5,10,13 | awk -F $'\t' 'BEGIN { OFS = FS } ($3 != "."){ print $1,$3,$2,$4; }'"""
    sh(cmd, infile=ab.fn, outfile=scoresfile)


def read_scores(scoresfile, opts=None, sort=False, trimsuffix=True):
    scores = {}
    _pid, _score, resolve = (
        (0.0, 0.0, "alignment")
        if opts is None
        else (opts.pid, opts.score, opts.resolve)
    )

    fp = must_open(scoresfile)
    logging.debug("Load scores file `{0}`".format(scoresfile))
    for row in fp:
        (new, old, identity, score) = row.strip().split("\t")
        if trimsuffix:
            old = re.sub(r"\.\d+$", "", old)
        if resolve == "alignment":
            match = re.search(r"\d+/\d+\s+\(\s*(\d+\.\d+)%\)", identity)
            pid = match.group(1)
            if float(pid) < _pid or float(score) < _score:
                continue
        else:
            pid = identity

        if new not in scores:
            scores[new] = []

        scores[new].append((new, old, float(pid), float(score)))

    if sort:
        for new in scores:
            scores[new].sort(key=lambda k: (-k[2], -k[3]))

    return scores


def annotate_chr(chr, chrbed, g, scores, abedline, opts, splits):
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
            print(accn, file=sys.stderr)
            for elem in scores[accn]:
                print("\t" + ", ".join([str(x) for x in elem[1:]]), file=sys.stderr)
                if opts.atg_name:
                    achr, arank = atg_name(elem[1])
                    if not achr or achr != current_chr:
                        continue

                accns.append(elem[1])
                if len(new) > 1:
                    if newgrp not in scores:
                        scores[newgrp] = []
                    scores[newgrp].append(elem)
                else:
                    accns[0:0] = [accn]
                    line.accn = ";".join([str(x) for x in accns])
                if len(scores[accn]) > 1:
                    break

        if len(new) > 1:
            splits.add(newgrp)
        else:
            abedline[line.accn] = line

    return abedline, splits


def process_splits(splits, scores, nbedline, abedline):
    for newgrp in splits:
        new = newgrp.split(";")
        print(new, file=sys.stderr)
        if newgrp in scores:
            best = {}
            scores[newgrp] = sorted(scores[newgrp], key=lambda x: (x[0], x[1]))
            scores[newgrp] = sorted(
                scores[newgrp], key=lambda x: float(x[3]), reverse=True
            )

            for elem in scores[newgrp]:
                if elem[1] not in best:
                    best[elem[1]] = elem[0]

            for n in new:
                line = nbedline[n]
                if n in scores:
                    accns = set()
                    scores[n] = sorted(scores[n], key=lambda x: x[1])
                    scores[n] = sorted(
                        scores[n], key=lambda x: float(x[3]), reverse=True
                    )
                    accns.add(n)
                    print("\t" + n, file=sys.stderr)
                    for elem in scores[n]:
                        if not elem[0] == n:
                            continue
                        print(
                            "\t\t" + ", ".join([str(x) for x in elem[1:]]),
                            file=sys.stderr,
                        )
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
                print("\t".join((a, b)), file=fp)

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
        print(b, file=fp)
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
    p.add_option(
        "-a",
        dest="gene_increment",
        default=10,
        type="int",
        help="Increment for continuous genes",
    )
    p.add_option(
        "-b",
        dest="gap_increment",
        default=1000,
        type="int",
        help="Increment for gaps",
    )
    p.add_option(
        "--pad0",
        default=6,
        type="int",
        help="Pad gene identifiers with 0",
    )
    p.add_option(
        "--spad0",
        default=4,
        type="int",
        help="Pad gene identifiers on small scaffolds",
    )
    p.add_option("--prefix", default="Bo", help="Genome prefix")
    p.add_option(
        "--jgi",
        default=False,
        action="store_true",
        help="Create JGI style identifier PREFIX.NN[G|TE]NNNNN.1",
    )
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
        isChr = chr[0].upper() == "C"
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
            print("\t".join((oldaccn, accn)), file=fw)
            r.accn = accn

    genes.print_to_file(newbedfile)
    logging.debug("Converted IDs written to `{0}`.".format(idsfile))
    logging.debug("Converted bed written to `{0}`.".format(newbedfile))


def parse_prefix(identifier):
    """
    Parse identifier such as a|c|le|d|li|re|or|AT4G00480.1 and return
    tuple of prefix string (separated at '|') and suffix (AGI identifier)
    """
    pf, id = (), identifier
    if "|" in identifier:
        pf, id = tuple(identifier.split("|")[:-1]), identifier.split("|")[-1]

    return pf, id


def reindex(args):
    """
    %prog reindex gffile pep.fasta ref.pep.fasta

    Reindex the splice isoforms (mRNA) in input GFF file, preferably
    generated after PASA annotation update

    In the input GFF file, there can be several types of mRNA within a locus:
    * CDS matches reference, UTR extended, inherits reference mRNA ID
    * CDS (slightly) different from reference, inherits reference mRNA ID
    * Novel isoform added by PASA, have IDs like "LOCUS.1.1", "LOCUS.1.2"
    * Multiple mRNA collapsed due to shared structure, have IDs like "LOCUS.1-LOCUS.1.1"

    In the case of multiple mRNA which have inherited the same reference mRNA ID,
    break ties by comparing the new protein with the reference protein using
    EMBOSS `needle` to decide which mRNA retains ID and which is assigned a new ID.

    All mRNA identifiers should follow the AGI naming conventions.

    When reindexing the isoform identifiers, order mRNA based on:
    * decreasing transcript length
    * decreasing support from multiple input datasets used to run pasa.consolidate()
    """
    from jcvi.formats.gff import make_index
    from jcvi.formats.fasta import Fasta
    from jcvi.apps.emboss import needle
    from jcvi.formats.base import FileShredder
    from tempfile import mkstemp

    p = OptionParser(reindex.__doc__)
    p.add_option(
        "--scores", type="str", help="read from existing EMBOSS `needle` scores file"
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    (
        gffile,
        pep,
        refpep,
    ) = args
    gffdb = make_index(gffile)
    reffasta = Fasta(refpep)

    if not opts.scores:
        fh, pairsfile = mkstemp(prefix="pairs", suffix=".txt", dir=".")
        fw = must_open(pairsfile, "w")

    conflict, novel = AutoVivification(), {}
    for gene in gffdb.features_of_type("gene", order_by=("seqid", "start")):
        geneid = atg_name(gene.id, retval="locus")
        novel[geneid] = []
        updated_mrna, hybrid_mrna = [], []
        for mrna in gffdb.children(
            gene, featuretype="mRNA", order_by=("seqid", "start")
        ):
            if re.match(atg_name_pat, mrna.id) is not None and "_" not in mrna.id:
                pf, mrnaid = parse_prefix(mrna.id)
                mlen = gffdb.children_bp(mrna, child_featuretype="exon")
                if "-" in mrna.id:
                    hybrid_mrna.append((mrna.id, mrna.start, mlen, len(pf)))
                else:
                    updated_mrna.append((mrna.id, mrna.start, mlen, len(pf)))

        for mrna in sorted(updated_mrna, key=lambda k: (k[1], -k[2], -k[3])):
            pf, mrnaid = parse_prefix(mrna[0])
            mstart, mlen = mrna[1], mrna[2]

            iso = atg_name(mrnaid, retval="iso")
            newiso = "{0}{1}".format(iso, re.sub(atg_name_pat, "", mrnaid))
            if iso == newiso:
                if iso not in conflict[geneid]:
                    conflict[geneid][iso] = []
                conflict[geneid][iso].append(
                    (mrna[0], iso, newiso, mstart, mlen, len(pf))
                )
            else:
                novel[geneid].append((mrna[0], None, newiso, mstart, mlen, len(pf)))

        for mrna in sorted(hybrid_mrna, key=lambda k: (k[1], -k[2], -k[3])):
            pf, mrnaid = parse_prefix(mrna[0])
            mstart, mlen = mrna[1], mrna[2]

            _iso, _newiso = [], []
            for id in sorted(mrnaid.split("-")):
                a = atg_name(id, retval="iso")
                b = "{0}{1}".format(a, re.sub(atg_name_pat, "", id))
                _iso.append(a)
                _newiso.append(b)

            _novel = None
            newiso = "-".join(str(x) for x in set(_newiso))
            for iso, niso in zip(_iso, _newiso):
                if iso == niso:
                    if iso not in conflict[geneid]:
                        conflict[geneid][iso] = [
                            (mrna[0], iso, newiso, mstart, mlen, len(pf))
                        ]
                        _novel = None
                        break

                _novel = True

            if _novel is not None:
                novel[geneid].append((mrna[0], None, newiso, mstart, mlen, len(pf)))

        if not opts.scores:
            for isoform in sorted(conflict[geneid]):
                mrnaid = "{0}.{1}".format(geneid, isoform)
                if mrnaid in reffasta.keys():
                    for mrna in conflict[geneid][isoform]:
                        print("\t".join(str(x) for x in (mrnaid, mrna[0])), file=fw)

    if not opts.scores:
        fw.close()
        needle([pairsfile, refpep, pep])
        FileShredder([pairsfile], verbose=False)
        scoresfile = "{0}.scores".format(pairsfile.rsplit(".")[0])
    else:
        scoresfile = opts.scores

    scores = read_scores(scoresfile, sort=True, trimsuffix=False)

    primary = {}
    for geneid in conflict:
        primary[geneid] = []
        for iso in sorted(conflict[geneid]):
            conflict[geneid][iso].sort(key=lambda k: (k[3], -k[4], -k[5]))
            _iso = "{0}.{1}".format(geneid, iso)
            if _iso not in scores:
                novel[geneid].extend(conflict[geneid][iso])
                continue
            top_score = scores[_iso][0][1]
            result = next(
                (i for i, v in enumerate(conflict[geneid][iso]) if v[0] == top_score),
                None,
            )
            if result is not None:
                primary[geneid].append(conflict[geneid][iso][result])
                del conflict[geneid][iso][result]
                if geneid not in novel:
                    novel[geneid] = []
                novel[geneid].extend(conflict[geneid][iso])
        novel[geneid].sort(key=lambda k: (k[3], -k[4], -k[5]))

    fw = must_open(opts.outfile, "w")
    for gene in gffdb.features_of_type("gene", order_by=("seqid", "start")):
        geneid = gene.id
        print(gene, file=fw)
        seen = []
        if geneid in primary:
            all_mrna = primary[geneid]
            all_mrna.extend(novel[geneid])
            for iso, mrna in enumerate(all_mrna):
                _mrna = gffdb[mrna[0]]
                _iso = mrna[1]
                if mrna not in novel[geneid]:
                    seen.append(int(mrna[1]))
                else:
                    mseen = 0 if len(seen) == 0 else max(seen)
                    _iso = (mseen + iso + 1) - len(seen)

                _mrnaid = "{0}.{1}".format(geneid, _iso)
                _mrna["ID"], _mrna["_old_ID"] = [_mrnaid], [_mrna.id]

                print(_mrna, file=fw)
                for c in gffdb.children(_mrna, order_by="start"):
                    c["Parent"] = [_mrnaid]
                    print(c, file=fw)
        else:
            for feat in gffdb.children(gene, order_by=("seqid", "start")):
                print(feat, file=fw)

    fw.close()


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
    p = OptionParser(publocus.__doc__)
    p.add_option("--locus_tag", default="MTR_", help="GenBank locus tag")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    locus_tag = opts.locus_tag

    index = AutoVivification()
    (idsfile,) = args
    fp = must_open(idsfile)
    for row in fp:
        locus, chrom, sep, rank, iso = atg_name(row, retval="locus,chr,sep,rank,iso")
        if None in (locus, chrom, sep, rank, iso):
            logging.warning("{0} is not a valid gene model identifier".format(row))
            continue
        if locus not in index.keys():
            pub_locus = gene_name(chrom, rank, prefix=locus_tag, sep=sep)
            index[locus]["pub_locus"] = pub_locus
            index[locus]["isos"] = set()

        index[locus]["isos"].add(int(iso))

    for locus in index:
        pub_locus = index[locus]["pub_locus"]
        index[locus]["isos"] = sorted(index[locus]["isos"])
        if len(index[locus]["isos"]) > 1:
            new = [chr(n + 64) for n in index[locus]["isos"] if n < 27]
            for i, ni in zip(index[locus]["isos"], new):
                print(
                    "\t".join(
                        x
                        for x in (
                            "{0}.{1}".format(locus, i),
                            "{0}{1}".format(pub_locus, ni),
                        )
                    )
                )
        else:
            print(
                "\t".join(
                    x
                    for x in (
                        "{0}.{1}".format(locus, index[locus]["isos"][0]),
                        pub_locus,
                    )
                )
            )


def augustus(args):
    """
    %prog augustus augustus.gff3 > reformatted.gff3

    AUGUSTUS does generate a gff3 (--gff3=on) but need some refinement.
    """
    from jcvi.formats.gff import Gff

    p = OptionParser(augustus.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (ingff3,) = args
    gff = Gff(ingff3)
    fw = must_open(opts.outfile, "w")
    seen = defaultdict(int)
    for g in gff:
        if g.type not in ("gene", "transcript", "CDS"):
            continue

        if g.type == "transcript":
            g.type = "mRNA"

        prefix = g.seqid + "_"
        pid = prefix + g.id
        newid = "{0}-{1}".format(pid, seen[pid]) if pid in seen else pid
        seen[pid] += 1
        g.attributes["ID"] = [newid]
        g.attributes["Parent"] = [(prefix + x) for x in g.attributes["Parent"]]
        g.update_attributes()
        print(g, file=fw)
    fw.close()


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

    (trnaout,) = args
    gffout = trnaout + ".gff3"
    fp = open(trnaout)
    fw = open(gffout, "w")

    next(fp)
    next(fp)
    row = next(fp)
    assert row.startswith("--------")

    for row in fp:
        atoms = [x.strip() for x in row.split("\t")]
        contig, trnanum, start, end, aa, codon, intron_start, intron_end, score = atoms

        start, end = int(start), int(end)
        orientation = "+"
        if start > end:
            start, end = end, start
            orientation = "-"

        source = "tRNAscan"
        type = "tRNA"
        if codon == "???":
            codon = "XXX"

        comment = "ID={0}.tRNA.{1};Name=tRNA-{2} (anticodon: {3})".format(
            contig, trnanum, aa, codon
        )

        print(
            "\t".join(
                str(x)
                for x in (
                    contig,
                    source,
                    type,
                    start,
                    end,
                    score,
                    orientation,
                    ".",
                    comment,
                )
            ),
            file=fw,
        )

    fw.close()
    sort([gffout, "-i"])


if __name__ == "__main__":
    main()
