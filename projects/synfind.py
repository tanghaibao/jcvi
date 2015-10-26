#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""

import os.path as op
import sys
import logging

from copy import deepcopy
from collections import defaultdict
from itertools import groupby

from jcvi.formats.base import get_number, must_open
from jcvi.utils.cbook import SummaryStats, gene_name
from jcvi.utils.grouper import Grouper
from jcvi.formats.blast import BlastLine
from jcvi.formats.bed import Bed
from jcvi.graphics.base import FancyArrow, plt, savefig, panel_labels, markup
from jcvi.graphics.glyph import CartoonRegion, RoundRect
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, symlink, sh


def main():

    actions = (
        ('cartoon', 'generate cartoon illustration of SynFind'),
        ('ecoli', 'gene presence absence analysis in ecoli'),
        ('athaliana', 'prepare pairs data for At alpha/beta/gamma'),
        ('grasses', 'validate SynFind pan-grass set against James'),
        # For benchmarking
        ('iadhore', 'wrap around iADHoRe'),
        ("mcscanx", 'wrap around MCScanX'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def write_lst(bedfile):
    pf = op.basename(bedfile).split(".")[0]
    mkdir(pf)
    bed = Bed(bedfile)
    stanza = []
    for seqid, bs in bed.sub_beds():
        fname = op.join(pf, "{0}.lst".format(seqid))
        fw = open(fname, "w")
        for b in bs:
            print >> fw, "{0}{1}".format(b.accn, b.strand)
        stanza.append((seqid, fname))
        fw.close()
    return pf, stanza


def iadhore(args):
    """
    %prog iadhore athaliana.athaliana.last athaliana.bed

    Wrap around iADHoRe.
    """
    p = OptionParser(iadhore.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    lastfile = args[0]
    bedfiles = args[1:]
    blast_table = "blast_table.txt"
    fp = open(lastfile)
    seen = set()
    for row in fp:
        c = BlastLine(row)
        a, b = c.query, c.subject
        a, b = gene_name(a), gene_name(b)
        if a > b:
            a, b = b, a
        seen.add((a, b))

    fw = open(blast_table, "w")
    for a, b in seen:
        print >> fw, "\t".join((a, b))
    fw.close()
    logging.debug("A total of {0} pairs written to `{1}`"\
            .format(len(seen), blast_table))

    fw = open("config.txt", "w")
    for bedfile in bedfiles:
        pf, stanza = write_lst(bedfile)
        print >> fw, "genome={0}".format(pf)
        for seqid, fname in stanza:
            print >> fw,  " ".join((seqid, fname))
        print >> fw

    print >> fw, "blast_table={0}".format(blast_table)
    print >> fw, "cluster_type=colinear"
    print >> fw, "tandem_gap=10"
    print >> fw, "prob_cutoff=0.001"
    print >> fw, "gap_size=20"
    print >> fw, "cluster_gap=20"
    print >> fw, "q_value=0.9"
    print >> fw, "anchor_points=4"
    print >> fw, "alignment_method=gg2"
    print >> fw, "max_gaps_in_alignment=20"
    print >> fw, "output_path=i-adhore_out"
    print >> fw, "number_of_threads=4"
    fw.close()


def extract_groups(g, txtfile):
    register = defaultdict(list)
    fp = open(txtfile)
    fp.next()
    for row in fp:
        if row[0] != '>':
            continue
        track, atg, myname, pairname = row.split()
        pairname = pairname.rstrip("ab").upper()
        register[pairname].append(atg)
    for pairname, genes in register.items():
        g.join(*genes)


def athaliana(args):
    """
    %prog athaliana J_a.txt J_bc.txt

    Prepare pairs data for At alpha/beta/gamma.
    """
    p = OptionParser(athaliana.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    atxt, bctxt = args
    g = Grouper()
    for txt in (atxt, bctxt):
        extract_groups(g, txt)
    for group in list(g):
        print ",".join(group)


def make_gff(bed, fw):
    prefix = op.basename(bed)[:2]
    bed = Bed(bed)
    nfeats = 0
    for b in bed:
        try:
            seqid = prefix + str(get_number(b.seqid))
        except:
            continue
        print >> fw, "\t".join(str(x) for x in \
            (seqid, b.accn, b.start, b.end))
        nfeats += 1
    fw.close()
    logging.debug("A total of {0} features converted to `{1}`".\
                    format(nfeats, fw.name))


def mcscanx(args):
    """
    %prog mcscanx athaliana.athaliana.last athaliana.bed

    Wrap around MCScanX.
    """
    p = OptionParser(mcscanx.__doc__)
    p.add_option("--path", default="../bin/MCScanX",
                 help="Path to MCScanX")
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    blastfile = args[0]
    bedfiles = args[1:]
    prefix = "_".join(op.basename(x)[:2] for x in bedfiles)
    symlink(blastfile, prefix + ".blast")
    allbedfile = prefix + ".gff"
    fw = open(allbedfile, "w")
    for bedfile in bedfiles:
        make_gff(bedfile, fw)
    fw.close()

    cmd = opts.path + " " + prefix
    sh(cmd)
    outfile = prefix + ".collinearity"
    fw = open(prefix + ".mcscanx", "w")
    fp = open(outfile)
    seen = set()
    for row in fp:
        if row[0] == "#" or row.strip() == "":
            continue
        pairs = row.split(":")[1].split()
        a, b = pairs[:2]
        a, b = gene_name(a), gene_name(b)
        seen.add((a, b))

    for a, b in sorted(seen):
        print >> fw, "\t".join((a, b))
    fw.close()
    logging.debug("Total pairs: {0}".format(len(seen)))


def grasses(args):
    """
    %prog grasses coge_master_table.txt james.txt

    Validate SynFind pan-grass set against James. This set can be generated:

    https://genomevolution.org/r/fhak
    """
    p = OptionParser(grasses.__doc__)
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    master, james = args

    fp = open(master)
    fp.next()
    master_store = defaultdict(set)
    for row in fp:
        atoms = row.split()
        s = set()
        for x in atoms[1:6]:
            m = x.split(",")
            s |= set(m)
        if '-' in s:
            s.remove('-')

        a = atoms[1]
        master_store[a] |= set(s)

    fp = open(james)
    fp.next()
    james_store = {}
    tandems = set()
    for row in fp:
        atoms = row.split()
        s = set()
        Os = set()
        for x in atoms[:-1]:
            m = x.split("||")
            if m[0].startswith("Os"):
                Os |= set(m)
            if m[0].startswith("http"):
                continue
            if m[0].startswith("chr"):
                m = ["proxy"]
            if "||" in x:
                tandems |= set(m)
            s |= set(m)

        for x in Os:
            james_store[x] = s

    jaccards = []
    corr_jaccards = []
    perfect_matches = 0
    corr_perfect_matches = 0
    for k, v in james_store.items():
        if k not in master_store:
            continue
        m = master_store[k]
        jaccard = len(v & m) * 100 / len(v | m)
        jaccards.append(jaccard)
        diff = (v ^ m ) - tandems
        corr_jaccard = 100 - len(diff) * 100 / len(v | m)
        corr_jaccards.append(corr_jaccard)
        if opts.verbose:
            print k
            print v
            print m
            print diff
            print jaccard
        if jaccard > 99:
            perfect_matches += 1
        if corr_jaccard > 99:
            corr_perfect_matches += 1

    logging.debug("Perfect matches: {0}".format(perfect_matches))
    logging.debug("Perfect matches (corrected): {0}".format(corr_perfect_matches))
    print "Jaccards:", SummaryStats(jaccards)
    print "Corrected Jaccards:", SummaryStats(corr_jaccards)


def ecoli(args):
    """
    %prog ecoli coge_master_table.txt query.bed

    Perform gene presence / absence analysis in Ecoli master spreadsheet. Ecoli
    spresheets can be downloaded below:

    Ecoli K12 MG1655 (K) as query
    Regenerate this analysis: https://genomevolution.org/r/fggo

    Ecoli O157:H7 EDL933 (O) as query
    Regenerate this analysis: https://genomevolution.org/r/fgt7

    Shigella flexneri 2a 301 (S) as query
    Regenerate this analysis: https://genomevolution.org/r/fgte

    Perform a similar analysis as in:
    Jin et al. (2002) Genome sequence of Shigella flexneri 2a: insights
    into pathogenicity through comparison with genomes of Escherichia
    coli K12 and O157. Nucleic Acid Research.
    """
    p = OptionParser(ecoli.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    master, querybed = args
    fp = open(master)
    header = fp.next()
    assert header[0] == '#'
    qorg = header.strip().split("\t")[1]
    qorg = qorg.split(":")[-1].strip()

    store = {}
    MISSING = ("proxy", "-")
    for row in fp:
        a, b, c = row.strip().split("\t")[1:4]
        store[a] = b in MISSING and c in MISSING

    bed = Bed(querybed)
    tags = []
    for i, b in enumerate(bed):
        accn = b.accn
        if accn not in store:
            logging.warn("missing {0}".format(accn))
            continue
        tags.append((store[accn], accn))

    large = 4  # large segments
    II = []
    II_large = []
    for missing, aa in groupby(tags, key=lambda x: x[0]):
        aa = list(aa)
        if not missing:
            continue
        glist = list(a for missing, a in aa)
        II.append(glist)
        size = len(glist)
        if size >= large:
            II_large.append(glist)

    fw = must_open(opts.outfile, "w")
    for a, t in zip((II, II_large), ("", ">=4 ")):
        nmissing = sum(len(x) for x in a)
        logging.debug("A total of {0} {1}-specific {2}islands found with {3} genes.".\
                        format(len(a), qorg, t, nmissing))

    for x in II:
        print >> fw, len(x), ",".join(x)


def plot_diagram(ax, x, y, A, B, tag, label):
    ax.text(x, y + .14, "{0}: {1}".format(tag, label), ha="center")
    strip = tag != 'G'
    A.draw(ax, x, y + .06, gene_len=.02, strip=strip)
    B.draw(ax, x, y, gene_len=.02, strip=strip)


def cartoon(args):
    """
    %prog synteny.py

    Generate cartoon illustration of SynFind.
    """
    p = OptionParser(cartoon.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x7")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Panel A
    A = CartoonRegion(41)
    A.draw(root, .35, .85, strip=False, color=False)
    x1, x2 = A.x1, A.x2
    lsg = "lightslategray"
    pad = .01
    xc, yc = .35, .88
    arrowlen = x2 - xc - pad
    arrowprops = dict(length_includes_head=True, width=.01, fc=lsg, lw=0,
                      head_length=arrowlen * .15, head_width=.03)
    p = FancyArrow(xc - pad, yc, -arrowlen, 0, shape="left", **arrowprops)
    root.add_patch(p)
    p = FancyArrow(xc + pad, yc, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    yt = yc + 4 * pad
    root.text((x1 + xc) / 2, yt, "20 genes upstream", ha="center")
    root.text((x2 + xc) / 2, yt, "20 genes downstream", ha="center")
    root.plot((xc,), (yc,), "o", mfc='w', mec=lsg, mew=2, lw=2, color=lsg)
    root.text(xc, yt, "Query gene", ha="center")

    # Panel B
    A.draw(root, .35, .7, strip=False)

    RoundRect(root, (.07, .49), .56, .14, fc='y', alpha=.2)
    a = deepcopy(A)
    a.evolve(mode='S', target=10)
    a.draw(root, .35, .6)
    b = deepcopy(A)
    b.evolve(mode='F', target=8)
    b.draw(root, .35, .56)
    c = deepcopy(A)
    c.evolve(mode='G', target=6)
    c.draw(root, .35, .52)

    for x in (a, b, c):
        root.text(.64, x.y, "Score={0}".format(x.nonwhites), va="center")

    # Panel C
    A.truncate_between_flankers()
    a.truncate_between_flankers()
    b.truncate_between_flankers()
    c.truncate_between_flankers(target=6)

    plot_diagram(root, .14, .2, A, a, "S", "syntenic")
    plot_diagram(root, .37, .2, A, b, "F", "missing, with both flankers")
    plot_diagram(root, .6, .2, A, c, "G", "missing, with one flanker")

    labels = ((.04, .95, 'A'), (.04, .75, 'B'), (.04, .4, 'C'))
    panel_labels(root, labels)

    # Descriptions
    xt = .85
    desc = ("Extract neighborhood",
            "of *window* size",
            "Count gene pairs within *window*",
            "Find regions above *score* cutoff",
            "Identify flankers",
            "Annotate syntelog class"
            )
    for yt, t in zip((.88, .84, .64, .6, .3, .26), desc):
        root.text(xt, yt, markup(t), ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
