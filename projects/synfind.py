#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""

import sys
import logging

from copy import deepcopy
from collections import defaultdict
from itertools import groupby

from jcvi.formats.base import must_open
from jcvi.utils.cbook import SummaryStats
from jcvi.formats.bed import Bed
from jcvi.graphics.base import FancyArrow, plt, savefig, panel_labels, markup
from jcvi.graphics.glyph import CartoonRegion, RoundRect
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('cartoon', 'generate cartoon illustration of SynFind'),
        ('islands', 'gene presence absence analysis in ecoli'),
        ('grasses', 'validate SynFind pan-grass set against James'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
            s |= set(m)

        for x in Os:
            james_store[x] = s

    jaccards = []
    for k, v in james_store.items():
        if k not in master_store:
            continue
        m = master_store[k]
        jaccard = len(v & m) * 100 / len(v | m)
        jaccards.append(jaccard)
        if opts.verbose:
            print k
            print v
            print m
            print jaccard

    j = SummaryStats(jaccards)
    print j


def islands(args):
    """
    %prog islands coge_master_table.txt query.bed

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
    p = OptionParser(islands.__doc__)
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
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

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
            "Generate anchors",
            "Group anchors within *window* apart",
            "Find regions above *score* cutoff",
            "Identify flankers",
            "Annotate syntelog class"
            )
    for yt, t in zip((.88, .84, .64, .6, .56, .3, .26), desc):
        root.text(xt, yt, markup(t), ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
