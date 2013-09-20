#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the various past manuscripts.
"""

import sys
import logging

import numpy as np

from jcvi.graphics.base import plt, _, Rectangle, Polygon, CirclePolygon, savefig
from jcvi.graphics.glyph import GeneGlyph, RoundLabel, RoundRect, \
        arrowprops, TextCircle, plot_cap
from jcvi.graphics.chromosome import Chromosome
from jcvi.utils.iter import pairwise
from jcvi.apps.base import OptionParser, ActionDispatcher, fname, debug
debug()


def main():

    actions = (
        # Brapa bites paper (Tang et al., 2012 Genetics)
        ('excision', 'show intra-chromosomal recombination'),
        ('bites', 'show the bites calling pipeline'),
        ('scenario', 'show step-wise genome merger events in brapa'),
        # Epoch paper (Woodhouse et al., 2012 Plant Cell)
        ('epoch', 'show the methods used in epoch paper'),
        # Cotton paper (Paterson et al., 2012 Nature)
        ('cotton', 'plot cotton macro- and micro-synteny (requires data)'),
        # Unpublished
        ('amborella', 'plot amborella macro- and micro-synteny (requires data)'),
        ('napusretention', 'plot retention rate along chr (requires data)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def napusretention(args):
    """
    %prog napusretention final.index scaffold.bed genes.bed chrX

    Plot retention rate along specified chromosome. File `final.index` contains
    three column with first (query), second (ortholog), third (homeolog).
    """
    from jcvi.algorithms.matrix import moving_sum
    from jcvi.formats.bed import Bed
    from jcvi.graphics.chromosome import ChromosomeMap

    p = OptionParser(napusretention.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="6x4")

    if len(args) != 4:
        sys.exit(not p.print_help())

    finalist, scaffoldbed, genesbed, chr = args
    bed = Bed(scaffoldbed)
    scaffolds = list(bed.sub_bed(chr))
    logging.debug("Imported {0} scaffolds from {1}.".format(len(scaffolds), chr))

    bed = Bed(genesbed)
    genes = list(x for x in bed.sub_bed(chr))
    ngenes = len(genes)
    logging.debug("Imported {0} genes from {1}.".format(ngenes, chr))

    ortho, homeo = {}, {}
    fp = open(finalist)
    for row in fp:
        query, ortholog, homeolog = row.split()
        ortho[query] = int(ortholog != '.')
        homeo[query] = int(homeolog != '.')

    genes = [x for x in genes if x.accn in ortho]
    ortholist, homeolist = [], []
    for x in genes:
        accn = x.accn
        ortholist.append(ortho[accn])
        homeolist.append(homeo[accn])

    window = 100
    halfw = window / 2
    ortholist = moving_sum(ortholist, window=window)[halfw : -halfw]
    homeolist = moving_sum(homeolist, window=window)[halfw : -halfw]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    xstart, xend = .2, .8
    ystart, yend = .2, .8
    pad = .02

    chrlen = max(x.end for x in scaffolds)
    ratio = (xend - xstart) / chrlen
    patchstart = [(xstart + x.start * ratio) for x in scaffolds]

    title = "Homolog retention in {0}-gene windows".format(window)
    subtitle = "{0}, {1} scaffolds, {2} genes".format(chr, len(scaffolds), ngenes)

    chr_map = ChromosomeMap(fig, root, xstart, xend, ystart, yend, pad, 0, 100, 2,
                    title, subtitle, patchstart=patchstart)

    # Legends
    y = ystart + pad * 4
    root.plot([.3, .33], [y, y], "r-", lw=2)
    root.text(.33 + pad, y, "Ortholog", va="center")
    y = ystart + pad * 2
    root.plot([.3, .33], [y, y], "g-", lw=2)
    root.text(.33 + pad, y, "Homeolog", va="center")

    ax = chr_map.axes
    starts = [x.start for x in genes][halfw : -halfw]

    ax.plot(starts, ortholist, "r-")
    ax.plot(starts, homeolist, "g-")

    ax.set_xlim(0, chrlen)
    ax.set_ylim(0, 100)
    ax.set_axis_off()

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + ".retention." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def amborella(args):
    """
    %prog amborella seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a composite figure that calls graphics.karyotype and graphic.synteny.
    """
    from jcvi.graphics.karyotype import Karyotype
    from jcvi.graphics.synteny import Synteny, draw_gene_legend

    p = OptionParser(amborella.__doc__)
    p.add_option("--tree",
                 help="Display trees on the bottom of the figure [default: %default]")
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 5:
        sys.exit(p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args
    switch = opts.switch
    tree = opts.tree

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)
    Synteny(fig, root, datafile, bedfile, slayout, switch=switch, tree=tree)

    # legend showing the orientation of the genes
    draw_gene_legend(root, .5, .68, .5)

    # annotate the WGD events
    fc = 'lightslategrey'
    x = .05
    radius = .012
    TextCircle(root, x, .86, '$\gamma$', radius=radius)
    TextCircle(root, x, .95, '$\epsilon$', radius=radius)
    root.plot([x, x], [.83, .9], ":", color=fc, lw=2)
    pts = plot_cap((x, .95), np.radians(range(-70, 250)), .02)
    x, y = zip(*pts)
    root.plot(x, y, ":", color=fc, lw=2)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "amborella"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def cotton(args):
    """
    %prog cotton seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a composite figure that calls graphics.karyotype and graphic.synteny.
    """
    from jcvi.graphics.karyotype import Karyotype
    from jcvi.graphics.synteny import Synteny, draw_gene_legend
    from jcvi.graphics.tree import draw_tree, read_trees

    p = OptionParser(cotton.__doc__)
    p.add_option("--depthfile",
                 help="Use depth info in this file [default: %default]")
    p.add_option("--tree",
                 help="Display trees on the bottom of the figure [default: %default]")
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 5:
        sys.exit(p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args
    switch = opts.switch
    tree = opts.tree
    depthfile = opts.depthfile

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    kt = Karyotype(fig, root, seqidsfile, klayout)
    st = Synteny(fig, root, datafile, bedfile, slayout, switch=switch)

    light = "lightslategrey"
    # Show the dup depth along the cotton chromosomes
    if depthfile:
        ymin, ymax = .9, .95
        root.text(.11, .96, "Cotton duplication level", color="gray", size=10)
        root.plot([.1, .95], [ymin, ymin], color="gray")
        root.text(.96, .9, "1x", color="gray", va="center")
        root.plot([.1, .95], [ymax, ymax], color="gray")
        root.text(.96, .95, "6x", color="gray", va="center")

        fp = open(depthfile)
        track = kt.tracks[0]  # Cotton
        depths = []
        for row in fp:
            a, b, depth = row.split()
            depth = int(depth)
            try:
                p = track.get_coords(a)
                depths.append((p, depth))
            except KeyError:
                pass

        depths.sort(key=lambda x: (x[0], -x[1]))
        xx, yy = zip(*depths)
        yy = [ymin + .01 * (x - 1) for x in yy]
        root.plot(xx, yy, "-", color=light)

    # legend showing the orientation of the genes
    draw_gene_legend(root, .5, .68, .5)

    # Zoom
    xpos = .835
    xmin, xmax = .18, .82
    ymin, ymax = ytop, .55
    lc = "k"
    kwargs = dict(lw=3, color=lc, mec=lc, mfc="w", zorder=3)
    root.plot((xpos, xpos), (ymax, .63), ":o", **kwargs)
    root.plot((xpos, xmin), (ymax, ymin), ":o", **kwargs)
    root.plot((xpos, xmax), (ymax, ymin), ":o", **kwargs)
    RoundRect(root, (.06, .17), .92, .35, fill=False, lw=2, ec=light)

    # Panels
    root.text(.05, .95, "a", size=20, fontweight="bold")
    root.text(.1, .45, "b", size=20, fontweight="bold")

    if tree:
        panel3 = fig.add_axes([.65, .05, .35, .35])
        trees = read_trees(tree)
        label, outgroup, tx = trees[0]
        draw_tree(panel3, tx, outgroup=outgroup, rmargin=.4, leafcolor="r")
        panel3.set_xlim(0, 1)
        panel3.set_ylim(0, 1)
        panel3.set_axis_off()

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cotton"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_diagram(ax, x, y, label="S", title="syntenic"):
    """
    Part of the diagrams that are re-used. (x, y) marks the center of the
    diagram. Label determines the modification to the "S" graph.
    """
    trackgap = .06
    tracklen = .12
    xa, xb = x - tracklen, x + tracklen
    ya, yb = y + trackgap, y - trackgap
    hsps = (((60, 150), (50, 130)),
           ((190, 225), (200, 240)),
           ((330, 280), (360, 310)))

    for yy in (ya, yb):
        ax.plot((xa, xb), (yy, yy), "-", color="gray", lw=2, zorder=1)

    ytip = .015
    mrange = 400
    m = lambda t: xa + t * 1. / mrange * tracklen * 2

    for i, ((a, b), (c, d)) in enumerate(hsps):
        fb = False
        if label == "FB" and i == 1:
            c, d = 270, 280
            fb = True
        if label == "G" and i == 0:
            c, d = 120, 65

        a, b, c, d = [m(t) for t in (a, b, c, d)]
        color = "g" if i == 1 else "r"
        GeneGlyph(ax, a, b, ya, 2 * ytip, fc=color)

        if i == 1 and label in ("F", "G", "FN"):
            pass
        else:
            if fb:
                GeneGlyph(ax, c, d, yb, 2 * ytip, fc='w', tip=0)
            else:
                GeneGlyph(ax, c, d, yb, 2 * ytip, fc='r')

        r = Polygon(((a, ya - ytip), (c, yb + ytip),
                      (d, yb + ytip), (b, ya - ytip)),
                      fc='r', alpha=.2)

        if i == 1 and label not in ("S", "FB"):
            pass
        elif i == 0 and label == "G":
            pass
        else:
            ax.add_patch(r)

    if label == "FN":
        ax.text(x + .005, yb, _("NNNNN"), ha="center", size=7)

    title = "{0}: {1}".format(label, title)
    ax.text(x, ya + 5 * ytip, _(title), size=8, ha="center")


def epoch(args):
    """
    %prog epoch

    Illustrate the methods used in Maggie's epoch paper, in particular, how to
    classifiy S/G/F/FB/FN for the genes.
    """
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    fig = plt.figure(1, (6, 4))
    root = fig.add_axes([0, 0, 1, 1])

    # Separators
    linestyle = dict(lw=2, color="b", alpha=.2, zorder=2)
    root.plot((0, 1), (.5, .5), "--", **linestyle)
    for i in (1./3, 2./3):
        root.plot((i, i), (.5, 1), "--", **linestyle)
    for i in (1./6, 3./6, 5./6):
        root.plot((i, i), (0, .5), "--", **linestyle)

    # Diagrams
    plot_diagram(root, 1./6, 3./4, "S", "syntenic")
    plot_diagram(root, 3./6, 3./4, "F", "missing, with both flankers")
    plot_diagram(root, 5./6, 3./4, "G", "missing, with one flanker")
    plot_diagram(root, 2./6, 1./4, "FB", "has non-coding matches")
    plot_diagram(root, 4./6, 1./4, "FN", "syntenic region has gap")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


def excision(args):
    """
    %prog excision

    Illustrate the mechanism of illegitimate recombination.
    """
    p = OptionParser(__doc__)
    opts, args = p.parse_args(args)

    fig = plt.figure(1, (5, 5))
    root = fig.add_axes([0, 0, 1, 1])

    plt.plot((.2, .8), (.6, .6), 'r-', lw=3)
    plt.plot((.4, .6), (.6, .6), 'b>-', mfc='g', mec='w', ms=12, lw=3)
    plt.plot((.3, .7), (.5, .5), 'r-', lw=3)
    plt.plot((.5, ), (.5, ), 'b>-', mfc='g', mec='w', ms=12, lw=3)

    # Circle excision
    plt.plot((.5, ), (.45, ), 'b>-', mfc='g', mec='w', ms=12, lw=3)
    circle = CirclePolygon((.5, .4), .05, fill=False, lw=3, ec="b")
    root.add_patch(circle)

    arrow_dist = .07
    ar_xpos, ar_ypos = .5, .52
    root.annotate(" ", (ar_xpos, ar_ypos),
            (ar_xpos, ar_ypos + arrow_dist),
            arrowprops=arrowprops)

    RoundLabel(root, .2, .64, "Gene")
    RoundLabel(root, .3, .54, "Excision")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


def bites(args):
    """
    %prog bites

    Illustrate the pipeline for automated bite discovery.
    """

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    fig = plt.figure(1, (6, 6))
    root = fig.add_axes([0, 0, 1, 1])

    # HSP pairs
    hsps = (((50, 150), (60, 180)),
           ((190, 250), (160, 235)),
           ((300, 360), (270, 330)),
           ((430, 470), (450, 490)),
           ((570, 620), (493, 543)),
           ((540, 555), (370, 385)),  # non-collinear hsps
          )

    titlepos = (.9, .65, .4)
    titles = ("Compare orthologous region",
              "Find collinear HSPs",
              "Scan paired gaps")
    ytip = .01
    mrange = 650.
    m = lambda x: x / mrange * .7 + .1
    for i, (ya, title) in enumerate(zip(titlepos, titles)):
        yb = ya - .1
        plt.plot((.1, .8), (ya, ya), "-", color="gray", lw=2, zorder=1)
        plt.plot((.1, .8), (yb, yb), "-", color="gray", lw=2, zorder=1)
        RoundLabel(root, .5, ya + 4 * ytip, title)
        root.text(.9, ya, _("A. thaliana"), ha="center", va="center")
        root.text(.9, yb, _("B. rapa"), ha="center", va="center")
        myhsps = hsps
        if i >= 1:
            myhsps = hsps[:-1]
        for (a, b), (c, d) in myhsps:
            a, b, c, d = [m(x) for x in (a, b, c, d)]
            r1 = Rectangle((a, ya - ytip), b - a, 2 * ytip, fc='r', lw=0, zorder=2)
            r2 = Rectangle((c, yb - ytip), d - c, 2 * ytip, fc='r', lw=0, zorder=2)
            r3 = Rectangle((a, ya - ytip), b - a, 2 * ytip, fill=False, zorder=3)
            r4 = Rectangle((c, yb - ytip), d - c, 2 * ytip, fill=False, zorder=3)
            r5 = Polygon(((a, ya - ytip), (c, yb + ytip),
                          (d, yb + ytip), (b, ya - ytip)),
                          fc='r', alpha=.2)
            rr = (r1, r2, r3, r4, r5)
            if i == 2:
                rr = rr[:-1]
            for r in rr:
                root.add_patch(r)

    # Gap pairs
    hspa, hspb = zip(*myhsps)
    gapa, gapb = [], []
    for (a, b), (c, d) in pairwise(hspa):
        gapa.append((b + 1, c - 1))
    for (a, b), (c, d) in pairwise(hspb):
        gapb.append((b + 1, c - 1))
    gaps = zip(gapa, gapb)
    tpos = titlepos[-1]

    yy = tpos - .05
    for i, ((a, b), (c, d)) in enumerate(gaps):
        i += 1
        a, b, c, d = [m(x) for x in (a, b, c, d)]
        xx = (a + b + c + d) / 4
        TextCircle(root, xx, yy, _(str(i)))

    # Bites
    ystart = .24
    ytip = .05
    bites = (("Bite(40=>-15)", True),
             ("Bite(50=>35)", False),
             ("Bite(70=>120)", False),
             ("Bite(100=>3)", True))
    for i, (bite, selected) in enumerate(bites):
        xx = .15 if (i % 2 == 0) else .55
        yy = ystart - i / 2 * ytip
        i += 1
        TextCircle(root, xx, yy, _(str(i)))
        color = "k" if selected else "gray"
        root.text(xx + ytip, yy, bite, size=10, color=color, va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


def scenario(args):
    """
    %prog scenario

    Illustration of the two-step genome merger process for B. rapa companion paper.
    """
    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    fig = plt.figure(1, (5, 5))
    root = fig.add_axes([0, 0, 1, 1])

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    # Layout format: (x, y, label, (chr lengths))
    anc = (.5, .9, "Ancestor", (1,))
    s1 = (.2, .6, "Genome I", (1,))
    s2 = (.5, .6, "Genome II", (1,))
    s3 = (.8, .6, "Genome III", (1,))
    tetra = (.35, .4, "Tetraploid I / II", (.5, .9))
    hexa = (.5, .1, "Hexaploid I / II / III", (.36, .46, .9))
    labels = (anc, s1, s2, s3, tetra, hexa)
    connections = ((anc, s1), (anc, s2), (anc, s3),\
            (s1, tetra), (s2, tetra),
            (tetra, hexa), (s3, hexa))

    xinterval = .02
    yratio = .05
    for xx, yy, label, chrl in labels:
        #RoundLabel(root, xx, yy, label)
        root.text(xx, yy, _(label), ha="center", va="center")
        offset = len(label) * .012
        for i, c in enumerate(chrl):
            ya = yy + yratio * c
            yb = yy - yratio * c
            Chromosome(root, xx - offset + i * xinterval, ya, yb, width=.01)

    # Comments
    comments = ((.15, .33, "II dominant"),
                (.25, .03, "III dominant"))

    for xx, yy, c in comments:
        root.text(xx, yy, _(c), size=9, ha="center", va="center")

    # Branches
    tip = .04
    for a, b in connections:
        xa, ya, la, chra = a
        xb, yb, lb, chrb = b
        plt.plot((xa, xb), (ya - tip, yb + 2 * tip), 'k-', lw=2, alpha=.5)

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


if __name__ == '__main__':
    main()
