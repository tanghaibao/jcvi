#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in various manuscripts.
"""

import os.path as op
import sys

import numpy as np

from jcvi.graphics.base import plt, Polygon, savefig
from jcvi.graphics.glyph import GeneGlyph, RoundRect, TextCircle, DoubleSquare, plot_cap
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
from jcvi.apps.base import OptionParser, ActionDispatcher, fname


def main():

    actions = (
        # Epoch paper (Woodhouse et al., 2012 Plant Cell)
        ('epoch', 'show the methods used in epoch paper'),
        # Cotton paper (Paterson et al., 2012 Nature)
        ('cotton', 'plot cotton macro- and micro-synteny (requires data)'),
        # Amborella paper (Albert et al., 2013 Science)
        ('amborella', 'plot amborella macro- and micro-synteny (requires data)'),
        # Mt4.0 paper (Tang et al., 2014 BMC Genomics)
        ('mtdotplots', 'plot Mt3.5 and Mt4.0 side-by-side'),
        # Unpublished
        ('litchi', 'plot litchi micro-synteny (requires data)'),
        ('birch', 'plot birch macro-synteny (requires data)'),
        ('oropetium', 'plot oropetium micro-synteny (requires data)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def join_nodes(root, coords, a, b, x, slope=2.4,
               fc="lightslategray", rectangle=True, circle=True):
    # Join node a and b to make an internal node
    ax, ay = coords[a]
    bx, by = coords[b]
    if ay < by:
        ax, ay, bx, by = bx, by, ax, ay
    if rectangle:
        nx, ny = x, (ay + by) / 2
        root.plot((nx, ax), (ay, ay), lw=2, color=fc)
        root.plot((nx, bx), (by, by), lw=2, color=fc)
        root.plot((nx, nx), (ay, by), lw=2, color=fc)
    else:
        dx = (abs(ay - by) / slope - abs(ax - bx)) / 2
        nx = max(ax, bx) + dx
        ny = by + (nx - bx) * slope
        root.plot((nx, ax), (ny, ay), lw=2, color=fc)
        root.plot((nx, bx), (ny, by), lw=2, color=fc)
    if circle:
        DoubleSquare(root, nx, ny, fc=fc)
    return nx, ny


def birch(args):
    """
    %prog birch seqids layout

    Plot birch macro-synteny, with an embedded phylogenetic tree to the right.
    """
    p = OptionParser(birch.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x6")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqids, layout = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    K = Karyotype(fig, root, seqids, layout)
    L = K.layout

    xs = .79
    dt = dict(rectangle=False, circle=False)
    # Embed a phylogenetic tree to the right
    coords = {}
    coords["Amborella"] = (xs, L[0].y)
    coords["Vitis"] = (xs, L[1].y)
    coords["Prunus"] = (xs, L[2].y)
    coords["Betula"] = (xs, L[3].y)
    coords["Populus"] = (xs, L[4].y)
    coords["Arabidopsis"] = (xs, L[5].y)
    coords["fabids"] = join_nodes(root, coords, "Prunus", "Betula", xs, **dt)
    coords["malvids"] = join_nodes(root, coords, \
                                   "Populus", "Arabidopsis", xs, **dt)
    coords["rosids"] = join_nodes(root, coords, "fabids", "malvids", xs, **dt)
    coords["eudicots"] = join_nodes(root, coords, "rosids", "Vitis", xs, **dt)
    coords["angiosperm"] = join_nodes(root, coords, \
                                      "eudicots", "Amborella", xs, **dt)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "birch"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def mtdotplots(args):
    """
    %prog mtdotplots Mt3.5 Mt4.0 medicago.medicago.lifted.1x1.anchors

    Plot Mt3.5 and Mt4.0 side-by-side. This is essentially combined from two
    graphics.dotplot() function calls as panel A and B.
    """
    from jcvi.graphics.dotplot import check_beds, dotplot

    p = OptionParser(mtdotplots.__doc__)
    p.set_beds()
    opts, args, iopts = p.set_image_options(args, figsize="16x8", dpi=90)

    if len(args) != 3:
        sys.exit(not p.print_help())

    a, b, ac = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    r1 = fig.add_axes([0, 0, .5, 1])
    r2 = fig.add_axes([.5, 0, .5, 1])
    a1 = fig.add_axes([.05, .1, .4, .8])
    a2 = fig.add_axes([.55, .1, .4, .8])

    anchorfile = op.join(a, ac)
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    dotplot(anchorfile, qbed, sbed, fig, r1, a1, is_self=is_self,
            genomenames="Mt3.5_Mt3.5")

    opts.qbed = opts.sbed = None
    anchorfile = op.join(b, ac)
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    dotplot(anchorfile, qbed, sbed, fig, r2, a2, is_self=is_self,
            genomenames="Mt4.0_Mt4.0")

    root.text(.03, .95, "A", ha="center", va="center", size=36)
    root.text(.53, .95, "B", ha="center", va="center", size=36)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "mtdotplots"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def oropetium(args):
    """
    %prog oropetium mcscan.out all.bed layout switch.ids

    Build a composite figure that calls graphis.synteny.
    """
    p = OptionParser(oropetium.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    if len(args) != 4:
        sys.exit(not p.print_help())

    datafile, bedfile, slayout, switch = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Synteny(fig, root, datafile, bedfile, slayout, switch=switch)

    # legend showing the orientation of the genes
    draw_gene_legend(root, .4, .57, .74)

    # On the left panel, make a species tree
    fc = 'lightslategrey'

    coords = {}
    xs, xp = .16, .03
    coords["oropetium"] = (xs, .7)
    coords["setaria"] = (xs, .6)
    coords["sorghum"] = (xs, .5)
    coords["rice"] = (xs, .4)
    coords["brachypodium"] = (xs, .3)
    xs -= xp
    coords["Panicoideae"] = join_nodes(root, coords, "setaria", "sorghum", xs)
    xs -= xp
    coords["BEP"] = join_nodes(root, coords, "rice", "brachypodium", xs)
    coords["PACMAD"] = join_nodes(root, coords, "oropetium", "Panicoideae", xs)
    xs -= xp
    coords["Poaceae"] = join_nodes(root, coords, "BEP", "PACMAD", xs)

    # Names of the internal nodes
    for tag in ("BEP", "Poaceae"):
        nx, ny = coords[tag]
        nx, ny = nx - .005, ny - .02
        root.text(nx, ny, tag, rotation=90, ha="right", va="top", color=fc)
    for tag in ("PACMAD",):
        nx, ny = coords[tag]
        nx, ny = nx - .005, ny + .02
        root.text(nx, ny, tag, rotation=90, ha="right", va="bottom", color=fc)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "oropetium"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def litchi(args):
    """
    %prog litchi mcscan.out all.bed layout switch.ids

    Build a composite figure that calls graphis.synteny.
    """
    p = OptionParser(litchi.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    if len(args) != 4:
        sys.exit(not p.print_help())

    datafile, bedfile, slayout, switch = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Synteny(fig, root, datafile, bedfile, slayout, switch=switch)

    # legend showing the orientation of the genes
    draw_gene_legend(root, .4, .7, .82)

    # On the left panel, make a species tree
    fc = 'lightslategrey'

    coords = {}
    xs, xp = .16, .03
    coords["lychee"] = (xs, .37)
    coords["clementine"] = (xs, .5)
    coords["cacao"] = (xs, .6)
    coords["strawberry"] = (xs, .7)
    coords["grape"] = (xs, .8)
    xs -= xp
    coords["Sapindales"] = join_nodes(root, coords, "clementine", "lychee", xs)
    xs -= xp
    coords["Rosid-II"] = join_nodes(root, coords, "cacao", "Sapindales", xs)
    xs -= xp
    coords["Rosid"] = join_nodes(root, coords, "strawberry", "Rosid-II", xs)
    xs -= xp
    coords["crown"] = join_nodes(root, coords, "grape", "Rosid", xs,
                                 circle=False)

    # Names of the internal nodes
    for tag in ("Rosid", "Rosid-II", "Sapindales"):
        nx, ny = coords[tag]
        nx, ny = nx - .01, ny - .02
        root.text(nx, ny, tag, rotation=90, ha="right", va="top", color=fc)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "litchi"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def amborella(args):
    """
    %prog amborella seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a composite figure that calls graphics.karyotype and graphics.synteny.
    """
    p = OptionParser(amborella.__doc__)
    p.add_option("--tree",
                 help="Display trees on the bottom of the figure [default: %default]")
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 5:
        sys.exit(not p.print_help())

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
    p = OptionParser(cotton.__doc__)
    p.add_option("--depthfile",
                 help="Use depth info in this file [default: %default]")
    p.add_option("--switch",
                 help="Rename the seqid with two-column file [default: %default]")
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 5:
        sys.exit(p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args
    switch = opts.switch
    depthfile = opts.depthfile

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    kt = Karyotype(fig, root, seqidsfile, klayout)
    Synteny(fig, root, datafile, bedfile, slayout, switch=switch)

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
    ytop = .9
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
        ax.text(x + .005, yb, "NNNNN", ha="center", size=7)

    title = "{0}: {1}".format(label, title)
    ax.text(x, ya + 5 * ytip, title, size=8, ha="center")


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


if __name__ == '__main__':
    main()
