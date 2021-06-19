#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in various manuscripts.
"""

import os.path as op
import sys
import logging

import numpy as np

from jcvi.graphics.base import Polygon, normalize_axes, panel_labels, plt, savefig
from jcvi.graphics.glyph import GeneGlyph, RoundRect, TextCircle, DoubleSquare, plot_cap
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
from jcvi.apps.base import OptionParser, ActionDispatcher, fname


def main():

    actions = (
        # Epoch paper (Woodhouse et al., 2012 Plant Cell)
        ("epoch", "show the methods used in epoch paper"),
        # Cotton paper (Paterson et al., 2012 Nature)
        ("cotton", "plot cotton macro- and micro-synteny (requires data)"),
        # Amborella paper (Albert et al., 2013 Science)
        ("amborella", "plot amborella macro- and micro-synteny (requires data)"),
        # Mt4.0 paper (Tang et al., 2014 BMC Genomics)
        ("mtdotplots", "plot Mt3.5 and Mt4.0 side-by-side"),
        # Oropetium paper (Vanburen et al., 2015 Nature)
        ("oropetium", "plot oropetium micro-synteny (requires data)"),
        # Pomegranate paper (Qin et al., 2017 Plant Journal)
        ("pomegranate", "plot pomegranate macro- and micro-synteny (requires data)"),
        # Unpublished
        ("birch", "plot birch macro-synteny (requires data)"),
        ("litchi", "plot litchi micro-synteny (requires data)"),
        ("utricularia", "plot utricularia micro-synteny (requires data)"),
        (
            "waterlilyGOM",
            "waterlily phylogeny and related infographics (requires data)",
        ),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def waterlilyGOM(args):
    """
    %prog mcmctree.tre table.csv

    Customized figure to plot phylogeny and related infographics.
    """
    from jcvi.graphics.tree import (
        LeafInfoFile,
        WGDInfoFile,
        draw_tree,
        parse_tree,
        draw_wgd_xy,
    )
    from jcvi.graphics.table import CsvTable, draw_table

    p = OptionParser(waterlilyGOM.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x9")

    if len(args) != 2:
        sys.exit(not p.print_help())

    (datafile, csvfile) = args
    outgroup = ["ginkgo"]

    logging.debug("Load tree file `{0}`".format(datafile))
    t, hpd = parse_tree(datafile)

    pf = datafile.rsplit(".", 1)[0]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    margin, rmargin = 0.15, 0.19  # Left and right margin
    leafinfo = LeafInfoFile("leafinfo.csv").cache
    wgdinfo = WGDInfoFile("wgdinfo.csv").cache
    groups = "Monocots,Eudicots,ANA-grade,Gymnosperms"

    draw_tree(
        root,
        t,
        hpd=hpd,
        margin=margin,
        rmargin=rmargin,
        supportcolor=None,
        internal=False,
        outgroup=outgroup,
        leafinfo=leafinfo,
        wgdinfo=wgdinfo,
        geoscale=True,
        groups=groups.split(","),
    )

    # Bottom right show legends for the WGD circles
    pad = 0.02
    ypad = 0.04
    xstart = 1 - rmargin + pad
    ystart = 0.2
    waterlily_wgdline = wgdinfo["waterlily"][0]
    ypos = ystart - 2 * ypad
    draw_wgd_xy(root, xstart, ypos, waterlily_wgdline)
    root.text(
        xstart + pad,
        ypos,
        "Nymphaealean WGD",
        color=waterlily_wgdline.color,
        va="center",
    )
    other_wgdline = wgdinfo["banana"][0]
    ypos = ystart - 3 * ypad
    draw_wgd_xy(root, xstart, ypos, other_wgdline)
    root.text(
        xstart + pad,
        ypos,
        "Other known WGDs",
        color=other_wgdline.color,
        va="center",
    )

    # Top left draw the comparison table
    csv_table = CsvTable(csvfile)
    draw_table(
        root,
        csv_table,
        extent=(0.02, 0.44, 0.55, 0.985),
        stripe_color="lavender",
        yinflation=iopts.w / iopts.h,
    )

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def pomegranate(args):
    """
    %prog cotton seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of WGD history of pineapple genome. The script calls both graphics.karyotype
    and graphic.synteny.
    """
    p = OptionParser(pomegranate.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x7")

    if len(args) != 5:
        sys.exit(not p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)
    Synteny(fig, root, datafile, bedfile, slayout)

    # legend showing the orientation of the genes
    draw_gene_legend(root, 0.42, 0.52, 0.48)

    labels = ((0.04, 0.96, "A"), (0.04, 0.52, "B"))
    panel_labels(root, labels)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "pomegranate-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def utricularia(args):
    from jcvi.graphics.synteny import main as synteny_main

    p = OptionParser(synteny_main.__doc__)
    p.add_option("--switch", help="Rename the seqid with two-column file")
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, bedfile, layoutfile = args
    switch = opts.switch

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    s = Synteny(
        fig, root, datafile, bedfile, layoutfile, loc_label=False, switch=switch
    )
    light = "lightslategrey"
    RoundRect(root, (0.02, 0.69), 0.96, 0.24, fill=False, lw=2, ec=light)
    RoundRect(root, (0.02, 0.09), 0.96, 0.48, fill=False, lw=2, ec=light)
    za, zb = s.layout[1].ratio, s.layout[-1].ratio  # zoom level
    if za != 1:
        root.text(
            0.96,
            0.89,
            "{}x zoom".format(za).replace(".0x", "x"),
            color=light,
            ha="right",
            va="center",
            size=14,
        )
    if zb != 1:
        root.text(
            0.96,
            0.12,
            "{}x zoom".format(zb).replace(".0x", "x"),
            color=light,
            ha="right",
            va="center",
            size=14,
        )

    # legend showing the orientation of the genes
    draw_gene_legend(root, 0.22, 0.3, 0.64, text=True)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def join_nodes(
    root, coords, a, b, x, slope=2.4, fc="lightslategray", rectangle=True, circle=True
):
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


def branch_length(ax, start, end, text, ha="left", va="bottom", color="r"):
    xs, ys = start
    xe, ye = end
    text = r"$\mathsf{" + text + "}$"
    ax.text((xs + xe) / 2, (ys + ye) / 2, text, ha=ha, va=va, color=color)


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

    xs = 0.79
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
    coords["malvids"] = join_nodes(root, coords, "Populus", "Arabidopsis", xs, **dt)
    coords["rosids"] = join_nodes(root, coords, "fabids", "malvids", xs, **dt)
    coords["eudicots"] = join_nodes(root, coords, "rosids", "Vitis", xs, **dt)
    coords["angiosperm"] = join_nodes(root, coords, "eudicots", "Amborella", xs, **dt)

    # Show branch length
    branch_length(root, coords["Amborella"], coords["angiosperm"], ">160.0")
    branch_length(root, coords["eudicots"], coords["angiosperm"], ">78.2", va="top")
    branch_length(root, coords["Vitis"], coords["eudicots"], "138.5")
    branch_length(root, coords["rosids"], coords["eudicots"], "19.8", va="top")
    branch_length(
        root, coords["Prunus"], coords["fabids"], "104.2", ha="right", va="top"
    )
    branch_length(root, coords["Arabidopsis"], coords["malvids"], "110.2", va="top")
    branch_length(
        root, coords["fabids"], coords["rosids"], "19.8", ha="right", va="top"
    )
    branch_length(root, coords["malvids"], coords["rosids"], "8.5", va="top")

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
    r1 = fig.add_axes([0, 0, 0.5, 1])
    r2 = fig.add_axes([0.5, 0, 0.5, 1])
    a1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])
    a2 = fig.add_axes([0.55, 0.1, 0.4, 0.8])

    anchorfile = op.join(a, ac)
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    dotplot(
        anchorfile, qbed, sbed, fig, r1, a1, is_self=is_self, genomenames="Mt3.5_Mt3.5"
    )

    opts.qbed = opts.sbed = None
    anchorfile = op.join(b, ac)
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    dotplot(
        anchorfile, qbed, sbed, fig, r2, a2, is_self=is_self, genomenames="Mt4.0_Mt4.0"
    )

    root.text(0.03, 0.95, "A", ha="center", va="center", size=36)
    root.text(0.53, 0.95, "B", ha="center", va="center", size=36)

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
    p.add_option("--extra", help="Extra features in BED format")
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    if len(args) != 4:
        sys.exit(not p.print_help())

    datafile, bedfile, slayout, switch = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Synteny(
        fig, root, datafile, bedfile, slayout, switch=switch, extra_features=opts.extra
    )

    # legend showing the orientation of the genes
    draw_gene_legend(root, 0.4, 0.57, 0.74, text=True, repeat=True)

    # On the left panel, make a species tree
    fc = "lightslategrey"

    coords = {}
    xs, xp = 0.16, 0.03
    coords["oropetium"] = (xs, 0.7)
    coords["setaria"] = (xs, 0.6)
    coords["sorghum"] = (xs, 0.5)
    coords["rice"] = (xs, 0.4)
    coords["brachypodium"] = (xs, 0.3)
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
        nx, ny = nx - 0.005, ny - 0.02
        root.text(nx, ny, tag, rotation=90, ha="right", va="top", color=fc)
    for tag in ("PACMAD",):
        nx, ny = coords[tag]
        nx, ny = nx - 0.005, ny + 0.02
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
    draw_gene_legend(root, 0.4, 0.7, 0.82)

    # On the left panel, make a species tree
    fc = "lightslategrey"

    coords = {}
    xs, xp = 0.16, 0.03
    coords["lychee"] = (xs, 0.37)
    coords["clementine"] = (xs, 0.5)
    coords["cacao"] = (xs, 0.6)
    coords["strawberry"] = (xs, 0.7)
    coords["grape"] = (xs, 0.8)
    xs -= xp
    coords["Sapindales"] = join_nodes(root, coords, "clementine", "lychee", xs)
    xs -= xp
    coords["Rosid-II"] = join_nodes(root, coords, "cacao", "Sapindales", xs)
    xs -= xp
    coords["Rosid"] = join_nodes(root, coords, "strawberry", "Rosid-II", xs)
    xs -= xp
    coords["crown"] = join_nodes(root, coords, "grape", "Rosid", xs, circle=False)

    # Names of the internal nodes
    for tag in ("Rosid", "Rosid-II", "Sapindales"):
        nx, ny = coords[tag]
        nx, ny = nx - 0.01, ny - 0.02
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
    p.add_option("--tree", help="Display trees on the bottom of the figure")
    p.add_option("--switch", help="Rename the seqid with two-column file")
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
    draw_gene_legend(root, 0.5, 0.68, 0.5)

    # annotate the WGD events
    fc = "lightslategrey"
    x = 0.05
    radius = 0.012
    TextCircle(root, x, 0.86, r"$\gamma$", radius=radius)
    TextCircle(root, x, 0.95, r"$\epsilon$", radius=radius)
    root.plot([x, x], [0.83, 0.9], ":", color=fc, lw=2)
    pts = plot_cap((x, 0.95), np.radians(range(-70, 250)), 0.02)
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
    p.add_option("--depthfile", help="Use depth info in this file")
    p.add_option("--switch", help="Rename the seqid with two-column file")
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
        ymin, ymax = 0.9, 0.95
        root.text(0.11, 0.96, "Cotton duplication level", color="gray", size=10)
        root.plot([0.1, 0.95], [ymin, ymin], color="gray")
        root.text(0.96, 0.9, "1x", color="gray", va="center")
        root.plot([0.1, 0.95], [ymax, ymax], color="gray")
        root.text(0.96, 0.95, "6x", color="gray", va="center")

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
        yy = [ymin + 0.01 * (x - 1) for x in yy]
        root.plot(xx, yy, "-", color=light)

    # legend showing the orientation of the genes
    draw_gene_legend(root, 0.5, 0.68, 0.5)

    # Zoom
    xpos = 0.835
    ytop = 0.9
    xmin, xmax = 0.18, 0.82
    ymin, ymax = ytop, 0.55
    lc = "k"
    kwargs = dict(lw=3, color=lc, mec=lc, mfc="w", zorder=3)
    root.plot((xpos, xpos), (ymax, 0.63), ":o", **kwargs)
    root.plot((xpos, xmin), (ymax, ymin), ":o", **kwargs)
    root.plot((xpos, xmax), (ymax, ymin), ":o", **kwargs)
    RoundRect(root, (0.06, 0.17), 0.92, 0.35, fill=False, lw=2, ec=light)

    # Panels
    root.text(0.05, 0.95, "a", size=20, fontweight="bold")
    root.text(0.1, 0.45, "b", size=20, fontweight="bold")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cotton"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_diagram(ax, x, y, label="S", title="syntenic", gradient=True):
    """
    Part of the diagrams that are re-used. (x, y) marks the center of the
    diagram. Label determines the modification to the "S" graph.
    """
    trackgap = 0.06
    tracklen = 0.12
    xa, xb = x - tracklen, x + tracklen
    ya, yb = y + trackgap, y - trackgap
    hsps = (((60, 150), (50, 130)), ((190, 225), (200, 240)), ((330, 280), (360, 310)))

    for yy in (ya, yb):
        ax.plot((xa, xb), (yy, yy), "-", color="gray", lw=2, zorder=1)

    ytip = 0.015
    mrange = 400
    m = lambda t: xa + t * 1.0 / mrange * tracklen * 2

    for i, ((a, b), (c, d)) in enumerate(hsps):
        fb = False
        if label == "FB" and i == 1:
            c, d = 270, 280
            fb = True
        if label == "G" and i == 0:
            c, d = 120, 65

        a, b, c, d = [m(t) for t in (a, b, c, d)]
        color = "g" if i == 1 else "r"
        GeneGlyph(ax, a, b, ya, 2 * ytip, fc=color, gradient=gradient, zorder=10)

        if i == 1 and label in ("F", "G", "FN"):
            pass
        else:
            if fb:
                GeneGlyph(
                    ax, c, d, yb, 2 * ytip, fc="w", tip=0, gradient=gradient, zorder=10
                )
            else:
                GeneGlyph(ax, c, d, yb, 2 * ytip, fc="r", gradient=gradient, zorder=10)

        r = Polygon(
            ((a, ya - ytip), (c, yb + ytip), (d, yb + ytip), (b, ya - ytip)),
            fc="r",
            alpha=0.2,
        )

        if i == 1 and label not in ("S", "FB"):
            pass
        elif i == 0 and label == "G":
            pass
        else:
            ax.add_patch(r)

    if label == "FN":
        ax.text(x + 0.005, yb, "NNNNN", ha="center", size=7)

    title = "{0}: {1}".format(label, title)
    ax.text(x, ya + 5 * ytip, title, size=8, ha="center")


def epoch(args):
    """
    %prog epoch

    Illustrate the methods used in Maggie's epoch paper, in particular, how to
    classifiy S/G/F/FB/FN for the genes.
    """
    p = OptionParser(__doc__)
    p.parse_args(args)

    fig = plt.figure(1, (6, 4))
    root = fig.add_axes([0, 0, 1, 1])

    # Separators
    linestyle = dict(lw=2, color="b", alpha=0.2, zorder=2)
    root.plot((0, 1), (0.5, 0.5), "--", **linestyle)
    for i in (1.0 / 3, 2.0 / 3):
        root.plot((i, i), (0.5, 1), "--", **linestyle)
    for i in (1.0 / 6, 3.0 / 6, 5.0 / 6):
        root.plot((i, i), (0, 0.5), "--", **linestyle)

    # Diagrams
    plot_diagram(root, 1.0 / 6, 3.0 / 4, "S", "syntenic")
    plot_diagram(root, 3.0 / 6, 3.0 / 4, "F", "missing, with both flankers")
    plot_diagram(root, 5.0 / 6, 3.0 / 4, "G", "missing, with one flanker")
    plot_diagram(root, 2.0 / 6, 1.0 / 4, "FB", "has non-coding matches")
    plot_diagram(root, 4.0 / 6, 1.0 / 4, "FN", "syntenic region has gap")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


if __name__ == "__main__":
    main()
