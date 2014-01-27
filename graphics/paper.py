#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the various past manuscripts.
"""

import sys
import logging

import numpy as np

from jcvi.graphics.base import plt, _, Rectangle, Polygon, CirclePolygon, \
        savefig, mpl, cm
from jcvi.graphics.glyph import GeneGlyph, RoundLabel, RoundRect, \
        arrowprops, TextCircle, plot_cap
from jcvi.graphics.chromosome import Chromosome
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
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
        # Amborella paper (Albert et al., 2013 Science)
        ('amborella', 'plot amborella macro- and micro-synteny (requires data)'),
        # Unpublished
        ('litchi', 'plot litchi micro-synteny (requires data)'),
        ('napus', 'plot napus macro-synteny (requires data)'),
        ('napusexp', 'plot expression values between homeologs (requires data)'),
        ('napuscov', 'plot coverage graphs between homeologs (requires data)'),
        ('napusdeletion', 'plot histogram for napus deletions (requires data)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


template = """# y, xstart, xend, rotation, color, label, va, bed
.56, {0}, {1}, 0, darkslategray, , top, AN.bed
.48, {2}, {3}, 0, darkslategray, , top, CN.bed
# edges
e, 0, 1, {4}
"""

def napuscov(args):
    """
    %prog napuscov chrA01 chrC01 chr.sizes data AN.CN.1x1.lifted.anchors.simple

    Plot coverage graphs between homeologs, the middle panel show the
    homeologous gene pairs.
    """
    from jcvi.graphics.coverage import Coverage, Sizes

    p = OptionParser(napuscov.__doc__)
    p.add_option("--order",
                default="swede,kale,h165,yudal,aviso,abu,bristol,bzh",
                help="The order to plot the tracks, comma-separated")
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 5:
        sys.exit(not p.print_help())

    chr1, chr2, sizes, datadir, simplefile = args
    order = opts.order
    if order:
        order = order.split(",")
    sizes = Sizes(sizes)
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    chr_size1 = sizes.get_size(chr1)
    chr_size2 = sizes.get_size(chr2)
    ratio = .8 / max(chr_size1, chr_size2)

    w1 = ratio * chr_size1
    w1_start = .5 - w1 / 2
    w1_end = .5 + w1 / 2
    canvas1 = (w1_start, .6, w1, .3)
    c = Coverage(fig, root, canvas1, chr1, (0, chr_size1), datadir,
                 order=order, gauge="top")

    w2 = ratio * chr_size2
    w2_start = .5 - w2 / 2
    w2_end = .5 + w2 / 2
    canvas2 = (w2_start, .15, w2, .3)
    c = Coverage(fig, root, canvas2, chr2, (0, chr_size2), datadir,
                 order=order, gauge="bottom")

    # Synteny panel
    seqidsfile = "seqids"
    fw = open(seqidsfile, "w")
    print >> fw, chr1
    print >> fw, chr2
    fw.close()
    logging.debug("File `{0}` written.".format(seqidsfile))

    klayout = "layout"
    fw = open(klayout, "w")
    print >> fw, template.format(w1_start, w1_end, w2_start, w2_end, simplefile)
    fw.close()
    logging.debug("File `{0}` written.".format(klayout))

    Karyotype(fig, root, seqidsfile, klayout)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr2 + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def napusdeletion(args):
    """
    %prog napusdeletion [deletion-genes|deletion-bases] C2-deletions boleracea.bed

    Plot histogram for napus deletions. Can plot deletion-genes or
    deletion-bases. The three largest segmental deletions will be highlighted
    along with a drawing of the C2 chromosome.
    """
    import math
    from jcvi.formats.bed import Bed
    from jcvi.graphics.glyph import TextCircle
    from jcvi.graphics.chromosome import HorizontalChromosome
    from jcvi.graphics.base import kb_formatter

    p = OptionParser(napusdeletion.__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    deletion_genes, deletions, bed = args
    dg = [int(x) for x in open(deletion_genes)]
    dsg, lsg = "darkslategray", "lightslategray"

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([.1, .1, .8, .8])
    minval = 2 if deletion_genes == "deleted-genes" else 2048
    bins = np.logspace(math.log(minval, 10), math.log(max(dg), 10), 16)
    n, bins, histpatches = ax.hist(dg, bins=bins, \
                                   fc=lsg, alpha=.75)
    ax.set_xscale('log', basex=2)
    if deletion_genes == "deleted-genes":
        ax.xaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
        ax.set_xlabel('No. of deleted genes in each segment')
    else:
        ax.xaxis.set_major_formatter(kb_formatter)
        ax.set_xlabel('No. of deleted bases in each segment')
    ax.yaxis.set_major_formatter(mpl.ticker.FormatStrFormatter('%d'))
    ax.set_ylabel('No. of segments')
    ax.patch.set_alpha(0.1)

    # Draw chromosome C2
    na, nb = .45, .85
    root.text((na + nb) / 2, .54, "ChrC02", ha="center")
    hc = HorizontalChromosome(root, na, nb, .5, height=.025,
                             fc=lsg, fill=True)

    order = Bed(bed).order
    fp = open(deletions)
    scale = lambda x: na + x * (nb - na) / 52886895
    for i, row in enumerate(fp):
        i += 1
        num, genes = row.split()
        genes = genes.split("|")
        ia, a = order[genes[0]]
        ib, b = order[genes[-1]]
        mi, mx = a.start, a.end
        mi, mx = scale(mi), scale(mx)
        root.add_patch(Rectangle((mi, .475), mx - mi, .05,
                       fc="red", ec="red"))
        if i == 1:   # offset between two adjacent regions for aesthetics
            mi -= .015
        elif i == 2:
            mi += .015
        TextCircle(root, mi, .44, str(i), fc="red")

    for i, mi in zip(range(1, 4), (.83, .78, .73)):
        TextCircle(root, mi, .2, str(i), fc="red")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = deletion_genes + ".pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def litchi(args):
    """
    %prog litchi mcscan.out all.bed layout switch.ids

    Build a composite figure that calls graphis.synteny.
    """
    from jcvi.graphics.glyph import DoubleSquare

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

    def join_nodes(root, coords, a, b, x, circle=True):
        # Join node a and b to make an internal node
        ax, ay = coords[a]
        bx, by = coords[b]
        nx, ny = x, (ay + by) / 2
        root.plot((nx, ax), (ay, ay), lw=2, color=fc)
        root.plot((nx, bx), (by, by), lw=2, color=fc)
        root.plot((nx, nx), (ay, by), lw=2, color=fc)
        if circle:
            DoubleSquare(root, nx, ny, fc=fc)
        return nx, ny

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


def napus(args):
    """
    %prog napus seqids layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of B. napus genome.
    """
    p = OptionParser(napus.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, klayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)

    fc = "lightslategrey"
    radius = .012
    ot = -.05  # use this to adjust vertical position of the left panel
    TextCircle(root, .1, .9 + ot, r'$\gamma$', radius=radius, fc=fc)
    root.text(.1, .88 + ot, r"$\times3$", ha="center", va="top", color=fc)
    TextCircle(root, .08, .79 + ot, r'$\alpha$', radius=radius, fc=fc)
    TextCircle(root, .12, .79 + ot, r'$\beta$', radius=radius, fc=fc)
    root.text(.1, .77 + ot, r"$\times3\times2\times2$", ha="center", va="top", color=fc)
    root.text(.1, .67 + ot, r"Brassica triplication", ha="center",
                va="top", color=fc, size=11)
    root.text(.1, .65 + ot, r"$\times3\times2\times2\times3$", ha="center", va="top", color=fc)
    root.text(.1, .42 + ot, r"Allo-tetraploidy", ha="center",
                va="top", color=fc, size=11)
    root.text(.1, .4 + ot, r"$\times3\times2\times2\times3\times2$", ha="center", va="top", color=fc)

    bb = dict(boxstyle="round", fc="w", ec="0.5", alpha=0.5)
    root.text(.5, .2 + ot, "Brassica napus", ha="center",
                size=18, color="g", bbox=bb)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "napus"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def napusexp(args):
    """
    %prog napusexp block exp layout napus.bed

    Plot a composite figure showing synteny and the expression level between
    homeologs in two tissues - total 4 lists of values. block file contains the
    gene pairs between AN and CN.
    """
    from matplotlib.colors import LogNorm
    from jcvi.graphics.base import red_purple as default_cm

    p = OptionParser(napusexp.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 4:
        sys.exit(not p.print_help())

    block, exp, layout, napusbed = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    s = Synteny(fig, root, block, napusbed, layout)
    draw_gene_legend(root, .22, .55, .57)

    # Import the expression values
    # Columns are: leaf-A, leaf-C, root-A, root-C
    pairs = [row.split() for row in open(block)]
    data = np.loadtxt(exp)
    nrows = len(pairs)
    assert data.shape[0] == nrows, "block and exp row counts mismatch"
    A = data[:, [2, 0]]
    C = data[:, [3, 1]]
    A = np.transpose(A)
    C = np.transpose(C)

    x, y, d, w, h = .18, .64, .008, .65, .08
    lsg = "lightslategrey"
    coords = s.gg  # Coordinates of the genes
    Ag = [p[0] for p in pairs]
    Cg = [p[1] for p in pairs]

    for y, Gg in ((.64, Ag), (.29, Cg)):
        root.add_patch(Rectangle((x - h, y - d), w + h + d, h + 2 * d, fill=False,
                                ec=lsg, lw=1))
        root.text(x - d, y + 3 * h / 4, "leaf", ha="right", va="center")
        root.text(x - d, y + h / 4, "root", ha="right", va="center")
        ty = y - 2 * d if y > .5 else y + h + 2 * d
        for i, g in enumerate(Gg):
            start, end = coords[g]
            sx, sy = start
            ex, ey = end
            assert sy == ey
            sy = sy + 2 * d if sy > .5 else sy - 2 * d
            root.plot(((sx + ex) / 2, x + w * (i + .5)/ nrows), (sy, ty),
                            lw=2, ls=":", color="k", alpha=.2)

    axA = fig.add_axes([x, .64, w, h])
    axC = fig.add_axes([x, .29, w, h])

    norm = LogNorm(1, 10000)
    p = axA.pcolormesh(A, cmap=default_cm, norm=norm)
    p = axC.pcolormesh(C, cmap=default_cm, norm=norm)
    axA.set_xlim(0, nrows)
    axC.set_xlim(0, nrows)

    x, y, w, h = .35, .17, .3, .03
    ax_colorbar = fig.add_axes([x, y, w, h])
    fig.colorbar(p, cax=ax_colorbar, orientation='horizontal')
    root.text(x - d, y + h / 2, "RPKM", ha="right", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    for x in (axA, axC, root):
        x.set_axis_off()

    pf = "napusexp"
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
