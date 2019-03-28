#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the Brassica napus genome manuscript (Chalhoub et al. Science 2014).
"""
from __future__ import print_function

import os.path as op
import sys
import logging

import numpy as np

from jcvi.graphics.base import plt, Rectangle, savefig, mpl, \
            adjust_spines, FancyArrowPatch, normalize_axes, panel_labels
from jcvi.graphics.glyph import TextCircle
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny
from jcvi.graphics.coverage import Coverage, Sizes, XYtrack, setup_gauge_ax
from jcvi.formats.base import LineFile
from jcvi.apps.base import OptionParser, ActionDispatcher


template_cov = """# y, xstart, xend, rotation, color, label, va, bed
.56, {0}, {1}, 0, darkslategray, , top, AN.bed
.48, {2}, {3}, 0, darkslategray, , top, CN.bed
# edges
e, 0, 1, AN.CN.1x1.lifted.simple
"""
template_f3a = r"""# y, xstart, xend, rotation, color, label, va, bed
.65, {0}, {1}, 0, gainsboro, \noindent\textit{{B. napus}} A$\mathsf{{_n}}$2\\(cv Darmor-\textit{{bzh}}), top, AN.bed
.55, {2}, {3}, 0, gainsboro, \textit{{B. rapa}} A$\mathsf{{_r}}$2, top, brapa.bed
.45, {4}, {5}, 0, gainsboro, \textit{{B. oleracea}} C$\mathsf{{_o}}$2, top, boleracea.bed
.35, {6}, {7}, 0, gainsboro, \noindent\textit{{B. napus}} C$\mathsf{{_n}}$2\\(cv Darmor-\textit{{bzh}}), top, CN.bed
# edges
e, 0, 1, AN.brapa.1x1.lifted.simple
e, 1, 2, brapa.boleracea.1x1.lifted.simple
e, 3, 2, CN.boleracea.1x1.lifted.simple"""

gap = .03


class F4ALayoutLine (object):

    def __init__(self, row, delimiter=",", datadir=None):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.region = args[0]
        self.seqid, se = self.region.split(":")
        start, end = se.split("-")
        self.start, self.end = int(start), int(end)
        self.center = (self.start + self.end) / 2
        self.span = self.end - self.start + 1
        self.box_region = args[1]
        self.y = float(args[2])
        self.i = int(args[3])


class F4ALayout (LineFile):

    def __init__(self, filename, delimiter=',', datadir=None):
        super(F4ALayout, self).__init__(filename)
        fp = open(filename)
        self.edges = []
        for row in fp:
            if row[0] == '#':
                continue
            self.append(F4ALayoutLine(row, delimiter=delimiter,
                                    datadir=datadir))


def main():

    actions = (
        ('ploidy', 'plot napus macro-synteny (requires data)'),
        ('expr', 'plot expression values between homeologs (requires data)'),
        ('cov', 'plot coverage graphs between homeologs (requires data)'),
        ('deletion', 'plot histogram for napus deletions (requires data)'),
        ('fig3', 'plot Figure-3'),
        ('fig4', 'plot Figure-4 (not in main text)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def calc_ratio(chrs, sizes):
    chr_sizes = [[sizes[x] for x in z] for z in chrs]
    chr_sum_sizes = [sum(x) for x in chr_sizes]
    ratio = .8 / max(chr_sum_sizes)
    return chr_sizes, chr_sum_sizes, ratio


def center_panel(chr, chr_size, ratio, gap=gap, shift=0):
    # Center two panels
    w = (ratio * chr_size + (len(chr) - 1) * gap) / 2
    return .5 - w + shift, .5 + w + shift


def make_seqids(chrs, seqidsfile="seqids"):
    seqidsfile = "seqids"
    fw = open(seqidsfile, "w")
    for chr in chrs:
        print(",".join(chr), file=fw)
    fw.close()
    logging.debug("File `{0}` written.".format(seqidsfile))
    return seqidsfile


def make_layout(chrs, chr_sizes, ratio, template, klayout="layout", shift=0):
    coords = []
    for chr, chr_size in zip(chrs, chr_sizes):
        coords.extend(center_panel(chr, chr_size, ratio, shift=shift))

    klayout = "layout"
    fw = open(klayout, "w")
    print(template.format(*coords), file=fw)
    fw.close()
    logging.debug("File `{0}` written.".format(klayout))

    return klayout


def cov(args):
    """
    %prog cov chrA01 chrC01 chr.sizes data AN.CN.1x1.lifted.anchors.simple

    Plot coverage graphs between homeologs, the middle panel show the
    homeologous gene pairs. Allow multiple chromosomes to multiple chromosomes.
    """
    p = OptionParser(cov.__doc__)
    p.add_option("--order",
                default="swede,kale,h165,yudal,aviso,abu,bristol,bzh",
                help="The order to plot the tracks, comma-separated")
    p.add_option("--reverse", default=False, action="store_true",
                help="Plot the order in reverse")
    p.add_option("--gauge_step", default=5000000, type="int",
                help="Step size for the base scale")
    p.add_option("--hlsuffix", default="regions.forhaibao",
                help="Suffix for the filename to be used to highlight regions")
    opts, args, iopts = p.set_image_options(args, figsize="11x8")

    if len(args) != 4:
        sys.exit(not p.print_help())

    chr1, chr2, sizesfile, datadir = args
    chr1 = chr1.split(",")
    chr2 = chr2.split(",")

    order = opts.order
    hlsuffix = opts.hlsuffix
    if order:
        order = order.split(",")
        if opts.reverse:
            order.reverse()
    sizes = Sizes(sizesfile).mapping
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    chrs = (chr1, chr2)
    chr_sizes, chr_sum_sizes, ratio = calc_ratio(chrs, sizes)
    chr_size1, chr_size2 = chr_sum_sizes
    chr_sizes1, chr_sizes2 = chr_sizes

    w1_start, w1_end = center_panel(chr1, chr_size1, ratio)
    w2_start, w2_end = center_panel(chr2, chr_size2, ratio)
    w1s = w1_start
    w2s = w2_start

    dsg = "gray"
    i = 0
    for c1, s1 in zip(chr1, chr_sizes1):
        w1 = ratio * s1
        plot_label = i == 0
        i += 1
        canvas1 = (w1s, .6, w1, .3)
        Coverage(fig, root, canvas1, c1, (0, s1), datadir,
                     order=order, gauge="top", plot_label=plot_label,
                     gauge_step=opts.gauge_step, palette=dsg,
                     cap=40, hlsuffix=hlsuffix)
        w1s += w1 + gap

    i = 0
    for c2, s2 in zip(chr2, chr_sizes2):
        w2 = ratio * s2
        plot_label = i == 0
        i += 1
        canvas2 = (w2s, .15, w2, .3)
        Coverage(fig, root, canvas2, c2, (0, s2), datadir,
                     order=order, gauge="bottom", plot_label=plot_label,
                     gauge_step=opts.gauge_step, palette=dsg,
                     cap=40, hlsuffix=hlsuffix)
        w2s += w2 + gap

    # Synteny panel
    seqidsfile = make_seqids(chrs)
    klayout = make_layout(chrs, chr_sum_sizes, ratio, template_cov)
    Karyotype(fig, root, seqidsfile, klayout, gap=gap,
              generank=False, sizes=sizes)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    chr2 = "_".join(chr2)
    if opts.reverse:
        chr2 += ".reverse"
    image_name = chr2 + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def conversion_track(order, filename, col, label, ax, color,
                     ypos=0, asterisk=False):
    ids = []
    fp = open(filename)
    for row in fp:
        if asterisk and row[0] != "*":
            continue
        if (not asterisk) and row[0] == "*":
            continue
        if asterisk:
            row = row[1:]
        atoms = row.split()
        gid = atoms[col].rsplit(".", 1)[0]
        gid = gid.replace('T', 'G')
        ids.append(gid)

    beds = [order[x][1] for x in ids if x in order]
    pts = [x.start for x in beds if x.seqid == label]
    if len(pts):
        logging.debug("A total of {0} converted loci imported.".format(len(pts)))
    else:
        logging.error("Array empty. Skipped scatterplot.")
        return

    ax.vlines(pts, [-1], [ypos], color=color)
    ax.set_axis_off()


def make_affix_axis(fig, t, yoffset, height=.001):
    x, y = t.xstart, t.y + yoffset
    w = t.xend - t.xstart
    ax = fig.add_axes([x, y, w, height])
    start, end = 0, t.total
    return ax


def fig3(args):
    """
    %prog fig3 chrA02,A02,C2,chrC02 chr.sizes all.bed data

    Napus Figure 3 displays alignments between quartet chromosomes, inset
    with read histograms.
    """
    from jcvi.formats.bed import Bed

    p = OptionParser(fig3.__doc__)
    p.add_option("--gauge_step", default=10000000, type="int",
                help="Step size for the base scale")
    opts, args, iopts = p.set_image_options(args, figsize="12x9")

    if len(args) != 4:
        sys.exit(not p.print_help())

    chrs, sizes, bedfile, datadir = args
    gauge_step = opts.gauge_step
    diverge = iopts.diverge
    rr, gg = diverge
    chrs = [[x] for x in chrs.split(",")]
    sizes = Sizes(sizes).mapping

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    chr_sizes, chr_sum_sizes, ratio = calc_ratio(chrs, sizes)

    # Synteny panel
    seqidsfile = make_seqids(chrs)
    klayout = make_layout(chrs, chr_sum_sizes, ratio, template_f3a, shift=.05)
    height = .07
    r = height / 4
    K = Karyotype(fig, root, seqidsfile, klayout, gap=gap,
                  height=height, lw=2, generank=False, sizes=sizes,
                  heightpad=r, roundrect=True, plot_label=False)

    # Chromosome labels
    for kl in K.layout:
        if kl.empty:
            continue
        lx, ly = kl.xstart, kl.y
        if lx < .11:
            lx += .1
            ly += .06
        label = kl.label
        root.text(lx - .015, ly, label, fontsize=15,
                  ha="right", va="center")

    # Inset with datafiles
    datafiles = ("chrA02.bzh.forxmgr", "parent.A02.per10kb.forxmgr",
                 "parent.C2.per10kb.forxmgr", "chrC02.bzh.forxmgr")
    datafiles = [op.join(datadir, x) for x in datafiles]
    tracks = K.tracks
    hlfile = op.join(datadir, "bzh.regions.forhaibao")
    xy_axes = []
    for t, datafile in zip(tracks, datafiles):
        ax = make_affix_axis(fig, t, -r, height=2 * r)
        xy_axes.append(ax)
        chr = t.seqids[0]
        xy = XYtrack(ax, datafile, color="lightslategray")
        start, end = 0, t.total
        xy.interpolate(end)
        xy.cap(ymax=40)
        xy.import_hlfile(hlfile, chr, diverge=diverge)
        xy.draw()
        ax.set_xlim(start, end)
        gauge_ax = make_affix_axis(fig, t, -r)
        adjust_spines(gauge_ax, ["bottom"])
        setup_gauge_ax(gauge_ax, start, end, gauge_step)

    # Converted gene tracks
    ax_Ar = make_affix_axis(fig, tracks[1], r, height=r/2)
    ax_Co = make_affix_axis(fig, tracks[2], r, height=r/2)

    order = Bed(bedfile).order
    for asterisk in (False, True):
        conversion_track(order, "data/Genes.Converted.seuil.0.6.AtoC.txt",
                         0, "A02", ax_Ar, rr, asterisk=asterisk)
        conversion_track(order, "data/Genes.Converted.seuil.0.6.AtoC.txt",
                         1, "C2", ax_Co, gg, asterisk=asterisk)
        conversion_track(order, "data/Genes.Converted.seuil.0.6.CtoA.txt",
                         0, "A02", ax_Ar, gg, ypos=1, asterisk=asterisk)
        conversion_track(order, "data/Genes.Converted.seuil.0.6.CtoA.txt",
                         1, "C2", ax_Co, rr, ypos=1, asterisk=asterisk)

    Ar, Co = xy_axes[1:3]
    annotations = ((Ar, "Bra028920 Bra028897", "center", "1DAn2+"),
                   (Ar, "Bra020081 Bra020171", "right", "2DAn2+"),
                   (Ar, "Bra020218 Bra020286", "left", "3DAn2+"),
                   (Ar, "Bra008143 Bra008167", "left", "4DAn2-"),
                   (Ar, "Bra029317 Bra029251", "right", "5DAn2+ (GSL)"),
                   (Co, "Bo2g001000 Bo2g001300", "left", "1DCn2-"),
                   (Co, "Bo2g018560 Bo2g023700", "right", "2DCn2-"),
                   (Co, "Bo2g024450 Bo2g025390", "left", "3DCn2-"),
                   (Co, "Bo2g081060 Bo2g082340", "left", "4DCn2+"),
                   (Co, "Bo2g161510 Bo2g164260", "right", "5DCn2-"))

    for ax, genes, ha, label in annotations:
        g1, g2 = genes.split()
        x1, x2 = order[g1][1].start, order[g2][1].start
        if ha == "center":
            x = (x1 + x2) / 2 * .8
        elif ha == "left":
            x = x2
        else:
            x = x1
        label = r"\textit{{{0}}}".format(label)
        color = rr if "+" in label else gg
        ax.text(x, 30, label, color=color, fontsize=9, ha=ha, va="center")

    ax_Ar.set_xlim(0, tracks[1].total)
    ax_Ar.set_ylim(-1, 1)
    ax_Co.set_xlim(0, tracks[2].total)
    ax_Co.set_ylim(-1, 1)

    # Plot coverage in resequencing lines
    gstep = 5000000
    order = "swede,kale,h165,yudal,aviso,abu,bristol".split(",")
    labels_dict = {"h165": "Resynthesized (H165)", "abu": "Aburamasari"}
    hlsuffix = "regions.forhaibao"
    chr1, chr2 = "chrA02", "chrC02"
    t1, t2 = tracks[0], tracks[-1]
    s1, s2 = sizes[chr1], sizes[chr2]

    canvas1 = (t1.xstart, .75, t1.xend - t1.xstart, .2)
    c = Coverage(fig, root, canvas1, chr1, (0, s1), datadir,
                 order=order, gauge=None, plot_chr_label=False,
                 gauge_step=gstep, palette="gray",
                 cap=40, hlsuffix=hlsuffix, labels_dict=labels_dict,
                 diverge=diverge)
    yys = c.yys
    x1, x2 = .37, .72
    tip = .02
    annotations = ((x1, yys[2] + .3 * tip, tip, tip / 2, "FLC"),
                   (x1, yys[3] + .6 * tip, tip, tip / 2, "FLC"),
                   (x1, yys[5] + .6 * tip, tip, tip / 2, "FLC"),
                   (x2, yys[0] + .9 * tip, -1.2 * tip, 0, "GSL"),
                   (x2, yys[4] + .9 * tip, -1.2 * tip, 0, "GSL"),
                   (x2, yys[6] + .9 * tip, -1.2 * tip, 0, "GSL"))

    arrowprops=dict(facecolor='black', shrink=.05, frac=.5,
                    width=1, headwidth=4)
    for x, y, dx, dy, label in annotations:
        label = r"\textit{{{0}}}".format(label)
        root.annotate(label, xy=(x, y), xytext=(x + dx, y + dy),
                      arrowprops=arrowprops, color=rr, fontsize=9,
                      ha="center", va="center")

    canvas2 = (t2.xstart, .05, t2.xend - t2.xstart, .2)
    Coverage(fig, root, canvas2, chr2, (0, s2), datadir,
                 order=order, gauge=None, plot_chr_label=False,
                 gauge_step=gstep, palette="gray",
                 cap=40, hlsuffix=hlsuffix, labels_dict=labels_dict,
                 diverge=diverge)

    pad = .03
    labels = ((.1, .67, "A"), (t1.xstart - 3 * pad, .95 + pad, "B"),
              (t2.xstart - 3 * pad, .25 + pad, "C"))
    panel_labels(root, labels)
    normalize_axes(root)

    image_name = "napus-fig3." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def fig4(args):
    """
    %prog fig4 layout data

    Napus Figure 4A displays an example deleted region for quartet chromosomes,
    showing read alignments from high GL and low GL lines.
    """
    p = OptionParser(fig4.__doc__)
    p.add_option("--gauge_step", default=200000, type="int",
                help="Step size for the base scale")
    opts, args, iopts = p.set_image_options(args, figsize="9x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    layout, datadir = args
    layout = F4ALayout(layout, datadir=datadir)

    gs = opts.gauge_step
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    block, napusbed, slayout = "r28.txt", "all.bed", "r28.layout"
    s = Synteny(fig, root, block, napusbed, slayout, chr_label=False)
    synteny_exts = [(x.xstart, x.xend) for x in s.rr]

    h = .1
    order = "bzh,yudal".split(",")
    labels = (r"\textit{B. napus} A$\mathsf{_n}$2",
              r"\textit{B. rapa} A$\mathsf{_r}$2",
              r"\textit{B. oleracea} C$\mathsf{_o}$2",
              r"\textit{B. napus} C$\mathsf{_n}$2")
    for t in layout:
        xstart, xend = synteny_exts[2 * t.i]
        canvas = [xstart, t.y, xend - xstart, h]
        root.text(xstart - h, t.y + h / 2, labels[t.i], ha="center", va="center")
        ch, ab = t.box_region.split(":")
        a, b = ab.split("-")
        vlines = [int(x) for x in (a, b)]
        Coverage(fig, root, canvas, t.seqid, (t.start, t.end), datadir,
                     order=order, gauge="top", plot_chr_label=False,
                     gauge_step=gs, palette="gray",
                     cap=40, hlsuffix="regions.forhaibao",
                     vlines=vlines)

    # Highlight GSL biosynthesis genes
    a, b = (3, "Bra029311"), (5, "Bo2g161590")
    for gid in (a, b):
        start, end = s.gg[gid]
        xstart, ystart = start
        xend, yend = end
        x = (xstart + xend) / 2
        arrow = FancyArrowPatch(posA=(x, ystart - .04),
                                posB=(x, ystart - .005),
                                arrowstyle="fancy,head_width=6,head_length=8",
                                lw=3, fc='k', ec='k', zorder=20)
        root.add_patch(arrow)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = "napus-fig4." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def deletion(args):
    """
    %prog deletion [deletion-genes|deletion-bases] C2-deletions boleracea.bed

    Plot histogram for napus deletions. Can plot deletion-genes or
    deletion-bases. The three largest segmental deletions will be highlighted
    along with a drawing of the C2 chromosome.
    """
    import math
    from jcvi.formats.bed import Bed
    from jcvi.graphics.chromosome import HorizontalChromosome
    from jcvi.graphics.base import kb_formatter

    p = OptionParser(deletion.__doc__)
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
    HorizontalChromosome(root, na, nb, .5, height=.025,
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


def ploidy(args):
    """
    %prog ploidy seqids layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of B. napus genome.
    """
    p = OptionParser(ploidy.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, klayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)

    fc = "darkslategrey"
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

    bb = dict(boxstyle="round,pad=.5", fc="w", ec="0.5", alpha=0.5)
    root.text(.5, .2 + ot, r"\noindent\textit{Brassica napus}\\"
                "(A$\mathsf{_n}$C$\mathsf{_n}$ genome)", ha="center",
                size=16, color="k", bbox=bb)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "napus"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def expr(args):
    """
    %prog expr block exp layout napus.bed

    Plot a composite figure showing synteny and the expression level between
    homeologs in two tissues - total 4 lists of values. block file contains the
    gene pairs between AN and CN.
    """
    from jcvi.graphics.base import red_purple as default_cm

    p = OptionParser(expr.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    if len(args) != 4:
        sys.exit(not p.print_help())

    block, exp, layout, napusbed = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    s = Synteny(fig, root, block, napusbed, layout)

    # Import the expression values
    # Columns are: leaf-A, leaf-C, root-A, root-C
    fp = open(exp)
    data = {}
    for row in fp:
        gid, lf, rt = row.split()
        lf, rt = float(lf), float(rt)
        data[gid] = (lf, rt)

    rA, rB = s.rr
    gA = [x.accn for x in rA.genes]
    gC = [x.accn for x in rB.genes]

    A = [data.get(x, (0, 0)) for x in gA]
    C = [data.get(x, (0, 0)) for x in gC]
    A = np.array(A)
    C = np.array(C)
    A = np.transpose(A)
    C = np.transpose(C)

    d, h = .01, .1
    lsg = "lightslategrey"
    coords = s.gg  # Coordinates of the genes
    axes = []
    for j, (y, gg) in enumerate(((.79, gA), (.24, gC))):
        r = s.rr[j]
        x = r.xstart
        w = r.xend - r.xstart
        ax = fig.add_axes([x, y, w, h])
        axes.append(ax)
        root.add_patch(Rectangle((x - h, y - d), w + h + d, h + 2 * d, fill=False,
                                ec=lsg, lw=1))
        root.text(x - d, y + 3 * h / 4, "root", ha="right", va="center")
        root.text(x - d, y + h / 4, "leaf", ha="right", va="center")
        ty = y - 2 * d if y > .5 else y + h + 2 * d
        nrows = len(gg)
        for i, g in enumerate(gg):
            start, end = coords[(j, g)]
            sx, sy = start
            ex, ey = end
            assert sy == ey
            sy = sy + 2 * d if sy > .5 else sy - 2 * d
            root.plot(((sx + ex) / 2, x + w * (i + .5) / nrows), (sy, ty),
                            lw=1, ls=":", color="k", alpha=.2)

    axA, axC = axes
    p = axA.pcolormesh(A, cmap=default_cm)
    p = axC.pcolormesh(C, cmap=default_cm)
    axA.set_xlim(0, len(gA))
    axC.set_xlim(0, len(gC))

    x, y, w, h = .35, .1, .3, .05
    ax_colorbar = fig.add_axes([x, y, w, h])
    fig.colorbar(p, cax=ax_colorbar, orientation='horizontal')
    root.text(x - d, y + h / 2, "RPKM", ha="right", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    for x in (axA, axC, root):
        x.set_axis_off()

    image_name = "napusf4b." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
