#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Plotting scripts for the vanilla genome paper.
"""
from __future__ import print_function

import logging
import sys

from jcvi.apps.base import ActionDispatcher, OptionParser
from jcvi.compara.synteny import AnchorFile, check_beds
from jcvi.formats.base import get_number
from jcvi.formats.bed import Bed, BedLine
from jcvi.graphics.base import normalize_axes, panel_labels, plt, savefig
from jcvi.graphics.glyph import TextCircle
from jcvi.graphics.synteny import Synteny, draw_gene_legend


def main():
    actions = (
        # Chromosome painting since WGD
        ("ancestral", "paint 14 chromosomes following alpha WGD (requires data)"),
        # main figures in text
        ("ploidy", "plot vanilla synteny (requires data)"),
        # Composite phylogeny - tree and ks
        ("phylogeny", "create a composite figure with tree and ks"),
        ("tree", "create a separate figure with tree"),
        ("ks", "create a separate figure with ks"),
        # Composite synteny - wgd and microsynteny
        ("synteny", "create a composite figure with wgd and microsynteny"),
        ("wgd", "create separate figures with wgd"),
        ("microsynteny", "create separate figures with microsynteny"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def phylogeny(args):
    """
    %prog phylogeny treefile ks.layout

    Create a composite figure with (A) tree and (B) ks.
    """
    from jcvi.graphics.tree import parse_tree, LeafInfoFile, WGDInfoFile, draw_tree

    p = OptionParser(phylogeny.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x12")

    (datafile, layoutfile) = args

    logging.debug("Load tree file `{0}`".format(datafile))
    t, hpd = parse_tree(datafile)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax1 = fig.add_axes([0, 0.4, 1, 0.6])
    ax2 = fig.add_axes([0.12, 0.065, 0.8, 0.3])

    supportcolor = "k"
    margin, rmargin = 0.1, 0.2  # Left and right margin
    leafinfo = LeafInfoFile("leafinfo.csv").cache
    wgdinfo = WGDInfoFile("wgdinfo.csv").cache
    outgroup = "ginkgo"

    # Panel A
    draw_tree(
        ax1,
        t,
        hpd=hpd,
        margin=margin,
        rmargin=rmargin,
        supportcolor=None,
        internal=False,
        outgroup=outgroup,
        reroot=False,
        leafinfo=leafinfo,
        wgdinfo=wgdinfo,
        geoscale=True,
    )

    from jcvi.apps.ks import Layout, KsPlot, KsFile

    # Panel B
    ks_min = 0.0
    ks_max = 3.0
    bins = 60
    fill = False
    layout = Layout(layoutfile)
    print(layout, file=sys.stderr)

    kp = KsPlot(ax2, ks_max, bins, legendp="upper right")
    for lo in layout:
        data = KsFile(lo.ksfile)
        data = [x.ng_ks for x in data]
        data = [x for x in data if ks_min <= x <= ks_max]
        kp.add_data(
            data,
            lo.components,
            label=lo.label,
            color=lo.color,
            marker=lo.marker,
            fill=fill,
            fitted=False,
            kde=True,
        )

    kp.draw(filename=None)

    normalize_axes([root, ax1])
    labels = ((0.05, 0.95, "A"), (0.05, 0.4, "B"))
    panel_labels(root, labels)

    image_name = "phylogeny.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def tree(args):
    """
    %prog tree treefile

    Create a tree figure.
    """
    from jcvi.graphics.tree import parse_tree, LeafInfoFile, WGDInfoFile, draw_tree

    p = OptionParser(tree.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x8")

    (datafile,) = args
    logging.debug("Load tree file `{0}`".format(datafile))
    t, hpd = parse_tree(datafile)

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax1 = fig.add_axes([0, 0, 1, 1])

    supportcolor = "k"
    margin, rmargin = 0.1, 0.2  # Left and right margin
    leafinfo = LeafInfoFile("leafinfo.csv").cache
    wgdinfo = WGDInfoFile("wgdinfo.csv").cache
    outgroup = "ginkgo"

    # Panel A
    draw_tree(
        ax1,
        t,
        hpd=hpd,
        margin=margin,
        rmargin=rmargin,
        supportcolor=None,
        internal=False,
        outgroup=outgroup,
        reroot=False,
        leafinfo=leafinfo,
        wgdinfo=wgdinfo,
        geoscale=True,
    )

    normalize_axes([ax1])
    image_name = "tree.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def ks(args):
    """
    %prog ks ks.layout

    Create a ks figure.
    """
    p = OptionParser(ks.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x4")

    (layoutfile,) = args

    from jcvi.apps.ks import Layout, KsPlot, KsFile

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax2 = fig.add_axes([0.12, 0.12, 0.8, 0.8])

    # Panel B
    ks_min = 0.0
    ks_max = 3.0
    bins = 60
    fill = False
    layout = Layout(layoutfile)
    print(layout, file=sys.stderr)

    kp = KsPlot(ax2, ks_max, bins, legendp="upper right")
    for lo in layout:
        data = KsFile(lo.ksfile)
        data = [x.ng_ks for x in data]
        data = [x for x in data if ks_min <= x <= ks_max]
        kp.add_data(
            data,
            lo.components,
            label=lo.label,
            color=lo.color,
            marker=lo.marker,
            fill=fill,
            fitted=False,
            kde=True,
        )

    kp.draw(filename=None)

    image_name = "ks.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def synteny(args):
    """
    %prog synteny vplanifoliaA_blocks.bed vplanifoliaA.sizes \
        b1.blocks all.bed b1.layout

    Create a composite figure with (A) wgd and (B) microsynteny.
    """
    from jcvi.graphics.chromosome import draw_chromosomes

    p = OptionParser(synteny.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x12")

    (bedfile, sizesfile, blocksfile, allbedfile, blockslayout) = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax1 = fig.add_axes([0, 0.5, 1, 0.5])
    ax2 = fig.add_axes([0.02, 0, 0.98, 0.5])

    # Panel A
    title = r"Genome duplication $\alpha^{O}$ event in $\textit{Vanilla}$"
    draw_chromosomes(
        ax1,
        bedfile,
        sizes=sizesfile,
        iopts=iopts,
        mergedist=200000,
        winsize=50000,
        imagemap=False,
        gauge=True,
        legend=False,
        title=title,
    )

    # Panel B
    draw_ploidy(fig, ax2, iopts, blocksfile, allbedfile, blockslayout)

    normalize_axes([root, ax1, ax2])
    labels = ((0.05, 0.95, "A"), (0.05, 0.5, "B"))
    panel_labels(root, labels)

    image_name = "synteny.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def wgd(args):
    """
    %prog wgd vplanifoliaA_blocks.bed vplanifoliaA.sizes

    Create a wgd figure.
    """
    from jcvi.graphics.chromosome import draw_chromosomes

    p = OptionParser(synteny.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x5")

    (bedfile, sizesfile) = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax1 = fig.add_axes([0, 0, 1, 1])

    title = r"Genome duplication $\alpha^{O}$ event in $\textit{Vanilla}$"
    draw_chromosomes(
        ax1,
        bedfile,
        sizes=sizesfile,
        iopts=iopts,
        mergedist=200000,
        winsize=50000,
        imagemap=False,
        gauge=True,
        legend=False,
        title=title,
    )

    normalize_axes([ax1])

    image_name = "wgd.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def microsynteny(args):
    """
    %prog microsynteny b1.blocks all.bed b1.layout

    Create a microsynteny figure.
    """
    p = OptionParser(synteny.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x6")

    (blocksfile, allbedfile, blockslayout) = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax2 = fig.add_axes([0, 0, 1, 1])

    draw_ploidy(fig, ax2, iopts, blocksfile, allbedfile, blockslayout)

    normalize_axes([ax2])

    image_name = "microsynteny.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def ancestral(args):
    """
    %prog ancestral vplanifoliaA.vplanifoliaA.anchors > vplanifoliaA_blocks.bed

    Paint 14 chromosomes following alpha WGD.
    """
    p = OptionParser(ancestral.__doc__)
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorsfile,) = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts)

    # We focus on the following chromosome pairs
    target_pairs = set(
        (
            (1, 1),
            (1, 6),
            (1, 8),
            (1, 13),
            (2, 4),
            (3, 12),
            (3, 14),
            (5, 6),
            (5, 8),
            (7, 9),
            (7, 11),
            (9, 10),
            (10, 11),
        )
    )

    def get_target(achr, bchr):
        if "chr" not in achr and "chr" not in bchr:
            return None
        achr, bchr = get_number(achr), get_number(bchr)
        if achr > bchr:
            achr, bchr = bchr, achr
        if (achr, bchr) in target_pairs:
            return achr, bchr
        return None

    def build_bedline(astart, aend, target_pair):
        # target_name = "{:02d}-{:02d}".format(*target_pair)
        target_name = [str(x) for x in target_pair if x in (1, 2, 3, 5, 7, 10)][0]
        return "\t".join(
            str(x) for x in (astart.seqid, astart.start, aend.end, target_name)
        )

    # Iterate through the blocks, store any regions that has hits to one of the
    # target_pairs
    ac = AnchorFile(anchorsfile)
    blocks = ac.blocks
    outbed = Bed()
    for i, block in enumerate(blocks):
        a, b, scores = zip(*block)
        a = [qorder[x] for x in a]
        b = [sorder[x] for x in b]
        astart, aend = min(a)[1], max(a)[1]
        bstart, bend = min(b)[1], max(b)[1]
        # Now convert to BED lines with new accn
        achr, bchr = astart.seqid, bstart.seqid
        target = get_target(achr, bchr)
        if target is None:
            continue
        outbed.add(build_bedline(astart, aend, target))
        outbed.add(build_bedline(bstart, bend, target))
    outbed.print_to_file(sorted=True)


def ploidy(args):
    """
    %prog ploidy b1.blocks all.bed b1.layout

    Build a figure that illustrates the WGD history of the vanilla genome.
    """
    p = OptionParser(ploidy.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="12x6")

    if len(args) != 3:
        sys.exit(not p.print_help())

    blocksfile, bedfile, blockslayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    draw_ploidy(fig, root, iopts, blocksfile, bedfile, blockslayout)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "vanilla-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def draw_ploidy(fig, root, iopts, blocksfile, bedfile, blockslayout):
    switchidsfile = "switch.ids"
    Synteny(
        fig,
        root,
        blocksfile,
        bedfile,
        blockslayout,
        scalebar=True,
        switch=switchidsfile,
    )

    # Legend showing the orientation of the genes
    draw_gene_legend(root, 0.2, 0.3, 0.53)

    # WGD labels
    radius = 0.025
    tau_color = "#bebada"
    alpha_color = "#bc80bd"
    label_color = "k"
    pad = 0.05
    for y in (0.74 + 1.5 * pad, 0.26 - 1.5 * pad):
        TextCircle(
            root,
            0.25,
            y,
            r"$\alpha^{O}$",
            radius=radius,
            fc=alpha_color,
            color=label_color,
            fontweight="bold",
        )
        TextCircle(
            root,
            0.75,
            y,
            r"$\alpha^{O}$",
            radius=radius,
            fc=alpha_color,
            color=label_color,
            fontweight="bold",
        )
    for y in (0.74 + 3 * pad, 0.26 - 3 * pad):
        TextCircle(
            root, 0.5, y, r"$\tau$", radius=radius, fc=tau_color, color=label_color
        )


if __name__ == "__main__":
    main()
