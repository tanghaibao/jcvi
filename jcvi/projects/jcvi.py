#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the JCVI manuscript.
"""

import sys

import networkx as nx

from ..apps.base import ActionDispatcher, OptionParser, logger
from ..assembly.geneticmap import draw_geneticmap_heatmap
from ..assembly.hic import draw_hic_heatmap
from ..assembly.kmer import draw_ks_histogram
from ..compara.pedigree import Pedigree, calculate_inbreeding
from ..compara.synteny import check_beds
from ..graphics.base import (
    cm,
    load_image,
    normalize_axes,
    panel_labels,
    plt,
    savefig,
    set1,
    setup_theme,
)
from ..graphics.chromosome import draw_chromosomes
from ..graphics.dotplot import dotplot
from ..graphics.karyotype import Karyotype
from ..graphics.landscape import draw_heatmaps, draw_multi_depth, draw_stacks
from ..graphics.synteny import Synteny, draw_gene_legend


def synteny(args):
    """
    %prog synteny grape.peach.anchors seqids layout blocks grape_peach.bed blocks.layout

    Plot synteny composite figure, including:
    A. Synteny dotplot
    B. Karyotype plot
    """
    p = OptionParser(synteny.__doc__)
    p.set_beds()
    opts, args, iopts = p.set_image_options(args, figsize="14x7")
    setup_theme(style="dark")

    if len(args) != 6:
        sys.exit(not p.print_help())

    anchorfile, seqidsfile, layoutfile, datafile, bedfile, blockslayoutfile = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    ax1_root = fig.add_axes((0, 0, 0.5, 1))
    ax1_canvas = fig.add_axes((0.05, 0.1, 0.4, 0.8))  # the dot plot
    ax2_root = fig.add_axes((0.5, 0.5, 0.5, 0.5))
    ax3_root = fig.add_axes((0.5, 0, 0.5, 0.5))

    # Panel A
    logger.info("Plotting synteny dotplot")
    qbed, sbed, _, _, is_self = check_beds(anchorfile, p, opts)
    dotplot(
        anchorfile,
        qbed,
        sbed,
        fig,
        ax1_root,
        ax1_canvas,
        is_self=is_self,
        chrlw=0.5,
        sepcolor=set1[3],
    )

    # Panel B
    logger.info("Plotting karyotype plot")
    Karyotype(ax2_root, seqidsfile, layoutfile)

    # Panel C
    logger.info("Plotting synteny blocks")
    Synteny(fig, ax3_root, datafile, bedfile, blockslayoutfile, pad=0.1, vpad=0.03)
    draw_gene_legend(root, 0.69, 0.8, 0.34)

    labels = ((0.02, 0.95, "A"), (0.52, 0.95, "B"), (0.52, 0.45, "C"))
    panel_labels(root, labels)
    normalize_axes(root, ax1_root, ax2_root, ax3_root)

    image_name = "synteny.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def diversity(args):
    """
    %prog diversity pedigree.ped VAR?_srtd.wgs.regions.bed.gz

    Plot diversity composite figure, including:
    A. Pedigree
    B. Depth distribution across genomes
    """
    p = OptionParser(diversity.__doc__)
    _, args, iopts = p.set_image_options(args, figsize="14x7")

    if len(args) < 2:
        sys.exit(not p.print_help())

    pedfile, bedfiles = args[0], args[1:]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    ax1_root = fig.add_axes((0, 0, 0.25, 1))
    ax2_root = fig.add_axes((0.25, 0, 0.75, 1))

    # Panel A
    logger.info("Plotting pedigree")
    ped = Pedigree(pedfile)
    pngfile = f"{pedfile}.png"
    inb = calculate_inbreeding(ped, ploidy=4, N=10000)

    G = ped.to_graph(inb, title="Pedigree of Variety1")
    A = nx.nx_agraph.to_agraph(G)
    dpi = 300
    A.draw(pngfile, prog="dot", args=f"-Gdpi={dpi}")
    logger.info("Pedigree graph written to `%s`", pngfile)

    # Show the image as is
    ax1_root.imshow(load_image(pngfile))
    ax1_root.set_axis_off()

    # Panel B
    logger.info("Plotting depth distribution across genomes")
    npanels = len(bedfiles)
    yinterval = 1.0 / npanels
    ypos = 1 - yinterval
    panel_roots, panel_axes = [], []
    for _ in range(npanels):
        panel_root = fig.add_axes((0.25, ypos, 0.75, yinterval))
        panel_ax = fig.add_axes(
            (0.25 + 0.1 * 0.75, ypos + 0.2 * yinterval, 0.8 * 0.75, 0.65 * yinterval)
        )
        panel_roots.append(panel_root)
        panel_axes.append(panel_ax)
        ypos -= yinterval

    draw_multi_depth(
        ax2_root,
        panel_roots,
        panel_axes,
        bedfiles,
        chrinfo_file="chrinfo.txt",
        titleinfo_file="titleinfo.txt",
        maxdepth=100,
        logscale=False,
    )

    labels = (
        (0.02, 0.95, "A"),
        (0.25 + 0.25 * 0.1, 0.95, "B"),
    )
    panel_labels(root, labels)
    normalize_axes(root, ax2_root)

    image_name = "diversity.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def landscape(args):
    """
    %prog landscape features.bed athaliana.sizes TAIR10_chr_all.fas Chr2

    Plot landscape composite figure, including:
    A. Example genomic features painted on Arabidopsis genome
    B. Landscape of genomic features across the genome
    """
    p = OptionParser(landscape.__doc__)
    _, args, iopts = p.set_image_options(args, figsize="12x8")

    if len(args) != 4:
        sys.exit(not p.print_help())

    bedfile, sizesfile, fastafile, ch = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))
    aspect_ratio = iopts.w / iopts.h
    ax1_root = fig.add_axes((0, 1 / 4, 0.4, 0.5 * aspect_ratio))
    ax2_root_extent = (0.4, 0.5, 0.6, 0.47)
    ax2_root = fig.add_axes(ax2_root_extent)
    ax3_root_extent = (0.41, 0, 0.6, 0.47)
    ax3_root = fig.add_axes(ax3_root_extent)

    # Panel A
    logger.info("Plotting example genomic features painted on Arabidopsis genome")
    draw_chromosomes(
        ax1_root,
        bedfile,
        sizesfile,
        iopts=iopts,
        mergedist=0,
        winsize=50000,
        gauge=True,
        legend=True,
        empty=False,
        title="*Arabidopsis* genome features",
    )

    # Panel B
    logger.info("Plotting landscape of genomic features across the genome")
    stacks = ["Repeats", "Exons"]
    heatmaps = ["Copia", "Gypsy", "Helitron", "hAT", "Exons"]
    window = 250000
    shift = 50000
    draw_stacks(
        fig,
        ax2_root,
        ax2_root_extent,
        stacks,
        fastafile,
        window,
        shift,
        top=5,
    )

    # Panel C
    draw_heatmaps(
        fig,
        ax3_root,
        ax3_root_extent,
        fastafile,
        "Chr2",
        stacks,
        heatmaps,
        window,
        shift,
        cmap=cm.viridis,
    )

    ax2_root.set_axis_off()
    ax3_root.set_axis_off()

    labels = ((0.02, 0.95, "A"), (0.42, 0.95, "B"), (0.42, 0.48, "C"))
    panel_labels(root, labels)
    normalize_axes(root, ax1_root)

    image_name = "landscape.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def genomebuild(args):
    """
    %prog genomebuild reads.histo geneticmap.matrix hic.resolution_500000.npy hic.resolution_500000.json

    Plot genome build composite figure, including:
    A. Read kmer histogram
    B. Genetic map concordance
    C. Hi-C contact map concordance
    """
    p = OptionParser(genomebuild.__doc__)
    _, args, iopts = p.set_image_options(args, figsize="21x7")

    if len(args) != 4:
        sys.exit(not p.print_help())

    reads_histo, mstmap, hic_matrix, hic_json = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))
    ax1_root = fig.add_axes((0, 0, 1 / 3, 1))
    ax2_root = fig.add_axes((1 / 3, 0, 1 / 3, 1))
    ax3_root = fig.add_axes((2 / 3, 0, 1 / 3, 1))
    ax1 = fig.add_axes((1 / 3 * 0.1, 0.1, 1 / 3 * 0.8, 0.8))
    ax2 = fig.add_axes((1 / 3 * 1.1, 0.1, 1 / 3 * 0.8, 0.8))
    ax3 = fig.add_axes((1 / 3 * 2.1, 0.1, 1 / 3 * 0.8, 0.8))

    # Panel A
    logger.info("Plotting read kmer histogram")
    _ = draw_ks_histogram(
        ax1,
        reads_histo,
        method="nbinom",
        coverage=0,
        vmin=2,
        vmax=200,
        species="*S. species* ‘Variety 1’",
        K=21,
        maxiter=100,
        peaks=False,
    )

    # Panel B
    logger.info("Plotting genetic map concordance")
    draw_geneticmap_heatmap(ax2_root, ax2, mstmap, 1000)

    # Panel C
    logger.info("Plotting Hi-C contact map concordance")
    draw_hic_heatmap(
        ax3_root,
        ax3,
        hic_matrix,
        hic_json,
        contig=None,
        groups_file="groups",
        title="*S. species* Hi-C contact map",
        vmin=1,
        vmax=6,
        plot_breaks=True,
    )

    labels = (
        (1 / 3 * 0.1, 0.95, "A"),
        (1 / 3 * 1.1, 0.95, "B"),
        (1 / 3 * 2.1, 0.95, "C"),
    )
    panel_labels(root, labels)
    normalize_axes(root, ax1_root, ax2_root, ax3_root)

    image_name = "genomebuild.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def main():

    actions = (
        ("synteny", "Plot synteny composite figure"),
        ("diversity", "Plot diversity composite figure"),
        ("genomebuild", "Plot genome build composite figure"),
        ("landscape", "Plot landscape composite figure"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
