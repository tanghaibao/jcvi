#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the JCVI manuscript.
"""

import sys

from ..apps.base import ActionDispatcher, OptionParser, logger
from ..assembly.geneticmap import draw_geneticmap_heatmap
from ..assembly.hic import draw_hic_heatmap
from ..assembly.kmer import draw_ks_histogram
from ..graphics.base import normalize_axes, panel_labels, plt, savefig
from ..graphics.landscape import draw_multi_depth


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
    normalize_axes([root, ax1_root, ax2_root, ax3_root])

    image_name = "genomebuild.pdf"
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

    # Panel B
    logger.info("Plotting depth distribution across genomes")
    npanels = len(bedfiles)
    yinterval = 1.0 / npanels
    ypos = 1 - yinterval
    panel_roots, panel_axes = [], []
    for _ in range(npanels):
        panel_root = root if npanels == 1 else fig.add_axes((0.25, ypos, 1, yinterval))
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
        (0.25 * 0.1, 0.95, "A"),
        (0.25 + 0.25 * 0.1, 0.95, "B"),
    )
    panel_labels(root, labels)
    normalize_axes([root, ax1_root, ax2_root])

    image_name = "diversity.pdf"
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def main():

    actions = (
        ("genomebuild", "Plot genome build composite figure"),
        ("diversity", "Plot diversity composite figure"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
