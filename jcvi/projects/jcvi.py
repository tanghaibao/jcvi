#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions in this script produce figures in the JCVI manuscript.
"""

from ..apps.base import ActionDispatcher, OptionParser, logger
from ..graphics.base import normalize_axes, panel_labels, plt, savefig


def genomebuild(args):
    """
    %prog genomebuild reads.histo geneticmap.matrix hic.resolution_500000.npy hic.resolution_500000.json

    Plot genome build composite figure, including:
    A. Read kmer histogram
    B. Genetic map concordance
    C. Hi-C contact map concordance
    """
    p = OptionParser(genomebuild.__doc__)
    _, args, iopts = p.set_image_options(args, figsize="15x5")

    if len(args) != 4:
        sys.exit(not p.print_help())

    reads_histo, geneticmap_matrix, hic_matrix, hic_json = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))
    ax1 = fig.add_axes((0.1, 0.1, 0.32, 0.8))
    ax2 = fig.add_axes((0.1, 0.1, 0.34, 0.8))
    ax3 = fig.add_axes((0.1, 0.1, 0.34, 0.8))

    # Panel A
    logger.info("Plotting read kmer histogram")

    # Panel B
    logger.info("Plotting genetic map concordance")

    # Panel C
    logger.info("Plotting Hi-C contact map concordance")

    labels = ((0.05, 0.95, "A"), (0.37, 0.95, "B"), (0.71, 0.95, "C"))
    panel_labels(root, labels)
    normalize_axes([root, ax1, ax2, ax3])

    image_name = "genomebuild.pdf"
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
