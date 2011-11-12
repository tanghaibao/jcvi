#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Assembly QC plots, including general statistics, base and mate coverages, and
scaffolding consistencies.
"""

import sys
import logging
import os.path as op

from optparse import OptionParser

from jcvi.formats.fasta import Fasta
from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.sizes import Sizes
from jcvi.assembly.base import calculate_A50
from jcvi.assembly.coverage import BedLine, Sizes, Coverage
from jcvi.algorithms.matrix import moving_average
from jcvi.graphics.base import plt, Rectangle, set_tex_axis, _
from jcvi.utils.cbook import thousands
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():
    actions = (
        ('A50', 'compare A50 graphics for a set of FASTA files'),
        ('coverage', 'performs QC graphics on given contig/scaffold'),
        ('scaffold', 'plot the alignment of the scaffold to other evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def scaffolding(ax, scaffoldID, blastf, qsizes, ssizes, qbed, sbed):

    from jcvi.graphics.blastplot import blastplot

    # qsizes, qbed are properties for the evidences
    # ssizes, sbed are properties for the current scaffoldID
    blastplot(ax, blastf, qsizes, ssizes, qbed, sbed, \
              style="circle", insetLabels=True, stripNames=True)

    # FPC_scf.bed => FPC
    fname = qbed.filename.split(".")[0].split("_")[0]
    xtitle = fname
    if xtitle == "FPC":
        ax.set_xticklabels([""] * len(ax.get_xticklabels()))
    ax.set_xlabel(_(xtitle), color="g")
    for x in ax.get_xticklines():
        x.set_visible(False)


def plot_one_scaffold(scaffoldID, ssizes, sbed, trios, imagename, iopts):
    ntrios = len(trios)
    fig = plt.figure(1, (14, 8))
    plt.cla()
    plt.clf()
    root = fig.add_axes([0, 0, 1, 1])
    axes = [fig.add_subplot(1, ntrios, x) for x in range(1, ntrios + 1)]
    scafsize = ssizes.get_size(scaffoldID)

    for trio, ax in zip(trios, axes):
        blastf, qsizes, qbed = trio
        scaffolding(ax, scaffoldID, blastf, qsizes, ssizes, qbed, sbed)

    root.text(.5, .95, _("{0}   (size={1})".\
            format(scaffoldID, thousands(scafsize))),
            size=18, ha="center", color='b')
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(imagename, dpi=iopts.dpi)
    logging.debug("Print image to `{0}` {1}".format(imagename, iopts))


def scaffold(args):
    """
    %prog scaffold scaffold.fasta synteny.blast synteny.sizes synteny.bed
                         physicalmap.blast physicalmap.sizes physicalmap.bed

    As evaluation of scaffolding, visualize external line of evidences:
    * Plot synteny to an external genome
    * Plot alignments to physical map
    * Plot alignments to genetic map (TODO)

    Each trio defines one panel to be plotted. blastfile defines the matchings
    between the evidences vs scaffolds. Then the evidence sizes, and evidence
    bed to plot dot plots.

    This script will plot a dot in the dot plot in the corresponding location
    the plots are one contig/scaffold per plot.
    """
    from jcvi.graphics.base import set_image_options
    from jcvi.utils.iter import grouper

    p = OptionParser(scaffold.__doc__)
    p.add_option("--cutoff", type="int", default=1000000,
            help="Plot scaffolds with size larger than [default: %default]")
    opts, args, iopts = set_image_options(p, args, figsize="14x8", dpi=150)

    if len(args) < 4 or len(args) % 3 != 1:
        sys.exit(not p.print_help())

    scafsizes = Sizes(args[0])
    trios = list(grouper(3, args[1:]))
    trios = [(a, Sizes(b), Bed(c)) for a, b, c in trios]

    for scaffoldID, scafsize in scafsizes.iter_sizes():
        if scafsize < opts.cutoff:
            continue
        logging.debug("Loading {0} (size={1})".format(scaffoldID,
            thousands(scafsize)))

        tmpname = scaffoldID + ".sizes"
        tmp = open(tmpname, "w")
        tmp.write("{0}\t{1}".format(scaffoldID, scafsize))
        tmp.close()

        tmpsizes = Sizes(tmpname)
        tmpsizes.close(clean=True)

        imagename = ".".join((scaffoldID, opts.format))
        plot_one_scaffold(scaffoldID, tmpsizes, None, trios, imagename, iopts)


def coverage(args):
    """
    %prog coverage prefix

    Expects data files including:
    1. `prefix.bedpe` draws Bezier curve between paired reads
    2. `prefix.sizes` draws length of the contig/scaffold
    3. `prefix.gaps.bed` mark the position of the gaps in sequence
    4. `prefix.bed.coverage` plots the base coverage
    5. `prefix.pairs.bed.coverage` plots the clone coverage

    See assembly.coverage.posmap() for the generation of these files.
    """
    p = OptionParser(coverage.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    prefix, = args
    scf = prefix

    # All these files *must* be present in the current folder
    bedpefile = prefix + ".bedpe"
    fastafile = prefix + ".fasta"
    sizesfile = prefix + ".sizes"
    gapsbedfile = prefix + ".gaps.bed"
    bedfile = prefix + ".bed"
    bedpefile = prefix + ".bedpe"
    pairsbedfile = prefix + ".pairs.bed"

    sizes = Sizes(fastafile).mapping
    size = sizes[scf]
    basecoverage = Coverage(bedfile, sizesfile)
    matecoverage = Coverage(pairsbedfile, sizesfile)

    import numpy as np
    from jcvi.graphics.glyph import Bezier

    fig = plt.figure(1, (8, 5))
    root = fig.add_axes([0, 0, 1, 1])

    # the scaffold
    root.add_patch(Rectangle((.1, .15), .8, .03, fc='k'))

    # basecoverage and matecoverage
    ax = fig.add_axes([.1, .45, .8, .45])
    bases = np.arange(1, size + 1)

    assert len(bases) == len(basecoverage)
    assert len(bases) == len(matecoverage)

    window = size / 200  # smooth the curve
    logging.debug("Coverage curve use window size of {0} bases.".format(window))
    basecoverage = moving_average(basecoverage, window=window)
    matecoverage = moving_average(matecoverage, window=window)

    baseline = ax.plot(bases, basecoverage, 'g-')
    mateline = ax.plot(bases, matecoverage, 'r-')
    legends = (_("Base coverage"), _("Mate coverage"))
    leg = ax.legend((baseline, mateline), legends, shadow=True, fancybox=True)
    leg.get_frame().set_alpha(.5)
    ax.set_xlim(0, size)

    # draw the read pairs
    fp = open(bedpefile)
    pairs = []
    for row in fp:
        scf, astart, aend, scf, bstart, bend, clonename = row.split()
        astart, bstart = int(astart), int(bstart)
        aend, bend = int(aend), int(bend)
        start = min(astart, bstart) + 1
        end = max(aend, bend)
        pairs.append((start, end))

    bpratio = .8 / size
    cutoff = 1000  # inserts smaller than this are not plotted
    # this convert from base => x-coordinate
    pos = lambda x: (.1 + x * bpratio)
    ypos = .15 + .03
    for start, end in pairs:
        dist = end - start

        if dist < cutoff:
            continue

        dist = min(dist, 10000)
        # 10Kb == .25 canvas height
        height = .25 * dist / 10000
        xstart = pos(start)
        xend = pos(end)
        p0 = (xstart, ypos)
        p1 = (xstart, ypos + height)
        p2 = (xend, ypos + height)
        p3 = (xend, ypos)
        Bezier(root, p0, p1, p2, p3)

    # gaps on the scaffold
    fp = open(gapsbedfile)
    for row in fp:
        b = BedLine(row)
        start, end = b.start, b.end
        xstart = pos(start)
        xend = pos(end)
        root.add_patch(Rectangle((xstart, .15), xend - xstart, .03, fc='w'))

    root.text(.5, .1, _(scf), color='b', ha="center")
    warn_msg = "Only the inserts > {0}bp are shown".format(cutoff)
    root.text(.5, .1, _(scf), color='b', ha="center")
    root.text(.5, .05, _(warn_msg), color='gray', ha="center")
    # clean up and output
    set_tex_axis(ax)
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = prefix + ".pdf"
    plt.savefig(figname, dpi=300)
    logging.debug("Figure saved to `{0}`".format(figname))


def generate_plot(filename, rplot="A50.rplot", rpdf="A50.pdf"):

    from jcvi.apps.r import RTemplate

    rplot_template = """
    library(ggplot2)

    data <- read.table("$rplot", header=T, sep="\t")
    g <- ggplot(data, aes(x=index, y=cumsize, group=fasta))
    g + geom_line(aes(colour=fasta)) +
    xlab("Contigs") + ylab("Cumulative size (Mb)") +
    opts(title="A50 plot", legend.position="top")

    ggsave(file="$rpdf")
    """

    rtemplate = RTemplate(rplot_template, locals())
    rtemplate.run()


def A50(args):
    """
    %prog A50 contigs_A.fasta contigs_B.fasta ...

    Plots A50 graphics, see blog post (http://blog.malde.org/index.php/a50/)
    """
    p = OptionParser(A50.__doc__)
    p.add_option("--overwrite", default=False, action="store_true",
            help="overwrite .rplot file if exists [default: %default]")
    p.add_option("--cutoff", default=0, type="int", dest="cutoff",
            help="use contigs above certain size [default: %default]")
    p.add_option("--stepsize", default=10, type="int", dest="stepsize",
            help="stepsize for the distribution [default: %default]")
    opts, args = p.parse_args(args)

    if not args:
        sys.exit(p.print_help())

    import numpy as np
    from jcvi.utils.table import loadtable

    stepsize = opts.stepsize  # use stepsize to speed up drawing
    rplot = "A50.rplot"
    if not op.exists(rplot) or opts.overwrite:
        fw = open(rplot, "w")
        header = "\t".join(("index", "cumsize", "fasta"))
        statsheader = ("Fasta", "L50", "N50", "Min", "Max", "Average", "Sum",
                "Counts")
        statsrows = []
        print >>fw, header
        for fastafile in args:
            f = Fasta(fastafile, index=False)
            ctgsizes = [length for k, length in f.itersizes()]
            ctgsizes = np.array(ctgsizes)

            a50, l50, n50 = calculate_A50(ctgsizes, cutoff=opts.cutoff)
            cmin, cmax, cmean = min(ctgsizes), max(ctgsizes), np.mean(ctgsizes)
            csum, counts = np.sum(ctgsizes), len(ctgsizes)
            cmean = int(round(cmean))
            statsrows.append((fastafile, l50, n50, cmin, cmax, cmean, csum,
                counts))

            logging.debug("`{0}` ctgsizes: {1}".format(fastafile, ctgsizes))

            tag = "{0} (L50={1})".format(\
                    op.basename(fastafile).rsplit(".", 1)[0], l50)
            logging.debug(tag)

            for i, s in zip(xrange(0, len(a50), stepsize), a50[::stepsize]):
                print >> fw, "\t".join((str(i), str(s / 1000000.), tag))
        fw.close()

        table = loadtable(statsheader, statsrows)
        print >> sys.stderr, table

    generate_plot(rplot)


if __name__ == '__main__':
    main()
