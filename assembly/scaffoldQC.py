#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging

from optparse import OptionParser

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, Rectangle, _
from jcvi.utils.cbook import thousands
from jcvi.apps.base import ActionDispatcher, debug
debug()


def scaffolding(ax, scaffoldID, blastf, qsizes, ssizes, qbed, sbed):

    from jcvi.graphics.blastplot import blastplot

    # qsizes, qbed are properties for the evidences
    # ssizes, sbed are properties for the current scaffoldID
    blastplot(ax, blastf, qsizes, ssizes, qbed, sbed, \
              style="circle", insetLabels=True)

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


def main():

    actions = (
        ('bed', 'make evidence alignments in the appropriate bed formats'),
        ('plot', 'plot the alignment of the scaffold to other evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed oldbedfile blastfile newbedfile

    Make evidence alignments in the appropriate bed formats.

    `newbedfile` contains lines that look like,

    scfxxxxxx 100 200 ctg3234 (physical map contig)
    scfxxxxxx 100 200 RAPAscfxxxxx (scaffold in a related organism)

    which really involves converting blast result, associating to a different
    scaffold or molecule through `oldbedfile`, and output `newbedfile`.
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, blastfile, newbedfile = args
    border = Bed(bedfile).order
    fw = open(newbedfile, "w")
    nbed = Bed()

    blast = Blast(blastfile)
    for b in blast.iter_line():
        query, subject = b.query, b.subject
        sstart, sstop = b.sstart, b.sstop
        if query not in border:
            continue
        qi, q = border[query]
        ctg = q.seqid

        bline = "\t".join(str(x) for x in (subject, sstart, sstop, ctg))
        nbed.append(BedLine(bline))

    nbed.sort(key=nbed.key)
    nbed.print_to_file(fw=fw)
    fw.close()


def plot(args):
    """
    %prog plot scaffold.fasta synteny.blast synteny.sizes syteny.bed
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

    p = OptionParser(plot.__doc__)
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


if __name__ == '__main__':
    main()
