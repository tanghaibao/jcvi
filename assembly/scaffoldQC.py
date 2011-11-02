#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging

from optparse import OptionParser

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, Rectangle, _, human_size_formatter
from jcvi.utils.cbook import thousands
from jcvi.apps.base import ActionDispatcher, debug
debug()


def scaffolding(ax, scaffoldID, scafsize, bedf):
    bed = list(bedf.sub_bed(scaffoldID))
    # order of the seqids of query, reference is always scaffoldID
    seen = set()
    for b in bed:
        accn = b.accn
        if accn in seen:
            continue
        seen.add(accn)

    ctgIDs = dict((x, i + 1) for i, x in enumerate(seen))
    data = []
    for b in bed:
        y = (b.start + b.end) / 2
        x = ctgIDs[b.accn]
        data.append((x, y))

    x, y = zip(*data)
    ax.plot(x, y, "k.", alpha=.75, ms=3)
    xlim = len(ctgIDs) + 1
    ylim = scafsize
    ax.set_xlim(0, xlim)
    ax.set_ylim(scafsize, 0)

    # FPC_scf.bed => FPC
    fname = bedf.filename.split(".")[0].split("_")[0]
    ctglabel = "contig" if fname == "FPC" else "scaffold"
    xtitle = "Each position is one {0} {1}".format(fname, ctglabel)
    ax.set_xlabel(_(xtitle), color="g")
    for x in ax.get_xticklines():
        x.set_visible(False)


def plot_one_scaffold(scaffoldID, scafsize, beds, imagename):
    nbeds = len(beds)
    fig = plt.figure(1, (9, 9))
    plt.cla()
    plt.clf()
    root = fig.add_axes([0, 0, 1, 1])
    axes = [fig.add_subplot(1, nbeds, x) for x in range(1, nbeds + 1)]

    formatter = human_size_formatter
    for bed, ax in zip(beds, axes):
        bed = Bed(bed)
        scaffolding(ax, scaffoldID, scafsize, bed)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_xticklabels([])

    root.text(.5, .95, _("{0}   (size={1})".\
            format(scaffoldID, thousands(scafsize))),
            size=18, ha="center", color='b')
    root.text(.5, .05, _("Look for vertical tracks"),
            color="gray", ha="center")
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(imagename, dpi=150)
    logging.debug("saving image to `{0}`".format(imagename))


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
    %prog plot scaffold.fasta synteny.bed physicalmap.bed

    As evaluation of scaffolding, visualize external line of evidences:
    * Plot synteny to an external genome
    * Plot alignments to physical map
    * Plot alignments to genetic map (TODO)

    For example, `physicalmap.bed` can look like
    scf7180014794474 104 104 ctg991

    and will plot a dot in the dot plot in the corresponding location
    the plots are one contig/scaffold per plot, but descend down the list in
    scaffolds.sizes, until the contig is smaller than the cutoff
    """
    from jcvi.graphics.base import set_format

    p = OptionParser(plot.__doc__)
    p.add_option("--cutoff", dest="cutoff", default=1000000,
            help="plot for contigs and scaffolds > [default: %default]")
    set_format(p)

    opts, args = p.parse_args(args)

    if len(args) < 3:
        sys.exit(not p.print_help())

    scafsizes = Sizes(args[0])
    beds = args[1:]

    for scaffoldID, scafsize in scafsizes.iter_sizes():
        if scafsize < opts.cutoff:
            continue
        logging.debug("Loading {0} (size={1})".format(scaffoldID,
            thousands(scafsize)))
        imagename = ".".join((scaffoldID, opts.format))
        plot_one_scaffold(scaffoldID, scafsize, beds, imagename)


if __name__ == '__main__':
    main()
