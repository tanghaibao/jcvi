#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog scaffold.sizes synteny.bed physicalmap.bed

As evaluation of scaffolding, visualize external line of evidences:
* Plot synteny to an external genome
* Plot alignments to physical map
* Plot alignments to genetic map (TODO)

For example, physicalmap.bed can look like
scf7180014794474 104 104 ctg991

and will plot a dot in the dot plot in the corresponding location
the plots are one contig/scaffold per plot, but descend down the list in
scaffolds.sizes, until the contig is smaller than the cutoff
"""

import sys
import logging

from optparse import OptionParser

from jcvi.graphics.base import plt, Rectangle, _, human_size_formatter
from jcvi.utils.cbook import thousands
from jcvi.formats.bed import Bed
from jcvi.formats.sizes import Sizes
from jcvi.apps.base import debug
debug()


def scaffolding(ax, scaffoldID, scafsize, bedf):
    bed = list(bedf.sub_bed(scaffoldID)) 
    # order of the seqids of query, reference is always scaffoldID
    seen = set() 
    for b in bed:
        accn = b.accn
        if accn in seen: continue
        seen.add(accn)

    ctgIDs = dict((x, i+1) for i, x in enumerate(seen))
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
    ctglabel = "contig" if fname=="FPC" else "scaffold"
    xtitle = "Each vertical position is {0} {1}".format(fname, ctglabel)
    ax.set_xlabel(_(xtitle), color="g")
    

def plot_one_scaffold(scaffoldID, scafsize, beds, imagename):
    nbeds = len(beds)
    fig = plt.figure(1, (10, 10))
    plt.cla()
    plt.clf()
    root = fig.add_axes([0, 0, 1, 1])
    axes = [fig.add_subplot(1, nbeds, x) for x in range(1, nbeds+1)] 

    formatter = human_size_formatter
    for bed, ax in zip(beds, axes):
        bed = Bed(bed)
        scaffolding(ax, scaffoldID, scafsize, bed)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_xticklabels([])

    root.text(.2, .95, _("{0}   (size={1})".format(scaffoldID, thousands(scafsize))), 
            size=18, color='b')
    root.text(.5, .05, _("Look for vertical tracks"), color="gray", ha="center")
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(imagename, dpi=150)
    logging.debug("saving image to `{0}`".format(imagename))


if __name__ == '__main__':
    
    p = OptionParser(__doc__)

    p.add_option("--cutoff", dest="cutoff", default=5*1e6,
            help="plot for contigs and scaffolds larger than [default: %default]")
    p.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    opts, args = p.parse_args()

    if len(args) < 3:
        sys.exit(p.print_help())

    scafsizes = Sizes(args[0])
    beds = args[1:]

    for scaffoldID, scafsize in scafsizes.iter_sizes():
        logging.debug("Loading {0} (size={1})".format(scaffoldID,
            thousands(scafsize)))
        if scafsize < opts.cutoff: break
        imagename = ".".join((scaffoldID, opts.format))
        plot_one_scaffold(scaffoldID, scafsize, beds, imagename)
