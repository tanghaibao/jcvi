#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog scaffoldID scaffold.sizes synteny.bed physicalmap.bed

As evaluation of scaffolding, visualize external line of evidences:
* Plot synteny to an external genome
* Plot alignments to physical map
* Plot alignments to genetic map (TODO)

For example, physicalmap.bed can look like
scf7180014794474 104 104 ctg991

and will plot a dot in the dot plot in the corresponding location
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


def scaffolding(ax, scaffoldID, scafsize, bed):
    bed = list(bed.sub_bed(scaffoldID)) 
    # order of the seqids of query, reference is always scaffoldID
    seen = {}
    for b in bed:
        accn = b.accn
        if accn in seen: continue
        seen[accn] = len(seen) + 1

    ctgIDs = [k for k, v in sorted(seen.items(), key=lambda x: x[1])]
    data = []
    for b in bed:
        y = (b.start + b.end) / 2
        x = seen[b.accn]
        data.append((x, y))

    x, y = zip(*data)
    ax.plot(x, y, "k.", ms=1)
    xlim = len(ctgIDs) + 1
    ylim = scafsize
    ax.set_xlim(0, xlim)
    ax.set_ylim(scafsize, 0)
    

if __name__ == '__main__':
    
    p = OptionParser(__doc__)

    p.add_option("--format", dest="format", default="png",
            help="generate image of format (png, pdf, ps, eps, svg, etc.)"
            "[default: %default]")

    opts, args = p.parse_args()

    if len(args) < 3:
        sys.exit(p.print_help())

    scaffoldID, scafsizes = args[:2]
    scafsize = Sizes(scafsizes).get_size(scaffoldID)
    logging.debug("Loading {0} (size={1})".format(scaffoldID,
        thousands(scafsize)))

    beds = args[2:]
    nbeds = len(beds)

    fig = plt.figure(1, (10, 10))
    root = fig.add_axes([0, 0, 1, 1])
    axes = [fig.add_subplot(1, nbeds, x) for x in range(1, nbeds+1)] 

    formatter = human_size_formatter
    for bed, ax in zip(beds, axes):
        bed = Bed(bed)
        scaffolding(ax, scaffoldID, scafsize, bed)
        ax.yaxis.set_major_formatter(formatter)
        ax.set_xticklabels([])

    root.text(.05, .95, _("{0} size={1}".format(scaffoldID, thousands(scafsize))), color='b')
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    imagename = ".".join((scaffoldID, opts.format))
    plt.savefig(imagename, dpi=150)
    logging.debug("saving image to `{0}`".format(imagename))

