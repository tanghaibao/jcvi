#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Process Hi-C output into AGP for chromosomal-scale scaffolding.
"""

import logging
import sys
import os.path as op

from jcvi.apps.base import OptionParser, ActionDispatcher, iglob
from jcvi.utils.natsort import natsorted
from jcvi.formats.agp import order_to_agp
from jcvi.formats.base import LineFile, must_open
from jcvi.formats.sizes import Sizes


class ContigOrderingLine(object):
    '''Stores one line in the ContigOrdering file
    '''
    def __init__(self, line, sep="|"):
        args = line.split()
        self.contig_id = args[0]
        self.contig_name = args[1].split(sep)[0]
        contig_rc = args[2]
        assert contig_rc in ('0', '1')
        self.strand = '+' if contig_rc == '0' else '-'
        self.orientation_score = args[3]
        self.gap_size_after_contig = args[4]


class ContigOrdering(LineFile):
    '''ContigOrdering file as created by LACHESIS, one per chromosome group.
    Header contains summary information per group, followed by list of contigs
    with given ordering.
    '''
    def __init__(self, filename):
        super(ContigOrdering, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row[0] == '#':
                continue
            orderline = ContigOrderingLine(row)
            self.append(orderline)

    def write_agp(self, obj, sizes, fw=sys.stdout, gapsize=100,
                  gaptype="contig", evidence="map"):
        '''Converts the ContigOrdering file into AGP format
        '''
        contigorder = [(x.contig_name, x.strand) for x in self]
        order_to_agp(obj, contigorder, sizes, fw,
                     gapsize=gapsize, gaptype=gaptype, evidence=evidence)


def main():

    actions = (
        ('agp', 'generate AGP file based on LACHESIS output'),
        ('heatmap', 'generate heatmap based on LACHESIS output'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def heatmap(args):
    """
    %prog heatmap main_results/ cached_data/ contigs.fasta

    Generate heatmap file based on LACHESISS output.
    """
    import math
    import numpy as np

    p = OptionParser(heatmap.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    mdir, cdir, contigsfasta = args
    orderingfiles = natsorted(iglob(mdir, "*.ordering"))
    sizes = Sizes(contigsfasta)
    start = 0
    bins = {}
    contig_ids = dict((name, str(i)) for (i, name) in \
                            enumerate(sizes.iter_names()))
    breaks = []
    for ofile in orderingfiles:
        co = ContigOrdering(ofile)
        for x in co:
            size = sizes.mapping[x.contig_name]
            contig_id = contig_ids[x.contig_name]
            end = start + int(math.ceil(size / 100000.))
            bins[contig_id] = (start, end)
            start = end
        breaks.append(start)

    totalbins = start
    M = np.zeros((totalbins, totalbins))
    glm = op.join(cdir, "all.GLM")
    for row in open(glm):
        if row[0] == '#':
            continue
        x, y, z = row.split()
        if x == 'X':
            continue
        if x not in bins or y not in bins:
            continue
        xstart, xend = bins[x]
        ystart, yend = bins[y]
        #z = float(z) / ((xend - xstart) * (yend - ystart))
        z = float(z)
        M[xstart:xend, ystart:yend] = z

    M = np.log10(M + 1)

    from jcvi.graphics.base import normalize_axes, plt, savefig

    fig = plt.figure(figsize=(iopts.w, iopts.h))
    plt.imshow(M, cmap=iopts.cmap, origin="lower", interpolation='none')
    ax = plt.gca()
    xlim = ax.get_xlim()
    for b in breaks[:-1]:
        plt.plot([b, b], xlim, 'w-')
        plt.plot(xlim, [b, b], 'w-')

    ax.set_xlim(xlim)
    ax.set_ylim(xlim)

    root = fig.add_axes([0, 0, 1, 1])
    normalize_axes(root)
    image_name = "heatmap." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def agp(args):
    """
    %prog agp main_results/ contigs.fasta

    Generate AGP file based on LACHESIS output.
    """
    p = OptionParser(agp.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    odir, contigsfasta = args
    fwagp = must_open(opts.outfile, 'w')
    orderingfiles = natsorted(iglob(odir, "*.ordering"))
    sizes = Sizes(contigsfasta).mapping
    contigs = set(sizes.keys())
    anchored = set()

    for ofile in orderingfiles:
        co = ContigOrdering(ofile)
        anchored |= set([x.contig_name for x in co])
        obj = op.basename(ofile).split('.')[0]
        co.write_agp(obj, sizes, fwagp)

    singletons = contigs - anchored
    logging.debug('Anchored: {}, Singletons: {}'.\
                  format(len(anchored), len(singletons)))

    for s in natsorted(singletons):
        order_to_agp(s, [(s, "?")], sizes, fwagp)


if __name__ == '__main__':
    main()
