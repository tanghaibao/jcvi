#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Process Hi-C output into AGP for chromosomal-scale scaffolding.
"""

import logging
import sys
import os.path as op
import numpy as np
import math

from functools import partial

from jcvi.apps.base import OptionParser, ActionDispatcher, iglob, mkdir
from jcvi.utils.natsort import natsorted
from jcvi.formats.agp import order_to_agp
from jcvi.formats.base import LineFile, must_open
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import normalize_axes, plt, savefig


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
        ('score', 'score the current LACHESIS CLM'),
        # Plotting
        ('heatmap', 'generate heatmap based on LACHESIS output'),
        ('animation', 'plot optimization history'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def animation(args):
    """
    %prog animation tour cached_data/ contigsfasta

    Plot optimization history.
    """
    p = OptionParser(animation.__doc__)
    p.add_option("--frames", default=500, type="int",
                 help="Only plot every N frames")
    opts, args, iopts = p.set_image_options(args, figsize="8x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    tourfile, cdir, contigsfasta = args
    fp = open(tourfile)
    odir = "animation"
    mkdir(odir)
    for row in fp:
        if row[0] == '>':
            label = row[1:].strip()
            pf, i, score = label.split("-")
            i = int(i)
            continue
        else:
            if i % opts.frames != 0:
                continue
            tour = ",".join(row.split())
            image_name = op.join(odir, label) + ".png"
            heatmap(["main_results", cdir, contigsfasta,
                     "--tour", tour, "--outfile", image_name])


def score(args):
    """
    %prog score main_results/ cached_data/ contigsfasta

    Score the current LACHESIS CLM.
    """
    from jcvi.algorithms.ec import GA_setup, GA_run

    p = OptionParser(score.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    mdir, cdir, contigsfasta = args
    orderingfiles = natsorted(iglob(mdir, "*.ordering"))
    sizes = Sizes(contigsfasta)
    contig_names = list(sizes.iter_names())
    contig_ids = dict((name, i) for (i, name) in enumerate(contig_names))

    oo = []
    # Load contact matrix
    glm = op.join(cdir, "all.GLM")
    N = len(contig_ids)
    M = np.zeros((N, N), dtype=int)
    fp = open(glm)
    for row in fp:
        if row[0] == '#':
            continue
        x, y, z = row.split()
        if x == 'X':
            continue
        M[int(x), int(y)] = int(z)

    fwtour = open("tour", "w")
    def callback(tour, gen, oo):
        fitness = tour.fitness if hasattr(tour, "fitness") else None
        label = "GA-{0}".format(gen)
        if fitness:
            fitness = "{0}".format(fitness).split(",")[0].replace("(", "")
            label += "-" + fitness
        print >> fwtour, ">" + label
        print >> fwtour, " ".join(contig_names[oo[x]] for x in tour)
        return tour

    for ofile in orderingfiles:
        co = ContigOrdering(ofile)
        for x in co:
            contig_id = contig_ids[x.contig_name]
            oo.append(contig_id)
        pf = op.basename(ofile).split(".")[0]
        print pf
        print oo

        tour, tour_sizes, tour_M = prepare_ec(oo, sizes, M)

        # Faster Cython version for evaluation
        from .chic import score_evaluate
        callbacki = partial(callback, oo=oo)
        toolbox = GA_setup(tour)
        toolbox.register("evaluate", score_evaluate,
                         tour_sizes=tour_sizes, tour_M=tour_M)
        tour, tour.fitness = GA_run(toolbox, npop=100, cpus=opts.cpus,
                                    callback=callbacki)
        print tour, tour.fitness
        break

    fwtour.close()


def prepare_ec(oo, sizes, M):
    """
    This prepares EC and converts from contig_id to an index.
    """
    tour = range(len(oo))
    tour_sizes = np.array([sizes.sizes[x] for x in oo])
    tour_M = M[oo, :][:, oo]
    return tour, tour_sizes, tour_M


def score_evaluate(tour, tour_sizes=None, tour_M=None):
    sizes_oo = np.array([tour_sizes[x] for x in tour])
    sizes_cum = np.cumsum(sizes_oo) - sizes_oo / 2
    s = 0
    size = len(tour)
    for ia in xrange(size):
        a = tour[ia]
        for ib in xrange(ia + 1, size):
            b = tour[ib]
            links = tour_M[a, b]
            dist = sizes_cum[ib] - sizes_cum[ia]
            if dist > 1e7:
                break
            s += links * 1. / dist
    return s,


def heatmap(args):
    """
    %prog heatmap main_results/ cached_data/ contigs.fasta

    Generate heatmap file based on LACHESISS output.
    """
    p = OptionParser(heatmap.__doc__)
    p.add_option("--tour", help="List of contigs separated by comma")
    p.set_outfile(outfile=None)
    opts, args, iopts = p.set_image_options(args, figsize="8x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    mdir, cdir, contigsfasta = args
    orderingfiles = natsorted(iglob(mdir, "*.ordering"))
    sizes = Sizes(contigsfasta)
    contig_names = list(sizes.iter_names())
    contig_ids = dict((name, i) for (i, name) in enumerate(contig_names))

    if opts.tour:
        tours = [opts.tour.split(",")]
    else:
        tours = []
        for ofile in orderingfiles:
            co = ContigOrdering(ofile)
            tour = [x.contig_name for x in co]
            tours.append(tour)

    totalbins, bins, breaks = make_bins(tours, sizes, contig_ids)

    glm = op.join(cdir, "all.GLM")
    M = read_glm(glm, totalbins, bins)

    image_name = opts.outfile or ("heatmap." + iopts.format)
    plot_heatmap(M, breaks, iopts, image_name)


def make_bins(tours, sizes, contig_ids):
    breaks = []
    start = 0
    bins = {}
    for tour in tours:
        for x in tour:
            size = sizes.mapping[x]
            contig_id = contig_ids[x]
            end = start + int(math.ceil(size / 100000.))
            bins[contig_id] = (start, end)
            start = end
        breaks.append(start)

    totalbins = start
    return totalbins, bins, breaks


def read_glm(glm, totalbins, bins):
    M = np.zeros((totalbins, totalbins))
    for row in open(glm):
        if row[0] == '#':
            continue
        x, y, z = row.split()
        if x == 'X':
            continue
        x, y = int(x), int(y)
        if x not in bins or y not in bins:
            continue
        xstart, xend = bins[x]
        ystart, yend = bins[y]
        #z = float(z) / ((xend - xstart) * (yend - ystart))
        z = float(z)
        M[xstart:xend, ystart:yend] = z

    M = np.log10(M + 1)
    return M


def plot_heatmap(M, breaks, iopts, image_name):
    fig = plt.figure(figsize=(iopts.w, iopts.h))
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.imshow(M, cmap=iopts.cmap, origin="lower", interpolation='none')
    xlim = ax.get_xlim()
    for b in breaks[:-1]:
        ax.plot([b, b], xlim, 'w-')
        ax.plot(xlim, [b, b], 'w-')

    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    label = op.basename(image_name).rsplit(".", 1)[0]
    ax.set_title(label)

    root = fig.add_axes([0, 0, 1, 1])
    normalize_axes(root)
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
