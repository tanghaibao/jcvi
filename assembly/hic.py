#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Process Hi-C output into AGP for chromosomal-scale scaffolding.
"""

import logging
import sys
import os
import os.path as op
import numpy as np
import math

from collections import defaultdict
from functools import partial

from jcvi.apps.base import OptionParser, ActionDispatcher, iglob, mkdir, symlink
from jcvi.assembly.allmaps import make_movie
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


class CLMFile:
    '''CLM file (modified) has the following format:

    tig00046211+ tig00063795+       1       53173
    tig00046211+ tig00063795-       1       116050
    tig00046211- tig00063795+       1       71155
    tig00046211- tig00063795-       1       134032
    tig00030676+ tig00077819+       7       136407 87625 87625 106905 102218 169660 169660
    tig00030676+ tig00077819-       7       126178 152952 152952 35680 118923 98367 98367
    tig00030676- tig00077819+       7       118651 91877 91877 209149 125906 146462 146462
    tig00030676- tig00077819-       7       108422 157204 157204 137924 142611 75169 75169
    '''
    def __init__(self, clmfile, idsfile, skiprecover=True):
        self.parse_ids(idsfile, skiprecover)
        self.parse_clm(clmfile)

    def parse_ids(self, idsfile, skiprecover):
        '''IDS file has a list of contigs that need to be ordered. 'recover',
        keyword, if available in the third column, is less confident.

        tig00015093     46912
        tig00035238     46779   recover
        tig00030900     119291
        '''
        logging.debug("Parse idsfile `{}`".format(idsfile))
        fp = open(idsfile)
        tigs = []
        for row in fp:
            atoms = row.split()
            tig, size = atoms[:2]
            if skiprecover and len(atoms) == 3 and atoms[2] == 'recover':
                continue
            tigs.append((tig, size))

        # Mapping tig names to their indices and sizes
        tig_to_idx = {}
        tig_to_size = {}
        for i, (tig, size) in enumerate(tigs):
            tig_to_idx[tig] = i
            tig_to_size[tig] = size

        self.tig_to_idx = tig_to_idx
        self.tig_to_size = tig_to_size
        self.ntigs = len(tig_to_idx)

    def parse_clm(self, clmfile):
        logging.debug("Parse clmfile `{}`".format(clmfile))
        fp = open(clmfile)
        contacts = defaultdict(list)
        for row in fp:
            atoms = row.strip().split('\t')
            assert len(atoms) == 3, "Malformed line `{}`".format(atoms)
            abtig, links, dists = atoms
            atig, btig = abtig.split()
            at, ao = atig[:-1], atig[-1]
            bt, bo = btig[:-1], btig[-1]
            if at not in self.tig_to_idx:
                continue
            if bt not in self.tig_to_idx:
                continue
            for dist in dists.split():
                contacts[(atig, btig)].append(int(dist))

        self.contacts = contacts
        print contacts


def main():

    actions = (
        ('agp', 'generate AGP file based on LACHESIS output'),
        ('score', 'score the current LACHESIS CLM'),
        ('optimize', 'optimize the contig order and orientation'),
        # Plotting
        ('heatmap', 'generate heatmap based on LACHESIS output'),
        ('heatmapmovie', 'plot heatmap optimization history'),
        ('syntenymovie', 'plot synteny optimization history'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def optimize(args):
    """
    %prog optimize test.clm test.ids

    Optimize the contig order and orientation, based on CLM file.
    """
    p = OptionParser(optimize.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    clmfile, idsfile = args
    # Load contact map
    clm = CLMFile(clmfile, idsfile)


def syntenymovie(args):
    """
    %prog syntenymovie tour contigs.ref.last

    Plot synteny optimization history.
    """
    from jcvi.compara.synteny import get_bed_filenames
    from jcvi.formats.bed import Bed
    from jcvi.formats.blast import Blast
    from jcvi.utils.cbook import gene_name
    from jcvi.graphics.dotplot import dotplot_main

    p = OptionParser(syntenymovie.__doc__)
    p.add_option("--frames", default=250, type="int",
                 help="Only plot every N frames")
    p.set_beds()
    opts, args, iopts = p.set_image_options(args, figsize="8x8")

    if len(args) != 2:
        sys.exit(not p.print_help())

    tourfile, lastfile = args
    qbedfile, sbedfile = get_bed_filenames(lastfile, p, opts)
    odir = "syntenymovie"
    mkdir(odir)

    tourfile = op.abspath(tourfile)
    lastfile = op.abspath(lastfile)
    qbedfile = op.abspath(qbedfile)
    sbedfile = op.abspath(sbedfile)
    cwd = os.getcwd()
    qbed = Bed(qbedfile, sorted=False)
    contig_to_beds = dict(qbed.sub_beds())
    os.chdir(odir)

    # Make anchorsfile
    anchorsfile = ".".join(op.basename(lastfile).split(".", 2)[:2]) + ".anchors"
    fw = open(anchorsfile, "w")
    for b in Blast(lastfile):
        print >> fw, "\t".join((gene_name(b.query), gene_name(b.subject),
                                str(int(b.score))))
    fw.close()

    # Symlink sbed
    symlink(sbedfile, op.basename(sbedfile))

    for i, label, tour in iter_tours(tourfile, frames=opts.frames):
        # Make BED file with new order
        qb = Bed()
        for contig in tour:
            for x in contig_to_beds[contig]:
                qb.append(x)
        qb.print_to_file(op.basename(qbedfile))
        # Plot dot plot, but do not sort contigs by name (otherwise losing
        # order)
        image_name = "{:06d}".format(i) + "." + iopts.format
        dotplot_main([anchorsfile, "--nosort", "--nosep",
                      "--title", label, "--outfile", image_name])

    os.chdir(cwd)
    make_movie(odir, odir)


def iter_tours(tourfile, frames=250):
    """
    Extract tours from tourfile. Tourfile contains a set of contig
    configurations, generated at each iteration of the genetic algorithm. Each
    configuration has two rows, first row contains iteration id and score,
    second row contains list of contigs, separated by comma.
    """
    fp = open(tourfile)

    for row in fp:
        if row[0] == '>':
            label = row[1:].strip()
            if label.count("-") == 2:
                pf, i, score = label.split("-")
                i = int(i)
            else:
                i = 0
            continue
        else:
            if i % frames != 0:
                continue
            yield i, label, row.split()

    fp.close()


def heatmapmovie(args):
    """
    %prog heatmapmovie tour cached_data/ contigsfasta

    Plot heatmap optimization history.
    """
    p = OptionParser(heatmapmovie.__doc__)
    p.add_option("--frames", default=250, type="int",
                 help="Only plot every N frames")
    opts, args, iopts = p.set_image_options(args, figsize="8x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    tourfile, cdir, contigsfasta = args
    odir = "heatmapmovie"
    mkdir(odir)

    for i, label, tour in iter_tours(tourfile):
        image_name = op.join(odir, "{:06d}".format(i)) + "." + iopts.format
        tour = ",".join(tour)
        heatmap(["main_results", cdir, contigsfasta,
                 "--tour", tour, "--outfile", image_name,
                 "--label", label])
    make_movie(odir, odir)


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
        print_tour(fwtour, tour, label, contig_names, oo)
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
        # Store INIT tour
        print_tour(fwtour, tour, "INIT", contig_names, oo)

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


def print_tour(fwtour, tour, label, contig_names, oo):
    print >> fwtour, ">" + label
    print >> fwtour, " ".join(contig_names[oo[x]] for x in tour)


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
    p.add_option("--label", help="Figure title")
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

    plot_heatmap(M, breaks, opts, iopts)


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


def plot_heatmap(M, breaks, opts, iopts):
    image_name = opts.outfile or ("heatmap." + iopts.format)

    fig = plt.figure(figsize=(iopts.w, iopts.h))
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.imshow(M, cmap=iopts.cmap, origin="lower", interpolation='none')
    xlim = ax.get_xlim()
    for b in breaks[:-1]:
        ax.plot([b, b], xlim, 'w-')
        ax.plot(xlim, [b, b], 'w-')

    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    label = opts.label or op.basename(image_name).rsplit(".", 1)[0]
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
