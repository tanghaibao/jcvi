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

from jcvi.algorithms.formula import outlier_cutoff
from jcvi.algorithms.ec import GA_setup, GA_run
from jcvi.algorithms.matrix import get_signs
from jcvi.apps.base import OptionParser, ActionDispatcher, backup, iglob, mkdir, symlink
from jcvi.apps.grid import Jobs
from jcvi.assembly.allmaps import make_movie
from jcvi.compara.synteny import check_beds, get_bed_filenames
from jcvi.formats.agp import order_to_agp
from jcvi.formats.base import LineFile, must_open
from jcvi.formats.bed import Bed
from jcvi.formats.sizes import Sizes
from jcvi.formats.blast import Blast
from jcvi.graphics.base import normalize_axes, plt, savefig
from jcvi.graphics.dotplot import dotplot
from jcvi.utils.cbook import gene_name
from jcvi.utils.natsort import natsorted


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
    def __init__(self, clmfile, skiprecover=True):
        self.clmfile = clmfile
        self.idsfile = clmfile.rsplit(".", 1)[0] + ".ids"
        self.parse_ids(skiprecover)
        self.parse_clm()

    def parse_ids(self, skiprecover):
        '''IDS file has a list of contigs that need to be ordered. 'recover',
        keyword, if available in the third column, is less confident.

        tig00015093     46912
        tig00035238     46779   recover
        tig00030900     119291
        '''
        idsfile = self.idsfile
        logging.debug("Parse idsfile `{}`".format(idsfile))
        fp = open(idsfile)
        tigs = []
        for row in fp:
            atoms = row.split()
            tig, size = atoms[:2]
            size = int(size)
            if skiprecover and len(atoms) == 3 and atoms[2] == 'recover':
                continue
            tigs.append((tig, size))

        # Arrange contig names and sizes
        _tigs, _sizes = zip(*tigs)
        self.contigs = set(_tigs)
        self.sizes = np.array(_sizes)
        self.tig_to_size = dict(tigs)

        # Initially all contigs are considered active
        self.active = set(_tigs)

    def parse_clm(self):
        clmfile = self.clmfile
        logging.debug("Parse clmfile `{}`".format(clmfile))
        fp = open(clmfile)
        contacts = {}
        orientations = defaultdict(list)
        for row in fp:
            atoms = row.strip().split('\t')
            assert len(atoms) == 3, "Malformed line `{}`".format(atoms)
            abtig, links, dists = atoms
            atig, btig = abtig.split()
            at, ao = atig[:-1], atig[-1]
            bt, bo = btig[:-1], btig[-1]
            if at not in self.tig_to_size:
                continue
            if bt not in self.tig_to_size:
                continue
            dists = [int(x) for x in dists.split()]
            contacts[(at, bt)] = dists
            orientations[(at, bt)].append((ao + bo, sum(dists)))

        self.contacts = contacts
        self.orientations = orientations

    def calculate_densities(self):
        """
        Calculate the density of inter-contig links per base. Strong contigs are
        considered to have high level of inter-contig links in the current
        partition.
        """
        active = self.active
        densities = defaultdict(int)
        for (at, bt), dists in self.contacts.items():
            if not (at in active and bt in active):
                continue
            densities[at] += len(dists)
            densities[bt] += len(dists)

        logdensities = {}
        for x, d in densities.items():
            s = self.tig_to_size[x]
            logd = np.log10(d * 1. / min(s, 500000))
            logdensities[x] = logd

        return logdensities

    def report_active(self):
        logging.debug("Active contigs: {} (length={})"\
                    .format(len(self.active), self.active_length))

    def activate(self, minsize=10000):
        """
        Select strong contigs in the current partition.
        """
        done = False
        self.report_active()
        while not done:
            logdensities = self.calculate_densities()
            lb, ub = outlier_cutoff(logdensities.values(), threshold=3)
            logging.debug("Lower bound: {}, Upper bound: {}".format(lb, ub))
            remove = set(x for x, d in logdensities.items() if (d < lb or d > ub))
            if remove:
                self.active -= remove
                self.report_active()
            else:
                done = True

        logging.debug("Remove contigs with size < {}".format(minsize))
        self.active = set(x for x in self.active if self.tig_to_size[x] >= minsize)
        self.report_active()

    @property
    def active_contigs(self):
        return list(self.active)

    @property
    def active_sizes(self):
        return np.array([self.tig_to_size[x] for x in self.active])

    @property
    def active_length(self):
        return self.active_sizes.sum()

    @property
    def N(self):
        return len(self.active)

    @property
    def tig_to_idx(self):
        return dict((x, i) for (i, x) in enumerate(self.active))

    @property
    def M(self):
        """
        Contact frequency matrix. Each cell contains how many inter-contig links
        between i-th and j-th contigs.
        """
        N = self.N
        tig_to_idx = self.tig_to_idx
        M = np.zeros((N, N), dtype=int)
        for (at, bt), dists in self.contacts.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            M[ai, bi] = M[bi, ai] = len(dists)
        return M

    @property
    def O(self):
        """
        Pairwise strandedness matrix. Each cell contains whether i-th and j-th
        contig are the same orientation +1, or opposite orientation -1.
        """
        N = self.N
        tig_to_idx = self.tig_to_idx
        O = np.zeros((N, N), dtype=int)
        for (at, bt), dists in self.orientations.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            mo, md = min(dists, key=lambda x: x[1])
            strandedness = 1 if mo[0] == mo[1] else -1
            O[ai, bi] = O[bi, ai] = strandedness
        return O


def main():

    actions = (
        # LACHESIS output processing
        ('agp', 'generate AGP file based on LACHESIS output'),
        ('score', 'score the current LACHESIS CLM'),
        # Scaffolding
        ('optimize', 'optimize the contig order and orientation'),
        ('density', 'estimate link density of contigs'),
        # Plotting
        ('movieframe', 'plot heatmap and synteny for a particular tour'),
        ('movie', 'plot heatmap optimization history in a tourfile'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def density(args):
    """
    %prog density test.clm

    Estimate link density of contigs.
    """
    p = OptionParser(density.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clmfile, = args
    clm = CLMFile(clmfile, skiprecover=False)
    pf = clmfile.rsplit(".", 1)[0]

    logdensities = clm.calculate_densities()
    densityfile = pf + ".density"
    fw = open(densityfile, "w")
    for name, logd in logdensities.items():
        s = clm.tig_to_size[name]
        print >> fw, "\t".join(str(x) for x in (name, s, logd))
    fw.close()
    logging.debug("Density written to `{}`".format(densityfile))

    clm.activate()


def optimize(args):
    """
    %prog optimize test.clm

    Optimize the contig order and orientation, based on CLM file.
    """
    p = OptionParser(optimize.__doc__)
    p.add_option("--startover", default=False, action="store_true",
                 help="Do not resume from existing tour file")
    p.set_outfile(outfile=None)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clmfile, = args
    startover = opts.startover
    # Load contact map
    clm = CLMFile(clmfile, skiprecover=False)
    clm.activate()
    N = clm.N

    tourfile = opts.outfile or clmfile.rsplit(".", 1)[0] + ".tour"
    if startover or (not op.exists(tourfile)):
        tour = range(N)  # Use starting (random) order otherwise
    else:
        logging.debug("File `{}` found".format(tourfile))
        tour = iter_last_tour(tourfile, clm)
        clm.active = set(tour)
        tig_to_idx = clm.tig_to_idx
        tour = [tig_to_idx[x] for x in tour]
        backup(tourfile)

    # Prepare input files
    tour_contigs = clm.active_contigs
    tour_sizes = clm.active_sizes
    tour_M = clm.M
    tour_O = clm.O
    oo = range(N)

    # Determine orientations
    signs = get_signs(tour_O, validate=False)

    fwtour = open(tourfile, "w")
    def callback(tour, gen, oo):
        fitness = tour.fitness if hasattr(tour, "fitness") else None
        label = "GA-{0}".format(gen)
        if fitness:
            fitness = "{0}".format(fitness).split(",")[0].replace("(", "")
            label += "-" + fitness
        print_tour(fwtour, tour, label, tour_contigs, oo, signs=signs)
        return tour

    # Store INIT tour
    print_tour(fwtour, tour, "INIT", tour_contigs, oo, signs=signs)

    # Faster Cython version for evaluation
    from .chic import score_evaluate
    callbacki = partial(callback, oo=oo)
    toolbox = GA_setup(tour)
    toolbox.register("evaluate", score_evaluate,
                     tour_sizes=tour_sizes, tour_M=tour_M)
    tour, tour.fitness = GA_run(toolbox, ngen=1000, npop=100, cpus=opts.cpus,
                                callback=callbacki)
    print tour, tour.fitness
    fwtour.close()


def prepare_synteny(tourfile, lastfile, odir, p, opts):
    """
    Prepare synteny plots for movie().
    """
    qbedfile, sbedfile = get_bed_filenames(lastfile, p, opts)
    qbedfile = op.abspath(qbedfile)
    sbedfile = op.abspath(sbedfile)

    qbed = Bed(qbedfile, sorted=False)
    contig_to_beds = dict(qbed.sub_beds())

    # Create a separate directory for the subplots and movie
    mkdir(odir)
    os.chdir(odir)
    logging.debug("Change into subdir `{}`".format(odir))

    # Make anchorsfile
    anchorsfile = ".".join(op.basename(lastfile).split(".", 2)[:2]) + ".anchors"
    fw = open(anchorsfile, "w")
    for b in Blast(lastfile):
        print >> fw, "\t".join((gene_name(b.query), gene_name(b.subject),
                                str(int(b.score))))
    fw.close()

    # Symlink sbed
    symlink(sbedfile, op.basename(sbedfile))

    return anchorsfile, qbedfile, contig_to_beds


def separate_tour_and_o(row):
    """
    The tour line typically contains contig list like:
    tig00044568+ tig00045748- tig00071055- tig00015093- tig00030900-

    This function separates the names from the orientations.
    """
    tour = []
    tour_o = []
    for contig in row.split():
        if contig[-1] in ('+', '-', '?'):
            tour.append(contig[:-1])
            tour_o.append(contig[-1])
        else:  # Unoriented
            tour.append(contig)
            tour_o.append('?')
    return tour, tour_o


def iter_last_tour(tourfile, clm):
    """
    Extract last tour from tourfile. The clm instance is also passed in to check
    if any contig is covered in the clm.
    """
    row = open(tourfile).readlines()[-1]
    _tour, tour_o = separate_tour_and_o(row)
    tour = []
    for tc in _tour:
        if tc not in clm.contigs:
            logging.debug("Contig `{}` in file `{}` not found in `{}`"\
                            .format(tc, tourfile, clm.idsfile))
            continue
        tour.append(tc)
    return tour


def iter_tours(tourfile, frames=1):
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
            tour, tour_o = separate_tour_and_o(row)
            yield i, label, tour, tour_o

    fp.close()


def movie(args):
    """
    %prog movie test.tour test.clm ref.contigs.last

    Plot optimization history.
    """
    p = OptionParser(movie.__doc__)
    p.add_option("--frames", default=500, type="int",
                 help="Only plot every N frames")
    p.set_beds()
    opts, args, iopts = p.set_image_options(args, figsize="16x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    tourfile, clmfile, lastfile = args
    tourfile = op.abspath(tourfile)
    clmfile = op.abspath(clmfile)
    lastfile = op.abspath(lastfile)
    cwd = os.getcwd()
    odir = op.basename(tourfile).rsplit(".", 1)[0] + "-movie"
    anchorsfile, qbedfile, contig_to_beds = \
                prepare_synteny(tourfile, lastfile, odir, p, opts)

    args = []
    for i, label, tour, tour_o in iter_tours(tourfile, frames=opts.frames):
        padi = "{:06d}".format(i)
        # Make sure the anchorsfile and bedfile has the serial number in,
        # otherwise parallelization may fail
        a, b = op.basename(anchorsfile).split(".", 1)
        ianchorsfile = a + "_" + padi + "." + b
        symlink(anchorsfile, ianchorsfile)

        # Make BED file with new order
        qb = Bed()
        for contig, o in zip(tour, tour_o):
            bedlines = contig_to_beds[contig]
            if o == '-':
                bedlines.reverse()
            for x in bedlines:
                qb.append(x)

        a, b = op.basename(qbedfile).split(".", 1)
        ibedfile = a + "_" + padi + "." + b
        qb.print_to_file(ibedfile)
        # Plot dot plot, but do not sort contigs by name (otherwise losing
        # order)
        image_name = padi + "." + iopts.format

        tour = ",".join(tour)
        args.append([[tour, clmfile, ianchorsfile,
                    "--outfile", image_name, "--label", label]])

    Jobs(movieframe, args).run()

    os.chdir(cwd)
    make_movie(odir, odir)


def score(args):
    """
    %prog score main_results/ cached_data/ contigsfasta

    Score the current LACHESIS CLM.
    """
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


def print_tour(fwtour, tour, label, contig_names, oo, signs=None):
    print >> fwtour, ">" + label
    if signs is not None:
        contig_o = []
        for x in tour:
            idx = oo[x]
            sign = {1: '+', 0: '?', -1: '-'}[signs[idx]]
            contig_o.append(contig_names[idx] + sign)
        print >> fwtour, " ".join(contig_o)
    else:
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


def movieframe(args):
    """
    %prog movieframe tour test.clm contigs.ref.anchors

    Draw heatmap and synteny in the same plot.
    """
    p = OptionParser(movieframe.__doc__)
    p.add_option("--label", help="Figure title")
    p.set_beds()
    p.set_outfile(outfile=None)
    opts, args, iopts = p.set_image_options(args, figsize="16x8",
                                            style="white", cmap="coolwarm")

    if len(args) != 3:
        sys.exit(not p.print_help())

    tour, clmfile, anchorsfile = args
    tour = tour.split(",")
    image_name = opts.outfile or ("movieframe." + iopts.format)
    label = opts.label or op.basename(image_name).rsplit(".", 1)[0]

    clm = CLMFile(clmfile, skiprecover=False)
    totalbins, bins, breaks = make_bins(tour, clm.tig_to_size)
    M = read_clm(clm, totalbins, bins)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])        # whole canvas
    ax1 = fig.add_axes([.05, .1, .4, .8])    # heatmap
    ax2 = fig.add_axes([.55, .1, .4, .8])    # dot plot
    ax2_root = fig.add_axes([.5, 0, .5, 1])  # dot plot canvas

    # Left axis: heatmap
    plot_heatmap(ax1, M, breaks, label, iopts)

    # Right axis: synteny
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts,
                sorted=False)
    dotplot(anchorsfile, qbed, sbed, fig, ax2_root, ax2, sep=False, title="")

    root.text(.5, .975, label, color="darkslategray", ha="center", va="center")
    normalize_axes(root)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def make_bins(tour, sizes):
    breaks = []
    start = 0
    bins = {}
    for x in tour:
        size = sizes[x]
        end = start + int(math.ceil(size / 100000.))
        bins[x] = (start, end)
        start = end
    breaks.append(start)

    totalbins = start
    return totalbins, bins, breaks


def read_clm(clm, totalbins, bins):
    M = np.zeros((totalbins, totalbins))
    for (x, y), dists in clm.contacts.items():
        if x not in bins or y not in bins:
            continue
        xstart, xend = bins[x]
        ystart, yend = bins[y]
        #z = float(z) / ((xend - xstart) * (yend - ystart))
        z = len(dists)
        M[xstart:xend, ystart:yend] = z
        M[ystart:yend, xstart:xend] = z

    M = np.log10(M + 1)
    return M


def plot_heatmap(ax, M, breaks, label, iopts):
    ax.imshow(M, cmap=iopts.cmap, origin="lower", interpolation='none')
    xlim = ax.get_xlim()
    for b in breaks[:-1]:
        ax.plot([b, b], xlim, 'w-')
        ax.plot(xlim, [b, b], 'w-')
    ax.set_xlim(xlim)
    ax.set_ylim(xlim)
    ax.set_xticklabels([int(x) for x in ax.get_xticks()],
                        family='Helvetica', color="gray")
    ax.set_yticklabels([int(x) for x in ax.get_yticks()],
                        family='Helvetica', color="gray")


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
