#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scaffold Ordering with Weighted Maps.
"""

import os.path as op
import sys
import logging

import numpy as np

from itertools import combinations
from collections import defaultdict
from scipy.stats import spearmanr

from jcvi.algorithms.formula import reject_outliers
from jcvi.algorithms.graph import merge_paths, longest_path_weighted_nodes
from jcvi.algorithms.lis import longest_monotonous_subseq_length
from jcvi.algorithms.matrix import determine_signs
from jcvi.formats.agp import AGP, order_to_agp, build as agp_build
from jcvi.formats.base import DictFile, FileMerger, must_open
from jcvi.formats.bed import Bed, BedLine, sort
from jcvi.formats.chain import fromagp
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import human_size
from jcvi.utils.counter import Counter
from jcvi.utils.grouper import Grouper
from jcvi.apps.base import OptionParser, ActionDispatcher, debug, sh, \
            need_update
debug()


START, END = "START", "END"
np.seterr(invalid="ignore")


class Scaffold (object):

    def __init__(self, seqid, mapc):
        self.markers = mapc.extract(seqid)
        self.seqid = seqid
        self.mapc = mapc

    @property
    def mlg_counts(self):
        return Counter([x.mlg for x in self.markers])

    def add_LG_pairs(self, G, mappair):
        # Computes co-occurrences of LG pairs
        cc = self.mlg_counts.items()
        mappair = sorted(mappair)
        for (ak, av), (bk, bv) in combinations(cc, 2):
            aks, bks = ak.split("-")[0], bk.split("-")[0]
            if sorted((aks, bks)) != mappair:
                continue
            weight = min(av, bv)
            G[ak, bk] += weight
            G[bk, ak] += weight


class ScaffoldOO (object):
    """
    This contains the routine to construct order and orientation for the
    scaffolds per partition.
    """
    def __init__(self, lgs, scaffolds, mapc, pivot, weights, sizes,
                 function=(lambda x: x.rank)):

        self.lgs = lgs
        self.lengths = mapc.compute_lengths(function)
        self.bins = mapc.get_bins(function)
        self.function = function
        self.sizes = sizes

        signs = self.assign_orientation(scaffolds, pivot, weights)
        scaffolds = zip(scaffolds, signs)

        tour = self.assign_order(scaffolds, pivot, weights)
        self.tour = tour

        for mlg in self.lgs:
            mapname, lg = mlg.rsplit("-", 1)
            if mapname == pivot:
                self.object = "chr{0}".format(lg)
                break

    def weighted_mean(self, a, weights):
        a, w = zip(*a)
        w = [weights[x] for x in w]
        return np.average(a, weights=w)

    def get_markers(self, lg, scaffold, orientation):
        xs = self.bins.get((lg, scaffold), [])
        if orientation < 0:
            xs = xs[::-1]
        return xs

    def get_rho(self, xy):
        if not xy:
            return 0
        x, y = zip(*xy)
        rho, p_value = spearmanr(x, y)
        if np.isnan(rho):
            rho = 0
        return rho

    def assign_order(self, scaffolds, pivot, weights):
        """
        The goal is to assign scaffold orders. To help order the scaffolds, two
        dummy node, START and END, mark the ends of the chromosome. We connect
        START to each scaffold (directed), and each scaffold to END.
        """
        f = self.function
        positions, paths, w = [], [], []
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            length = self.lengths[lg]
            position = {}
            for a, ao in scaffolds:
                xa = self.get_markers(lg, a, ao)
                if not xa:
                    continue
                d = np.median([f(x) for x in xa])
                assert 0 <= d <= length
                position[a] = d

            if mapname == pivot:
                pivot_position = position

            positions.append(position)
            w.append(weights[mapname])

        for position in positions:
            # Sort the order based on median distance
            path = sorted((v, k) for k, v in position.items())
            vv, path = zip(*path)

            # Flip order if path goes in the opposite direction to the pivot
            common = []
            for a, ap in position.items():
                if a not in pivot_position:
                    continue
                pp = pivot_position[a]
                common.append((ap, pp))

            rho = self.get_rho(common)
            if rho < 0:
                path = path[::-1]
            paths.append(path)

        G = merge_paths(paths, weights=w)
        for s, so in scaffolds:  # Connect everything to DUMMY ends
            G.add_edge(START, s)
            G.add_edge(s, END)
        tour, total_size = longest_path_weighted_nodes(G, START, END, self.sizes)

        # Remove dummy nodes
        assert tour[0] == START and tour[-1] == END
        tour = tour[1:-1]
        logging.debug("Best order contains {0} scaffolds (Score={1})".\
                        format(len(tour), total_size))

        scaffolds_oo = dict(scaffolds)
        recode = {0: '?', 1: '+', -1: '-'}
        tour = [(x, recode[scaffolds_oo[x]]) for x in tour]
        return tour

    def get_orientation(self, si, sj):
        '''
        si, sj are two number series. To compute whether these two series have
        same orientation or not. We combine them in the two orientation
        configurations and compute length of the longest monotonous series.
        '''
        if not si or not sj:
            return 0
        # Same orientation configuration
        a = longest_monotonous_subseq_length(si + sj)
        b = longest_monotonous_subseq_length(sj + si)
        # Opposite orientation configuration
        c = longest_monotonous_subseq_length(si + sj[::-1])
        d = longest_monotonous_subseq_length(sj[::-1] + si)
        return max(a, b)[0] - max(c, d)[0]

    def assign_orientation(self, scaffolds, pivot, weights):
        bins = self.bins
        f = self.function
        signs = defaultdict(list)
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            oo = []
            nmarkers = {}
            series = []
            for s in scaffolds:
                xs = bins.get((lg, s), [])
                nmarkers[s] = len(xs)
                if not xs:
                    series.append([])
                    oo.append(0)
                    continue
                physical_to_cm = [(x.pos, f(x)) for x in xs]
                rho = self.get_rho(physical_to_cm)
                oo.append(rho)
                series.append([f(x) for x in xs])

            if mapname == pivot:
                pivot_oo = oo

            for i, j in combinations(range(len(scaffolds)), 2):
                si, sj = series[i], series[j]
                d = self.get_orientation(si, sj)

                if not d:
                    continue
                signs[(i, j)].append((d, mapname))

        for e, v in signs.items():
            signs[e] = self.weighted_mean(v, weights)

        signs_edges = sorted((a, b, w) for (a, b), w in signs.items())
        signs = determine_signs(scaffolds, signs_edges)

        # Finally flip this according to pivot map, then weight by #_markers
        nmarkers = [nmarkers[x] for x in scaffolds]
        flipr = signs * np.sign(np.array(pivot_oo)) * nmarkers
        if sum(flipr) < 0:
            signs = - signs
        return signs


class CSVMapLine (object):

    def __init__(self, row, sep=",", mapname=None):
        # ScaffoldID,ScaffoldPosition,LinkageGroup,GeneticPosition
        args = row.strip().split(sep)
        self.seqid = args[0]
        self.pos = int(args[1])
        self.lg = args[2]
        self.cm = float(args[3])
        self.mapname = mapname

    @property
    def bedline(self):
        marker = "{0}-{1}:{2:.6f}".format(self.mapname, self.lg, self.cm)
        track = "{0}:{1}".format(self.seqid, self.pos)
        return "\t".join(str(x) for x in \
                (self.seqid, self.pos - 1, self.pos, marker, track))


class Marker (object):

    def __init__(self, b):
        self.seqid = b.seqid
        self.pos = b.start
        self.mlg, cm = b.accn.split(":")
        self.mapname, self.lg = b.accn.split("-")
        self.cm = float(cm)
        self.accn = b.accn
        self.rank = -1

    def __str__(self):
        return "\t".join(str(x) for x in
                    (self.seqid, self.pos - 1, self.pos,
                     self.accn, self.rank))

    __repr__ = __str__


class Map (list):

    def __init__(self, filename):
        bed = Bed(filename)
        for b in bed:
            self.append(Marker(b))

        self.nmarkers = len(self)
        self.mlg = set(x.mlg for x in self)
        logging.debug("Map contains {0} markers in {1} linkage groups.".\
                      format(self.nmarkers, len(self.mlg)))
        self.ranks = self.compute_ranks()

    def extract(self, seqid):
        r = [x for x in self if x.seqid == seqid]
        r.sort(key=lambda x: x.pos)
        return r

    def extract_mlg(self, mlg):
        r = [x for x in self if x.mlg == mlg]
        return r

    def compute_ranks(self):
        ranks = {}  # Store the length for each linkage group
        for mlg in self.mlg:
            mlg_set = [x for x in self if x.mlg == mlg]
            mlg_set.sort(key=lambda x: x.cm)
            for rank, marker in enumerate(mlg_set):
                marker.rank = rank
            ranks[mlg] = mlg_set
        return ranks

    def compute_lengths(self, function):
        lengths = {}
        for mlg, v in self.ranks.items():
            lengths[mlg] = max(function(x) for x in v)
        return lengths

    def get_bins(self, function, remove_outliers=True):
        s = defaultdict(list)
        for m in self:
            s[(m.mlg, m.seqid)].append(m)

        if remove_outliers:
            for pair, markers in s.items():
                s[pair] = self.remove_outliers(markers, function)
        return s

    def remove_outliers(self, markers, function):
        data = [function(x) for x in markers]
        reject = reject_outliers(data)
        clean_markers = [m for m, r in zip(markers, reject) if not r]
        return clean_markers

    @property
    def seqids(self):
        return sorted(set(x.seqid for x in self))

    @property
    def mapnames(self):
        return sorted(set(x.mapname for x in self))

    @property
    def mlgs(self):
        return sorted(set(x.mlg for x in self))


class Weights (DictFile):

    def __init__(self, filename, mapnames, cast=int):
        super(Weights, self).__init__(filename, cast=cast)
        self.update_maps(mapnames)

    def update_maps(self, mapnames, default=1):
        for m in mapnames:
            if m in self:
                continue
            self[m] = default
            logging.debug("Weight for `{0}` set to {1}.".format(m, default))

    def get_pivot(self, mapnames):
        return max((w, m) for m, w in self.items() if m in mapnames)


class Layout (object):

    def __init__(self, mlgsizes):

        self.mlgsizes = mlgsizes
        self.partition()
        self.calculate_coords()

    def partition(self, N=2):
        # Partition LGs into two sides with approximately similar sum of sizes
        endtime = [0] * N
        parts = []
        for i in xrange(N):
            parts.append([])
        # LPT greedy algorithm, sort by LG size decreasing
        for mlg, mlgsize in sorted(self.mlgsizes.items(), key=lambda x: - x[-1]):
            mt, mi = min((x, i) for (i, x) in enumerate(endtime))
            endtime[mi] += mlgsize
            parts[mi].append((mlg, mlgsize))
        self.parts = parts

    def calculate_coords(self, r=.8, gapsize=.1):
        # Find the larger partition
        part_sizes = []
        for p in self.parts:
            ps = sum(ms for m, ms in p)
            part_sizes.append((ps, len(p) - 1))
        max_part_size, ngaps = max(part_sizes)
        gaps = gapsize * ngaps
        ratio = (r - gaps) / max_part_size
        self.ratio = ratio

        coords = {}
        for x, p, (ps, ngaps) in zip((.25, .75), self.parts, part_sizes):
            gaps = gapsize * ngaps
            ystart = (1 + ratio * ps + gaps) / 2
            for m, ms in p:
                mlen = ratio * ms
                coords[m] = (x, ystart - mlen, ystart)
                ystart -= mlen + gapsize
        self.coords = coords


def main():

    actions = (
        ('merge', 'merge csv maps and convert to bed format'),
        ('path', 'construct golden path given a set of genetic maps'),
        ('build', 'build associated FASTA and CHAIN file'),
        ('plot', 'plot matches between goldenpath and maps for single object'),
        ('plotall', 'plot matches between goldenpath and maps for all objects'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def merge(args):
    """
    %prog merge map1 map2 map3 ...

    Convert csv maps to bed format.

    Each input map is csv formatted, for example:

    ScaffoldID,ScaffoldPosition,LinkageGroup,GeneticPosition
    scaffold_2707,11508,1,0
    scaffold_2707,11525,1,1.00000000000001e-05
    scaffold_759,81336,1,49.7317510625759
    """
    p = OptionParser(merge.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    maps = args
    fp = must_open(maps)
    b = Bed()
    mapnames = set()
    for row in fp:
        mapname = fp.filename().split(".")[0]
        mapnames.add(mapname)
        try:
            m = CSVMapLine(row, mapname=mapname)
            b.append(BedLine(m.bedline))
        except ValueError:  # header
            continue

    b.print_to_file(filename=opts.outfile, sorted=True)

    assert len(maps) == len(mapnames), "You have a collision in map names"
    weightsfile = "weights.txt"
    if op.exists(weightsfile):
        logging.debug("Weights file `{0}` found. Will not overwrite.".\
                        format(weightsfile))
        return

    fw = open(weightsfile, "w")
    for mapname in sorted(mapnames):
        weight = 1
        print >> fw, mapname, weight
    logging.debug("Weights file written to `{0}`.".format(weightsfile))


def best_no_ambiguous(d, label):
    best, best_value = max(d.items(), key=lambda x: x[1])
    if d.values().count(best_value) > 1:  # tie
        print >> sys.stderr, "AMBIGUOUS", label, d
        return None, None
    return best, best_value


def path(args):
    """
    %prog path map.bed weights.txt scaffolds.fasta

    Construct golden path given a set of genetic maps. The respective weight for
    each map is given in file `weights.txt`. The map with the highest weight is
    considered the pivot map. The final output is an AGP file that contain the
    ordered scaffolds.
    """
    distance_choices = ("cM", "rank")

    p = OptionParser(path.__doc__)
    p.add_option("--distance", default="rank", choices=distance_choices,
                 help="Distance function when building consensus")
    p.add_option("--gapsize", default=100, type="int",
                 help="Insert gaps of size")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, weightsfile, fastafile = args

    cc = Map(bedfile)
    mapnames = cc.mapnames
    allseqids = cc.seqids
    weights = Weights(weightsfile, mapnames)
    pivot_weight, pivot = weights.get_pivot(mapnames)

    logging.debug("Pivot map: `{0}` (weight={1}).".format(pivot, pivot_weight))
    gapsize = opts.gapsize
    function = opts.distance
    function = (lambda x: x.cm) if function == "cM" else \
               (lambda x: x.rank)

    # Partition the linkage groups into consensus clusters
    C = Grouper()
    # Initialize the partitions
    for mlg in cc.mlgs:
        C.join(mlg)
    for mapname in mapnames:
        if mapname == pivot:
            continue
        # Compute co-occurrence between LG pairs
        G = defaultdict(int)
        for s in allseqids:
            s = Scaffold(s, cc)
            s.add_LG_pairs(G, (pivot, mapname))
        # Convert edge list to adj list
        nodes = defaultdict(list)
        for (a, b), w in G.items():
            nodes[a].append((b, w))
        # Find the best pivot LG every non-pivot LG matches to
        for n, neighbors in nodes.items():
            if n.split("-")[0] == pivot:
                continue
            neighbors = dict(neighbors)
            best_neighbor, best_value = best_no_ambiguous(neighbors, n)
            if best_neighbor is None:
                continue
            C.join(n, best_neighbor)

    partitions = defaultdict(list)
    # Partition the scaffolds and assign them to one consensus
    for s in allseqids:
        s = Scaffold(s, cc)
        seqid = s.seqid
        counts = {}
        for mlg, count in s.mlg_counts.items():
            consensus = C[mlg]
            mapname = mlg.split("-")[0]
            mw = weights[mapname]
            if consensus not in counts:
                counts[consensus] = 0
            counts[consensus] += count * mw
        best_consensus, best_value = best_no_ambiguous(counts, seqid)
        if best_consensus is None:
            continue
        partitions[best_consensus].append(seqid)

    # Perform OO within each partition
    pf = bedfile.rsplit(".", 1)[0]
    agpfile = pf + ".chr.agp"
    tourfile = pf + ".tour"
    sizes = Sizes(fastafile).mapping
    fwagp = must_open(agpfile, "w")
    fwtour = must_open(tourfile, "w")
    solutions = []
    for lgs, scaffolds in sorted(partitions.items()):
        tag = "|".join(lgs)
        s = ScaffoldOO(lgs, scaffolds, cc, pivot, weights, sizes,
                       function=function)

        for fw in (sys.stderr, fwtour):
            print >> fw, ">{0} ({1})".format(s.object, tag)
            print >> fw, " ".join("".join(x) for x in s.tour)
        solutions.append(s)
    fwtour.close()

    for s in sorted(solutions, key=lambda x: x.object):
        order_to_agp(s.object, s.tour, sizes, fwagp, gapsize=gapsize,
                     gaptype="map")
    fwagp.close()

    logging.debug("AGP file written to `{0}`.".format(agpfile))
    logging.debug("Tour file written to `{0}`.".format(tourfile))


def write_unplaced_agp(agpfile, scaffolds, unplaced_agp):
    agp = AGP(agpfile)
    scaffolds_seen = set(x.component_id for x in agp)
    sizes = Sizes(scaffolds).mapping
    fwagp = must_open(unplaced_agp, "w")
    for s in sorted(sizes.keys()):
        if s in scaffolds_seen:
            continue
        order_to_agp(s, [(s, "?")], sizes, fwagp)
    logging.debug("Write unplaced AGP to `{0}`.".format(unplaced_agp))


def build(args):
    """
    %prog build agpfile scaffolds.fasta map.bed

    Build associated genome FASTA file and CHAIN file that can be used to lift
    old coordinates to new coordinates. The CHAIN file will be used to lift the
    original marker positions to new positions in the reconstructed genome. The
    new positions of the markers will be reported in *.lifted.bed.
    """
    p = OptionParser(build.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    chr_agp, scaffolds, mapbed = args
    pf = chr_agp.split(".")[0]
    chr_fasta = pf + ".chr.fasta"
    if need_update((chr_agp, scaffolds), chr_fasta):
        agp_build([chr_agp, scaffolds, chr_fasta])

    unplaced_agp = pf + ".unplaced.agp"
    if need_update((chr_agp, scaffolds), unplaced_agp):
        write_unplaced_agp(chr_agp, scaffolds, unplaced_agp)

    unplaced_fasta = pf + ".unplaced.fasta"
    if need_update((unplaced_agp, scaffolds), unplaced_fasta):
        agp_build([unplaced_agp, scaffolds, unplaced_fasta])

    combined_agp = pf + ".agp"
    if need_update((chr_agp, unplaced_agp), combined_agp):
        FileMerger((chr_agp, unplaced_agp), combined_agp).merge()

    combined_fasta = pf + ".fasta"
    if need_update((chr_fasta, unplaced_fasta), combined_fasta):
        FileMerger((chr_fasta, unplaced_fasta), combined_fasta).merge()

    chainfile = pf + ".chain"
    if need_update((combined_agp, scaffolds, combined_fasta), chainfile):
        fromagp([combined_agp, scaffolds, combined_fasta])

    liftedbed = mapbed.rsplit(".", 1)[0] + ".lifted.bed"
    if need_update((mapbed, chainfile), liftedbed):
        cmd = "liftOver -minMatch=1 {0} {1} {2} unmapped".\
                format(mapbed, chainfile, liftedbed)
        sh(cmd)

    sort([liftedbed, "-i"])  # Sort bed in place


def plot(args):
    """
    %prog plot seqid map.lifted.bed agpfile weightsfile

    Plot the matchings between the reconstructed pseudomolecules and the maps.
    Two types of visualizations are available in one canvas:

    1. Parallel axes, and matching markers are shown in connecting lines;
    2. Scatter plot.
    """
    from jcvi.graphics.base import plt, savefig, normalize_axes, set2
    from jcvi.graphics.chromosome import Chromosome, GeneticMap, \
                HorizontalChromosome

    p = OptionParser(plot.__doc__)
    p.add_option("--links", default=10, type="int",
                 help="Only plot matchings more than")
    opts, args, iopts = p.set_image_options(args, figsize="10x6")

    if len(args) != 4:
        sys.exit(not p.print_help())

    seqid, bedfile, agpfile, weightsfile = args
    links = opts.links

    cc = Map(bedfile)
    allseqids = cc.seqids
    mapnames = cc.mapnames
    weights = Weights(weightsfile, mapnames)
    assert seqid in allseqids, "{0} not in {1}".format(seqid, allseqids)

    s = Scaffold(seqid, cc)
    mlgs = [k for k, v in s.mlg_counts.items() if v >= links]
    mlgsizes = {}
    for mlg in mlgs:
        mm = cc.extract_mlg(mlg)
        mlgsize = max(x.cm for x in mm)
        mlgsizes[mlg] = mlgsize

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax1 = fig.add_axes([0, 0, .5, 1])
    ax2 = fig.add_axes([.5, 0, .5, 1])

    # Find the layout first
    ystart, ystop = .9, .1
    L = Layout(mlgsizes)
    coords = L.coords

    tip = .02
    marker_pos = {}
    # Palette
    colors = dict((mapname, set2[i]) for i, mapname in enumerate(mapnames))
    colors = dict((mlg, colors[mlg.split("-")[0]]) for mlg in mlgs)
    # Parallel coordinates
    for mlg, (x, y1, y2) in coords.items():
        mm = cc.extract_mlg(mlg)
        markers = [(m.accn, m.cm) for m in mm]  # exhaustive marker list
        xy = [(m.pos, m.cm) for m in mm if m.seqid == seqid]
        mx, my = zip(*xy)
        rho, p_value = spearmanr(mx, my)
        flip = rho < 0

        g = GeneticMap(ax1, x, y1, y2, markers, tip=tip, flip=flip)
        extra = -2 * tip if x < .5 else 2 * tip
        ha = "right" if x < .5 else "left"
        mapname = mlg.split("-")[0]
        label = "{0} (weight={1})".format(mlg, weights[mapname])
        ax1.text(x + extra, (y1 + y2) / 2, label, color=colors[mlg],
                 ha=ha, va="center", rotation=90)
        marker_pos.update(g.marker_pos)

    agp = AGP(agpfile)
    agp = [x for x in agp if x.object == seqid]
    chrsize = max(x.object_end for x in agp)

    # Pseudomolecules in the center
    r = ystart - ystop
    ratio = r / chrsize
    f = lambda x: (ystart - ratio * x)
    patchstart = [f(x.object_beg) for x in agp if not x.is_gap]
    Chromosome(ax1, .5, ystart, ystop, width=2 * tip, patch=patchstart, lw=2)

    label = "{0} ({1})".format(seqid, human_size(chrsize, precision=0))
    ax1.text(.5, ystart + tip, label, ha="center")

    scatter_data = defaultdict(list)
    # Connecting lines
    for b in s.markers:
        marker_name = b.accn
        if marker_name not in marker_pos:
            continue

        cx = .5
        cy = f(b.pos)
        mx = coords[b.mlg][0]
        my = marker_pos[marker_name]

        extra = -tip if mx < cx else tip
        extra *= 1.25  # leave boundaries for aesthetic reasons
        cx += extra
        mx -= extra
        ax1.plot((cx, mx), (cy, my), "-", color=colors[b.mlg])
        scatter_data[b.mlg].append((b.pos, b.cm))

    # Scatter plot, same data as parallel coordinates
    xstart, xstop = sorted((ystart, ystop))
    f = lambda x: (xstart + ratio * x)
    patchstart = [f(x.object_beg) for x in agp if not x.is_gap]
    HorizontalChromosome(ax2, xstart, xstop, ystop,
                         height=2 * tip, patch=patchstart, lw=2)

    gap = .03
    ratio = (r - gap * len(mlgs) - tip) / sum(mlgsizes.values())

    for mlg, mlgsize in sorted(mlgsizes.items(), key=lambda x: - x[-1]):
        height = ratio * mlgsize
        ystart -= height
        xx = .5 + xstart / 2
        width = r / 2
        color = colors[mlg]
        ax = fig.add_axes([xx, ystart, width, height])
        ystart -= gap
        sd = scatter_data[mlg]
        xx, yy = zip(*sd)
        ax.plot(xx, yy, ".", color=color)
        ax.set_xlim(0, chrsize)
        ax.set_ylim(0, mlgsize)
        ax.set_xticks([])
        while height / len(ax.get_yticks()) < .03:
            ax.set_yticks(ax.get_yticks()[::2])  # Sparsify the ticks
        yticklabels = [int(x) for x in ax.get_yticks()]
        ax.set_yticklabels(yticklabels, family='Helvetica')
        ax.set_ylabel(mlg, color=color)

    normalize_axes((ax1, ax2, root))
    image_name = seqid + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    plt.close(fig)


def plotall(args):
    """
    %prog plot map.lifted.bed agpfile weightsfile

    Plot the matchings between the reconstructed pseudomolecules and the maps.
    This command will plot each reconstructed object (non-singleton).
    """
    p = OptionParser(plotall.__doc__)
    p.add_option("--links", default=10, type="int",
                 help="Only plot matchings more than")
    opts, args, iopts = p.set_image_options(args, figsize="10x6")

    if len(args) != 3:
        sys.exit(not p.print_help())

    mapsbed, agpfile, weightsfile = args
    agp = AGP(agpfile)
    objects = [ob for ob, lines in agp.iter_object() if len(lines) > 1]
    for seqid in sorted(objects):
        plot([seqid, mapsbed, agpfile, weightsfile,
              "--links={0}".format(opts.links)])


if __name__ == '__main__':
    main()
