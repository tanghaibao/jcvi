#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scaffold Ordering with Weighted Maps.
"""

import os.path as op
import sys
import logging

import numpy as np
import networkx as nx

from itertools import combinations
from collections import defaultdict
from scipy.stats import spearmanr

from jcvi.algorithms.formula import reject_outliers
from jcvi.algorithms.lis import longest_monotonous_subseq_length
from jcvi.algorithms.tsp import hamiltonian
from jcvi.algorithms.matrix import determine_signs
from jcvi.formats.agp import AGP, order_to_agp, build as agp_build
from jcvi.formats.base import DictFile, FileMerger, must_open
from jcvi.formats.bed import Bed, BedLine, sort
from jcvi.formats.chain import fromagp
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import human_size
from jcvi.utils.counter import Counter
from jcvi.utils.grouper import Grouper
from jcvi.utils.iter import flatten
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


class LinkageGroup (object):

    def __init__(self, lg, markers, function=(lambda x: x.rank)):
        # Three dicts that are keyed by scaffold id
        self.lg = lg
        self.markers = markers
        self.series = dict((k, [function(x) for x in v]) \
                            for k, v in markers.items())
        self.position = dict((k, np.median(v)) \
                            for k, v in self.series.items())
        self.guide = dict((k, np.median([x.cm for x in v])) \
                            for k, v in markers.items())


class ScaffoldOO (object):
    """
    This contains the routine to construct order and orientation for the
    scaffolds per partition.
    """
    def __init__(self, lgs, scaffolds, mapc, pivot, weights, sizes,
                 function=(lambda x: x.rank), cutoff=.01):

        self.lgs = lgs
        self.lengths = mapc.compute_lengths(function)
        self.bins = mapc.get_bins(function)
        self.function = function
        self.sizes = sizes

        signs = self.assign_orientation(scaffolds, pivot, weights)
        scaffolds = zip(scaffolds, signs)

        tour, linkage_groups = self.assign_order(scaffolds, pivot, weights,
                                 cutoff=cutoff)
        tour = self.fix_orientation(tour, linkage_groups, weights)

        recode = {0: '?', 1: '+', -1: '-'}
        tour = [(x, recode[o]) for x, o in tour]
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

    def get_markers(self, lg, scaffold, orientation=0):
        xs = self.bins.get((lg, scaffold), [])
        if orientation < 0:
            xs = xs[::-1]
        return xs

    def get_series(self, lg, scaffold, orientation=0):
        xs = self.get_markers(lg, scaffold, orientation=orientation)
        return [self.function(x) for x in xs]

    def get_rho(self, xy):
        if not xy:
            return 0
        x, y = zip(*xy)
        rho, p_value = spearmanr(x, y)
        if np.isnan(rho):
            rho = 0
        return rho

    def assign_order(self, scaffolds, pivot, weights, cutoff=.01):
        """
        The goal is to assign scaffold orders. To help order the scaffolds, two
        dummy node, START and END, mark the ends of the chromosome. We connect
        START to each scaffold (directed), and each scaffold to END.
        """
        f = self.function
        linkage_groups = []
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            length = self.lengths[lg]
            markers = {}
            position = {}
            guide = {}  # cM guide to ward against 'close binning'
            for a, ao in scaffolds:
                xa = self.get_markers(lg, a, orientation=ao)
                if not xa:
                    continue
                d = np.median([f(x) for x in xa])
                assert 0 <= d <= length
                markers[a] = xa
                position[a] = d

            if mapname == pivot:
                pivot_position = position

            LG = LinkageGroup(lg, markers, function=f)
            linkage_groups.append(LG)

        paths = []
        rhos = []
        for lg in linkage_groups:
            position = lg.position
            guide = lg.guide
            # Sort the order based on median distance
            path = sorted((v, guide[k], k) for k, v in position.items())
            vv, gg, path = zip(*path)

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
            rhos.append(rho)

            print lg.lg
            print lg.position
            print lg.guide
            print path

        # Preparation of TSP
        distances = defaultdict(list)
        for path, rho, lg in zip(paths, rhos, linkage_groups):
            mapname = lg.lg.split("-")[0]
            position = lg.position
            length = self.lengths[lg.lg]
            for a, b in combinations(path, 2):
                d = abs(position[a] - position[b])
                distances[a, b].append((d, mapname))
            for p in path:
                adist, bdist = position[p], length - position[p]
                if rho < 0:
                    adist, bdist = bdist, adist
                distances[START, p].append((adist, mapname))
                distances[p, END].append((bdist, mapname))

        G = nx.DiGraph()
        for (a, b), v in distances.items():
            d = self.weighted_mean(v, weights)
            G.add_edge(a, b, weight=d)
            if a == START or b == END:
                continue
            G.add_edge(b, a, weight=d)

        logging.debug("Graph size: |V|={0}, |E|={1}.".format(len(G), G.size()))

        L = nx.all_pairs_dijkstra_path_length(G)
        for (a, ao), (b, bo) in combinations(scaffolds, 2):
            if G.has_edge(a, b):
                continue
            l = L[a][b]
            G.add_edge(a, b, weight=l)
            G.add_edge(b, a, weight=l)

        edges = []
        for a, b, d in G.edges(data=True):
            edges.append((a, b, d['weight']))

        tour = hamiltonian(edges, symmetric=False, precision=2)
        print tour
        assert tour[0] == START and tour[-1] == END
        tour = tour[1:-1]

        scaffolds_oo = dict(scaffolds)
        tour = [(x, scaffolds_oo[x]) for x in tour]
        return tour, linkage_groups

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
        f = self.function
        signs = defaultdict(list)
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            oo = []
            nmarkers = {}
            series = []
            for s in scaffolds:
                xs = self.get_markers(lg, s)
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
                pivot_nmarkers = nmarkers

            for i, j in combinations(range(len(scaffolds)), 2):
                si, sj = series[i], series[j]
                d = self.get_orientation(si, sj)

                if not d:
                    continue
                signs[i, j].append((d, mapname))

        for e, v in signs.items():
            signs[e] = self.weighted_mean(v, weights)

        signs_edges = sorted((a, b, w) for (a, b), w in signs.items())
        signs = determine_signs(scaffolds, signs_edges)

        # Finally flip this according to pivot map, then weight by #_markers
        nmarkers = [pivot_nmarkers[x] for x in scaffolds]
        flipr = signs * np.sign(np.array(pivot_oo)) * nmarkers
        print "scaffolds", scaffolds
        print "nmarkers", nmarkers
        print "pivot_oo", pivot_oo
        print "signs", signs
        print "flipr", flipr
        print "sum", sum(flipr)
        if sum(flipr) < 0:
            signs = - signs
        return signs

    def fix_orientation(self, tour, linkage_groups, weights):
        """
        Test each scaffold if flipping will increass longest monotonous chain
        length.
        """
        orientations = dict(tour)  # old configuration here
        scaffold_oo = defaultdict(list)
        scaffolds, oos = zip(*tour)
        for mlg in linkage_groups:
            lg = mlg.lg
            mapname = lg.split("-")[0]
            for s, o in tour:
                i = scaffolds.index(s)
                L = [self.get_series(lg, x, xo) for x, xo in tour[:i]]
                U = [self.get_series(lg, x, xo) for x, xo in tour[i + 1:]]
                L, U = list(flatten(L)), list(flatten(U))
                M = self.get_series(lg, s)
                plus = longest_monotonous_subseq_length(L + M + U)
                minus = longest_monotonous_subseq_length(L + M[::-1] + U)
                print lg, s, o, M, plus, minus
                d = plus[0] - minus[0]
                scaffold_oo[s].append((d, mapname))  # reset orientation

        for s, v in scaffold_oo.items():
            d = self.weighted_mean(v, weights)
            if abs(d) > .5:  # only update when there is good evidence
                orientations[s] = np.sign(d)
            print s, v, d

        tour = [(x, orientations[x]) for x in scaffolds]
        return tour


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
        self.report()

    def report(self):
        self.nmarkers = len(self)
        self.seqids = sorted(set(x.seqid for x in self))
        self.mapnames = sorted(set(x.mapname for x in self))
        self.mlgs = sorted(set(x.mlg for x in self))
        logging.debug("Map contains {0} markers in {1} linkage groups.".\
                      format(self.nmarkers, len(self.mlgs)))
        self.ranks = self.compute_ranks()

    def extract(self, seqid):
        r = [x for x in self if x.seqid == seqid]
        return sorted(r, key=lambda x: x.pos)

    def extract_mlg(self, mlg):
        r = [x for x in self if x.mlg == mlg]
        return sorted(r, key=lambda x: x.cm)

    def compute_ranks(self):
        ranks = {}  # Store the length for each linkage group
        for mlg in self.mlgs:
            mlg_set = self.extract_mlg(mlg)
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


class Weights (DictFile):

    def __init__(self, filename, mapnames, cast=int):
        super(Weights, self).__init__(filename, cast=cast)
        self.maps = [x.split()[0] for x in must_open(filename)]
        self.update_maps(mapnames)

    def update_maps(self, mapnames, default=1):
        for m in mapnames:
            if m in self:
                continue
            self[m] = default
            logging.debug("Weight for `{0}` set to {1}.".format(m, default))

    def get_pivot(self, mapnames):
        # Break ties by occurence in file
        return max((w, -self.maps.index(m), m) \
                    for m, w in self.items() if m in mapnames)


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
    p.add_option("--cutoff", default=.01, type="float",
                 help="Down-weight distance <= cM apart, 0 to disable")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, weightsfile, fastafile = args
    gapsize = opts.gapsize
    function = opts.distance
    cutoff = opts.cutoff

    cc = Map(bedfile)
    mapnames = cc.mapnames
    allseqids = cc.seqids
    weights = Weights(weightsfile, mapnames)
    pivot_weight, o, pivot = weights.get_pivot(mapnames)

    logging.debug("Pivot map: `{0}` (weight={1}).".format(pivot, pivot_weight))
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
                       function=function, cutoff=cutoff)

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

    rhos = {}
    # Parallel coordinates
    for mlg, (x, y1, y2) in coords.items():
        mm = cc.extract_mlg(mlg)
        markers = [(m.accn, m.cm) for m in mm]  # exhaustive marker list
        xy = [(m.pos, m.cm) for m in mm if m.seqid == seqid]
        mx, my = zip(*xy)
        rho, p_value = spearmanr(mx, my)
        rhos[mlg] = rho
        flip = rho < 0

        g = GeneticMap(ax1, x, y1, y2, markers, tip=tip, flip=flip)
        extra = -3 * tip if x < .5 else 3 * tip
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
    pp = [x.object_beg for x in agp if not x.is_gap]
    patchstart = [f(x) for x in pp]
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
        ax.vlines(pp, [0], [mlgsize], colors="snow")
        ax.plot(xx, yy, ".", color=color)
        ax.text(.5, 1 - .5 * gap / height, r"$\rho$={0:.3f}".format(rhos[mlg]),
                    ha="center", va="top", transform=ax.transAxes, color="gray")
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
