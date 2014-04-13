#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scaffold Ordering with Weighted Maps.
"""

import sys
import logging

import numpy as np

from itertools import combinations
from collections import defaultdict
from scipy.stats import spearmanr

from jcvi.algorithms.matrix import determine_signs
from jcvi.algorithms.formula import reject_outliers
from jcvi.algorithms.tsp import hamiltonian
from jcvi.formats.agp import order_to_agp, build as agp_build
from jcvi.formats.base import must_open
from jcvi.formats.bed import Bed, BedLine, sort
from jcvi.formats.chain import fromagp
from jcvi.formats.sizes import Sizes
from jcvi.utils.grouper import Grouper
from jcvi.utils.counter import Counter
from jcvi.apps.base import OptionParser, ActionDispatcher, debug, sh, \
            need_update
debug()


START, END = "START", "END"


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
    def __init__(self, lgs, scaffolds, mapc, pivot, weights,
                 function=(lambda x: x.rank), precision=0):

        self.lgs = lgs
        self.lengths = mapc.compute_lengths(function)
        self.bins = mapc.get_bins(function)
        self.function = function
        self.precision = precision

        signs, flip = self.assign_orientation(scaffolds, pivot, weights)
        if flip:
            signs = - signs
        scaffolds = zip(scaffolds, signs)

        tour = self.assign_order(scaffolds, pivot, weights)
        self.tour = tour

        for lg in self.lgs:
            if lg.split("-")[0] == pivot:
                self.object = lg
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

    def assign_order(self, scaffolds, pivot, weights):
        """
        The goal is to assign scaffold orders. To help order the scaffolds, two
        dummy node, START and END, mark the ends of the chromosome. We connect
        START to each scaffold (directed), and each scaffold to END.
        """
        distances = defaultdict(list)
        f = self.function
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            length = self.lengths[lg]
            positions = {}
            for a, ao in scaffolds:
                xa = self.get_markers(lg, a, ao)
                if not xa:
                    continue
                d = np.median([f(x) for x in xa])
                assert 0 <= d <= length
                positions[a] = d

            for (a, ao), (b, bo) in combinations(scaffolds, 2):
                if a not in positions or b not in positions:
                    continue
                d = abs(positions[a] - positions[b])
                distances[a, b].append((d, mapname))
                distances[b, a].append((d, mapname))

            # Only connect ends for the pivot map, this ensure that the starting
            # and ending scaffold is always part of pivot map. The terminal
            # scaffolds for the less weighted maps are 'folded' in.
            if mapname != pivot:
                continue

            for a, ao in scaffolds:
                if a not in positions:
                    continue
                d = positions[a]
                distances[START, a].append((d, mapname))
                distances[a, END].append((length - d, mapname))

        for e, v in distances.items():
            distances[e] = self.weighted_mean(v, weights)

        distance_edges = sorted((a, b, w) for (a, b), w in distances.items())
        tour = hamiltonian(distance_edges, symmetric=False, precision=self.precision)
        # Remove dummy nodes
        assert tour[0] == START and tour[-1] == END
        tour = tour[1:-1]
        assert len(tour) == len(scaffolds), \
                "Tour ({0}) != Scaffolds ({1})".format(len(tour), len(scaffolds))

        scaffolds_oo = dict(scaffolds)
        recode = {0: '?', 1: '+', -1: '-'}
        tour = [(x, recode[scaffolds_oo[x]]) for x in tour]
        return tour

    def assign_orientation(self, scaffolds, pivot, weights):
        bins = self.bins
        f = self.function
        signs = defaultdict(list)
        for lg in self.lgs:
            mapname = lg.split("-")[0]
            oo = []
            nmarkers = {}
            for s in scaffolds:
                xs = bins.get((lg, s), [])
                nmarkers[s] = len(xs)
                if not xs:
                    oo.append(0)
                    continue
                physical = [x.pos for x in xs]
                cm = [f(x) for x in xs]
                rho, p_value = spearmanr(physical, cm)
                if np.isnan(rho):
                    rho = 0
                oo.append(rho)

            if mapname == pivot:
                pivot_oo = oo

            for i, j in combinations(range(len(scaffolds)), 2):
                si, sj = scaffolds[i], scaffolds[j]
                ni, nj = nmarkers[si], nmarkers[sj]
                orientation = oo[i] * oo[j]
                if not orientation:
                    continue
                d = np.sign(orientation) * ni * nj
                signs[(i, j)].append((d, mapname))

        for e, v in signs.items():
            signs[e] = self.weighted_mean(v, weights)

        signs_edges = sorted((a, b, w) for (a, b), w in signs.items())
        signs = determine_signs(scaffolds, signs_edges)

        # Finally flip this according to pivot map, then weight by #_markers
        nmarkers = [nmarkers[x] for x in scaffolds]
        flipr = signs * np.sign(np.array(pivot_oo)) * nmarkers
        flip = sum(flipr) < 0
        return signs, flip


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

    def __init__(self, bedfile):
        bed = Bed(bedfile)
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


def main():

    actions = (
        ('merge', 'merge csv maps and convert to bed format'),
        ('path', 'construct golden path given a set of genetic maps'),
        ('build', 'build associated FASTA and CHAIN file'),
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
    fw = open(weightsfile, "w")
    for mapname in sorted(mapnames):
        weight = 1
        print >> fw, mapname, weight
    logging.debug("Weights file written to `{0}`.".format(weightsfile))


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
    # Read in the weights
    weights = {}
    fp = open(weightsfile)
    for row in fp:
        mapname, w = row.split()
        weights[mapname] = int(w)

    pivot_weight, pivot = max((w, m) for m, w in weights.items())
    logging.debug("Pivot map: `{0}` (weight={1}).".format(pivot, pivot_weight))
    gapsize = opts.gapsize
    function = opts.distance
    precision = 3 if function == "cM" else 0
    function = (lambda x: x.cm) if function == "cM" else \
               (lambda x: x.rank)

    cc = Map(bedfile)
    mapnames = cc.mapnames
    allseqids = cc.seqids
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
            best_neighbor = max(neighbors, key=lambda x: x[-1])
            C.join(n, best_neighbor[0])

    partitions = defaultdict(list)
    # Partition the scaffolds and assign them to one consensus
    for s in allseqids:
        s = Scaffold(s, cc)
        seqid = s.seqid
        counts = {}
        for mlg, count in s.mlg_counts.items():
            consensus = C[mlg]
            if consensus not in counts:
                counts[consensus] = 0
            counts[consensus] += count
        best_assignment = max(counts.items(), key=lambda x: x[1])
        best_consensus, best_value = best_assignment
        if counts.values().count(best_value) > 1:  # tie
            print >> sys.stderr, "AMBIGUOUS", seqid, counts
            continue
        partitions[best_consensus].append(seqid)

    # Perform OO within each partition
    agpfile = bedfile.rsplit(".", 1)[0] + ".agp"
    sizes = Sizes(fastafile).mapping
    fwagp = must_open(agpfile, "w")
    for lgs, scaffolds in sorted(partitions.items()):
        print >> sys.stderr, lgs
        s = ScaffoldOO(lgs, scaffolds, cc, pivot, weights,
                       function=function, precision=precision)
        print >> sys.stderr, s.tour
        order_to_agp(s.object, s.tour, sizes, fwagp, gapsize=gapsize,
                     gaptype="map")
    fwagp.close()

    logging.debug("AGP file written to `{0}`.".format(agpfile))


def build(args):
    """
    %prog build agpfile scaffolds.fasta genome.fasta map.bed

    Build associated genome FASTA file and CHAIN file that can be used to lift
    old coordinates to new coordinates. The CHAIN file will be used to lift the
    original marker positions to new positions in the reconstructed genome.
    """
    p = OptionParser(build.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    agpfile, scaffolds, genome, mapbed = args
    if need_update((agpfile, scaffolds), genome):
        agp_build([agpfile, scaffolds, genome])

    chainfile = agpfile.rsplit(".", 1)[0] + ".chain"
    if need_update((agpfile, scaffolds, genome), chainfile):
        fromagp([agpfile, scaffolds, genome])

    liftedbed = mapbed.rsplit(".", 1)[0] + ".lifted.bed"
    if need_update((mapbed, chainfile), liftedbed):
        cmd = "liftOver -minMatch=1 {0} {1} {2} unmapped".\
                format(mapbed, chainfile, liftedbed)
        sh(cmd)

    sort([liftedbed, "-i"])  # Sort bed in place


if __name__ == '__main__':
    main()
