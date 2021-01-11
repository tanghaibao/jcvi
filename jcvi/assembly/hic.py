#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Process Hi-C output into AGP for chromosomal-scale scaffolding.
"""
from __future__ import print_function

import array
import json
import logging
import math
import os
import os.path as op
import sys
from collections import defaultdict
from functools import partial
from multiprocessing import Pool

import numpy as np

from natsort import natsorted

from jcvi.algorithms.ec import GA_run, GA_setup
from jcvi.algorithms.formula import outlier_cutoff
from jcvi.algorithms.matrix import get_signs
from jcvi.apps.base import ActionDispatcher, OptionParser, backup, iglob, mkdir, symlink
from jcvi.apps.grid import Jobs
from jcvi.assembly.allmaps import make_movie
from jcvi.compara.synteny import check_beds, get_bed_filenames
from jcvi.formats.agp import order_to_agp
from jcvi.formats.base import LineFile, must_open
from jcvi.formats.bed import Bed
from jcvi.formats.blast import Blast
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import (
    markup,
    normalize_axes,
    plt,
    savefig,
    ticker,
    human_readable,
)
from jcvi.graphics.dotplot import dotplot
from jcvi.utils.cbook import gene_name, human_size

# Map orientations to ints
FF = {"+": 1, "-": -1, "?": 1}
RR = {"+": -1, "-": 1, "?": -1}
LB = 18  # Lower bound for golden_array()
UB = 29  # Upper bound for golden_array()
BB = UB - LB + 1  # Span for golden_array()
ACCEPT = "[green]ACCEPT"
REJECT = "[red]REJECT"
BINSIZE = 50000


class ContigOrderingLine(object):
    """Stores one line in the ContigOrdering file"""

    def __init__(self, line, sep="|"):
        args = line.split()
        self.contig_id = args[0]
        self.contig_name = args[1].split(sep)[0]
        contig_rc = args[2]
        assert contig_rc in ("0", "1")
        self.strand = "+" if contig_rc == "0" else "-"
        self.orientation_score = args[3]
        self.gap_size_after_contig = args[4]


class ContigOrdering(LineFile):
    """ContigOrdering file as created by LACHESIS, one per chromosome group.
    Header contains summary information per group, followed by list of contigs
    with given ordering.
    """

    def __init__(self, filename):
        super(ContigOrdering, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            if row[0] == "#":
                continue
            orderline = ContigOrderingLine(row)
            self.append(orderline)

    def write_agp(
        self, obj, sizes, fw=sys.stdout, gapsize=100, gaptype="contig", evidence="map"
    ):
        """Converts the ContigOrdering file into AGP format"""
        contigorder = [(x.contig_name, x.strand) for x in self]
        order_to_agp(
            obj,
            contigorder,
            sizes,
            fw,
            gapsize=gapsize,
            gaptype=gaptype,
            evidence=evidence,
        )


class CLMFile:
    """CLM file (modified) has the following format:

    tig00046211+ tig00063795+       1       53173
    tig00046211+ tig00063795-       1       116050
    tig00046211- tig00063795+       1       71155
    tig00046211- tig00063795-       1       134032
    tig00030676+ tig00077819+       5       136407 87625 87625 106905 102218
    tig00030676+ tig00077819-       5       126178 152952 152952 35680 118923
    tig00030676- tig00077819+       5       118651 91877 91877 209149 125906
    tig00030676- tig00077819-       5       108422 157204 157204 137924 142611
    """

    def __init__(self, clmfile, skiprecover=False):
        self.name = op.basename(clmfile).rsplit(".", 1)[0]
        self.clmfile = clmfile
        self.idsfile = clmfile.rsplit(".", 1)[0] + ".ids"
        self.parse_ids(skiprecover)
        self.parse_clm()
        self.signs = None

    def parse_ids(self, skiprecover):
        """IDS file has a list of contigs that need to be ordered. 'recover',
        keyword, if available in the third column, is less confident.

        tig00015093     46912
        tig00035238     46779   recover
        tig00030900     119291
        """
        idsfile = self.idsfile
        logging.debug("Parse idsfile `{}`".format(idsfile))
        fp = open(idsfile)
        tigs = []
        for row in fp:
            if row[0] == "#":  # Header
                continue
            atoms = row.split()
            tig, _, size = atoms
            size = int(size)
            if skiprecover and len(atoms) == 3 and atoms[2] == "recover":
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
        contacts_oriented = defaultdict(dict)
        orientations = defaultdict(list)
        for row in fp:
            atoms = row.strip().split("\t")
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
            contacts[(at, bt)] = len(dists)
            gdists = golden_array(dists)
            contacts_oriented[(at, bt)][(FF[ao], FF[bo])] = gdists
            contacts_oriented[(bt, at)][(RR[bo], RR[ao])] = gdists
            strandedness = 1 if ao == bo else -1
            orientations[(at, bt)].append((strandedness, dists))

        self.contacts = contacts
        self.contacts_oriented = contacts_oriented
        # Preprocess the orientations dict
        for (at, bt), dists in orientations.items():
            dists = [(s, d, hmean_int(d)) for (s, d) in dists]
            strandedness, md, mh = min(dists, key=lambda x: x[-1])
            orientations[(at, bt)] = (strandedness, len(md), mh)
        self.orientations = orientations

    def calculate_densities(self):
        """
        Calculate the density of inter-contig links per base. Strong contigs
        considered to have high level of inter-contig links in the current
        partition.
        """
        active = self.active
        densities = defaultdict(int)
        for (at, bt), links in self.contacts.items():
            if not (at in active and bt in active):
                continue
            densities[at] += links
            densities[bt] += links

        logdensities = {}
        for x, d in densities.items():
            s = self.tig_to_size[x]
            logd = np.log10(d * 1.0 / min(s, 500000))
            logdensities[x] = logd

        return logdensities

    def report_active(self):
        logging.debug(
            "Active contigs: {} (length={})".format(self.N, self.active_sizes.sum())
        )

    def activate(self, tourfile=None, minsize=10000, backuptour=True):
        """
        Select contigs in the current partition. This is the setup phase of the
        algorithm, and supports two modes:

        - "de novo": This is useful at the start of a new run where no tours
          available. We select the strong contigs that have significant number
          of links to other contigs in the partition. We build a histogram of
          link density (# links per bp) and remove the contigs that appear as
          outliers. The orientations are derived from the matrix decomposition
          of the pairwise strandedness matrix O.

        - "hotstart": This is useful when there was a past run, with a given
          tourfile. In this case, the active contig list and orientations are
          derived from the last tour in the file.
        """
        if tourfile and (not op.exists(tourfile)):
            logging.debug("Tourfile `{}` not found".format(tourfile))
            tourfile = None

        if tourfile:
            logging.debug("Importing tourfile `{}`".format(tourfile))
            tour, tour_o = iter_last_tour(tourfile, self)
            self.active = set(tour)
            tig_to_idx = self.tig_to_idx
            tour = [tig_to_idx[x] for x in tour]
            signs = sorted([(x, FF[o]) for (x, o) in zip(tour, tour_o)])
            _, signs = zip(*signs)
            self.signs = np.array(signs, dtype=int)
            if backuptour:
                backup(tourfile)
            tour = array.array("i", tour)
        else:
            self.report_active()
            while True:
                logdensities = self.calculate_densities()
                lb, ub = outlier_cutoff(list(logdensities.values()))
                logging.debug("Log10(link_densities) ~ [{}, {}]".format(lb, ub))
                remove = set(
                    x
                    for x, d in logdensities.items()
                    if (d < lb and self.tig_to_size[x] < minsize * 10)
                )
                if remove:
                    self.active -= remove
                    self.report_active()
                else:
                    break

            logging.debug("Remove contigs with size < {}".format(minsize))
            self.active = set(x for x in self.active if self.tig_to_size[x] >= minsize)
            tour = range(self.N)  # Use starting (random) order otherwise
            tour = array.array("i", tour)

            # Determine orientations
            self.flip_all(tour)

        self.report_active()
        self.tour = tour

        return tour

    def evaluate_tour_M(self, tour):
        """Use Cythonized version to evaluate the score of a current tour"""
        from .chic import score_evaluate_M

        return score_evaluate_M(tour, self.active_sizes, self.M)

    def evaluate_tour_P(self, tour):
        """Use Cythonized version to evaluate the score of a current tour,
        with better precision on the distance of the contigs.
        """
        from .chic import score_evaluate_P

        return score_evaluate_P(tour, self.active_sizes, self.P)

    def evaluate_tour_Q(self, tour):
        """Use Cythonized version to evaluate the score of a current tour,
        taking orientation into consideration. This may be the most accurate
        evaluation under the right condition.
        """
        from .chic import score_evaluate_Q

        return score_evaluate_Q(tour, self.active_sizes, self.Q)

    def flip_log(self, method, score, score_flipped, tag):
        logging.debug("{}: {} => {} {}".format(method, score, score_flipped, tag))

    def flip_all(self, tour):
        """Initialize the orientations based on pairwise O matrix."""
        if self.signs is None:  # First run
            score = 0
        else:
            old_signs = self.signs[: self.N]
            (score,) = self.evaluate_tour_Q(tour)

        # Remember we cannot have ambiguous orientation code (0 or '?') here
        self.signs = get_signs(self.O, validate=False, ambiguous=False)
        (score_flipped,) = self.evaluate_tour_Q(tour)
        if score_flipped >= score:
            tag = ACCEPT
        else:
            self.signs = old_signs[:]
            tag = REJECT
        self.flip_log("FLIPALL", score, score_flipped, tag)
        return tag

    def flip_whole(self, tour):
        """Test flipping all contigs at the same time to see if score improves."""
        (score,) = self.evaluate_tour_Q(tour)
        self.signs = -self.signs
        (score_flipped,) = self.evaluate_tour_Q(tour)
        if score_flipped > score:
            tag = ACCEPT
        else:
            self.signs = -self.signs
            tag = REJECT
        self.flip_log("FLIPWHOLE", score, score_flipped, tag)
        return tag

    def flip_one(self, tour):
        """Test flipping every single contig sequentially to see if score
        improves.
        """
        n_accepts = n_rejects = 0
        any_tag_ACCEPT = False
        for i, t in enumerate(tour):
            if i == 0:
                (score,) = self.evaluate_tour_Q(tour)
            self.signs[t] = -self.signs[t]
            (score_flipped,) = self.evaluate_tour_Q(tour)
            if score_flipped > score:
                n_accepts += 1
                tag = ACCEPT
            else:
                self.signs[t] = -self.signs[t]
                n_rejects += 1
                tag = REJECT
            self.flip_log(
                "FLIPONE ({}/{})".format(i + 1, len(self.signs)),
                score,
                score_flipped,
                tag,
            )
            if tag == ACCEPT:
                any_tag_ACCEPT = True
                score = score_flipped
        logging.debug("FLIPONE: N_accepts={} N_rejects={}".format(n_accepts, n_rejects))
        return ACCEPT if any_tag_ACCEPT else REJECT

    def prune_tour(self, tour, cpus):
        """Test deleting each contig and check the delta_score; tour here must
        be an array of ints.
        """
        while True:
            (tour_score,) = self.evaluate_tour_M(tour)
            logging.debug("Starting score: {}".format(tour_score))
            active_sizes = self.active_sizes
            M = self.M
            args = []
            for i, t in enumerate(tour):
                stour = tour[:i] + tour[i + 1 :]
                args.append((t, stour, tour_score, active_sizes, M))

            # Parallel run
            p = Pool(processes=cpus)
            results = list(p.imap(prune_tour_worker, args))
            assert len(tour) == len(
                results
            ), "Array size mismatch, tour({}) != results({})".format(
                len(tour), len(results)
            )

            # Identify outliers
            active_contigs = self.active_contigs
            idx, log10deltas = zip(*results)
            lb, ub = outlier_cutoff(log10deltas)
            logging.debug("Log10(delta_score) ~ [{}, {}]".format(lb, ub))

            remove = set(active_contigs[x] for (x, d) in results if d < lb)
            self.active -= remove
            self.report_active()

            tig_to_idx = self.tig_to_idx
            tour = [active_contigs[x] for x in tour]
            tour = array.array("i", [tig_to_idx[x] for x in tour if x not in remove])
            if not remove:
                break

        self.tour = tour
        self.flip_all(tour)

        return tour

    @property
    def active_contigs(self):
        return list(self.active)

    @property
    def active_sizes(self):
        return np.array([self.tig_to_size[x] for x in self.active])

    @property
    def N(self):
        return len(self.active)

    @property
    def oo(self):
        return range(self.N)

    @property
    def tig_to_idx(self):
        return dict((x, i) for (i, x) in enumerate(self.active))

    @property
    def M(self):
        """
        Contact frequency matrix. Each cell contains how many inter-contig
        links between i-th and j-th contigs.
        """
        N = self.N
        tig_to_idx = self.tig_to_idx
        M = np.zeros((N, N), dtype=int)
        for (at, bt), links in self.contacts.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            M[ai, bi] = M[bi, ai] = links
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
        for (at, bt), (strandedness, md, mh) in self.orientations.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            score = strandedness * md
            O[ai, bi] = O[bi, ai] = score
        return O

    @property
    def P(self):
        """
        Contact frequency matrix with better precision on distance between
        contigs. In the matrix M, the distance is assumed to be the distance
        between mid-points of two contigs. In matrix Q, however, we compute
        harmonic mean of the links for the orientation configuration that is
        shortest. This offers better precision for the distance between big
        contigs.
        """
        N = self.N
        tig_to_idx = self.tig_to_idx
        P = np.zeros((N, N, 2), dtype=int)
        for (at, bt), (strandedness, md, mh) in self.orientations.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            P[ai, bi, 0] = P[bi, ai, 0] = md
            P[ai, bi, 1] = P[bi, ai, 1] = mh
        return P

    @property
    def Q(self):
        """
        Contact frequency matrix when contigs are already oriented. This is s a
        similar matrix as M, but rather than having the number of links in the
        cell, it points to an array that has the actual distances.
        """
        N = self.N
        tig_to_idx = self.tig_to_idx
        signs = self.signs
        Q = np.ones((N, N, BB), dtype=int) * -1  # Use -1 as the sentinel
        for (at, bt), k in self.contacts_oriented.items():
            if not (at in tig_to_idx and bt in tig_to_idx):
                continue
            ai = tig_to_idx[at]
            bi = tig_to_idx[bt]
            ao = signs[ai]
            bo = signs[bi]
            Q[ai, bi] = k[(ao, bo)]
        return Q


def hmean_int(a, a_min=5778, a_max=1149851):
    """Harmonic mean of an array, returns the closest int"""
    from scipy.stats import hmean

    return int(round(hmean(np.clip(a, a_min, a_max))))


def golden_array(a, phi=1.61803398875, lb=LB, ub=UB):
    """Given list of ints, we aggregate similar values so that it becomes an
    array of multiples of phi, where phi is the golden ratio.

    phi ^ 14 = 843
    phi ^ 33 = 7881196

    So the array of counts go between 843 to 788196. One triva is that the
    exponents of phi gets closer to integers as N grows. See interesting
    discussion here:
    <https://www.johndcook.com/blog/2017/03/22/golden-powers-are-nearly-integers/>
    """
    counts = np.zeros(BB, dtype=int)
    for x in a:
        c = int(round(math.log(x, phi)))
        if c < lb:
            c = lb
        if c > ub:
            c = ub
        counts[c - lb] += 1
    return counts


def prune_tour_worker(arg):
    """Worker thread for CLMFile.prune_tour()"""
    from .chic import score_evaluate_M

    t, stour, tour_score, active_sizes, M = arg
    (stour_score,) = score_evaluate_M(stour, active_sizes, M)
    delta_score = tour_score - stour_score
    log10d = np.log10(delta_score) if delta_score > 1e-9 else -9
    return t, log10d


def main():

    actions = (
        # LACHESIS output processing
        ("agp", "generate AGP file based on LACHESIS output"),
        ("score", "score the current LACHESIS CLM"),
        # Simulation
        ("simulate", "simulate CLM data"),
        # Scaffolding
        ("optimize", "optimize the contig order and orientation"),
        ("density", "estimate link density of contigs"),
        # Plotting
        ("movieframe", "plot heatmap and synteny for a particular tour"),
        ("movie", "plot heatmap optimization history in a tourfile"),
        # Reference-based analytics
        ("bam2mat", "convert bam file to .npy format used in plotting"),
        ("mergemat", "combine counts from multiple .npy data files"),
        ("heatmap", "plot heatmap based on .npy file"),
        ("dist", "plot distance distribution based on .dist.npy file"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fit_power_law(xs, ys):
    """Fit power law distribution.

    See reference:
    http://mathworld.wolfram.com/LeastSquaresFittingPowerLaw.html
    Assumes the form Y = A * X^B, returns

    Args:
        xs ([int]): X vector
        ys ([float64]): Y vector

    Returns:
        (A, B), the coefficients
    """
    import math

    sum_logXlogY, sum_logXlogX, sum_logX, sum_logY = 0, 0, 0, 0
    N = len(xs)
    for i in range(N):
        if not xs[i] or not ys[i]:
            continue
        logXs, logYs = math.log(xs[i]), math.log(ys[i])
        sum_logXlogY += logXs * logYs
        sum_logXlogX += logXs * logXs
        sum_logX += logXs
        sum_logY += logYs

    B = (N * sum_logXlogY - sum_logX * sum_logY) / (
        N * sum_logXlogX - sum_logX * sum_logX
    )
    A = math.exp((sum_logY - B * sum_logX) / N)
    logging.debug("Power law Y = {:.1f} * X ^ {:.4f}".format(A, B))
    label = "$Y={:.1f} \\times X^{{ {:.4f} }}$".format(A, B)
    return A, B, label


def dist(args):
    """
    %prog dist input.dist.npy genome.json

    Plot histogram based on .dist.npy data file. The .npy file stores an array
    with link counts per dist bin, with the bin starts stored in the genome.json.
    """
    import seaborn as sns
    import pandas as pd
    from jcvi.graphics.base import human_base_formatter, markup

    p = OptionParser(dist.__doc__)
    p.add_option("--title", help="Title of the histogram")
    p.add_option("--xmin", default=300, help="Minimum distance")
    p.add_option("--xmax", default=6000000, help="Maximum distance")
    opts, args, iopts = p.set_image_options(args, figsize="6x6")

    if len(args) != 2:
        sys.exit(not p.print_help())

    npyfile, jsonfile = args
    pf = npyfile.rsplit(".", 1)[0]
    header = json.loads(open(jsonfile).read())
    distbin_starts = np.array(header["distbinstarts"], dtype="float64")
    distbin_sizes = np.array(header["distbinsizes"], dtype="float64")
    a = np.load(npyfile)

    xmin, xmax = opts.xmin, opts.xmax
    (size,) = min(distbin_sizes.shape, distbin_starts.shape, a.shape)
    df = pd.DataFrame()
    xstart, xend = (
        np.searchsorted(distbin_starts, xmin),
        np.searchsorted(distbin_starts, xmax),
    )
    df["BinStart"] = distbin_starts[xstart:xend]
    df["LinkDensity"] = a[xstart:xend] / distbin_sizes[xstart:xend]
    ax = sns.lineplot(
        x="BinStart", y="LinkDensity", data=df, lw=3, color="lightslategray"
    )
    tx = df["BinStart"]
    A, B, label = fit_power_law(tx, df["LinkDensity"])
    ty = A * tx ** B
    ax.plot(tx, ty, "r:", lw=3, label=label)
    ax.legend()
    if opts.title:
        ax.set_title(markup(opts.title))
    ax.set_xlabel("Link size (bp)")
    ax.set_ylabel("Density (\# of links per bp)")
    ax.set_xscale("log", nonposx="clip")
    ax.set_yscale("log", nonposy="clip")
    ax.xaxis.set_major_formatter(human_base_formatter)

    image_name = pf + "." + opts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def generate_groups(groupsfile):
    """Parse 'groups' file. The 'groups' file has the following format,
    for example:

    seq1,seq2 b
    seq1 g
    seq2 g

    Args:
        groupsfile (str): Path to the groups file
    """
    data = []
    with open(groupsfile) as fp:
        for row in fp:
            seqids, color = row.split()
            yield seqids, color


def heatmap(args):
    """
    %prog heatmap input.npy genome.json

    Plot heatmap based on .npy data file. The .npy stores a square matrix with
    bins of genome, and cells inside the matrix represent number of links
    between bin i and bin j. The `genome.json` contains the offsets of each
    contig/chr so that we know where to draw boundary lines, or extract per
    contig/chromosome heatmap.

    If a 'groups' file is given (with --groups), we will draw squares on the
    heatmap. The 'groups' file has the following format, for example:

    seq1,seq2 b
    seq1 g
    seq2 g

    This will first draw a square around seq1+seq2 with blue color, then seq1
    and seq2 individually with green color.
    """
    p = OptionParser(heatmap.__doc__)
    p.add_option("--title", help="Title of the heatmap")
    p.add_option("--groups", help="Groups file, see doc")
    p.add_option("--vmin", default=1, type="int", help="Minimum value in the heatmap")
    p.add_option("--vmax", default=6, type="int", help="Maximum value in the heatmap")
    p.add_option("--chr", help="Plot this contig/chr only")
    p.add_option(
        "--nobreaks",
        default=False,
        action="store_true",
        help="Do not plot breaks (esp. if contigs are small)",
    )
    opts, args, iopts = p.set_image_options(
        args, figsize="11x11", style="white", cmap="coolwarm", format="png", dpi=120
    )

    if len(args) != 2:
        sys.exit(not p.print_help())

    npyfile, jsonfile = args
    contig = opts.chr
    groups = list(generate_groups(opts.groups)) if opts.groups else []

    # Load contig/chromosome starts and sizes
    header = json.loads(open(jsonfile).read())
    resolution = header.get("resolution")
    assert resolution is not None, "`resolution` not found in `{}`".format(jsonfile)
    logging.debug("Resolution set to {}".format(resolution))
    # Load the matrix
    A = np.load(npyfile)

    # Select specific submatrix
    if contig:
        contig_start = header["starts"][contig]
        contig_size = header["sizes"][contig]
        contig_end = contig_start + contig_size
        A = A[contig_start:contig_end, contig_start:contig_end]

    # Convert seqids to positions for each group
    new_groups = []
    for seqids, color in groups:
        seqids = seqids.split(",")
        assert all(
            x in header["starts"] for x in seqids
        ), f"{seqids} contain ids not found in starts"
        assert all(
            x in header["sizes"] for x in seqids
        ), f"{seqids} contain ids not found in sizes"
        start = min(header["starts"][x] for x in seqids)
        end = max(header["starts"][x] + header["sizes"][x] for x in seqids)
        position_seqids = []
        for seqid in seqids:
            seqid_start = header["starts"][seqid]
            seqid_size = header["sizes"][seqid]
            position_seqids.append((seqid_start + seqid_size / 2, seqid))
        new_groups.append((start, end, position_seqids, color))

    # Several concerns in practice:
    # The diagonal counts may be too strong, this can either be resolved by
    # masking them. Or perform a log transform on the entire heatmap.
    B = A.astype("float64")
    B += 1.0
    B = np.log(B)
    vmin, vmax = opts.vmin, opts.vmax
    B[B < vmin] = vmin
    B[B > vmax] = vmax
    print(B)
    logging.debug(
        "Matrix log-transformation and thresholding ({}-{}) done".format(vmin, vmax)
    )

    # Canvas
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # whole canvas
    ax = fig.add_axes([0.05, 0.05, 0.9, 0.9])  # just the heatmap

    breaks = list(header["starts"].values())
    breaks += [header["total_bins"]]  # This is actually discarded
    breaks = sorted(breaks)[1:]
    if contig or opts.nobreaks:
        breaks = []
    plot_heatmap(ax, B, breaks, iopts, groups=new_groups, binsize=resolution)

    # Title
    pf = npyfile.rsplit(".", 1)[0]
    title = opts.title
    if contig:
        title += "-{}".format(contig)
    root.text(
        0.5,
        0.98,
        markup(title),
        color="darkslategray",
        size=18,
        ha="center",
        va="center",
    )

    normalize_axes(root)
    image_name = pf + "." + iopts.format
    # macOS sometimes has way too verbose output
    logging.getLogger().setLevel(logging.CRITICAL)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def mergemat(args):
    """
    %prog mergemat *.npy

    Combine counts from multiple .npy data files.
    """
    p = OptionParser(mergemat.__doc__)
    p.set_outfile(outfile="out")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    npyfiles = args
    A = np.load(npyfiles[0])
    logging.debug(
        "Load `{}`: matrix of shape {}:; sum={}".format(npyfiles[0], A.shape, A.sum())
    )
    for npyfile in npyfiles[1:]:
        B = np.load(npyfile)
        A += B
        logging.debug("Load `{}`: sum={}".format(npyfiles[0], A.sum()))

    pf = opts.outfile
    np.save(pf, A)
    logging.debug("Combined {} files into `{}.npy`".format(len(npyfiles), pf))


def get_seqstarts(bamfile, N, seqids=None):
    """Go through the SQ headers and pull out all sequences with size
    greater than the resolution settings, i.e. contains at least a few cells
    """
    import pysam

    bamfile = pysam.AlignmentFile(bamfile, "rb")
    seqsize = {}
    for kv in bamfile.header["SQ"]:
        if kv["LN"] < 10 * N:
            continue
        seqsize[kv["SN"]] = kv["LN"] // N + 1

    allseqs = seqids or natsorted(seqsize.keys())
    allseqsizes = np.array([seqsize[x] for x in allseqs])
    seqstarts = np.cumsum(allseqsizes)
    seqstarts = np.roll(seqstarts, 1)
    total_bins = seqstarts[0]
    seqstarts[0] = 0
    seqstarts = dict(zip(allseqs, seqstarts))
    seqid_sizes = dict((x, seqsize[x]) for x in allseqs)

    return seqstarts, seqid_sizes, total_bins


def get_distbins(start=100, bins=2000, ratio=1.01):
    """Get exponentially sized bins for link length"""
    b = np.ones(bins, dtype="float64")
    b[0] = 100
    for i in range(1, bins):
        b[i] = b[i - 1] * ratio
    bins = np.around(b).astype(dtype="int")
    binsizes = np.diff(bins)
    return bins, binsizes


def bam2mat(args):
    """
    %prog bam2mat input.bam

    Convert bam file to .mat format, which is simply numpy 2D array. Important
    parameter is the resolution, which is the cell size. Small cell size lead
    to more fine-grained heatmap, but leads to large .mat size and slower
    plotting.
    """
    import pysam
    from jcvi.utils.cbook import percentage

    p = OptionParser(bam2mat.__doc__)
    p.add_option(
        "--resolution",
        default=500000,
        type="int",
        help="Resolution when counting the links",
    )
    p.add_option(
        "--seqids",
        default=None,
        help="Use a given seqids file, a single line with seqids joined by comma",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bamfilename,) = args
    pf = bamfilename.rsplit(".", 1)[0]
    N = opts.resolution
    pf += f".resolution_{N}"
    bins = 1500  # Distance distribution bins
    minsize = 100  # Record distance if it is at least minsize
    seqids = (
        open(opts.seqids).readline().strip().split(",")
        if op.exists(opts.seqids)
        else None
    )

    seqstarts, seqsize, total_bins = get_seqstarts(bamfilename, N, seqids=seqids)
    distbinstarts, distbinsizes = get_distbins(start=minsize, bins=bins)

    # Store the starts and sizes into a JSON file
    jsonfile = pf + ".json"
    fwjson = open(jsonfile, "w")
    header = {
        "starts": seqstarts,
        "sizes": seqsize,
        "total_bins": total_bins,
        "distbinstarts": list(distbinstarts),
        "distbinsizes": list(distbinsizes),
        "resolution": N,
    }

    # int64 will not be able to deserialize with Python 3
    # Here is a workaround:
    # https://stackoverflow.com/questions/11942364/typeerror-integer-is-not-json-serializable-when-serializing-json-in-python
    def default(o):
        if isinstance(o, np.int64):
            return int(o)
        raise TypeError

    json.dump(header, fwjson, sort_keys=True, indent=4, default=default)
    fwjson.close()
    logging.debug("Contig bin starts written to `{}`".format(jsonfile))

    print(sorted(seqstarts.items(), key=lambda x: x[-1]))
    logging.debug("Initialize matrix of size {}x{}".format(total_bins, total_bins))
    A = np.zeros((total_bins, total_bins), dtype="int")
    B = np.zeros(bins, dtype="int")

    # Find the bin ID of each read
    def bin_number(chr, pos):
        return seqstarts[chr] + pos // N

    def distbin_number(dist, start=minsize, ratio=1.01):
        return int(round(math.log(dist * 1.0 / start, ratio)))

    bamfile = pysam.AlignmentFile(bamfilename, "rb")
    # Check all reads, rules borrowed from LACHESIS
    # https://github.com/shendurelab/LACHESIS/blob/master/src/GenomeLinkMatrix.cc#L1476
    j = k = 0
    for c in bamfile:
        j += 1
        if j % 100000 == 0:
            print("{} reads counted".format(j), file=sys.stderr)

        if c.is_qcfail and c.is_duplicate:
            continue
        if c.is_secondary and c.is_supplementary:
            continue
        if c.mapping_quality == 0:
            continue
        if not c.is_paired:
            continue
        if c.is_read2:  # Take only one read
            continue

        # pysam v0.8.3 does not support keyword reference_name
        achr = bamfile.getrname(c.reference_id)
        apos = c.reference_start
        bchr = bamfile.getrname(c.next_reference_id)
        bpos = c.next_reference_start
        if achr not in seqstarts or bchr not in seqstarts:
            continue
        if achr == bchr:
            dist = abs(apos - bpos)
            if dist < minsize:
                continue
            db = distbin_number(dist)
            B[db] += 1

        abin, bbin = bin_number(achr, apos), bin_number(bchr, bpos)
        A[abin, bbin] += 1
        if abin != bbin:
            A[bbin, abin] += 1

        k += 1

    logging.debug("Total reads counted: {}".format(percentage(2 * k, j)))
    bamfile.close()
    np.save(pf, A)
    logging.debug("Link counts written to `{}.npy`".format(pf))
    np.save(pf + ".dist", B)
    logging.debug("Link dists written to `{}.dist.npy`".format(pf))


def simulate(args):
    """
    %prog simulate test

    Simulate CLM and IDS files with given names.

    The simulator assumes several distributions:
    - Links are distributed uniformly across genome
    - Log10(link_size) are distributed normally
    - Genes are distributed uniformly
    """
    p = OptionParser(simulate.__doc__)
    p.add_option("--genomesize", default=10000000, type="int", help="Genome size")
    p.add_option("--genes", default=1000, type="int", help="Number of genes")
    p.add_option("--contigs", default=100, type="int", help="Number of contigs")
    p.add_option("--coverage", default=10, type="int", help="Link coverage")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (pf,) = args
    GenomeSize = opts.genomesize
    Genes = opts.genes
    Contigs = opts.contigs
    Coverage = opts.coverage
    PE = 500
    Links = int(GenomeSize * Coverage / PE)

    # Simulate the contig sizes that sum to GenomeSize
    # See also:
    # <https://en.wikipedia.org/wiki/User:Skinnerd/Simplex_Point_Picking>
    (ContigSizes,) = np.random.dirichlet([1] * Contigs, 1) * GenomeSize
    ContigSizes = np.array(np.round_(ContigSizes, decimals=0), dtype=int)
    ContigStarts = np.zeros(Contigs, dtype=int)
    ContigStarts[1:] = np.cumsum(ContigSizes)[:-1]

    # Write IDS file
    idsfile = pf + ".ids"
    fw = open(idsfile, "w")
    print("#Contig\tRECounts\tLength", file=fw)
    for i, s in enumerate(ContigSizes):
        print("tig{:04d}\t{}\t{}".format(i, s // (4 ** 4), s), file=fw)
    fw.close()

    # Simulate the gene positions
    GenePositions = np.sort(np.random.randint(0, GenomeSize, size=Genes))
    write_last_and_beds(pf, GenePositions, ContigStarts)

    # Simulate links, uniform start, with link distances following 1/x, where x
    # is the distance between the links. As an approximation, we have links
    # between [1e3, 1e7], so we map from uniform [1e-7, 1e-3]
    LinkStarts = np.sort(np.random.randint(1, GenomeSize, size=Links))
    a, b = 1e-7, 1e-3
    LinkSizes = np.array(
        np.round_(1 / ((b - a) * np.random.rand(Links) + a), decimals=0), dtype="int"
    )
    LinkEnds = LinkStarts + LinkSizes

    # Find link to contig membership
    LinkStartContigs = np.searchsorted(ContigStarts, LinkStarts) - 1
    LinkEndContigs = np.searchsorted(ContigStarts, LinkEnds) - 1

    # Extract inter-contig links
    InterContigLinks = (LinkStartContigs != LinkEndContigs) & (
        LinkEndContigs != Contigs
    )
    ICLinkStartContigs = LinkStartContigs[InterContigLinks]
    ICLinkEndContigs = LinkEndContigs[InterContigLinks]
    ICLinkStarts = LinkStarts[InterContigLinks]
    ICLinkEnds = LinkEnds[InterContigLinks]

    # Write CLM file
    write_clm(
        pf,
        ICLinkStartContigs,
        ICLinkEndContigs,
        ICLinkStarts,
        ICLinkEnds,
        ContigStarts,
        ContigSizes,
    )


def write_last_and_beds(pf, GenePositions, ContigStarts):
    """
    Write LAST file, query and subject BED files.
    """
    qbedfile = pf + "tigs.bed"
    sbedfile = pf + "chr.bed"
    lastfile = "{}tigs.{}chr.last".format(pf, pf)
    qbedfw = open(qbedfile, "w")
    sbedfw = open(sbedfile, "w")
    lastfw = open(lastfile, "w")

    GeneContigs = np.searchsorted(ContigStarts, GenePositions) - 1
    for i, (c, gstart) in enumerate(zip(GeneContigs, GenePositions)):
        gene = "gene{:05d}".format(i)
        tig = "tig{:04d}".format(c)
        start = ContigStarts[c]
        cstart = gstart - start
        print("\t".join(str(x) for x in (tig, cstart, cstart + 1, gene)), file=qbedfw)
        print(
            "\t".join(str(x) for x in ("chr1", gstart, gstart + 1, gene)), file=sbedfw
        )
        lastatoms = [gene, gene, 100] + [0] * 8 + [100]
        print("\t".join(str(x) for x in lastatoms), file=lastfw)

    qbedfw.close()
    sbedfw.close()
    lastfw.close()


def write_clm(
    pf,
    ICLinkStartContigs,
    ICLinkEndContigs,
    ICLinkStarts,
    ICLinkEnds,
    ContigStarts,
    ContigSizes,
):
    """
    Write CLM file from simulated data.
    """
    clm = defaultdict(list)
    for start, end, linkstart, linkend in zip(
        ICLinkStartContigs, ICLinkEndContigs, ICLinkStarts, ICLinkEnds
    ):
        start_a = ContigStarts[start]
        start_b = start_a + ContigSizes[start]
        end_a = ContigStarts[end]
        end_b = end_a + ContigSizes[end]
        if linkend >= end_b:
            continue
        clm[(start, end)].append(
            (linkstart - start_a, start_b - linkstart, linkend - end_a, end_b - linkend)
        )

    clmfile = pf + ".clm"
    fw = open(clmfile, "w")

    def format_array(a):
        return [str(x) for x in sorted(a) if x > 0]

    for (start, end), links in sorted(clm.items()):
        start = "tig{:04d}".format(start)
        end = "tig{:04d}".format(end)
        nlinks = len(links)
        if not nlinks:
            continue
        ff = format_array([(b + c) for a, b, c, d in links])
        fr = format_array([(b + d) for a, b, c, d in links])
        rf = format_array([(a + c) for a, b, c, d in links])
        rr = format_array([(a + d) for a, b, c, d in links])
        print("{}+ {}+\t{}\t{}".format(start, end, nlinks, " ".join(ff)), file=fw)
        print("{}+ {}-\t{}\t{}".format(start, end, nlinks, " ".join(fr)), file=fw)
        print("{}- {}+\t{}\t{}".format(start, end, nlinks, " ".join(rf)), file=fw)
        print("{}- {}-\t{}\t{}".format(start, end, nlinks, " ".join(rr)), file=fw)
    fw.close()


def density(args):
    """
    %prog density test.clm

    Estimate link density of contigs.
    """
    p = OptionParser(density.__doc__)
    p.add_option(
        "--save",
        default=False,
        action="store_true",
        help="Write log densitites of contigs to file",
    )
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (clmfile,) = args
    clm = CLMFile(clmfile)
    pf = clmfile.rsplit(".", 1)[0]

    if opts.save:
        logdensities = clm.calculate_densities()
        densityfile = pf + ".density"
        fw = open(densityfile, "w")
        for name, logd in logdensities.items():
            s = clm.tig_to_size[name]
            print("\t".join(str(x) for x in (name, s, logd)), file=fw)
        fw.close()
        logging.debug("Density written to `{}`".format(densityfile))

    tourfile = clmfile.rsplit(".", 1)[0] + ".tour"
    tour = clm.activate(tourfile=tourfile, backuptour=False)
    clm.flip_all(tour)
    clm.flip_whole(tour)
    clm.flip_one(tour)


def optimize(args):
    """
    %prog optimize test.clm

    Optimize the contig order and orientation, based on CLM file.
    """
    p = OptionParser(optimize.__doc__)
    p.add_option(
        "--skiprecover",
        default=False,
        action="store_true",
        help="Do not import 'recover' contigs",
    )
    p.add_option(
        "--startover",
        default=False,
        action="store_true",
        help="Do not resume from existing tour file",
    )
    p.add_option("--skipGA", default=False, action="store_true", help="Skip GA step")
    p.set_outfile(outfile=None)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (clmfile,) = args
    startover = opts.startover
    runGA = not opts.skipGA
    cpus = opts.cpus

    # Load contact map
    clm = CLMFile(clmfile, skiprecover=opts.skiprecover)

    tourfile = opts.outfile or clmfile.rsplit(".", 1)[0] + ".tour"
    if startover:
        tourfile = None
    tour = clm.activate(tourfile=tourfile)

    fwtour = open(tourfile, "w")
    # Store INIT tour
    print_tour(fwtour, clm.tour, "INIT", clm.active_contigs, clm.oo, signs=clm.signs)

    if runGA:
        for phase in range(1, 3):
            tour = optimize_ordering(fwtour, clm, phase, cpus)
            tour = clm.prune_tour(tour, cpus)

    # Flip orientations
    phase = 1
    while True:
        tag1, tag2 = optimize_orientations(fwtour, clm, phase, cpus)
        if tag1 == REJECT and tag2 == REJECT:
            logging.debug("Terminating ... no more {}".format(ACCEPT))
            break
        phase += 1

    fwtour.close()


def optimize_ordering(fwtour, clm, phase, cpus):
    """
    Optimize the ordering of contigs by Genetic Algorithm (GA).
    """
    from .chic import score_evaluate_M

    # Prepare input files
    tour_contigs = clm.active_contigs
    tour_sizes = clm.active_sizes
    tour_M = clm.M
    tour = clm.tour
    signs = clm.signs
    oo = clm.oo

    def callback(tour, gen, phase, oo):
        fitness = tour.fitness if hasattr(tour, "fitness") else None
        label = "GA{}-{}".format(phase, gen)
        if fitness:
            fitness = "{0}".format(fitness).split(",")[0].replace("(", "")
            label += "-" + fitness
        if gen % 20 == 0:
            print_tour(fwtour, tour, label, tour_contigs, oo, signs=signs)
        return tour

    callbacki = partial(callback, phase=phase, oo=oo)
    toolbox = GA_setup(tour)
    toolbox.register("evaluate", score_evaluate_M, tour_sizes=tour_sizes, tour_M=tour_M)
    tour, tour_fitness = GA_run(
        toolbox, ngen=1000, npop=100, cpus=cpus, callback=callbacki
    )
    clm.tour = tour

    return tour


def optimize_orientations(fwtour, clm, phase, cpus):
    """
    Optimize the orientations of contigs by using heuristic flipping.
    """
    # Prepare input files
    tour_contigs = clm.active_contigs
    tour = clm.tour
    oo = clm.oo

    print_tour(
        fwtour, tour, "FLIPALL{}".format(phase), tour_contigs, oo, signs=clm.signs
    )
    tag1 = clm.flip_whole(tour)
    print_tour(
        fwtour, tour, "FLIPWHOLE{}".format(phase), tour_contigs, oo, signs=clm.signs
    )
    tag2 = clm.flip_one(tour)
    print_tour(
        fwtour, tour, "FLIPONE{}".format(phase), tour_contigs, oo, signs=clm.signs
    )

    return tag1, tag2


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
    mkdir(odir, overwrite=True)
    os.chdir(odir)
    logging.debug("Change into subdir `{}`".format(odir))

    # Make anchorsfile
    anchorsfile = ".".join(op.basename(lastfile).split(".", 2)[:2]) + ".anchors"
    fw = open(anchorsfile, "w")
    for b in Blast(lastfile):
        print(
            "\t".join((gene_name(b.query), gene_name(b.subject), str(int(b.score)))),
            file=fw,
        )
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
        if contig[-1] in ("+", "-", "?"):
            tour.append(contig[:-1])
            tour_o.append(contig[-1])
        else:  # Unoriented
            tour.append(contig)
            tour_o.append("?")
    return tour, tour_o


def iter_last_tour(tourfile, clm):
    """
    Extract last tour from tourfile. The clm instance is also passed in to see
    if any contig is covered in the clm.
    """
    row = open(tourfile).readlines()[-1]
    _tour, _tour_o = separate_tour_and_o(row)
    tour = []
    tour_o = []
    for tc, to in zip(_tour, _tour_o):
        if tc not in clm.contigs:
            logging.debug(
                "Contig `{}` in file `{}` not found in `{}`".format(
                    tc, tourfile, clm.idsfile
                )
            )
            continue
        tour.append(tc)
        tour_o.append(to)
    return tour, tour_o


def iter_tours(tourfile, frames=1):
    """
    Extract tours from tourfile. Tourfile contains a set of contig
    configurations, generated at each iteration of the genetic algorithm. Each
    configuration has two rows, first row contains iteration id and score,
    second row contains list of contigs, separated by comma.
    """
    fp = open(tourfile)

    i = 0
    for row in fp:
        if row[0] == ">":
            label = row[1:].strip()
            if label.startswith("GA"):
                pf, j, score = label.split("-", 2)
                j = int(j)
            else:
                j = 0
            i += 1
        else:
            if j % frames != 0:
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
    p.add_option("--frames", default=500, type="int", help="Only plot every N frames")
    p.add_option(
        "--engine",
        default="ffmpeg",
        choices=("ffmpeg", "gifsicle"),
        help="Movie engine, output MP4 or GIF",
    )
    p.set_beds()
    opts, args, iopts = p.set_image_options(
        args, figsize="16x8", style="white", cmap="coolwarm", format="png", dpi=300
    )

    if len(args) != 3:
        sys.exit(not p.print_help())

    tourfile, clmfile, lastfile = args
    tourfile = op.abspath(tourfile)
    clmfile = op.abspath(clmfile)
    lastfile = op.abspath(lastfile)
    cwd = os.getcwd()
    odir = op.basename(tourfile).rsplit(".", 1)[0] + "-movie"
    anchorsfile, qbedfile, contig_to_beds = prepare_synteny(
        tourfile, lastfile, odir, p, opts
    )

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
            if contig not in contig_to_beds:
                continue
            bedlines = contig_to_beds[contig][:]
            if o == "-":
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
        args.append(
            [[tour, clmfile, ianchorsfile, "--outfile", image_name, "--label", label]]
        )

    Jobs(movieframe, args).run()

    os.chdir(cwd)
    make_movie(odir, odir, engine=opts.engine, format=iopts.format)


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
        if row[0] == "#":
            continue
        x, y, z = row.split()
        if x == "X":
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
        print(pf)
        print(oo)

        tour, tour_sizes, tour_M = prepare_ec(oo, sizes, M)
        # Store INIT tour
        print_tour(fwtour, tour, "INIT", contig_names, oo)

        # Faster Cython version for evaluation
        from .chic import score_evaluate_M

        callbacki = partial(callback, oo=oo)
        toolbox = GA_setup(tour)
        toolbox.register(
            "evaluate", score_evaluate_M, tour_sizes=tour_sizes, tour_M=tour_M
        )
        tour, tour.fitness = GA_run(
            toolbox, npop=100, cpus=opts.cpus, callback=callbacki
        )
        print(tour, tour.fitness)
        break

    fwtour.close()


def print_tour(fwtour, tour, label, contig_names, oo, signs=None):
    print(">" + label, file=fwtour)
    if signs is not None:
        contig_o = []
        for x in tour:
            idx = oo[x]
            sign = {1: "+", 0: "?", -1: "-"}[signs[idx]]
            contig_o.append(contig_names[idx] + sign)
        print(" ".join(contig_o), file=fwtour)
    else:
        print(" ".join(contig_names[oo[x]] for x in tour), file=fwtour)


def prepare_ec(oo, sizes, M):
    """
    This prepares EC and converts from contig_id to an index.
    """
    tour = range(len(oo))
    tour_sizes = np.array([sizes.sizes[x] for x in oo])
    tour_M = M[oo, :][:, oo]
    return tour, tour_sizes, tour_M


def score_evaluate(tour, tour_sizes=None, tour_M=None):
    """SLOW python version of the evaluation function. For benchmarking
    purposes only. Do not use in production.
    """
    sizes_oo = np.array([tour_sizes[x] for x in tour])
    sizes_cum = np.cumsum(sizes_oo) - sizes_oo / 2
    s = 0
    size = len(tour)
    for ia in range(size):
        a = tour[ia]
        for ib in range(ia + 1, size):
            b = tour[ib]
            links = tour_M[a, b]
            dist = sizes_cum[ib] - sizes_cum[ia]
            if dist > 1e7:
                break
            s += links * 1.0 / dist
    return (s,)


def movieframe(args):
    """
    %prog movieframe tour test.clm contigs.ref.anchors

    Draw heatmap and synteny in the same plot.
    """
    p = OptionParser(movieframe.__doc__)
    p.add_option("--label", help="Figure title")
    p.set_beds()
    p.set_outfile(outfile=None)
    opts, args, iopts = p.set_image_options(
        args, figsize="16x8", style="white", cmap="coolwarm", format="png", dpi=120
    )

    if len(args) != 3:
        sys.exit(not p.print_help())

    tour, clmfile, anchorsfile = args
    tour = tour.split(",")
    image_name = opts.outfile or ("movieframe." + iopts.format)
    label = opts.label or op.basename(image_name).rsplit(".", 1)[0]

    clm = CLMFile(clmfile)
    totalbins, bins, breaks = make_bins(tour, clm.tig_to_size)
    M = read_clm(clm, totalbins, bins)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # whole canvas
    ax1 = fig.add_axes([0.05, 0.1, 0.4, 0.8])  # heatmap
    ax2 = fig.add_axes([0.55, 0.1, 0.4, 0.8])  # dot plot
    ax2_root = fig.add_axes([0.5, 0, 0.5, 1])  # dot plot canvas

    # Left axis: heatmap
    plot_heatmap(ax1, M, breaks, iopts)

    # Right axis: synteny
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts, sorted=False)
    dotplot(anchorsfile, qbed, sbed, fig, ax2_root, ax2, sep=False, title="")

    root.text(0.5, 0.98, clm.name, color="g", ha="center", va="center")
    root.text(0.5, 0.95, label, color="darkslategray", ha="center", va="center")
    normalize_axes(root)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def make_bins(tour, sizes):
    breaks = []
    start = 0
    bins = {}
    for x in tour:
        size = sizes[x]
        end = start + int(round(size * 1.0 / BINSIZE))
        bins[x] = (start, end)
        start = end
    breaks.append(start)

    totalbins = start
    return totalbins, bins, breaks


def read_clm(clm, totalbins, bins):
    M = np.zeros((totalbins, totalbins))
    for (x, y), z in clm.contacts.items():
        if x not in bins or y not in bins:
            continue
        xstart, xend = bins[x]
        ystart, yend = bins[y]
        M[xstart:xend, ystart:yend] = z
        M[ystart:yend, xstart:xend] = z

    M = np.log10(M + 1)
    return M


def plot_heatmap(ax, M, breaks, iopts, groups=[], plot_breaks=False, binsize=BINSIZE):
    """Plot heatmap illustrating the contact probabilities in Hi-C data.

    Args:
        ax (pyplot.axes): Matplotlib axis
        M (np.array): 2D numpy-array
        breaks (List[int]): Positions of chromosome starts. Can be None.
        iopts (OptionParser options): Graphical options passed in from commandline
        groups (List, optional): [(start, end, [(position, seqid)], color)]. Defaults to [].
        plot_breaks (bool): Whether to plot white breaks. Defaults to False.
        binsize (int, optional): Resolution of the heatmap. Defaults to BINSIZE.
    """
    import seaborn as sns

    cmap = sns.cubehelix_palette(rot=0.5, as_cmap=True)
    ax.imshow(M, cmap=cmap, interpolation="none")
    _, xmax = ax.get_xlim()
    xlim = (0, xmax)
    if plot_breaks:
        for b in breaks[:-1]:
            ax.plot([b, b], xlim, "w-")
            ax.plot(xlim, [b, b], "w-")

    def simplify_seqid(seqid):
        seqid = seqid.replace("_", "")
        if seqid[:3].lower() == "chr":
            seqid = seqid[3:]
        return seqid.lstrip("0")

    for start, end, position_seqids, color in groups:
        # Plot a square
        ax.plot([start, start], [start, end], "-", color=color)
        ax.plot([start, end], [start, start], "-", color=color)
        ax.plot([start, end], [end, end], "-", color=color)
        ax.plot([end, end], [start, end], "-", color=color)
        for position, seqid in position_seqids:
            seqid = simplify_seqid(seqid)
            ax.text(position, end, seqid, ha="center", va="top")

    ax.set_xlim(xlim)
    ax.set_ylim((xlim[1], xlim[0]))  # Flip the y-axis so the origin is at the top
    ax.set_xticklabels(ax.get_xticks(), family="Helvetica", color="gray")
    ax.set_yticklabels(ax.get_yticks(), family="Helvetica", color="gray", rotation=90)
    ax.tick_params(left=True, bottom=True, labelleft=True, labelbottom=True)
    formatter = ticker.FuncFormatter(
        lambda x, pos: human_readable(int(x) * binsize, pos, base=True)
    )
    ax.xaxis.set_major_formatter(formatter)
    ax.yaxis.set_major_formatter(formatter)
    binlabel = "Resolution = {} per bin".format(human_size(binsize, precision=0))
    ax.set_xlabel(binlabel)


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
    fwagp = must_open(opts.outfile, "w")
    orderingfiles = natsorted(iglob(odir, "*.ordering"))
    sizes = Sizes(contigsfasta).mapping
    contigs = set(sizes.keys())
    anchored = set()

    for ofile in orderingfiles:
        co = ContigOrdering(ofile)
        anchored |= set([x.contig_name for x in co])
        obj = op.basename(ofile).split(".")[0]
        co.write_agp(obj, sizes, fwagp)

    singletons = contigs - anchored
    logging.debug("Anchored: {}, Singletons: {}".format(len(anchored), len(singletons)))

    for s in natsorted(singletons):
        order_to_agp(s, [(s, "?")], sizes, fwagp)


if __name__ == "__main__":
    main()
