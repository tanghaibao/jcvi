#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# sugarcane.py
# projects
#
# Created by Haibao Tang on 12/02/19
# Copyright © 2019 Haibao Tang. All rights reserved.
#
import os.path as op
import sys

from collections import Counter, defaultdict
from itertools import combinations, groupby, product
from random import random, sample
from typing import Dict

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

from ..apps.base import ActionDispatcher, OptionParser, logger, mkdir
from ..formats.blast import Blast
from ..graphics.base import adjust_spines, markup, normalize_axes, savefig
from ..utils.validator import validate_in_choices

SoColor = "#7436a4"  # Purple
SsColor = "#5a8340"  # Green

# Computed using prepare(), corrected with real sizes
ChrSizes = {
    "SO-chr01": 148750011,
    "SO-chr02": 119865146,
    "SO-chr03": 103845728,
    "SO-chr04": 104559946,
    "SO-chr05": 93134056,
    "SO-chr06": 74422021,
    "SO-chr07": 81308893,
    "SO-chr08": 71010813,
    "SO-chr09": 86380266,
    "SO-chr10": 73923121,
    "SS-chr01": 114519418,
    "SS-chr02": 119157314,
    "SS-chr03": 85009228,
    "SS-chr04": 79762909,
    "SS-chr05": 90584537,
    "SS-chr06": 95848354,
    "SS-chr07": 83589369,
    "SS-chr08": 64028871,
}


# Simulate genome composition
class Genome:
    def __init__(self, name, prefix, ploidy, haploid_chromosome_count):
        """
        Simulate a genome with given ploidy and haploid_chromosome_count. Example:

        >>> print(Genome("t", "pf", 2, 3))
        test: pf-chr01_a,pf-chr01_b,pf-chr02_a,pf-chr02_b,pf-chr03_a,pf-chr03_b
        """
        self.name = name
        chromosomes = []
        for i in range(haploid_chromosome_count):
            chromosomes += [
                f"{prefix}-chr{i + 1:02d}_{chr(ord('a') + j)}" for j in range(ploidy)
            ]
        self.chromosomes = chromosomes

    def __len__(self):
        return len(self.chromosomes)

    @classmethod
    def make(cls, name, chromosomes):
        genome = Genome(name, "", 0, 0)
        genome.chromosomes = chromosomes
        return genome

    @property
    def gamete(self):
        """Randomly generate a gamete from current genome that"""
        self.chromosomes.sort()
        gamete_chromosomes = []

        # Check for any chromosome that have 2 identical copies, if so, we will assume disomic
        # inheritance for that chromosome and always keep one and only copy
        duplicate_chromosomes = []
        singleton_chromosomes = []
        for chromosome, chromosomes in groupby(self.chromosomes):
            chromosomes = list(chromosomes)
            ncopies = len(chromosomes)
            duplicate_chromosomes += [chromosome] * (ncopies // 2)
            if ncopies % 2 == 1:
                singleton_chromosomes.append(chromosome)

        # Get one copy of each duplicate chromosome first
        gamete_chromosomes += duplicate_chromosomes

        def prefix(x):
            return x.split("_", 1)[0]

        # Randomly assign the rest, singleton chromosomes
        for group, chromosomes in groupby(singleton_chromosomes, key=prefix):
            chromosomes = list(chromosomes)
            halfn = len(chromosomes) // 2
            # Odd number, e.g. 5, equal chance to be 2 or 3
            if len(chromosomes) % 2 != 0 and random() < 0.5:
                halfn += 1
            gamete_chromosomes += sorted(sample(chromosomes, halfn))
        return Genome.make(self.name + " gamete", gamete_chromosomes)

    def mate_nplusn(self, name, other_genome, verbose=True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (n+n)", file=sys.stderr
            )
        f1_chromosomes = sorted(
            self.gamete.chromosomes + other_genome.gamete.chromosomes
        )
        return Genome.make(name, f1_chromosomes)

    def mate_nx2plusn(self, name, other_genome, verbose=True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (2xn+n)",
                file=sys.stderr,
            )
        f1_chromosomes = sorted(
            2 * self.gamete.chromosomes + other_genome.gamete.chromosomes
        )
        return Genome.make(name, f1_chromosomes)

    def mate_2nplusn(self, name, other_genome, verbose=True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (2n+n)",
                file=sys.stderr,
            )
        f1_chromosomes = sorted(self.chromosomes + other_genome.gamete.chromosomes)
        return Genome.make(name, f1_chromosomes)

    def __str__(self):
        return self.name + ": " + ",".join(self.chromosomes)

    @property
    def summary(self):
        def prefix(x, sep="-"):
            return x.split(sep, 1)[0]

        def size(chromosomes):
            return sum(ChrSizes[prefix(x, sep="_")] for x in chromosomes)

        # Chromosome count
        total_count = 0
        total_unique = 0
        total_size = 0
        total_so_size = 0
        ans = []
        for group, chromosomes in groupby(self.chromosomes, prefix):
            chromosomes = list(chromosomes)
            uniq_chromosomes = set(chromosomes)
            group_count = len(chromosomes)
            group_unique = len(uniq_chromosomes)
            group_so_size = size({x for x in uniq_chromosomes if x[:2] == "SO"})
            group_size = size(uniq_chromosomes)
            total_count += group_count
            total_unique += group_unique
            total_so_size += group_so_size
            total_size += group_size
            ans.append((group, group_count, group_unique, group_so_size, group_size))
        ans.append(("Total", total_count, total_unique, total_so_size, total_size))
        return ans

    def print_summary(self):
        print("[SUMMARY]")
        for group, group_count, group_unique in self.summary:
            print(f"{group}: count={group_count}, unique={group_unique}")


class GenomeSummary:
    def __init__(self, SO_data, SS_data, percent_SO_data):
        self.SO_data = SO_data
        self.SS_data = SS_data
        self.percent_SO_data = percent_SO_data
        self.percent_SS_data = [100 - x for x in percent_SO_data]

    def _summary(self, a, tag, precision=0):
        mean, min, max = (
            round(np.mean(a), precision),
            round(np.min(a), precision),
            round(np.max(a), precision),
        )
        s = f"*{tag}* chr: {mean:.0f}"
        if min == mean and max == mean:
            return s
        return s + f" ({min:.0f}-{max:.0f})"

    def _percent_summary(self, a, tag, precision=1):
        mean, min, max = (
            round(np.mean(a), precision),
            round(np.min(a), precision),
            round(np.max(a), precision),
        )
        s = f"*{tag}*%: {mean:.1f}%"
        print(s)
        if min == mean and max == mean:
            return s
        return s + f" ({min:.1f}-{max:.1f}%)"

    @property
    def percent_SO_summary(self):
        return self._percent_summary(self.percent_SO_data, "So")

    @property
    def percent_SS_summary(self):
        return self._percent_summary(self.percent_SS_data, "Ss")

    @property
    def SO_summary(self):
        return self._summary(self.SO_data, "So")

    @property
    def SS_summary(self):
        return self._summary(self.SS_data, "Ss")


def simulate_F1(SO, SS, mode="nx2+n", verbose=False):
    SO_SS_F1 = (
        SO.mate_nx2plusn("SOxSS F1", SS, verbose=verbose)
        if mode == "nx2+n"
        else SO.mate_2nplusn("SOxSS F1", SS, verbose=verbose)
    )
    if verbose:
        SO_SS_F1.print_summary()
    return SO_SS_F1


def simulate_F2(SO, SS, mode="nx2+n", verbose=False):
    SO_SS_F1 = simulate_F1(SO, SS, mode=mode, verbose=verbose)
    SO_SS_F2_nplusn = SO_SS_F1.mate_nplusn("SOxSS F2", SO_SS_F1, verbose=verbose)
    if verbose:
        SO_SS_F2_nplusn.print_summary()
    return SO_SS_F2_nplusn


def simulate_F1intercross(SO, SS, mode="nx2+n", verbose=False):
    SO_SS_F1_1 = simulate_F1(SO, SS, mode=mode, verbose=verbose)
    SO_SS_F1_2 = simulate_F1(SO, SS, mode=mode, verbose=verbose)
    SO_SS_F1intercross_nplusn = SO_SS_F1_1.mate_nplusn(
        "SOxSS F1 intercross", SO_SS_F1_2, verbose=verbose
    )
    return SO_SS_F1intercross_nplusn


def simulate_BCn(n, SO, SS, mode="nx2+n", verbose=False):
    SS_SO_F1 = simulate_F1(SO, SS, mode=mode, verbose=verbose)
    SS_SO_BC1, SS_SO_BC2_nplusn, SS_SO_BC3_nplusn, SS_SO_BC4_nplusn = (
        None,
        None,
        None,
        None,
    )
    # BC1
    if n >= 1:
        SS_SO_BC1 = (
            # Expecting one more round of female restitution in BC1
            SO.mate_nx2plusn("SSxSO BC1", SS_SO_F1, verbose=verbose)
            if mode == "nx2+n"
            else SO.mate_2nplusn("SSxSO BC1", SS_SO_F1, verbose=verbose)
        )
    # BC2
    if n >= 2:
        SS_SO_BC2_nplusn = SO.mate_nplusn("SSxSO BC2", SS_SO_BC1, verbose=verbose)
    # BC3
    if n >= 3:
        SS_SO_BC3_nplusn = SO.mate_nplusn(
            "SSxSO BC3", SS_SO_BC2_nplusn, verbose=verbose
        )
    # BC4
    if n >= 4:
        SS_SO_BC4_nplusn = SO.mate_nplusn(
            "SSxSO BC4", SS_SO_BC3_nplusn, verbose=verbose
        )
    return [
        None,
        SS_SO_BC1,
        SS_SO_BC2_nplusn,
        SS_SO_BC3_nplusn,
        SS_SO_BC4_nplusn,
    ][n]


def plot_summary(ax, samples: list[Genome]) -> GenomeSummary:
    """Plot the distribution of chromosome numbers given simulated samples.

    Args:
        ax (Axes): Matplotlib axes.
        samples (list[Genome]): Summarized genomes.

    Returns:
        GenomeSummary: Summary statistics of simulated genomes.
    """
    SO_data = []
    SS_data = []
    percent_SO_data = []
    for sample in samples:
        summary = sample.summary
        try:
            _, _, group_unique, _, _ = [x for x in summary if x[0] == "SO"][0]
        except:
            group_unique = 0
        SO_data.append(group_unique)
        try:
            _, _, group_unique, _, _ = [x for x in summary if x[0] == "SS"][0]
        except:
            group_unique = 0
        SS_data.append(group_unique)
        total_tag, _, _, total_so_size, total_size = summary[-1]
        assert total_tag == "Total"
        percent_SO = total_so_size * 100.0 / total_size
        percent_SO_data.append(percent_SO)
    # Avoid overlapping bars
    SS_counter, SO_counter = Counter(SS_data), Counter(SO_data)
    overlaps = SS_counter.keys() & SO_counter.keys()
    shift = 0.5  # used to offset bars a bit to avoid cluttering
    if overlaps:
        for overlap in overlaps:
            logger.debug(f"Modify bar offsets at {overlap} due to SS and SO overlaps")
            SS_counter[overlap - shift] = SS_counter[overlap]
            del SS_counter[overlap]
            SO_counter[overlap + shift] = SO_counter[overlap]
            del SO_counter[overlap]

    def modify_range_end(d: dict, value: int):
        if value not in d:
            return
        # Has data at the range end, but no adjacent data points (i.e. isolated bar)
        if value in d and (value - 1 in d or value + 1 in d):
            return
        logger.debug(f"Modify bar offsets at {value} due to end of range ends")
        d[value - shift if value else value + shift] = d[80]
        del d[value]

    modify_range_end(SS_counter, 0)
    modify_range_end(SS_counter, 80)
    modify_range_end(SO_counter, 0)
    modify_range_end(SO_counter, 80)

    x, y = zip(*sorted(SS_counter.items()))
    ax.bar(np.array(x), y, color=SsColor, ec=SsColor)
    x, y = zip(*sorted(SO_counter.items()))
    ax.bar(np.array(x), y, color=SoColor, ec=SoColor)
    ax.set_xlim(80, 0)
    ax.set_ylim(0, len(samples) / 2)
    ax.set_yticks([])
    summary = GenomeSummary(SO_data, SS_data, percent_SO_data)

    # Write the stats summary within the plot
    summary_style = dict(
        size=9,
        ha="center",
        va="center",
        transform=ax.transAxes,
    )
    ax.text(0.75, 0.85, markup(summary.SS_summary), color=SsColor, **summary_style)
    ax.text(
        0.75, 0.65, markup(summary.percent_SS_summary), color=SsColor, **summary_style
    )
    ax.text(0.25, 0.85, markup(summary.SO_summary), color=SoColor, **summary_style)
    ax.text(
        0.25, 0.65, markup(summary.percent_SO_summary), color=SoColor, **summary_style
    )

    return summary


def write_chromosomes(genomes: list[Genome], filename: str):
    """Write simulated chromosomes to file

    Args:
        genomes (list[Genome]): List of simulated genomes.
        filename (str): File path to write to.
    """
    print(f"Write chromosomes to `{filename}`", file=sys.stderr)
    with open(filename, "w") as fw:
        for genome in genomes:
            print(genome, file=fw)


def write_SO_percent(summary: GenomeSummary, filename: str):
    """Write SO % to file

    Args:
        summary (GenomeSummary): List of simulated genomes.
        filename (str): File path to write to.
    """
    print(f"Write SO percent to `{filename}`", file=sys.stderr)
    with open(filename, "w") as fw:
        print("\n".join(str(x) for x in sorted(summary.percent_SO_data)), file=fw)


def simulate(args):
    """
    %prog simulate [2n+n|nx2+n]

    Run simulation on female restitution. There are two modes:
    - 2n+n: merger between a somatic and a germline
    - nx2+n: merger between a doubled germline and a germline

    These two modes would impact the sequence diversity in the progeny
    genome in F1, F2, BCn ... the goal of this simulation, is thus to
    understand the mode and the spread of such diversity in the hybrid
    progenies.
    """
    sns.set_style("darkgrid")

    p = OptionParser(simulate.__doc__)
    p.add_argument(
        "--verbose",
        default=False,
        action="store_true",
        help="Verbose logging during simulation",
    )
    p.add_argument("-N", default=10000, type=int, help="Number of simulated samples")
    opts, args, iopts = p.set_image_options(args, figsize="6x6")
    if len(args) != 1:
        sys.exit(not p.print_help())

    (mode,) = args
    validate_in_choices(mode, ["2n+n", "nx2+n"], "Mode")
    logger.info(f"Transmission: {mode}")

    # Construct a composite figure with 6 tracks
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    rows = 6
    ypad = 0.05
    yinterval = (1 - 2 * ypad) / (rows + 1)
    yy = 1 - ypad
    xpad = 0.2
    xwidth = 0.7

    # Axes are vertically stacked, and share x-axis
    axes = []
    yy_positions = []  # Save yy positions so we can show details to the right laterr
    for idx in range(rows):
        yy_positions.append(yy)
        yy -= yinterval
        ax = fig.add_axes([xpad, yy, xwidth, yinterval * 0.85])
        if idx != rows - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        axes.append(ax)
    ax1, ax2, ax3, ax4, ax5, ax6 = axes

    # Prepare the simulated data
    # Simulate two parents
    SS = Genome("SS", "SS", 10, 8)
    SO = Genome("SO", "SO", 8, 10)

    verbose = opts.verbose
    N = opts.N
    all_F1s = [simulate_F1(SO, SS, mode=mode, verbose=verbose) for _ in range(N)]
    all_F2s = [simulate_F2(SO, SS, mode=mode, verbose=verbose) for _ in range(N)]
    all_BC1s = [simulate_BCn(1, SO, SS, mode=mode, verbose=verbose) for _ in range(N)]
    all_BC2s = [simulate_BCn(2, SO, SS, mode=mode, verbose=verbose) for _ in range(N)]
    all_BC3s = [simulate_BCn(3, SO, SS, mode=mode, verbose=verbose) for _ in range(N)]
    all_BC4s = [simulate_BCn(4, SO, SS, mode=mode, verbose=verbose) for _ in range(N)]

    # Plotting
    all_F1s_summary = plot_summary(ax1, all_F1s)
    all_F2s_summary = plot_summary(ax2, all_F2s)
    plot_summary(ax3, all_BC1s)
    plot_summary(ax4, all_BC2s)
    plot_summary(ax5, all_BC3s)
    plot_summary(ax6, all_BC4s)

    # Show title to the left
    xx = xpad / 2
    for (title, subtitle), yy in zip(
        (
            (r"$\mathrm{F_1}$", None),
            (r"$\mathrm{F_2}$", None),
            (r"$\mathrm{BC_1}$", None),
            (r"$\mathrm{BC_2}$", None),
            (r"$\mathrm{BC_3}$", None),
            (r"$\mathrm{BC_4}$", None),
        ),
        yy_positions,
    ):
        if subtitle:
            yy -= 0.06
        else:
            yy -= 0.07
        root.text(
            xx,
            yy,
            title,
            color="darkslategray",
            ha="center",
            va="center",
            fontweight="semibold",
        )
        if subtitle:
            yy -= 0.02
            root.text(
                xx, yy, subtitle, color="lightslategray", ha="center", va="center"
            )

    axes[-1].set_xlabel("Number of unique chromosomes")
    adjust_spines(axes[-1], ["bottom"], outward=True)
    normalize_axes(root)

    # Title
    mode_title = r"$n_1\times2 + n_2$" if mode == "nx2+n" else r"$2n + n$"
    root.text(0.5, 0.95, f"Transmission: {mode_title}", ha="center")

    savefig(f"{mode}.pdf", dpi=120)

    outdir = f"simulations_{mode}"
    mkdir(outdir)
    # Write chromosomes to disk
    for genomes, filename in (
        (all_F1s, "all_F1s"),
        (all_F2s, "all_F2s"),
        (all_BC1s, "all_BC1s"),
        (all_BC2s, "all_BC2s"),
        (all_BC3s, "all_BC3s"),
        (all_BC4s, "all_BC4s"),
    ):
        write_chromosomes(genomes, op.join(outdir, filename))

    # Write the SO percent in simulated samples so that we can compute P-value
    for summary, SO_percent_filename in (
        (all_F1s_summary, "all_F1s_SO_percent"),
        (all_F2s_summary, "all_F2s_SO_percent"),
    ):
        write_SO_percent(summary, op.join(outdir, SO_percent_filename))


def _get_sizes(filename, prefix_length, tag, target_size=None):
    """Returns a dictionary of chromome lengths from a given file.

    Args:
        filename ([str]): Path to the input file. Input file is 2-column file
        with rows `seqid length`.
        prefix_length (int): Extract first N characters.
        tag (str): Prepend `tag-` to the seqid.
        target_size (int): Expected genome size. Defaults to None.
    """
    from collections import defaultdict

    sizes_list = defaultdict(list)
    with open(filename) as fp:
        for row in fp:
            if not row.startswith("Chr"):
                continue
            name, size = row.split()
            idx = int(name[3:prefix_length])
            size = int(size)
            name = f"{tag}-chr{idx:02d}"
            sizes_list[name].append(size)

    # Get the average length
    sizes = dict(
        (name, int(round(np.mean(size_list)))) for name, size_list in sizes_list.items()
    )
    print(sizes)
    if target_size is None:
        return sizes

    total_size = sum(sizes.values())
    correction_factor = target_size / total_size
    print(
        f"{tag} total:{total_size} target:{target_size} correction:{correction_factor:.2f}x"
    )
    return dict(
        (name, int(round(correction_factor * size))) for name, size in sizes.items()
    )


def prepare(args):
    """
    %prog SoChrLen.txt SsChrLen.txt

    Calculate lengths from real sugarcane data.
    """
    p = OptionParser(prepare.__doc__)
    _, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    solist, sslist = args
    # The haploid set of LA Purple is 957.2 Mb and haploid set of US56-14-4 is 732.5 Mb
    sizes = _get_sizes(solist, 5, "SO", target_size=int(957.2 * 1e6))
    sizes.update(_get_sizes(sslist, 4, "SS", target_size=int(732.5 * 1e6)))
    print(sizes)


def get_genome_wide_pct(summary: str) -> Dict[tuple, list]:
    """Collect genome-wide ungapped percent identity.
    Specifically, from file `SS_SR_SO.summary.txt`.

    Args:
        summary (str): File that contains per chromosome pct identity info,
        collected via `formats.blast.summary()`.

    Returns:
        Dict[tuple, list]: Genome pair to list of pct identities.
    """
    COLUMNS = "filename, identicals, qry_gapless, qry_gapless_pct, ref_gapless, ref_gapless_pct, qryspan, pct_qryspan, refspan, pct_refspan".split(
        ", "
    )
    df = pd.read_csv(summary, sep="\t", names=COLUMNS)
    data_by_genomes = defaultdict(list)
    for _, row in df.iterrows():
        filename = row["filename"]
        # e.g. SO_Chr01A.SO_Chr01B.1-1.blast
        chr1, chr2 = filename.split(".")[:2]
        genome1, chr1 = chr1.split("_")
        genome2, chr2 = chr2.split("_")
        chr1, chr2 = chr1[:5], chr2[:5]
        if (  # Special casing for SS certain chromosomes that are non-collinear with SO/SR
            genome1 != "SS"
            and genome2 == "SS"
            and chr2 not in ("Chr01", "Chr03", "Chr04")
        ):
            continue
        qry_gapless_pct, ref_gapless_pct = (
            row["qry_gapless_pct"],
            row["ref_gapless_pct"],
        )
        data_by_genomes[(genome1, genome2)] += [qry_gapless_pct, ref_gapless_pct]
    return data_by_genomes


def get_anchors_pct(filename: str, min_pct: int = 94) -> list:
    """Collect CDS-wide ungapped percent identity.

    Args:
        filename (str): Input file name, which is a LAST file.

    Returns:
        list: List of pct identities from this LAST file.
    """
    blast = Blast(filename)
    pct = []
    for c in blast:
        if c.pctid < min_pct:
            continue
        identicals = c.hitlen - c.nmismatch - c.ngaps
        qstart, qstop = c.qstart, c.qstop
        sstart, sstop = c.sstart, c.sstop
        qrycovered = qstop - qstart + 1
        refcovered = sstop - sstart + 1
        pct.append(identicals * 100 / qrycovered)
        pct.append(identicals * 100 / refcovered)
    return pct


def divergence(args):
    """
    %prog divergence SS_SR_SO.summary.txt

    Plot divergence between and within SS/SR/SO genomes.
    """
    sns.set_style("white")

    p = OptionParser(divergence.__doc__)
    p.add_argument("--title", default="Gapless", help="Plot title")
    p.add_argument(
        "--xmin",
        default=94,
        type=int,
        help="Minimum percent identity in the histogram",
    )
    opts, args, iopts = p.set_image_options(args, figsize="8x8")
    if len(args) != 1:
        sys.exit(not p.print_help())

    (summary,) = args
    data_by_genomes = get_genome_wide_pct(summary)
    # Print summary statistics
    print("Genome-wide ungapped percent identity:")
    for (genome1, genome2), pct in sorted(data_by_genomes.items()):
        print(genome1, genome2, np.mean(pct), np.std(pct))

    # Plotting genome-wide divergence
    fig = plt.figure(figsize=(iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    SPECIES_CONFIG = {
        "SS": {"label": "S. spontaneum", "pos": (0.5, 0.67)},
        "SR": {"label": "S. robustum", "pos": (0.2, 0.3)},
        "SO": {"label": "S. officinarum", "pos": (0.8, 0.3)},
    }
    # Get median for each genome pair
    medians = {}
    for g1, g2 in product(SPECIES_CONFIG.keys(), repeat=2):
        g1, g2 = sorted((g1, g2))
        d = data_by_genomes[(g1, g2)]
        medians[(g1, g2)] = np.median(d)
    for g, config in SPECIES_CONFIG.items():
        x, y = config["pos"]
        text = f'*{config["label"]}*' + f"\n{medians[(g, g)]:.1f} %"
        text = markup(text)
        root.text(x, y, text, color="darkslategray", ha="center", va="center")

    # Connect lines
    PAD, YPAD = 0.09, 0.045
    for g1, g2 in combinations(SPECIES_CONFIG.keys(), 2):
        g1, g2 = sorted((g1, g2))
        x1, y1 = SPECIES_CONFIG[g1]["pos"]
        x2, y2 = SPECIES_CONFIG[g2]["pos"]
        x1, x2 = (x1 + PAD, x2 - PAD) if x1 < x2 else (x1 - PAD, x2 + PAD)
        if y1 != y2:
            y1, y2 = (y1 + YPAD, y2 - YPAD) if y1 < y2 else (y1 - YPAD, y2 + YPAD)
        xmid, ymid = (x1 + x2) / 2, (y1 + y2) / 2
        text = f"{medians[(g1, g2)]:.1f} %"
        text = markup(text)
        root.text(xmid, ymid, text, ha="center", va="center", backgroundcolor="w")
        root.plot([x1, x2], [y1, y2], "-", lw=4, color="darkslategray")

    # Pct identity histograms
    PCT_CONFIG = {
        ("SS", "SS"): {"pos": (0.5, 0.82)},
        ("SR", "SR"): {"pos": (0.1, 0.2)},
        ("SO", "SO"): {"pos": (0.9, 0.2)},
        ("SR", "SS"): {"pos": (0.3 - PAD, 0.55)},
        ("SO", "SS"): {"pos": (0.7 + PAD, 0.55)},
        ("SO", "SR"): {"pos": (0.5, 0.18)},
    }
    HIST_WIDTH = 0.15
    xmin = opts.xmin
    for genome_pair, config in PCT_CONFIG.items():
        x, y = config["pos"]
        ax = fig.add_axes(
            [x - HIST_WIDTH / 2, y - HIST_WIDTH / 2, HIST_WIDTH, HIST_WIDTH]
        )
        d = data_by_genomes[genome_pair]
        binwidth = (100 - xmin) / 20
        sns.histplot(d, ax=ax, binwidth=binwidth, kde=False)
        ax.set_xlim(xmin, 100)
        ax.get_yaxis().set_visible(False)
        ax.set_xticks([xmin, 100])
        adjust_spines(ax, ["bottom"], outward=True)
        ax.spines["bottom"].set_color("lightslategray")

    title = opts.title
    italic_title = markup(f"*{title}*")
    root.text(
        0.5,
        0.95,
        f"{italic_title} identities between and within SS/SR/SO genomes",
        size=14,
        ha="center",
        va="center",
    )
    normalize_axes(root)
    image_name = f"SO_SR_SS.{title}." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def main():

    actions = (
        ("prepare", "Calculate lengths from real sugarcane data"),
        ("simulate", "Run simulation on female restitution"),
        # Plotting scripts to illustrate divergence between and within genomes
        ("divergence", "Plot divergence between and within SS/SR/SO genomes"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
