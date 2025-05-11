#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# sugarcane.py
# projects
#
# Created by Haibao Tang on 12/02/19
# Copyright © 2019 Haibao Tang. All rights reserved.
#
"""
Simulate sugarcane genomes and analyze the diversity in the progeny genomes.
"""

from collections import Counter, defaultdict
from dataclasses import dataclass
from enum import Enum
from itertools import combinations, groupby, product
from multiprocessing import Pool, cpu_count
import os.path as op
from random import randint, random
import sys
from typing import Dict, List, Tuple
from uuid import uuid4

import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns

from ..apps.base import ActionDispatcher, OptionParser, flatten, logger, mkdir
from ..formats.blast import Blast
from ..graphics.base import (
    Rectangle,
    adjust_spines,
    get_shades,
    markup,
    normalize_axes,
    savefig,
)
from ..graphics.chromosome import Chromosome as ChromosomePlot
from ..utils.cbook import short_float

SoColor = "#7436a4"  # Purple
SsColor = "#5a8340"  # Green
HAPLOID_GENE_COUNT = 160
SO_PLOIDY = 8
SS_PLOIDY = 16
SO_GENE_COUNT = SO_PLOIDY * HAPLOID_GENE_COUNT
SS_GENE_COUNT = SS_PLOIDY * HAPLOID_GENE_COUNT
SO_CHROM_COUNT = 10
SS_CHROM_COUNT = 8
CROSSES = ("F1", "BC1", "BC2", "BC3", "BC4")
TITLES = (
    (r"$\mathrm{F_1}$", None),
    (r"$\mathrm{BC_1}$", None),
    (r"$\mathrm{BC_2}$", None),
    (r"$\mathrm{BC_3}$", None),
    (r"$\mathrm{BC_4}$", None),
)


class CrossMode(Enum):
    """
    How the F1 is generated.
    """

    nplusn = "n+n"
    nx2plusn = "nx2+n"
    twoplusnFDR = "2n+n_FDR"
    twoplusnSDR = "2n+n_SDR"


@dataclass(frozen=True)
class Gene:
    """
    A gene in a chromosome.
    """

    chrom: str
    haplotype: str
    idx: int


@dataclass
class Chromosome(list):
    """
    A chromosome with genes.
    """

    subgenome: str
    chrom: str
    genes: List[Gene]
    uuid: str

    def __init__(self, subgenome: str, chrom: str, haplotype: str, gene_count: int):
        self.subgenome = subgenome
        self.chrom = chrom
        self.genes = [Gene(chrom, haplotype, i + 1) for i in range(gene_count)]
        self.uuid = str(uuid4())

    def __str__(self) -> str:
        """Merge genes into consecutive ranges."""
        ranges = []
        for haplotype, genes in groupby(self.genes, key=lambda x: x.haplotype):
            genes = list(genes)
            start, end = genes[0].idx, genes[-1].idx
            ranges.append(f"{haplotype}{start}-{end}")
        return f"{self.subgenome}-{self.chrom}|{','.join(ranges)}"

    def __len__(self) -> int:
        return len(self.genes)

    def duplicate(self) -> "Chromosome":
        """Duplicate the chromosome."""
        return Chromosome.make(self.subgenome, self.chrom, self.genes)

    @classmethod
    def make(cls, subgenome: str, chrom: str, genes: List[Gene]) -> "Chromosome":
        chrom = Chromosome(subgenome, chrom, "", 0)
        chrom.genes = genes
        return chrom

    def num_matching_genes(self, other: "Chromosome") -> int:
        """Count the number of matching genes between two chromosomes"""
        return sum(1 for x, y in zip(self.genes, other.genes) if x == y)

    def truncate_n(self, n: int):
        """Truncate the chromosome to n genes."""
        self.genes = self.genes[:n]

    def truncate(self, ratio: float):
        """Truncate the chromosome to a ratio of genes."""
        before = str(self)
        n = round(len(self) * ratio)
        self.truncate_n(n)
        logger.info("Truncate %s to %s (%d)", before, str(self), n)


# Simulate genome composition
class Genome:
    """
    A genome with a collection of chromosomes.
    """

    name: str
    chromosomes: List[Chromosome]

    def __init__(
        self,
        name: str,
        subgenome: str,
        ploidy: int,
        haploid_chromosome_count: int,
        haploid_gene_count: int,
    ):
        """
        Simulate a genome with given ploidy and haploid_chromosome_count. Example:

        >>> print(Genome("test", "pf", 2, 3, 90))
        test: pf-chr01|a1-30;pf-chr01|b1-30;pf-chr02|a1-30;pf-chr02|b1-30;pf-chr03|a1-30;pf-chr03|b1-30
        """
        self.name = name
        chrom_gene_count = (
            haploid_gene_count // haploid_chromosome_count
            if haploid_chromosome_count
            else 0
        )
        chromosomes = []
        for i in range(haploid_chromosome_count):
            chromosomes += [
                Chromosome(
                    subgenome, f"chr{i + 1:02d}", chr(ord("a") + j), chrom_gene_count
                )
                for j in range(ploidy)
            ]
        self.chromosomes = chromosomes

    def __str__(self):
        return self.name + ": " + ";".join(str(_) for _ in self.chromosomes)

    def ploidy(self, chrom: str) -> int:
        """
        Return the ploidy of a chromosome.
        """
        return sum(1 for x in self.chromosomes if x.chrom == chrom)

    @classmethod
    def from_str(cls, s: str) -> "Genome":
        """
        Parse a string representation of a genome.

        >>> s = "test: pf-chr01|a1-30;pf-chr01|b1-30;pf-chr02|a1-30;pf-chr02|b1-30;pf-chr03|a1-30;pf-chr03|b1-30"
        >>> str(Genome.from_str(s)) == s
        True
        """
        name, chroms = s.split(": ")
        genome = Genome(name, "", 0, 0, 0)
        chromosomes = []
        for x in chroms.split(";"):
            chrom, genes = x.split("|")
            subgenome, chrom = chrom.split("-")
            cgenes = []
            for genes in genes.split(","):
                haplotype = genes[0]
                start, end = [int(x) for x in genes[1:].split("-")]
                cgenes += [Gene(chrom, haplotype, i) for i in range(start, end + 1)]
            chrom = Chromosome.make(subgenome, chrom, cgenes)
            chromosomes.append(chrom)
        genome.chromosomes = chromosomes
        return genome

    @classmethod
    def make(cls, name: str, chromosomes: List[Chromosome]) -> "Genome":
        genome = Genome(name, "", 0, 0, 0)
        genome.chromosomes = sorted(chromosomes, key=genome._sort_key)
        return genome

    def _sort_key(self, x: Chromosome):
        return x.subgenome, x.chrom

    def pair_with_weights(self) -> List[List[str]]:
        """
        Pair chromosomes by similarity, results are sorted in descending similarity.
        """
        G = nx.Graph()
        scores = {}
        self.chromosomes.sort(key=self._sort_key)
        for _, chromosomes in groupby(self.chromosomes, key=self._sort_key):
            for a, b in combinations(chromosomes, 2):
                weight = a.num_matching_genes(b) + 1  # +1 so we don't have zero weight
                G.add_edge(a.uuid, b.uuid, weight=weight)
                scores[tuple(sorted((a.uuid, b.uuid)))] = weight
        # Find the maximum matching
        matching = nx.max_weight_matching(G)
        return sorted(matching, key=lambda x: scores[tuple(sorted(x))], reverse=True)

    def pair_chromosomes(
        self, inplace=False
    ) -> Tuple[List[List[Chromosome]], List[Chromosome]]:
        """
        Pair chromosomes by similarity.
        """
        matching = self.pair_with_weights()
        # Partition the chromosomes into paired and singleton
        paired = set()
        pairs = []
        for a, b in matching:
            paired.add(a)
            paired.add(b)
            pair = [x for x in self.chromosomes if x.uuid in (a, b)]
            pairs.append(pair)
        singletons = [x for x in self.chromosomes if x.uuid not in paired]
        if inplace:
            reordered = []
            for a, b in pairs:
                reordered += [a, b]
            reordered += singletons
            reordered.sort(key=self._sort_key)
            self.chromosomes = reordered
        return pairs, singletons

    def _crossover_chromosomes(
        self, a: Chromosome, b: Chromosome, sdr: bool
    ) -> List[Chromosome]:
        """
        Crossover two chromosomes.
        """
        if a.subgenome != b.subgenome or a.chrom != b.chrom or len(a) != len(b):
            raise ValueError("Incompatible chromosomes. Crossover failed.")
        subgenome, chrom = a.subgenome, a.chrom
        n = len(a)
        if random() < 0.5:
            a, b = b, a
        if random() < 0.5:
            crossover_point = randint(0, n // 2 - 1)
            recombinant = a.genes[:crossover_point] + b.genes[crossover_point:]
            non_recombinant = b.genes
        else:
            crossover_point = randint(n // 2, n - 1)
            recombinant = a.genes[:crossover_point] + b.genes[crossover_point:]
            non_recombinant = a.genes
        recombinant = Chromosome.make(subgenome, chrom, recombinant)
        non_recombinant = Chromosome.make(subgenome, chrom, non_recombinant)
        return [recombinant, non_recombinant] if sdr else [recombinant]

    def _gamete(self, sdr: bool):
        """Randomly generate a gamete from current genome."""
        gamete_chromosomes = []
        paired_chromosomes, singleton_chromosomes = self.pair_chromosomes()
        for a, b in paired_chromosomes:
            gamete_chromosomes += self._crossover_chromosomes(a, b, sdr=sdr)

        # Randomly assign the rest, singleton chromosomes
        for a in singleton_chromosomes:
            if random() < 0.5:
                gamete_chromosomes.append(a)

        tag = "sdr gamete" if sdr else "gamete"
        return Genome.make(f"{self.name} {tag}", gamete_chromosomes)

    @property
    def gamete(self):
        return self._gamete(sdr=False)

    @property
    def gamete_fdr(self):
        return Genome.make(f"{self.name} fdr gamete", self.chromosomes)

    @property
    def gamete_sdr(self):
        return self._gamete(sdr=True)

    def mate_nplusn(self, name: str, other_genome: "Genome", verbose: bool = True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (n+n)", file=sys.stderr
            )
        f1_chromosomes = self.gamete.chromosomes + other_genome.gamete.chromosomes
        return Genome.make(name, f1_chromosomes)

    def mate_nx2plusn(self, name: str, other_genome: "Genome", verbose: bool = True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (nx2+n)",
                file=sys.stderr,
            )
        gamete = self.gamete.chromosomes
        duplicate = [x.duplicate() for x in gamete]
        f1_chromosomes = gamete + duplicate + other_genome.gamete.chromosomes
        return Genome.make(name, f1_chromosomes)

    def mate_2nplusn_FDR(self, name: str, other_genome: "Genome", verbose: bool = True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (2n+n_FDR)",
                file=sys.stderr,
            )
        f1_chromosomes = self.gamete_fdr.chromosomes + other_genome.gamete.chromosomes
        return Genome.make(name, f1_chromosomes)

    def mate_2nplusn_SDR(self, name: str, other_genome: "Genome", verbose: bool = True):
        if verbose:
            print(
                f"Crossing '{self.name}' x '{other_genome.name}' (2n+n_SDR)",
                file=sys.stderr,
            )
        f1_chromosomes = sorted(
            self.gamete_sdr.chromosomes + other_genome.gamete.chromosomes
        )
        return Genome.make(name, f1_chromosomes)

    @property
    def summary(self) -> List[Tuple[str, int, int, int]]:
        ans = []
        total_chrom_count = 0
        total_so_size = 0
        total_ss_size = 0

        self.chromosomes.sort(key=lambda x: x.subgenome)
        for subgenome, chromosomes in groupby(
            self.chromosomes, key=lambda x: x.subgenome
        ):
            chromosomes = list(chromosomes)
            uniq_genes = set(flatten(x.genes for x in chromosomes))
            group_count = len(uniq_genes)
            group_chrom_count = group_count / len(chromosomes[0])
            group_so_size = group_count if subgenome == "SO" else 0
            group_ss_size = group_count if subgenome == "SS" else 0
            group_size = group_so_size + group_ss_size
            ans.append(
                (
                    subgenome,
                    group_chrom_count,
                    group_so_size / group_size,
                    group_ss_size / group_size,
                )
            )
            total_chrom_count += group_chrom_count
            total_so_size += group_so_size
            total_ss_size += group_ss_size
        total_size = total_so_size + total_ss_size
        ans.append(
            (
                "Total",
                total_chrom_count,
                total_so_size / total_size,
                total_ss_size / total_size,
            )
        )
        return ans

    def print_summary(self):
        print("[SUMMARY]")
        for group, group_chrom_count, group_so_size, group_ss_size in self.summary:
            print(
                f"{group}: chrom_count={group_chrom_count}, SO={group_so_size}, SS={group_ss_size}"
            )

    def truncate_chromosome(self, subgenome: str, chrom: str, ratio: float):
        """
        Truncate the specified chromosome.
        """
        for c in self.chromosomes:
            if c.subgenome == subgenome and c.chrom == chrom:
                c.truncate(ratio)


class GenomeSummary:
    def __init__(self, SO_data, SS_data, percent_SO_data, percent_SS_data):
        self.SO_data = SO_data
        self.SS_data = SS_data
        self.percent_SO_data = percent_SO_data
        self.percent_SS_data = percent_SS_data

    def _summary(self, a, tag, precision=0) -> Tuple[str, str]:
        mean, mn, mx = (
            round(np.mean(a), precision),
            round(np.min(a), precision),
            round(np.max(a), precision),
        )
        s = f"*{tag}* chr: {mean:.0f}"
        if mn == mean and mx == mean:
            return s, ""
        return s, f" ({mn:.0f}-{mx:.0f})"

    def _percent_summary(self, a, tag, precision=1) -> Tuple[str, str]:
        mean, mn, mx = (
            round(np.mean(a), precision),
            round(np.min(a), precision),
            round(np.max(a), precision),
        )
        s = f"*{tag}*%: {short_float(mean, precision)}%"
        print(s)
        if mn == mean and mx == mean:
            return s, ""
        return s, f"({short_float(mn, precision)}-{short_float(mx, precision)}%)"

    @property
    def percent_SO_summary(self) -> Tuple[str, str]:
        return self._percent_summary(self.percent_SO_data, "So")

    @property
    def percent_SS_summary(self) -> Tuple[str, str]:
        return self._percent_summary(self.percent_SS_data, "Ss")

    @property
    def SO_summary(self) -> Tuple[str, str]:
        return self._summary(self.SO_data, "So")

    @property
    def SS_summary(self) -> Tuple[str, str]:
        return self._summary(self.SS_data, "Ss")


def simulate_F1(SO: Genome, SS: Genome, mode: CrossMode, verbose: bool = False):
    SO_SS_F1 = None
    if mode == CrossMode.nx2plusn:
        SO_SS_F1 = SO.mate_nx2plusn("SOxSS F1", SS, verbose=verbose)
    elif mode == CrossMode.twoplusnFDR:
        SO_SS_F1 = SO.mate_2nplusn_FDR("SOxSS F1", SS, verbose=verbose)
    elif mode == CrossMode.twoplusnSDR:
        SO_SS_F1 = SO.mate_2nplusn_SDR("SOxSS F1", SS, verbose=verbose)
    if verbose and SO_SS_F1:
        SO_SS_F1.print_summary()
    return SO_SS_F1


def simulate_BC1(SO: Genome, SS_SO_F1: Genome, mode: CrossMode, verbose=False):
    SS_SO_BC1 = None
    if mode == CrossMode.nx2plusn:
        SS_SO_BC1 = SO.mate_nx2plusn("SOxSS BC1", SS_SO_F1, verbose=verbose)
    elif mode == CrossMode.twoplusnFDR:
        SS_SO_BC1 = SO.mate_2nplusn_FDR("SOxSS BC1", SS_SO_F1, verbose=verbose)
    elif mode == CrossMode.twoplusnSDR:
        SS_SO_BC1 = SO.mate_2nplusn_SDR("SOxSS BC1", SS_SO_F1, verbose=verbose)
    return SS_SO_BC1


def simulate_F1_to_BC4(
    SO: Genome, SS: Genome, mode: CrossMode, verbose=False
) -> List[Genome]:
    SS_SO_F1 = simulate_F1(SO, SS, mode=mode, verbose=verbose)
    SS_SO_BC1 = simulate_BC1(SO, SS_SO_F1, mode=mode, verbose=verbose)
    SS_SO_BC2_nplusn = SO.mate_nplusn("SOxSS BC2", SS_SO_BC1, verbose=verbose)
    SS_SO_BC3_nplusn = SO.mate_nplusn("SOxSS BC3", SS_SO_BC2_nplusn, verbose=verbose)
    SS_SO_BC4_nplusn = SO.mate_nplusn("SOxSS BC4", SS_SO_BC3_nplusn, verbose=verbose)
    return [
        SS_SO_F1,
        SS_SO_BC1,
        SS_SO_BC2_nplusn,
        SS_SO_BC3_nplusn,
        SS_SO_BC4_nplusn,
    ]


def plot_summary(ax, samples: List[Genome]) -> GenomeSummary:
    """Plot the distribution of chromosome numbers given simulated samples.

    Args:
        ax (Axes): Matplotlib axes.
        samples (List[Genome]): Summarized genomes.

    Returns:
        GenomeSummary: Summary statistics of simulated genomes.
    """
    SO_data = []
    SS_data = []
    percent_SO_data = []
    percent_SS_data = []
    for s in samples:
        summary = s.summary
        SO_summary = [x for x in summary if x[0] == "SO"]
        group_unique = 0
        if SO_summary:
            _, group_unique, _, _ = SO_summary[0]
        SO_data.append(round(group_unique))
        SS_summary = [x for x in summary if x[0] == "SS"]
        group_unique = 0
        if SS_summary:
            _, group_unique, _, _ = SS_summary[0]
        SS_data.append(round(group_unique))
        total_tag, _, total_so_size, total_ss_size = summary[-1]
        assert total_tag == "Total"
        percent_SO_data.append(total_so_size * 100)
        percent_SS_data.append(total_ss_size * 100)
    # Avoid overlapping bars
    SS_counter, SO_counter = Counter(SS_data), Counter(SO_data)
    print(SS_counter, sum(SS_counter.values()))
    print(SO_counter, sum(SO_counter.values()))
    overlaps = SS_counter.keys() & SO_counter.keys()
    shift = 0.5  # used to offset bars a bit to avoid cluttering
    if overlaps:
        for overlap in overlaps:
            logger.debug("Modify bar offsets at %s due to SS and SO overlaps", overlap)
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
        logger.debug("Modify bar offsets at %d due to end of range ends", value)
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
    ymax = len(samples) * 0.8
    ax.set_xlim(80, 0)
    ax.set_ylim(0, ymax)
    ax.set_yticks([])
    summary = GenomeSummary(SO_data, SS_data, percent_SO_data, percent_SS_data)

    # Write the stats summary within the plot
    summary_style = dict(size=8, va="center")
    SO_peak = SO_counter.most_common(1)[0][0]
    SS_peak = SS_counter.most_common(1)[0][0]
    SO_single = len(SO_counter) == 1
    SS_single = len(SS_counter) == 1

    for xpos, ypos, single, text, color, ha in zip(
        [SO_peak] * 4 + [SS_peak] * 4,
        ([ymax * 0.85] * 2 + [ymax * 0.65] * 2) * 2,
        [SO_single] * 4 + [SS_single] * 4,
        summary.SO_summary
        + summary.percent_SO_summary
        + summary.SS_summary
        + summary.percent_SS_summary,
        [SoColor] * 4 + [SsColor] * 4,
        ["right", "left"] * 4,
    ):
        # Offset some text to avoid overlapping
        if abs(SS_peak - SO_peak) < 20 and xpos == SO_peak:
            xpos += 12 if SO_peak > SS_peak else -12
        PAD = 1 if single else 0.25
        if ha == "left":
            xpos -= PAD
        else:
            xpos += PAD
        ax.text(xpos, ypos, markup(text), color=color, ha=ha, **summary_style)
    return summary


def write_chromosomes(genomes: List[Genome], filename: str):
    """Write simulated chromosomes to file

    Args:
        genomes (List[Genome]): List of simulated genomes.
        filename (str): File path to write to.
    """
    print(f"Write chromosomes to `{filename}`", file=sys.stderr)
    with open(filename, "w", encoding="utf-8") as fw:
        for genome in genomes:
            print(genome, file=fw)


def get_mode_title(mode: CrossMode) -> str:
    """
    Get the title of the mode.
    """
    if mode == CrossMode.nx2plusn:
        mode_title = r"$n_1\times2 + n_2$"
    elif mode == CrossMode.twoplusnFDR:
        mode_title = r"$2n + n$ (FDR)"
    elif mode == CrossMode.twoplusnSDR:
        mode_title = r"$n_1^*\times2 + n$ (SDR)"
    else:
        mode_title = "Unknown"
    return mode_title


def simulate(args):
    """
    %prog simulate [2n+n_FDR|2n+n_SDR|nx2+n]

    Run simulation on female restitution. There are several modes:
    - 2n+n_FDR: merger between a somatic and a germline
    - 2n+n_SDR: merger between a recombined germline and a germline (not yet supported)
    - nx2+n: merger between a doubled germline and a germline

    These modes would impact the sequence diversity in the progeny
    genome in F1, BCn ... the goal of this simulation, is thus to
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
    p.add_argument("-N", default=1000, type=int, help="Number of simulated samples")
    p.add_argument("--ss-ploidy", default=16, type=int, help="SS ploidy")
    opts, args, iopts = p.set_image_options(args, figsize="6x6")
    if len(args) != 1:
        sys.exit(not p.print_help())

    (mode,) = args
    SS_PLOIDY = opts.ss_ploidy
    SS_GENE_COUNT = SS_PLOIDY * HAPLOID_GENE_COUNT
    mode = CrossMode(mode)
    logger.info("Transmission: %s", mode)

    # Construct a composite figure with 6 tracks
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    rows = 5
    ypad = 0.05
    yinterval = (1 - 2 * ypad) / (rows + 1)
    yy = 1 - ypad
    xpad = 0.2
    xwidth = 0.7

    # Axes are vertically stacked, and share x-axis
    axes = []
    yy_positions = []  # Save yy positions so we can show details to the right later
    for idx in range(rows):
        yy_positions.append(yy)
        yy -= yinterval
        ax = fig.add_axes([xpad, yy, xwidth, yinterval * 0.85])
        if idx != rows - 1:
            plt.setp(ax.get_xticklabels(), visible=False)
        axes.append(ax)
    ax1, ax2, ax3, ax4, ax5 = axes

    # Prepare the simulated data
    # Simulate two parents
    SS = Genome("SS", "SS", SS_PLOIDY, SS_CHROM_COUNT, HAPLOID_GENE_COUNT)
    SO = Genome("SO", "SO", SO_PLOIDY, SO_CHROM_COUNT, HAPLOID_GENE_COUNT)

    verbose = opts.verbose
    N = opts.N
    with Pool(cpu_count()) as pool:
        results = pool.starmap(simulate_F1_to_BC4, [(SO, SS, mode, verbose)] * N)
    all_F1s, all_BC1s, all_BC2s, all_BC3s, all_BC4s = zip(*results)

    # Plotting
    plot_summary(ax1, all_F1s)
    plot_summary(ax2, all_BC1s)
    plot_summary(ax3, all_BC2s)
    plot_summary(ax4, all_BC3s)
    plot_summary(ax5, all_BC4s)

    # Show title to the left
    xx = xpad / 2
    for (title, subtitle), yy in zip(TITLES, yy_positions):
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
    mode_title = get_mode_title(mode)
    root.text(0.5, 0.95, f"Transmission: {mode_title}", ha="center")

    savefig(f"{mode}.ss_ploidy_{SS_PLOIDY}.pdf", dpi=120)

    outdir = f"simulations_{mode}"
    mkdir(outdir)
    # Write chromosomes to disk
    for genomes, cross in zip(
        [all_F1s, all_BC1s, all_BC2s, all_BC3s, all_BC4s], CROSSES
    ):
        write_chromosomes(genomes, op.join(outdir, f"all_{cross}"))


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


def plot_genome(
    ax,
    x: float,
    y: float,
    genome: Genome,
    haplotype_gene_height: Dict[str, float],
    haplotype_colors: Dict[str, str],
    chrom_width: float = 0.012,
    gap_width: float = 0.008,
    pair: bool = False,
):
    """
    Plot the genome in the axes, centered around (x, y).
    """
    target = "chr01"  # Arbitrary target chromosome
    ploidy = genome.ploidy(target)
    total_width = ploidy * (chrom_width + gap_width) - gap_width
    xx = x - total_width / 2
    if pair:
        genome.pair_chromosomes(inplace=True)
        print(
            genome.name, ":", [str(x) for x in genome.chromosomes if x.chrom == target]
        )
    for chrom in genome.chromosomes:
        if chrom.chrom != target:
            continue
        gene_height = haplotype_gene_height[chrom.subgenome]
        height = gene_height * len(chrom)
        ChromosomePlot(ax, xx, y, y - height, ec="lightslategray")
        yy = y
        subgenome = chrom.subgenome
        for haplotype, genes in groupby(chrom.genes, key=lambda x: x.haplotype):
            genes = list(genes)
            g1, g2 = genes[0].idx - 1, genes[-1].idx
            patch_height = gene_height * (g2 - g1)
            color = haplotype_colors[subgenome, haplotype]
            ax.add_patch(
                Rectangle(
                    (xx - chrom_width / 2, yy - patch_height),
                    chrom_width,
                    patch_height,
                    fc=color,
                    lw=0,
                )
            )
            yy -= patch_height
        xx += chrom_width + gap_width


def truncate_SS_chromosome(genome: Genome, cross: str, ss_ploidy: int) -> Genome:
    """
    Truncate the last SS chromosome to illustrate half/quarter ploidy.
    """
    subgenome, chrom = "SS", "chr01"
    if ss_ploidy == 8:
        if cross == "BC3":
            genome.truncate_chromosome(subgenome, chrom, 0.5)
        elif cross == "BC4":
            genome.truncate_chromosome(subgenome, chrom, 0.25)
    elif ss_ploidy == 16:
        if cross == "BC4":
            genome.truncate_chromosome(subgenome, chrom, 0.5)
    else:
        raise ValueError("Unsupported SS ploidy for custom mode")
    return genome


def chromosome(args):
    """
    %prog chromosome [2n+n_FDR|2n+n_SDR|nx2+n]
    """
    p = OptionParser(simulate.__doc__)
    p.add_argument("-k", default=0, type=int, help="Plot k-th simulated genomes")
    p.add_argument("--ss-ploidy", default=16, type=int, help="SS ploidy")
    p.add_argument(
        "--custom",
        default=False,
        action="store_true",
        help="Custom mode to truncate BC3/4 SS chromosome",
    )
    opts, args, iopts = p.set_image_options(args, figsize="6x6")
    if len(args) != 1:
        sys.exit(not p.print_help())

    (mode,) = args
    SS_PLOIDY = opts.ss_ploidy
    SS_GENE_COUNT = SS_PLOIDY * HAPLOID_GENE_COUNT
    mode = CrossMode(mode)
    logger.info("Transmission: %s", mode)
    if opts.custom:
        logger.info("Custom mode: truncate SS chromosome")

    # Construct a composite figure with 6 tracks
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    rows = 5
    ypad = 0.05
    yinterval = (1 - 2 * ypad) / (rows + 1)
    yy = 1 - ypad
    xpad = 0.2

    SS = Genome("SS", "SS", SS_PLOIDY, SS_CHROM_COUNT, HAPLOID_GENE_COUNT)
    SO = Genome("SO", "SO", SO_PLOIDY, SO_CHROM_COUNT, HAPLOID_GENE_COUNT)

    indir = f"simulations_{mode}"
    genomes = []
    for cross in CROSSES:
        filename = op.join(indir, f"all_{cross}")
        with open(filename, encoding="utf-8") as fp:
            for i, row in enumerate(fp):
                if i != opts.k:
                    continue
                genome = Genome.from_str(row)
                if opts.custom:
                    truncate_SS_chromosome(genome, cross, opts.ss_ploidy)
                break
        genomes.append((cross, genome))

    yy_positions = []  # Save yy positions so we can show details to the right later
    for idx in range(rows + 1):
        yy_positions.append(yy)
        yy -= yinterval
        if idx != rows:
            root.plot([0, 1], [yy, yy], lw=0.5, color="gray")

    # Show title to the left
    xx = xpad / 2
    titles = (("Progenitors", None),) + TITLES
    for (title, _), yy in zip(titles, yy_positions):
        yy -= 0.07
        root.text(
            xx,
            yy,
            markup(title),
            color="darkslategray",
            ha="center",
            va="center",
            fontweight="semibold",
        )

    SO_colors = get_shades(SoColor, SO_PLOIDY)
    SS_colors = get_shades(SsColor, SS_PLOIDY)
    SO_haplotypes = [("SO", chr(ord("a") + i)) for i in range(SO_PLOIDY)]
    SS_haplotypes = [("SS", chr(ord("a") + i)) for i in range(SS_PLOIDY)]
    SO_haplotype_colors = dict(zip(SO_haplotypes, SO_colors))
    SS_haplotype_colors = dict(zip(SS_haplotypes, SS_colors))
    haplotype_colors = {**SO_haplotype_colors, **SS_haplotype_colors}

    # Plotting
    chrom_height = 0.1
    SO_GENE_HEIGHT = chrom_height / (HAPLOID_GENE_COUNT // SO_CHROM_COUNT)
    SS_GENE_HEIGHT = chrom_height / (HAPLOID_GENE_COUNT // SS_CHROM_COUNT)
    haplotype_gene_height = {"SO": SO_GENE_HEIGHT, "SS": SS_GENE_HEIGHT}
    yy = 0.92
    plot_genome(root, 0.35, yy, SO, haplotype_gene_height, haplotype_colors)
    plot_genome(
        root,
        0.75 if SS_PLOIDY == 16 else 0.66,
        yy,
        SS,
        haplotype_gene_height,
        haplotype_colors,
    )
    # Plot big cross sign
    root.text(
        0.5, yy - chrom_height / 2, r"$\times$", ha="center", va="center", fontsize=36
    )
    # Genome labels
    root.text(
        0.215,
        yy - chrom_height / 2,
        markup(f"*So*\n({SO_PLOIDY}x)"),
        ha="center",
        va="center",
        color=SoColor,
    )
    root.text(
        0.945 if SS_PLOIDY == 16 else 1 - 0.215,
        yy - chrom_height / 2,
        markup(f"*Ss*\n({SS_PLOIDY}x)"),
        ha="center",
        va="center",
        color=SsColor,
    )

    for _, genome in genomes:
        yy -= yinterval
        plot_genome(
            root, 0.5, yy, genome, haplotype_gene_height, haplotype_colors, pair=True
        )

    # Title
    mode_title = get_mode_title(mode)
    root.text(0.5, 0.95, f"Transmission: {mode_title}", ha="center")
    normalize_axes(root)

    savefig(f"{mode}.chromosome.ss_ploidy_{SS_PLOIDY}.pdf", dpi=120)


def main():

    actions = (
        ("simulate", "Run simulation on female restitution"),
        # Plot the simulated chromosomes
        ("chromosome", "Plot the chromosomes of the simulated genomes"),
        # Plotting scripts to illustrate divergence between and within genomes
        ("divergence", "Plot divergence between and within SS/SR/SO genomes"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
