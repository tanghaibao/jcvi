#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# standalone_plotter.py
# notebooks
#
# Created by Haibao Tang on 12/02/19
# Copyright © 2019 Haibao Tang. All rights reserved.
#

import os.path as op
import sys
import seaborn as sns
from random import choice, random, sample
from itertools import groupby
from collections import Counter
import numpy as np
import matplotlib.pyplot as plt

from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir
from jcvi.graphics.base import savefig


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
        """ Randomly generate a gamete from current genome that
        """
        self.chromosomes.sort()
        gamete_chromosomes = []

        duplicate_chromosomes = []
        singleton_chromosomes = []
        for chromosome, chromosomes in groupby(self.chromosomes):
            chromosomes = list(chromosomes)
            ncopies = len(chromosomes)
            duplicate_chromosomes += [chromosome] * (ncopies // 2)
            if ncopies % 2 == 1:
                singleton_chromosomes.append(chromosome)

        # 1. Check for any chromosome that have 2 identical copies, if so, we will assume disomic
        # inheritance for that chromosome and always keep one and only copy
        gamete_chromosomes += duplicate_chromosomes

        def prefix(x):
            return x.split("_", 1)[0]

        for group, chromosomes in groupby(singleton_chromosomes, key=prefix):
            chromosomes = list(chromosomes)
            halfN = len(chromosomes) // 2
            # Odd number, e.g. 5, equal chance to be 2 or 3
            if len(chromosomes) % 2 != 0 and random() < 0.5:
                halfN += 1
            gamete_chromosomes += sorted(sample(chromosomes, halfN))
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

    def mate_2xnplusn(self, name, other_genome, verbose=True):
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
        def prefix(x):
            return x.split("-", 1)[0]

        total_count = 0
        total_unique = 0
        ans = []
        for group, chromosomes in groupby(self.chromosomes, prefix):
            chromosomes = list(chromosomes)
            group_count = len(chromosomes)
            group_unique = len(set(chromosomes))
            total_count += group_count
            total_unique += group_unique
            ans.append((group, group_count, group_unique))
        ans.append(("Total", total_count, total_unique))
        return ans

    def print_summary(self):
        print("[SUMMARY]")
        for group, group_count, group_unique in self.summary:
            print(f"{group}: count={group_count}, unique={group_unique}")


# Simulate two parents
SS = Genome("SS", "SS", 10, 8)
SO = Genome("SO", "SO", 8, 10)


def simulate_F1(verbose=True):
    SO_SS_F1_2xnplusn = SO.mate_2xnplusn("SOxSS F1", SS, verbose=verbose)
    if verbose:
        SO_SS_F1_2xnplusn.print_summary()
    return SO_SS_F1_2xnplusn


# Start F2 simulation (mode: 2xn + n)
def simulate_F2(verbose=True):
    SO_SS_F1_2xnplusn = simulate_F1(verbose=verbose)
    SO_SS_F2_nplusn = SO_SS_F1_2xnplusn.mate_nplusn(
        "SOxSS F2", SO_SS_F1_2xnplusn, verbose=verbose
    )
    if verbose:
        SO_SS_F2_nplusn.print_summary()
    return SO_SS_F2_nplusn


def simulate_BCn(n, verbose=True):
    SS_SO_F1_2xnplusn = simulate_F1(verbose=verbose)
    SS_SO_BC1_2xnplusn, SS_SO_BC2_nplusn, SS_SO_BC3_nplusn, SS_SO_BC4_nplusn = (
        None,
        None,
        None,
        None,
    )
    # BC1
    if n >= 1:
        SS_SO_BC1_2xnplusn = SO.mate_2xnplusn(
            "SSxSO BC1", SS_SO_F1_2xnplusn, verbose=verbose
        )
    # BC2
    if n >= 2:
        SS_SO_BC2_nplusn = SO.mate_nplusn(
            "SSxSO BC2", SS_SO_BC1_2xnplusn, verbose=verbose
        )
    # BC3
    if n >= 3:
        SS_SO_BC3_nplusn = SO.mate_nplusn(
            "SSxSO BC3", SS_SO_BC2_nplusn, verbose=verbose
        )
    if n >= 4:
        SS_SO_BC4_nplusn = SO.mate_nplusn(
            "SSxSO BC4", SS_SO_BC3_nplusn, verbose=verbose
        )
    return [
        None,
        SS_SO_BC1_2xnplusn,
        SS_SO_BC2_nplusn,
        SS_SO_BC3_nplusn,
        SS_SO_BC4_nplusn,
    ][n]


def plot_summary(ax, samples, title):
    SO_data = []
    SS_data = []
    for sample in samples:
        try:
            group, group_count, group_unique = [
                x for x in sample.summary if x[0] == "SO"
            ][0]
        except:
            group_unique = 0
        SO_data.append(group_unique)
        try:
            group, group_count, group_unique = [
                x for x in sample.summary if x[0] == "SS"
            ][0]
        except:
            group_unique = 0
        SS_data.append(group_unique)
    x, y = zip(*sorted(Counter(SS_data).items()))
    ss = ax.bar(x, y, color="purple", alpha=0.5, ec="purple")
    x, y = zip(*sorted(Counter(SO_data).items()))
    so = ax.bar(x, y, color="orange", alpha=0.5, ec="orange")
    ax.set_xlim(0, 80)
    ax.set_yticks([])
    ax.set_ylabel(
        title, rotation=0, labelpad=20, fontdict=dict(fontweight="bold", color="blue")
    )
    SO_mean, SO_std = np.mean(SO_data), np.std(SO_data)
    SS_mean, SS_std = np.mean(SS_data), np.std(SS_data)
    ax.legend(
        [so, ss],
        [
            f"So chrs ({SO_mean:.1f}±{SO_std:.1f})",
            f"Ss chrs ({SS_mean:.1f}±{SS_std:.1f})",
        ],
        loc="upper right",
    )


def write_chromosomes(genomes, filename):
    """ Write simulated chromosomes to file

    Args:
        genomes (List[Genome]): List of simulated genomes.
        filename: File path to write to.
    """
    print(f"Write chromosomes to `{filename}`", file=sys.stderr)
    with open(filename, "w") as fw:
        for genome in genomes:
            print(genome, file=fw)


def simulate(args):
    """
    %prog simulate

    Run simulation on female restitution.
    """
    p = OptionParser(simulate.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="6x12")
    if len(args) != 0:
        sys.exit(not p.print_help())

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax1 = plt.subplot(611)
    ax2 = plt.subplot(612, sharex=ax1)
    ax3 = plt.subplot(613, sharex=ax1)
    ax4 = plt.subplot(614, sharex=ax1)
    ax5 = plt.subplot(615, sharex=ax1)
    ax6 = plt.subplot(616, sharex=ax1)
    for ax in (ax1, ax2, ax3, ax4, ax5):
        plt.setp(ax.get_xticklabels(), visible=False)

    # Prepare the simulated data
    all_F1s = [simulate_F1(verbose=False) for x in range(1000)]
    all_F2s = [simulate_F2(verbose=False) for x in range(1000)]
    all_BC1s = [simulate_BCn(1, verbose=False) for x in range(1000)]
    all_BC2s = [simulate_BCn(2, verbose=False) for x in range(1000)]
    all_BC3s = [simulate_BCn(3, verbose=False) for x in range(1000)]
    all_BC4s = [simulate_BCn(4, verbose=False) for x in range(1000)]

    plot_summary(ax1, all_F1s, "$F_1$")
    plot_summary(ax2, all_F2s, "$F_2$")
    plot_summary(ax3, all_BC1s, "$BC_1$")
    plot_summary(ax4, all_BC2s, "$BC_2$")
    plot_summary(ax5, all_BC3s, "$BC_3$")
    plot_summary(ax6, all_BC4s, "$BC_4$")
    ax6.set_xlabel("Number of unique chromosomes")
    savefig("plotter.png", dpi=120)
    savefig("plotter.pdf", dpi=120)

    outdir = "simulations"
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


def main():

    actions = (("simulate", "Run simulation on female restitution"),)
    p = ActionDispatcher(actions)
    p.dispatch(globals())


if __name__ == "__main__":
    main()
