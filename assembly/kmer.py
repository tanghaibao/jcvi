#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deals with K-mers and K-mer distribution from reads or genome
"""

import sys

from optparse import OptionParser

from jcvi.utils.iter import pairwise
from jcvi.graphics.base import plt, _
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('histogram', 'plot the histogram based on meryl K-mer distribution'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def histogram(args):
    """
    %prog histogram meryl.histogram species N totalKmers

    Plot the histogram based on meryl K-mer distribution, species and N are
    only used to annotate the graphic
    """
    p = OptionParser(histogram.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(p.print_help())

    histfile, species, N, totalKmers = args
    fp = open(histfile)
    hist = {}
    for row in fp:
        K, counts = row.split()[:2]
        K, counts = int(K), int(counts)
        Kcounts = K * counts
        hist[K] = Kcounts

    history = ["drop"]
    for a, b in pairwise(sorted(hist.items())):
        Ka, ca = a
        Kb, cb = b
        if ca <= cb:
            status = "rise"
        else:
            status = "drop"
        if history[-1] != status:
            history.append(status)
        #print Ka, ca, status
        if history == ["drop", "rise", "drop"]:
            break

    Total_Kmers = int(totalKmers)
    Kmer_coverage = Ka
    Genome_size = Total_Kmers * 1. / Ka / 1e6

    Total_Kmers_msg = "Total {0}-mers: {1}".format(N, Total_Kmers)
    Kmer_coverage_msg = "{0}-mer coverage: {1}".format(N, Kmer_coverage)
    Genome_size_msg = "Estimated genome size: {0:.1f}Mb".format(Genome_size)

    for msg in (Total_Kmers_msg, Kmer_coverage_msg, Genome_size_msg):
        print >> sys.stderr, msg

    fig = plt.figure(1, (6, 6))
    counts = sorted((a, b) for a, b in hist.items() if a <= 100)
    x, y = zip(*counts)
    plt.plot(x, y, 'g-', lw=2, alpha=.5)

    ax = plt.gca()
    ax.text(.5, .9, _(Total_Kmers_msg),
            ha="center", color='b', transform=ax.transAxes)
    ax.text(.5, .8, _(Kmer_coverage_msg),
            ha="center", color='b', transform=ax.transAxes)
    ax.text(.5, .7, _(Genome_size_msg),
            ha="center", color='b', transform=ax.transAxes)

    title = "{0} genome {1}-mer histogram".format(species, N)
    ax.set_title(_(title), color='r')
    xlabel, ylabel = "Coverage (X)", "Counts"
    ax.set_xlabel(_(xlabel), color='r')
    ax.set_ylabel(_(ylabel), color='r')

    imagename = histfile.split(".")[0] + ".pdf"
    plt.savefig(imagename, dpi=100)
    print >>sys.stderr, "Image saved to {0}.".format(imagename)


if __name__ == '__main__':
    main()
