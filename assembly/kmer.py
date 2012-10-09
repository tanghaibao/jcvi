#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Deals with K-mers and K-mer distribution from reads or genome
"""

import os.path as op
import sys
import logging

from optparse import OptionParser

from jcvi.utils.iter import pairwise
from jcvi.graphics.base import plt, asciiplot, _, set_human_axis, savefig
from jcvi.apps.base import ActionDispatcher, sh, debug, need_update
debug()


def main():

    actions = (
        ('jellyfish', 'dump histogram using `jellyfish`'),
        ('meryl', 'dump histogram using `meryl`'),
        ('histogram', 'plot the histogram based on meryl K-mer distribution'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def jellyfish(args):
    """
    %prog jellyfish *.fastq

    Run jellyfish to dump histogram to be used in kmer.histogram().
    """
    from jcvi.apps.base import getfilesize
    from jcvi.utils.cbook import human_size
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(jellyfish.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastqfiles = args
    K = opts.K

    totalfilesize = sum(getfilesize(x) for x in fastqfiles)
    hashsize = int(totalfilesize / 46)
    #hashsize = max(hashsize, 4000000000)  # based on msr-ca

    logging.debug("Total file size: {0}, hashsize (-s): {1}".\
                    format(human_size(totalfilesize,
                           a_kilobyte_is_1024_bytes=True), hashsize))

    offset = guessoffset([fastqfiles[0]])
    assert all(guessoffset([x]) == offset for x in fastqfiles[1:])

    jfpf = "jf-{0}".format(K)
    jfdb = jfpf + "_0"

    cmd = "jellyfish count -t 32 -p 126 -C -r -o {0}".format(jfpf)
    cmd += " -s {0} -m {1} --min-quality 5".format(hashsize, K)
    cmd += " --quality-start {0}".format(offset)
    cmd += " " + " ".join(fastqfiles)

    if need_update(fastqfiles, jfdb):
        sh(cmd)

    jfhisto = jfpf + ".histogram"
    cmd = "jellyfish histo -t 3 {0} -o {1}".format(jfdb, jfhisto)

    if need_update(jfdb, jfhisto):
        sh(cmd)


def meryl(args):
    """
    %prog meryl merylfile

    Run meryl to dump histogram to be used in kmer.histogram(). The merylfile
    are the files ending in .mcidx or .mcdat.
    """
    p = OptionParser(meryl.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    merylfile, = args
    pf, sf = op.splitext(merylfile)
    outfile = pf + ".histogram"
    cmd = "meryl -Dh -s {0}".format(pf)
    sh(cmd, outfile=outfile)

    return outfile


def histogram(args):
    """
    %prog histogram meryl.histogram species K

    Plot the histogram based on meryl K-mer distribution, species and N are
    only used to annotate the graphic. Find out totalKmers when running
    kmer.meryl().
    """
    p = OptionParser(histogram.__doc__)
    p.add_option("--maskone", default=False, action="store_true",
            help="Reduce the error peak [default: %default]")
    p.add_option("--pdf", default=False, action="store_true",
            help="Print PDF instead of ASCII plot [default: %default]")
    p.add_option("--coverage", default=0, type="int",
            help="Kmer coverage [default: auto]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    histfile, species, N = args
    KMERYL, KSOAP, KALLPATHS = range(3)
    kformats = ("Meryl", "Soap", "AllPaths")
    kformat = KMERYL

    ascii = not opts.pdf
    fp = open(histfile)
    hist = {}
    totalKmers = 0

    # Guess the format of the Kmer histogram
    soap = False
    for row in fp:
        if row.startswith("# 1:"):
            kformat = KALLPATHS
            break
        if len(row.split()) == 1:
            kformat = KSOAP
            break
    fp.seek(0)

    logging.debug("Guessed format: {0}".format(kformats[kformat]))

    for rowno, row in enumerate(fp):
        if row[0] == '#':
            continue
        if kformat == KSOAP:
            K = rowno + 1
            counts = int(row.strip())
        else:  # meryl histogram
            K, counts = row.split()[:2]
            K, counts = int(K), int(counts)

        Kcounts = K * counts
        totalKmers += Kcounts
        hist[K] = Kcounts

    if opts.maskone:
        hist[1] = 0

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
        if history == ["drop", "rise", "drop"]:
            break

    Total_Kmers = int(totalKmers)
    coverage = opts.coverage
    Kmer_coverage = Ka if not coverage else coverage
    Genome_size = Total_Kmers * 1. / Kmer_coverage / 1e6

    Total_Kmers_msg = "Total {0}-mers: {1}".format(N, Total_Kmers)
    Kmer_coverage_msg = "{0}-mer coverage: {1}".format(N, Kmer_coverage)
    Genome_size_msg = "Estimated genome size: {0:.1f}Mb".format(Genome_size)

    for msg in (Total_Kmers_msg, Kmer_coverage_msg, Genome_size_msg):
        print >> sys.stderr, msg

    counts = sorted((a, b) for a, b in hist.items() if a <= 100)
    x, y = zip(*counts)
    title = "{0} genome {1}-mer histogram".format(species, N)

    if ascii:
        return asciiplot(x, y, title=title)

    fig = plt.figure(1, (6, 6))
    plt.plot(x, y, 'g-', lw=2, alpha=.5)

    ax = plt.gca()
    ax.text(.5, .9, _(Total_Kmers_msg),
            ha="center", color='b', transform=ax.transAxes)
    ax.text(.5, .8, _(Kmer_coverage_msg),
            ha="center", color='b', transform=ax.transAxes)
    ax.text(.5, .7, _(Genome_size_msg),
            ha="center", color='b', transform=ax.transAxes)

    ax.set_title(_(title), color='r')
    xlabel, ylabel = "Coverage (X)", "Counts"
    ax.set_xlabel(_(xlabel), color='r')
    ax.set_ylabel(_(ylabel), color='r')
    set_human_axis(ax)

    imagename = histfile.split(".")[0] + ".pdf"
    savefig(imagename, dpi=100)


if __name__ == '__main__':
    main()
