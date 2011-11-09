#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, Rectangle, _
from jcvi.utils.cbook import thousands
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('rnaseq', 'evaluate completeness, contiguity, and chimer of RNAseq'),
        ('plot', 'plot the alignment of the scaffold to other evidences'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rnaseq(args):
    """
    %prog rnaseq blastfile ref.fasta

    Evaluate de-novo RNA-seq assembly against a reference gene set (same or
    closely related organism). Ideally blatfile needs to be supermap'd.

    Following metric is used (Martin et al. 2010, Rnnotator paper):
    Accuracy: % of contigs share >=95% identity with ref genome (TODO)
    Completeness: % of ref genes covered by contigs to >=80% of their lengths
    Contiguity: % of ref genes covered by a *single* contig >=80% of lengths
    Chimer: % of contigs that contain two or more annotated genes >= 50bp
    """
    from jcvi.algorithms.supermap import supermap

    p = OptionParser(rnaseq.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    blastfile, reffasta = args
    sizes = Sizes(reffasta).mapping
    known_genes = len(sizes)

    querysupermap = blastfile + ".query.supermap"
    refsupermap = blastfile + ".ref.supermap"

    if not op.exists(querysupermap):
        supermap(blastfile, filter="query")
    if not op.exists(refsupermap):
        supermap(blastfile, filter="ref")

    blast = Blast(querysupermap)
    chimers = 0
    goodctg80 = set()
    goodctg50 = set()
    for ctg, hits in blast.iter_hits():
        bps = defaultdict(int)
        for x in hits:
            bps[x.subject] += abs(x.sstop - x.sstart) + 1

        valid_hits = bps.items()
        for vh, length in valid_hits:
            rsize = sizes[vh]
            ratio = length * 100. / rsize
            if ratio >= 80:
                goodctg80.add(ctg)
            if ratio >= 50:
                goodctg50.add(ctg)

        # Chimer
        if len(valid_hits) > 1:
            chimers += 1

    blast = Blast(refsupermap)
    goodref80 = set()
    goodref50 = set()
    bps = defaultdict(int)
    for x in blast.iter_line():
        bps[x.subject] += abs(x.sstop - x.sstart) + 1

    for vh, length in bps.items():
        rsize = sizes[vh]
        ratio = length * 100. / rsize
        if ratio >= 80:
            goodref80.add(vh)
        if ratio >= 50:
            goodref50.add(vh)

    print >> sys.stderr, "Reference set: `{0}`,  # of transcripts {1}".\
            format(reffasta, known_genes)
    print >> sys.stderr, "A total of {0} contigs map to 80% of a reference"\
            " transcript".format(len(goodctg80))
    print >> sys.stderr, "A total of {0} contigs map to 50% of a reference"\
            " transcript".format(len(goodctg50))
    print >> sys.stderr, "A total of {0} reference transcripts ({1:.1f}%) have 80% covered" \
            .format(len(goodref80), len(goodref80) * 100. / known_genes)


def scaffolding(ax, scaffoldID, blastf, qsizes, ssizes, qbed, sbed):

    from jcvi.graphics.blastplot import blastplot

    # qsizes, qbed are properties for the evidences
    # ssizes, sbed are properties for the current scaffoldID
    blastplot(ax, blastf, qsizes, ssizes, qbed, sbed, \
              style="circle", insetLabels=True, stripNames=True)

    # FPC_scf.bed => FPC
    fname = qbed.filename.split(".")[0].split("_")[0]
    xtitle = fname
    if xtitle == "FPC":
        ax.set_xticklabels([""] * len(ax.get_xticklabels()))
    ax.set_xlabel(_(xtitle), color="g")
    for x in ax.get_xticklines():
        x.set_visible(False)


def plot_one_scaffold(scaffoldID, ssizes, sbed, trios, imagename, iopts):
    ntrios = len(trios)
    fig = plt.figure(1, (14, 8))
    plt.cla()
    plt.clf()
    root = fig.add_axes([0, 0, 1, 1])
    axes = [fig.add_subplot(1, ntrios, x) for x in range(1, ntrios + 1)]
    scafsize = ssizes.get_size(scaffoldID)

    for trio, ax in zip(trios, axes):
        blastf, qsizes, qbed = trio
        scaffolding(ax, scaffoldID, blastf, qsizes, ssizes, qbed, sbed)

    root.text(.5, .95, _("{0}   (size={1})".\
            format(scaffoldID, thousands(scafsize))),
            size=18, ha="center", color='b')
    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    plt.savefig(imagename, dpi=iopts.dpi)
    logging.debug("Print image to `{0}` {1}".format(imagename, iopts))


def plot(args):
    """
    %prog plot scaffold.fasta synteny.blast synteny.sizes synteny.bed
                         physicalmap.blast physicalmap.sizes physicalmap.bed

    As evaluation of scaffolding, visualize external line of evidences:
    * Plot synteny to an external genome
    * Plot alignments to physical map
    * Plot alignments to genetic map (TODO)

    Each trio defines one panel to be plotted. blastfile defines the matchings
    between the evidences vs scaffolds. Then the evidence sizes, and evidence
    bed to plot dot plots.

    This script will plot a dot in the dot plot in the corresponding location
    the plots are one contig/scaffold per plot.
    """
    from jcvi.graphics.base import set_image_options
    from jcvi.utils.iter import grouper

    p = OptionParser(plot.__doc__)
    p.add_option("--cutoff", type="int", default=1000000,
            help="Plot scaffolds with size larger than [default: %default]")
    opts, args, iopts = set_image_options(p, args, figsize="14x8", dpi=150)

    if len(args) < 4 or len(args) % 3 != 1:
        sys.exit(not p.print_help())

    scafsizes = Sizes(args[0])
    trios = list(grouper(3, args[1:]))
    trios = [(a, Sizes(b), Bed(c)) for a, b, c in trios]

    for scaffoldID, scafsize in scafsizes.iter_sizes():
        if scafsize < opts.cutoff:
            continue
        logging.debug("Loading {0} (size={1})".format(scaffoldID,
            thousands(scafsize)))

        tmpname = scaffoldID + ".sizes"
        tmp = open(tmpname, "w")
        tmp.write("{0}\t{1}".format(scaffoldID, scafsize))
        tmp.close()

        tmpsizes = Sizes(tmpname)
        tmpsizes.close(clean=True)

        imagename = ".".join((scaffoldID, opts.format))
        plot_one_scaffold(scaffoldID, tmpsizes, None, trios, imagename, iopts)


if __name__ == '__main__':
    main()
