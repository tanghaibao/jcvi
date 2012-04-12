#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Create chromosome landscape plots that are similar to the ones used in soybean
and sorghum paper.
"""


import sys
import logging

import numpy as np

from math import ceil
from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.sizes import Sizes
from jcvi.formats.base import BaseFile
from jcvi.formats.bed import Bed, bins
from jcvi.algorithms.matrix import moving_sum
from jcvi.graphics.base import plt, _, set_image_options, \
        Rectangle, CirclePolygon
from jcvi.utils.cbook import human_size, autoscale
from jcvi.apps.base import ActionDispatcher, debug
debug()


# Colors picked from Schmutz soybean genome paper using ColorPic
palette = ["#ACABD5","#DBF0F5", "#3EA77A", "#FBF5AB", "#C162A6"] + \
          list("rgbymck")
gray = "#CCCCCB"
Registration = {"Gypsy": "LTR-RT/Gypsy",
                "Copia": "LTR-RT/Copia",
                "hAT": "DNA-TE/hAT",
                "Helitron": "DNA-TE/Helitron",
                "Tourist": "DNA-TE/Tourist",
                "Introns": "Genes (introns)",
                "Exons": "Genes (exons)"}


class BinFile (BaseFile):

    def __init__(self, filename):
        super(BinFile, self).__init__(filename)
        self.mapping = defaultdict(list)

        fp = open(filename)
        for row in fp:
            chr, len, binlen = row.split()
            len, binlen = int(len), int(binlen)
            self.mapping[chr].append((len, binlen))
        fp.close()


def main():

    actions = (
        ('stack', 'create landscape plote with genic/te composition'),
        ('heatmap', 'similar to stack but adding heatmap'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add_window_options(p):
    p.add_option("--window", default=500000, type="int",
                 help="Size of window [default: %default]")
    p.add_option("--shift", default=100000, type="int",
                 help="Size of shift [default: %default]")


def check_window_options(opts):
    window = opts.window
    shift = opts.shift
    assert window % shift == 0, "--window must be divisible by --shift"

    return window, shift


def heatmap(args):
    """
    %prog heatmap fastafile chr1

    Combine stack plot with heatmap to show abundance of various tracks along
    given chromosome. Need to give multiple beds to --stacks and --heatmaps
    """
    p = OptionParser(heatmap.__doc__)
    p.add_option("--stacks",
                 default="Exons,Introns,DNA_transposons,Retrotransposons",
                 help="Features to plot in stackplot [default: %default]")
    p.add_option("--heatmaps",
                 default="Copia,Gypsy,hAT,Helitron,Introns,Exons",
                 help="Features to plot in heatmaps [default: %default]")
    p.add_option("--meres", default=None,
                 help="Extra centromere / telomere features [default: %default]")
    add_window_options(p)
    opts, args, iopts = set_image_options(p, args, figsize="8x5")

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, chr = args
    window, shift = check_window_options(opts)

    stacks = opts.stacks.split(",")
    heatmaps = opts.heatmaps.split(",")
    stackbeds = [x + ".bed" for x in stacks]
    heatmapbeds = [x + ".bed" for x in heatmaps]
    stackbins = get_binfiles(stackbeds, fastafile, shift)
    heatmapbins = get_binfiles(heatmapbeds, fastafile, shift)

    window, shift = check_window_options(opts)
    margin = .06
    inner = .015
    clen = Sizes(fastafile).mapping[chr]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Gauge
    ratio = draw_gauge(root, margin, clen, rightmargin=4 * margin)
    yinterval = .3
    xx = margin
    yy = 1 - margin
    yy -= yinterval
    xlen = clen / ratio
    if "_" in chr:
        ca, cb = chr.split("_")
        cc = ca[0] + cb

    root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner, color=gray))
    ax = fig.add_axes([xx, yy, xlen, yinterval - inner])

    nbins = clen / shift
    if clen % shift:
        nbins += 1

    owindow = clen / 100
    if owindow > window:
        window = owindow / shift * shift

    stackplot(ax, stackbins, nbins, palette, chr, window, shift)
    root.text(xx + inner, yy + yinterval - 2 * inner, cc, va="top")

    # Legends
    xx += xlen + .01
    yspace = (yinterval - inner) / (len(stackbins) + 1)
    yy = 1 - margin - yinterval
    for s, p in zip(stacks, palette):
        s = s.replace("_", " ")
        s = Registration.get(s, s)

        yy += yspace
        root.add_patch(Rectangle((xx, yy), inner, inner, color=p, lw=0))
        root.text(xx + 1.5 * inner, yy, s, size=10)

    yh = .05  # Heatmap height
    # Heatmaps
    xx = margin
    yy = 1 - margin - yinterval - inner
    for s, p in zip(heatmaps, heatmapbins):
        s = s.replace("_", " ")
        s = Registration.get(s, s)

        yy -= yh
        m = stackarray(p, chr, window, shift)

        Y = np.array([m, m])
        root.imshow(Y, extent=(xx, xx + xlen, yy, yy + yh - inner),
                    interpolation="nearest", aspect="auto")
        root.text(xx + xlen + .01, yy, s, size=10)

    yy -= yh

    meres = opts.meres
    if meres:
        bed = Bed(meres)
        for b in bed:
            if b.seqid != chr:
                continue
            pos = (b.start + b.end) / 2
            cpos = pos / ratio
            xx = margin + cpos
            accn = b.accn.capitalize()
            root.add_patch(CirclePolygon((xx, yy), radius=.01, fc="m", ec="m"))
            root.text(xx + .014, yy, _(accn), va="center", color="m")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = chr + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


def draw_gauge(ax, margin, maxl, rightmargin=None, optimal=7):
    # Draw a gauge on the top of the canvas
    rightmargin = rightmargin or margin
    ax.plot([margin, 1 - rightmargin], [1 - margin, 1 - margin], "k-", lw=2)

    best_stride = autoscale(maxl)
    nintervals = int(round(maxl * 1. / best_stride))
    newl = nintervals * best_stride

    xx, yy = margin, 1 - margin
    tip = .005
    xinterval = (1 - margin - rightmargin) / nintervals
    l = human_size(best_stride)
    if l[-1] == 'b':
        suffix = target = l[-2:]

    for i in xrange(0, newl + 1, best_stride):
        l = human_size(i, precision=0, target=target)
        if l[-1] == 'b':
            l, suffix = l[:-2], l[-2:]

        ax.plot([xx, xx], [yy, yy + tip], "k-", lw=2)
        ax.text(xx, yy + 2 * tip, _(l), ha="center")
        xx += xinterval

    xx += 4 * tip - xinterval
    ax.text(xx, yy + 2 * tip, _(suffix))

    return best_stride / xinterval


def get_binfiles(bedfiles, fastafile, shift):
    binfiles = [bins([x, "--binsize={0}".format(shift), fastafile]) for x in bedfiles]
    binfiles = [BinFile(x) for x in binfiles]
    return binfiles


def stackarray(binfile, chr, window, shift):
    mn = binfile.mapping[chr]
    m, n = zip(*mn)

    m = np.array(m, dtype="float")
    n = np.array(n, dtype="float")

    w = window / shift
    m = moving_sum(m, window=w)
    n = moving_sum(n, window=w)
    m /= n

    return m


def stackplot(ax, binfiles, nbins, palette, chr, window, shift):
    t = np.arange(nbins, dtype="float") + .5
    m = np.zeros(nbins, dtype="float")
    zorders = range(10)[::-1]
    for binfile, p, z in zip(binfiles, palette, zorders):
        s = stackarray(binfile, chr, window, shift)
        m += s
        ax.fill_between(t, m, color=p, lw=0, zorder=z)

    ax.set_xlim(0, nbins)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


def stack(args):
    """
    %prog stack fastafile

    Create landscape plots that show the amounts of genic sequences, and repetitive
    sequences along the chromosomes.
    """
    p = OptionParser(stack.__doc__)
    p.add_option("--top", default=17, type="int",
                 help="Draw the first N chromosomes [default: %default]")
    p.add_option("--stacks",
                 default="Exons,Introns,DNA_transposons,Retrotransposons",
                 help="Features to plot in stackplot [default: %default]")
    add_window_options(p)
    opts, args, iopts = set_image_options(p, args, figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    top = opts.top
    window, shift = check_window_options(opts)

    bedfiles = [x + ".bed" for x in opts.stacks.split(",")]
    binfiles = get_binfiles(bedfiles, fastafile, shift)

    sizes = Sizes(fastafile)
    s = list(sizes.iter_sizes())[:top]
    maxl = max(x[1] for x in s)
    margin = .06
    inner = .015   # y distance between tracks

    pf = fastafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    max_len = s
    # Gauge
    ratio = draw_gauge(root, margin, maxl)

    # Per chromosome
    yinterval = (1 - 2 * margin) / (top + 1)
    xx = margin
    yy = 1 - margin
    for chr, clen in s:
        yy -= yinterval
        xlen = clen / ratio
        if "_" in chr:
            ca, cb = chr.split("_")
            cc = ca[0] + cb

        root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner, color=gray))
        ax = fig.add_axes([xx, yy, xlen, yinterval - inner])

        nbins = clen / shift
        if clen % shift:
            nbins += 1

        stackplot(ax, binfiles, nbins, palette, chr, window, shift)
        root.text(xx - .03, yy + .5 * (yinterval - inner), cc, ha="center", va="center")

        ax.set_xlim(0, nbins)
        ax.set_ylim(0, 1)
        ax.set_axis_off()

    # Legends
    yy -= yinterval
    xx = margin
    for b, p in zip(bedfiles, palette):
        b = b.rsplit(".", 1)[0].replace("_", " ")
        b = Registration.get(b, b)

        root.add_patch(Rectangle((xx, yy), inner, inner, color=p, lw=0))
        xx += 2 * inner
        root.text(xx, yy, _(b))
        xx += len(b) * .01 + inner

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
