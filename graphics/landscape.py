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
from jcvi.algorithms.matrix import moving_sum
from jcvi.graphics.base import plt, _, set_image_options, Rectangle
from jcvi.utils.cbook import human_size
from jcvi.apps.base import ActionDispatcher, debug
debug()


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


def draw_gauge(ax, margin, maxl, optimal=7):
    # Draw a gauge on the top of the canvas
    ax.plot([margin, 1 - margin], [1 - margin, 1 - margin], "k-", lw=2)
    digits = len(str(maxl))
    p = 10 ** (digits - 2)
    two_digits = int(ceil(maxl * 1. / p))
    strides = [(x, ceil(two_digits / x)) for x in (1, 2, 5)]

    best_stride = min(strides, key=lambda x: abs(x[1] - optimal))[0]
    interval = best_stride * p
    #nintervals = int(ceil(maxl * 1. / interval))
    nintervals = maxl / interval
    newl = nintervals * interval

    xx, yy = margin, 1 - margin
    tip = .005
    xinterval = (1 - 2 * margin) / nintervals
    for i in xrange(0, newl + 1, interval):
        l = human_size(i, precision=0)
        if l[-1] == 'b':
            l, suffix = l[:-2], l[-2:]

        ax.plot([xx, xx], [yy, yy + tip], "k-", lw=2)
        ax.text(xx, yy + 2 * tip, _(l), ha="center")
        xx += xinterval

    xx += 4 * tip - xinterval
    ax.text(xx, yy + 2 * tip, _(suffix))
    return interval / xinterval


def stackplot(binfile, chr, t, window, shift, palette="rgbymck"):
    mn = binfile.mapping[chr]
    m, n = zip(*mn)

    m = np.array(m, dtype="float")
    n = np.array(n, dtype="float")

    w = window / shift
    m = moving_sum(m, window=w)
    n = moving_sum(n, window=w)
    m /= n

    m[0] = m[-1] = 0

    return m


def stack(args):
    """
    %prog fastafile bedfile1 bedfile2 bedfile3 ...

    Create landscape plots that show the amounts of genic sequences, and repetitive
    sequences along the chromosomes.
    """
    from jcvi.formats.bed import bins

    p = OptionParser(__doc__)
    p.add_option("--top", default=17, type="int",
                 help="Draw the first N chromosomes [default: %default]")
    p.add_option("--window", default=500000, type="int",
                 help="Size of window [default: %default]")
    p.add_option("--shift", default=100000, type="int",
                 help="Size of shift [default: %default]")
    opts, args, iopts = set_image_options(p, args, figsize="8x8")

    if len(args) < 2:
        sys.exit(not p.print_help())

    fastafile = args[0]
    bedfiles = args[1:]

    top = opts.top
    window = opts.window
    shift = opts.shift
    assert window % shift == 0, "--window must be divisible by --shift"
    binfiles = [bins([x, "--binsize={0}".format(shift), fastafile]) for x in bedfiles]
    binfiles = [BinFile(x) for x in binfiles]

    sizes = Sizes(fastafile)
    s = list(sizes.iter_sizes())[:top]
    maxl = max(x[1] for x in s)
    margin = .06
    inner = .02   # y distance between tracks

    pf = fastafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    max_len = s
    # Gauge
    ratio = draw_gauge(root, margin, maxl)

    # Per chromosome
    yinterval = (1 - 2 * margin) / (top + 1)
    palette = "rgbymck"
    xx = margin
    yy = 1 - margin
    for chr, clen in s:
        yy -= yinterval
        xlen = clen / ratio
        if "_" in chr:
            ca, cb = chr.split("_")
            cc = ca[0] + cb

        root.add_patch(Rectangle((xx, yy), xlen, yinterval - inner,
                     color="pink", alpha=.5))
        ax = fig.add_axes([xx, yy, xlen, yinterval - inner])

        ax.set_axis_off()
        root.text(xx, yy + yinterval - inner, _(cc), va="top")

        nbins = clen / shift
        if clen % shift:
            nbins += 1

        t = np.arange(nbins)
        m = np.zeros(nbins, dtype="float")
        zorders = range(10)[::-1]
        for binfile, p, z in zip(binfiles, palette, zorders):
            s = stackplot(binfile, chr, t, window, shift)
            m += s
            ax.fill(t, m, p, lw=0, zorder=z)

        ax.set_xlim(0, len(t))
        ax.set_ylim(0, 1)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    logging.debug("Print image to `{0}` {1}".format(image_name, iopts))
    plt.savefig(image_name, dpi=iopts.dpi)
    plt.rcdefaults()


if __name__ == '__main__':
    main()
