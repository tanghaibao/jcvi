#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the ALLMAPS manuscript (un-published)
"""

import sys

from jcvi.graphics.base import plt, savefig, normalize_axes, panel_labels
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('lms', 'ALLMAPS cartoon to illustrate LMS metric'),
        ('allmapsQC', 'plot ALLMAPS accuracy across a range of simulated data'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def normalize_lms_axis(ax):
    ax.set_xlim(0, 110)
    ax.set_ylim(0, 110)
    yticklabels = [int(x) for x in ax.get_yticks()]
    ax.set_yticklabels(yticklabels, family='Helvetica')
    ax.set_xticks([])
    ax.set_ylabel("Map (cM)")


def lms(args):
    """
    %prog lms

    ALLMAPS cartoon to illustrate LMS metric.
    """
    from random import randint
    from jcvi.graphics.chromosome import HorizontalChromosome

    p = OptionParser(lms.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="6x6", dpi=300)

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Panel A
    w, h = .7, .35
    ax = fig.add_axes([.15, .6, w, h])

    xdata = [x + randint(-3, 3) for x in range(10, 110, 10)]
    ydata = [x + randint(-3, 3) for x in range(10, 110, 10)]
    ydata[3:7] = ydata[3:7][::-1]
    xydata = zip(xdata, ydata)
    lis = xydata[:3] + [xydata[4]] + xydata[7:]
    lds = xydata[3:7]
    xlis, ylis = zip(*lis)
    xlds, ylds = zip(*lds)
    ax.plot(xlis, ylis, "r-", lw=12, alpha=.3,
                              solid_capstyle="round", solid_joinstyle="round")
    ax.plot(xlds, ylds, "g-", lw=12, alpha=.3,
                              solid_capstyle="round", solid_joinstyle="round")
    ax.plot(xdata, ydata, "k.", mec="k", mfc="w", mew=3, ms=12)
    HorizontalChromosome(root, .15, .15 + w, .57, height=.02, lw=2)
    root.text(.15 + w / 2, .55, "Chromosome location (bp)", ha="center", va="top")

    ax.text(80, 30, "LIS = 7", color="r", ha="center", va="center")
    ax.text(80, 20, "LDS = 4", color="g", ha="center", va="center")
    ax.text(80, 10, "LMS = $max$(LIS, LDS) = 7", ha="center", va="center")
    normalize_lms_axis(ax)

    # Panel B
    w = .37
    p = (0, 45, 75, 110)
    ax = fig.add_axes([.1, .12, w, h])
    xdata = [x for x in range(10, 110, 10)]
    ydata = ydata_orig = [x for x in range(10, 110, 10)]
    ydata = ydata[:4] + ydata[7:] + ydata[4:7][::-1]
    xydata = zip(xdata, ydata)
    lis = xydata[:7]
    xlis, ylis = zip(*lis)
    ax.plot(xlis, ylis, "r-", lw=12, alpha=.3,
                              solid_capstyle="round", solid_joinstyle="round")
    ax.plot(xdata, ydata, "k.", mec="k", mfc="w", mew=3, ms=12)
    ax.vlines(p, 0, 110, colors="beige", lw=3)
    normalize_lms_axis(ax)
    patch = [.1 + w * x / 110. for x in p]
    HorizontalChromosome(root, .1, .1 + w, .09, patch=patch,
                         height=.02, lw=2)
    scaffolds = ("a", "b", "c")
    for i, s in enumerate(scaffolds):
        xx = (patch[i] + patch[i + 1]) / 2
        root.text(xx, .09, s, va="center", ha="center")
    root.text(.1 + w / 2, .04, "LMS($a||b||c$) = 7", ha="center")

    # Panel C
    ax = fig.add_axes([.6, .12, w, h])
    patch = [.6 + w * x / 110. for x in p]
    ydata = ydata_orig
    ax.plot(xdata, ydata, "r-", lw=12, alpha=.3,
                              solid_capstyle="round", solid_joinstyle="round")
    ax.plot(xdata, ydata, "k.", mec="k", mfc="w", mew=3, ms=12)
    ax.vlines(p, [0], [110], colors="beige", lw=3)
    normalize_lms_axis(ax)
    HorizontalChromosome(root, .6, .6 + w, .09, patch=patch,
                         height=.02, lw=2)
    scaffolds = ("a", "-c", "b")
    for i, s in enumerate(scaffolds):
        xx = (patch[i] + patch[i + 1]) / 2
        root.text(xx, .09, s, va="center", ha="center")
    root.text(.6 + w / 2, .04, "LMS($a||-c||b$) = 10", ha="center")

    labels = ((.05, .95, 'A'), (.05, .48, 'B'), (.55, .48, 'C'))
    panel_labels(root, labels)

    normalize_axes(root)

    pf = "lms"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def import_data(datafile):
    data = []
    fp = open(datafile)
    fp.readline()
    for row in fp:
        x, y = row.split()[:2]
        x, y = float(x), float(y)
        data.append((x, y))
    return data


def subplot(ax, data, xlabel, ylabel, xlim=None, ylim=1.1,
                      xcast=float, ycast=float):
    x, y = zip(*data)
    ax.plot(x, y, "ko:", mec="k", mfc="w", ms=4)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    if xlim:
        ax.set_xlim(0, xlim)
    if ylim:
        ax.set_ylim(0, ylim)
    xticklabels = [xcast(x) for x in ax.get_xticks()]
    yticklabels = [ycast(x) for x in ax.get_yticks()]
    ax.set_xticklabels(xticklabels, family='Helvetica')
    ax.set_yticklabels(yticklabels, family='Helvetica')


def allmapsQC(args):
    """
    %prog allmapsQC inversion.txt translocation.txt maps.txt

    Plot ALLMAPS accuracy across a range of simulated datasets.
    """
    p = OptionParser(allmapsQC.__doc__)
    opts, args, iopts = p.set_image_options(args, dpi=300)

    if len(args) != 4:
        sys.exit(not p.print_help())

    dataA, dataB, dataC, dataD = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    A = fig.add_axes([.12, .62, .35, .35])
    B = fig.add_axes([.62, .62, .35, .35])
    C = fig.add_axes([.12, .12, .35, .35])
    D = fig.add_axes([.62, .12, .35, .35])
    dataA = import_data(dataA)
    dataB = import_data(dataB)
    dataC = import_data(dataC)
    dataD = import_data(dataD)
    subplot(A, dataA, "Inversion error rate", "Accuracy", xlim=.5)
    subplot(B, dataB, "Translocation error rate", "Accuracy", xlim=.5)
    subplot(C, dataC, "Number of input maps", "Accuracy", xcast=int)
    subplot(D, dataD, "Number of input maps", "Accuracy", xcast=int)

    labels = ((.03, .97, "A"), (.53, .97, "B"),
              (.03, .47, "C"), (.53, .47, "D"))
    panel_labels(root, labels)

    normalize_axes(root)
    image_name = "simulation." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
