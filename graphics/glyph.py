#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Gradient gene features
"""

import sys
import logging

from optparse import OptionParser

import numpy as np
from jcvi.apps.base import ActionDispatcher, debug
from jcvi.graphics.base import plt, Rectangle, CirclePolygon, _
debug()

tstep = .05
Timing = np.arange(0, 1 + tstep, tstep)

class Bezier (object):
    """
    Cubic bezier curve, see the math:
    <http://www.moshplant.com/direct-or/bezier/math.html>
    p0 : origin, p1, p2 :control, p3: destination
    """
    def __init__(self, ax, p0, p1, p2, p3, color='m', alpha=.2):
        pts = (p0, p1, p2, p3)
        px, py = zip(*pts)
        xt = self.get_array(px)
        yt = self.get_array(py)

        ax.plot(xt, yt, "-", color=color, alpha=alpha)

    def get_array(self, pts, t=Timing):
        p0, p1, p2, p3 = pts

        # Get the coeffiencients
        c = 3*(p1 - p0)
        b = 3*(p2 - p1) -c
        a = p3 - p0 - c - b

        tsquared = t ** 2
        tcubic = tsquared * t
        return a * tcubic + b * tsquared + c * t + p0


class RoundLabel (object):
    """Round rectangle around the text label
    """
    def __init__(self, ax, x1, x2, t, **kwargs):

        ax.text(x1, x2, _(t), ha="center",
            bbox=dict(boxstyle="round",fill=False))


class DCircle:
    """Circle with a double-line margin
    """
    def __init__(self, ax, x, y, radius=.01, **kwargs):

      ax.add_patch(CirclePolygon((x, y), radius * 1.4,
          resolution=50, fc="w", ec="k"))
      ax.add_patch(CirclePolygon((x, y), radius,
          resolution=50, **kwargs))


class Glyph (object):
    """Draws gradient rectangle
    """
    def __init__(self, ax, x1, x2, y, height=.03, 
            fc="gray", **kwargs):

        width = x2 - x1
        # Frame around the gradient rectangle
        p1 = (x1, y - .5 * height)
        ax.add_patch(Rectangle(p1, width, height, fc=fc, 
            lw=0, **kwargs))
        # Several overlaying patches
        for cascade in np.arange(.1, .55, .05):
            p1 = (x1, y - height * cascade)
            ax.add_patch(Rectangle(p1, width, 2 * cascade * height,
                fc='w', lw=0, alpha=.1))


def main():
    """
    %prog demo

    Draw sample gene features to illustrate the various fates of duplicate
    genes - to be used in a book chapter.
    """
    p = OptionParser(main.__doc__)
    opts, args = p.parse_args()

    fig = plt.figure(1, (8, 5))
    root = fig.add_axes([0, 0, 1, 1])

    panel_space = .23
    dup_space = .025
    # Draw a gene and two regulatory elements at these arbitrary locations
    locs = [(.5, .9), # ancestral gene
            (.5, .9 - panel_space + dup_space), # identical copies
            (.5, .9 - panel_space - dup_space),
            (.5, .9 - 2 * panel_space + dup_space), # degenerate copies
            (.5, .9 - 2 * panel_space - dup_space), 
            (.2, .9 - 3 * panel_space + dup_space), # sub-functionalization
            (.2, .9 - 3 * panel_space - dup_space), 
            (.5, .9 - 3 * panel_space + dup_space), # neo-functionalization
            (.5, .9 - 3 * panel_space - dup_space), 
            (.8, .9 - 3 * panel_space + dup_space), # non-functionalization
            (.8, .9 - 3 * panel_space - dup_space), 
            ]

    default_regulator = "gm"
    regulators = [default_regulator,
            default_regulator, default_regulator,
            "wm", default_regulator,
            "wm", "gw",
            "wb", default_regulator,
            "ww", default_regulator,
            ]
    
    width = .24
    for i, (xx, yy) in enumerate(locs):
        regulator = regulators[i]
        x1, x2 = xx - .5 * width, xx + .5 * width
        Glyph(root, x1, x2, yy)
        if i == 9:  # upper copy for non-functionalization
            continue

        # coding region
        x1, x2 = xx - .16 * width, xx + .45 * width
        Glyph(root, x1, x2, yy, fc="k")

        # two regulatory elements
        x1, x2 = xx - .4 * width, xx - .28 * width
        for xx, fc in zip((x1, x2), regulator):
            if fc == 'w':
                continue

            DCircle(root, xx, yy, fc=fc)

        rotation = 30
        tip = .02
        if i == 0:
            ya = yy + tip
            root.text(x1, ya, _("Flower"), rotation=rotation, va="bottom")
            root.text(x2, ya, _("Root"), rotation=rotation, va="bottom")
        elif i == 7:
            ya = yy + tip
            root.text(x2, ya, _("Leaf"), rotation=rotation, va="bottom")

    # Draw arrows between panels (center)
    arrow_dist = .08
    ar_xpos = .5
    arrowprops = dict(arrowstyle="fancy", fc="k", alpha=.5, 
            connectionstyle="arc3,rad=-0.05")
    for ar_ypos in (.3, .53, .76):
        root.annotate(" ", (ar_xpos, ar_ypos), 
                (ar_xpos, ar_ypos + arrow_dist),
                arrowprops=arrowprops)

    ar_ypos = .3
    for ar_xpos in (.2, .8):
        root.annotate(" ", (ar_xpos, ar_ypos),
                (.5, ar_ypos + arrow_dist), 
                arrowprops=arrowprops)

    # Duplication, Degeneration
    xx = .6
    ys = (.76, .53)
    processes = ("Duplication", "Degeneration")
    for yy, process in zip(ys, processes):
        root.text(xx, yy + .02, process, fontweight="bold")

    # Label of fates
    xs = (.2, .5, .8)
    fates = ("Subfunctionalization", "Neofunctionalization",
            "Nonfunctionalization")
    yy = .05
    for xx, fate in zip(xs, fates):
        RoundLabel(root, xx, yy, fate)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = "demo.pdf"
    plt.savefig(figname, dpi=300)
    logging.debug("Figure saved to `{0}`".format(figname))


if __name__ == '__main__':
    main()
