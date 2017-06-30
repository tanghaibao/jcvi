#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wheel plot that shows continuous data in radial axes.
"""

import numpy as np
import sys


from jcvi.graphics.base import plt, savefig, normalize_axes, set2
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('wheel', 'wheel plot that shows continuous data in radial axes'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def closed_plot(ax, theta, r, *args, **kwargs):
    theta = list(theta) + [theta[0]]
    r = list(r) + [r[0]]
    ax.plot(theta, r, *args, **kwargs)


def sector(ax, theta_min, theta_max, theta_pad, r, *args, **kwargs):
    theta = np.linspace(theta_min - theta_pad, theta_max + theta_pad, num=100)
    r = len(theta) * [r]
    theta = list(theta) + [0]
    r = list(r) + [0]
    closed_plot(ax, theta, r, *args, **kwargs)


def wheel(args):
    """
    %prog wheel datafile

    Wheel plot that shows continous data in radial axes.
    """
    p = OptionParser(wheel.__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    categories = 20
    ax = plt.subplot(111, projection='polar')
    R = 2

    # Baseline
    r = [R / 2] * categories
    theta = np.linspace(0, 2 * np.pi, endpoint=False, num=categories)
    closed_plot(ax, theta, r, "k:", lw=2)

    # Grid
    for t in theta:
        ax.plot([t, t], [0, R], color="gainsboro")
    ax.axis('off')

    # Sectors (groupings)
    groups = {1: [0, 1, 2, 3],
              2: [4, 5, 6, 7, 8],
              3: [9, 10, 11],
              4: [12, 13, 14, 15, 16],
              5: [17],
              6: [18, 19]}

    theta_interval = 2 * np.pi / categories
    theta_pad = theta_interval / 2 * .9
    for color, group in zip(set2, groups.values()):
        tmin, tmax = min(group), max(group)
        sector(ax, theta[tmin], theta[tmax], theta_pad, R * .95,
                   "-", color=color)

    # Data
    r = np.random.uniform(0, R, size=categories)
    closed_plot(ax, theta, r, color="lightslategray")
    ax.set_rmax(R)

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
