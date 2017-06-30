#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wheel plot that shows continuous data in radial axes.
"""

import numpy as np
import sys

from collections import OrderedDict
from itertools import groupby

from jcvi.graphics.base import plt, savefig, normalize_axes, set2
from jcvi.apps.base import OptionParser, ActionDispatcher


R = 30

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
    r = list(r) + [-R]
    closed_plot(ax, theta, r, *args, **kwargs)


def parse_data(datafile):
    data = {}
    fp = open(datafile)
    for row in fp:
        label, score, percentile = row.split(",")
        label = label.strip()
        score = float(score.strip())
        data[label] = score
    return data


def parse_groups(groupsfile):
    groups = OrderedDict()
    fp = open(groupsfile)
    for row in fp:
        group, label = row.split(",")
        group = group.strip()
        label = label.strip()
        groups[label] = group
    return groups


def wheel(args):
    """
    %prog wheel datafile.csv groups.csv

    Wheel plot that shows continous data in radial axes.
    """
    p = OptionParser(wheel.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 2:
        sys.exit(not p.print_help())

    datafile, groupsfile = args
    df = parse_data(datafile)
    groups = parse_groups(groupsfile)
    labels = [g for g in groups if g in df]
    print labels
    df = [df[g] for g in labels]
    print df
    groups = [groups[g] for g in labels]
    print groups

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    categories = len(df)
    ax = plt.subplot(111, projection='polar')

    brewer = ["#a6cee3", "#1f78b4", "#b2df8a", "#33a02c",
              "#fb9a99", "#e31a1c", "#fdbf6f", "#ff7f00",
              "#cab2d6", "#6a3d9a", "#ffff99", "#b15928"]
    # Baseline
    theta = np.linspace(0, 2 * np.pi, endpoint=False, num=categories)
    _theta = np.linspace(0, 2 * np.pi)
    _r = len(_theta) * [0]
    closed_plot(ax, _theta, _r, "k:", lw=2)

    # Grid
    for t in theta:
        ax.plot([t, t], [-R, R], color="gainsboro")
    ax.axis('off')

    # Sectors (groupings)
    gg = {}
    for group, c in groupby(enumerate(groups), lambda x: x[1]):
        c = [x[0] for x in list(c)]
        gg[group] = c
    print gg

    theta_interval = 2 * np.pi / categories
    theta_pad = theta_interval / 2 * .9
    for color, group in zip(brewer, gg.values()):
        tmin, tmax = min(group), max(group)
        sector(ax, theta[tmin], theta[tmax], theta_pad, R * .95,
                   "-", color=color)

    # Data
    r = df
    closed_plot(ax, theta, r, color="lightslategray")
    for color, group in zip(brewer, gg.values()):
        color_theta = [theta[x] for x in group]
        color_r = [r[x] for x in group]
        ax.plot(color_theta, color_r, "o")

    ax.set_rmin(-R)
    ax.set_rmax(R)

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
