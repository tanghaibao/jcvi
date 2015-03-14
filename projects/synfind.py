#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""

import sys

from jcvi.graphics.base import plt, savefig, panel_labels
from jcvi.graphics.glyph import CartoonRegion
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('cartoon', 'generate cartoon illustration of SynFind'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def cartoon(args):
    """
    %prog synteny.py

    Generate cartoon illustration of SynFind.
    """
    p = OptionParser(cartoon.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    CartoonRegion(root, .1, .9, .5, 41)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
