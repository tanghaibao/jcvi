#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the pineapple genome manuscript (unpublished).
"""

import sys

from jcvi.graphics.base import plt, savefig, mpl, normalize_axes, panel_labels
from jcvi.graphics.glyph import TextCircle
from jcvi.graphics.karyotype import Karyotype
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('ploidy', 'plot pineapple macro-synteny (requires data)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ploidy(args):
    """
    %prog ploidy seqids layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of WGD history of pineapple genome.
    """
    p = OptionParser(ploidy.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x6")

    if len(args) != 2:
        sys.exit(not p.print_help())

    seqidsfile, klayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "pineapple-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
