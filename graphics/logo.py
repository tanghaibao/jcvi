#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog text

Generate logo with certain color and font. Font must be present in:
jcvi/graphics/fonts

You can download a number of free TTF fonts from:
<http://www.urbanfonts.com>
"""


import sys

from jcvi.graphics.base import plt, savefig, fontprop, available_fonts
from jcvi.apps.base import OptionParser


def main():
    p = OptionParser(__doc__)
    p.add_option("--customfont", default="Airswing.ttf", choices=available_fonts,
                 help="Custom font name")
    p.add_option("--color", default="limegreen",
                 help="Font color [default: %default]")
    p.add_option("--size", default=36, type="int",
                 help="Font size [default: %default]")
    opts, args, iopts = p.set_image_options(figsize='2x1', dpi=60, format='png')

    if len(args) != 1:
        sys.exit(not p.print_help())

    text, = args

    plt.rcdefaults()
    fig = plt.figure(1, (iopts.w, iopts.h))
    ax = fig.add_axes([0, 0, 1, 1])

    ax.text(.5, .5, text, color=opts.color, ha="center", va="center")
    fontprop(ax, opts.customfont, size=opts.size)

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()

    image_name = text + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
