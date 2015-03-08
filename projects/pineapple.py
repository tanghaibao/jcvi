#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the pineapple genome manuscript (unpublished).
"""

import sys

from jcvi.graphics.base import plt, savefig, panel_labels
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
from jcvi.graphics.glyph import TextCircle
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('ploidy', 'plot pineapple macro-synteny (requires data)'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ploidy(args):
    """
    %prog cotton seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of WGD history of pineapple genome. The script calls both graphics.karyotype
    and graphic.synteny.
    """
    p = OptionParser(ploidy.__doc__)
    p.add_option("--switch", help="Rename the seqid with two-column file")
    opts, args, iopts = p.set_image_options(args, figsize="9x7")

    if len(args) != 5:
        sys.exit(not p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)
    Synteny(fig, root, datafile, bedfile, slayout, switch=opts.switch)

    # legend showing the orientation of the genes
    draw_gene_legend(root, .27, .37, .52)

    # annotate the WGD events
    fc = 'lightslategrey'
    x = .09
    radius = .012
    TextCircle(root, x, .825, r'$\tau$', radius=radius)
    TextCircle(root, x, .8, r'$\sigma$', radius=radius)
    TextCircle(root, x, .72, r'$\rho$', radius=radius)
    for ypos in (.825, .8, .72):
        root.text(.12, ypos, r"$\times2$", color=fc, ha="center", va="center")
    root.plot([x, x], [.85, .77], ":", color=fc, lw=2)
    root.plot([x, x], [.75, .67], ":", color=fc, lw=2)

    labels = ((.04, .96, 'A'), (.04, .54, 'B'))
    panel_labels(root, labels)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "pineapple-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
