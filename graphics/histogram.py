#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use R ggplot2 library to plot histogram
"""

import sys

from optparse import OptionParser

from jcvi.apps.R import RTemplate
from jcvi.apps.base import ActionDispatcher, debug
debug()

histogram_template = """
library(ggplot2)
data <- read.table('$numbersfile')
data <- data[data > $vmin]
data <- data[data < $vmax]
data <- data.frame($xlabel=data)
qplot($xlabel, data=data, geom='histogram', main='$title')
ggsave('$outfile')
"""


def histogram(numbersfile, vmin, vmax, xlabel, title):
    """
    Generate histogram using number from numberfile, and use only numbers in the
    range of (vmin, vmax)
    """
    outfile = numbersfile + '.pdf'

    rtemplate = RTemplate(histogram_template, locals())
    rtemplate.run()


def main():
    """
    %prog numbers.txt

    """
    p = OptionParser(main.__doc__)
    p.add_option("--vmin", dest="vmin", default=-1e9, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--vmax", dest="vmax", default=1e9, type="int",
            help="numbers larger than this is ignored [default: %default]")
    p.add_option("--xlabel", dest="xlabel", default="value",
            help="label on the X-axis")
    p.add_option("--title", dest="title", default="title",
            help="title of the plot")
    opts, args = p.parse_args()

    if len(args) != 1:
        sys.exit(p.print_help())

    numbersfile = args[0]
    vmin, vmax = opts.vmin, opts.vmax
    xlabel, title = opts.xlabel, opts.title

    histogram(numbersfile, vmin, vmax, xlabel, title)


if __name__ == '__main__':
    main()
