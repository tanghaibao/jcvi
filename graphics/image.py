#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import sys
import logging

import matplotlib.image as mpimg
from jcvi.graphics.base import plt, savefig
from jcvi.apps.base import OptionParser, ActionDispatcher, fname


def main():

    actions = (
        ('seeds', 'extract seed color from images'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def seeds(args):
    """
    %prog seeds pngfile

    Extract seed color from pngfile.
    """
    p = OptionParser(seeds.__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    img = mpimg.imread(pngfile)
    w, h, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(pngfile, w, h))

    plt.figure(1, (iopts.w, iopts.h))
    plt.imshow(img)

    image_name = fname() + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
