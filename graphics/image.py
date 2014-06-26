#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import os.path as op
import sys
import logging

import matplotlib.image as mpimg
from jcvi.graphics.base import plt, savefig

from skimage import filter, color, img_as_float
from skimage.transform import hough_ellipse
from skimage.draw import ellipse_perimeter
from skimage import exposure
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('seeds', 'extract seed color from images'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def display_histogram(ax_hist, img, bins=256):
    img = img_as_float(img)
    ax_cdf = ax_hist.twinx()
    ax_hist.hist(img.ravel(), bins=bins, histtype='step', color='black')
    ax_hist.ticklabel_format(axis='y', style='scientific', scilimits=(0, 0))
    ax_hist.set_xlabel('Pixel intensity')
    ax_hist.set_xlim(0, 1)

     # Display cumulative distribution
    img_cdf, bins = exposure.cumulative_distribution(img, bins)
    ax_cdf.plot(bins, img_cdf, 'r')
    ax_cdf.set_yticks([])


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
    pf = op.basename(pngfile).split(".")[0]
    img = mpimg.imread(pngfile)
    w, h, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(pngfile, w, h))
    img = img[200:500, 100:400]

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

    img_gray = color.rgb2gray(img)
    #img_gray[img_gray < .2] = 1
    #edges = filter.sobel(img_gray)
    edges = filter.canny(img_gray)

    #result = hough_ellipse(edges, min_size=40, max_size=40)
    #print result
    ax1.set_title('Original picture')
    ax1.imshow(img)

    ax2.set_title('Edge detection')
    display_histogram(ax2, img_gray)
    edges = color.gray2rgb(edges)
    #img = color.gray2rgb(img_gray)
    #ax2.imshow(edges)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
