#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import os.path as op
import sys
import logging
from math import sqrt

import numpy as np
import matplotlib.image as mpimg
from matplotlib.colors import rgb2hex
from jcvi.graphics.base import plt, savefig, normalize_axes, Rectangle

from scipy import ndimage
from skimage import filter, color, img_as_float
from skimage.transform import hough_ellipse, probabilistic_hough_line
from skimage.draw import ellipse_perimeter
from skimage import exposure
from skimage.morphology import disk, closing
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('seeds', 'extract seed color from images'),
        ('rgb', 'extract rgb from image'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rgb2ints(rgbx):
    r, g, b = rgbx
    r *= 255
    g *= 255
    b *= 255
    return int(round(r)), int(round(g)), int(round(b))


def rgb(args):
    """
    %prog rgb pngfile

    Extract RGB from image.
    """
    p = OptionParser(rgb.__doc__)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    pf = op.basename(pngfile).split(".")[0]
    img = load_image(pngfile)

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(6, 3))
    ax1.set_title('Original picture')
    ax1.imshow(img)

    img_rgbi = []
    intensity = lambda x: sqrt((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 3)
    for row in img:
        for r, g, b, a in row:
            its = intensity((r, g, b))
            if its > .95:
                continue
            img_rgbi.append((r, g, b))
    img_rgbi.sort(key=intensity)
    npixels = len(img_rgbi)
    logging.debug("A total of {0} pixels imported".format(npixels))
    rgb_1q = img_rgbi[npixels / 4 * 3]
    rgb_m = img_rgbi[npixels / 2]
    rgb_3q = img_rgbi[npixels / 4]

    yy = .6
    ax2.text(.5, .8, pf.replace("_", "-"), color="g", ha="center", va="center")
    for t, rgbx in zip(("1st quartile", "Median", "3rd quartile"),
                       (rgb_1q, rgb_m, rgb_3q)):
        ax2.add_patch(Rectangle((.55, yy - .04), .08, .08, lw=0,
                      fc=rgb2hex(rgbx)))
        hashtag = rgb2ints(rgbx)
        hashtag = ",".join(str(x) for x in hashtag)
        print >> sys.stderr, pf, t, hashtag
        ax2.text(.5, yy, t, ha="right", va="center")
        ax2.text(.65, yy, hashtag, va="center")
        yy -= .1
    normalize_axes(ax2)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


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


def load_image(pngfile):
    img = mpimg.imread(pngfile)
    w, h, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(pngfile, w, h))
    return img


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
    img = load_image(pngfile)
    img = img[200:500, 100:400]

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(8, 4))

    img_gray = color.rgb2gray(img)
    #img_gray[img_gray < .2] = 1
    #edges = filter.sobel(img_gray)
    edges = filter.canny(img_gray)
    #filled = ndimage.binary_fill_holes(edges)
    #selem = disk(1)
    #closed = closing(edges, selem)
    lines = probabilistic_hough_line(edges, threshold=30,
                                    line_length=30, line_gap=1)

    #result = hough_ellipse(edges, min_size=40, max_size=40)
    #print result
    ax1.set_title('Original picture')
    ax1.imshow(img)

    ax2.set_title('Edge detection')
    #display_histogram(ax2, img_gray)
    edges = color.gray2rgb(edges)
    w, h, c = edges.shape
    #closed = color.gray2rgb(closed)
    ax2.imshow(edges)
    for p0, p1 in lines:
        ax, ay = p0[0:2]
        bx, by = p1[0:2]
        plt.plot((ax, bx), (ay, by), "r-")
    ax2.set_xlim(0, w)
    ax2.set_ylim(h, 0)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
