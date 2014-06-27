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
from skimage import exposure
from skimage.color import label2rgb
from skimage.measure import regionprops
from skimage.morphology import disk, closing
from jcvi.algorithms.formula import reject_outliers
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('seeds', 'extract seed color from images'),
        ('rgb', 'extract rgb from image'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


intensity = lambda x: sqrt((x[0] * x[0] + x[1] * x[1] + x[2] * x[2]) / 3)


def rgb2ints(rgbx):
    r, g, b = rgbx
    r *= 255
    g *= 255
    b *= 255
    return int(round(r)), int(round(g)), int(round(b))


def pixel_stats(img):
    img.sort(key=intensity)
    npixels = len(img)
    logging.debug("A total of {0} pixels imported".format(npixels))
    rgb_1q = img[npixels / 4 * 3]
    rgb_m = img[npixels / 2]
    rgb_3q = img[npixels / 4]
    return npixels, rgb_1q, rgb_m, rgb_3q


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
    for row in img:
        for r, g, b, a in row:
            its = intensity((r, g, b))
            if its > .99:
                continue
            img_rgbi.append((r, g, b))

    npixels, rgb_1q, rgb_m, rgb_3q = pixel_stats(img_rgbi)

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
    fig, (ax1, ax2, ax3) = plt.subplots(ncols=3, nrows=1, figsize=(8, 6))

    img_gray = color.rgb2gray(img)
    edges = filter.canny(img_gray)
    selem = disk(1)
    closed = closing(edges, selem)
    filled = ndimage.binary_fill_holes(closed)

    w, h = img_gray.shape
    max_size = w * h / 2
    label_objects, nb_labels = ndimage.label(filled)
    print nb_labels
    sizes = np.bincount(label_objects.ravel())
    print sizes
    mask_sizes = np.logical_and(sizes > 1000, sizes < max_size)
    print mask_sizes
    cleaned = mask_sizes[label_objects]

    #label_objects, nb_labels = ndimage.label(cleaned)
    #print nb_labels
    #mask_sizes = np.invert(reject_outliers(sizes))
    #cleaned = mask_sizes[label_objects]

    label_objects, nb_labels = ndimage.label(cleaned)
    print nb_labels
    #img_label_overlay = label2rgb(label_objects, image=img)

    ax1.set_title('Original picture')
    ax1.imshow(img)

    ax2.set_title('Object detection')
    edges = color.gray2rgb(edges)
    w, h, c = edges.shape
    ax2.imshow(img)
    img_ravel = []
    for row in img:
        for r, g, b in row:
            img_ravel.append((r, g, b))
    label_ravel = label_objects.ravel()
    assert len(img_ravel) == len(label_ravel)
    data = []
    for i, region in enumerate(regionprops(label_objects)):
        i += 1
        pixels = [pix for pix, label in zip(img_ravel, label_ravel) \
                        if label == i]
        npixels, rgb_1q, rgb_m, rgb_3q = pixel_stats(pixels)
        data.append((i, npixels, rgb_m))
        # draw rectangle around segmented coins
        minr, minc, maxr, maxc = region.bbox
        rect = Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, ec='w', lw=1)
        ax2.add_patch(rect)
        mc, mr = (minc + maxc) / 2, (minr + maxr) / 2
        ax2.text(mc, mr, "\#{0}".format(i), color='w',
                    ha="center", va="center")
    ax2.set_xlim(0, h)
    ax2.set_ylim(w, 0)

    yy = .7
    for i, npixels, rgbx in data:
        itag =  "\#{0}:".format(i)
        pixeltag = "{0} pixels".format(npixels)
        rgbx = [x / 255. for x in rgbx]
        hashtag = ",".join(str(x) for x in rgb2ints(rgbx))
        print >> sys.stderr, itag, pixeltag, hashtag
        ax3.text(.05, yy, itag, va="center")
        ax3.text(.2, yy, pixeltag, va="center")
        yy -= .04
        ax3.add_patch(Rectangle((.2, yy - .025), .25, .05, lw=0,
                      fc=rgb2hex(rgbx)))
        ax3.text(.5, yy, hashtag, va="center")
        yy -= .06
    normalize_axes(ax3)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
