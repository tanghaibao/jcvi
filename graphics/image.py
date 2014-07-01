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
from matplotlib.colors import rgb2hex
from jcvi.graphics.base import plt, savefig, normalize_axes, Rectangle

from wand.image import Image
from scipy import ndimage
from skimage import filter, color, img_as_float
from skimage import exposure
from skimage.measure import regionprops
from skimage.morphology import disk, closing
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, need_update


np.seterr(all="ignore")


class Seed (object):

    def __init__(self, img):
        pass


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


def load_image(pngfile, resize=1000, format="jpeg", rotate=0):
    resizefile = pngfile.rsplit(".", 1)[0] + ".resize.jpg"
    if need_update(pngfile, resizefile):
        img = Image(filename=pngfile)
        if rotate:
            img.rotate(rotate)
        w, h = img.size
        nw, nh = resize, resize * h / w
        img.resize(nw, nh)
        img.format = format
        img.save(filename=resizefile)
        logging.debug("Image resized from ({0} x {1}px) to ({2}px x {3}px)".\
                        format(w, h, nw, nh))

    img = plt.imread(resizefile)
    img = np.flipud(img)
    w, h, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(resizefile, w, h))
    return img


def seeds(args):
    """
    %prog seeds [pngfile|jpgfile]

    Extract seed color from [pngfile|jpgfile]. Use --rows and --cols to crop image.
    """
    p = OptionParser(seeds.__doc__)
    g1 = OptionGroup(p, "Image manipulation")
    g1.add_option("--rows", default=':',
                help="Crop rows e.g. `:800` takes first 800 rows")
    g1.add_option("--cols", default=':',
                help="Crop cols e.g. `800:` takes last 800 cols")
    g1.add_option("--rotate", default=0, type="int",
                help="Rotate degrees clockwise")
    p.add_option_group(g1)

    g2 = OptionGroup(p, "Object recognition")
    g2.add_option("--minsize", default=.0005, type="float",
                help="Min ratio of object to image")
    g2.add_option("--maxsize", default=.5, type="float",
                help="Max ratio of object to image")
    g2.add_option("--count", default=10, type="int",
                help="Report max number of objects")
    p.add_option_group(g2)

    g3 = OptionGroup(p, "De-noise")
    g3.add_option("--sigma", default=3, type="int",
                help="Canny edge detection, higher for noisy image")
    g3.add_option("--kernel", default=3, type="int",
                help="Edge closure, higher for noisy image")
    p.add_option_group(g3)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    pf = op.basename(pngfile).split(".")[0]
    sigma, kernel = opts.sigma, opts.kernel

    img = load_image(pngfile, rotate=opts.rotate)
    w, h, c = img.shape
    # Crop image
    ra, rb = opts.rows.split(":")
    ca, cb = opts.cols.split(":")
    ra = 0 if ra == '' else int(ra)
    rb = w if rb == '' else int(rb)
    ca = 0 if ca == '' else int(ca)
    cb = h if cb == '' else int(cb)
    if opts.rows != ':' or opts.cols != ':':
        img = img[ra:rb, ca:cb]
        logging.debug("Crop image to {0}:{1} {2}:{3}".format(ra, rb, ca, cb))

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, figsize=(12, 6))

    img_gray = color.rgb2gray(img)
    edges = filter.canny(img_gray, sigma=opts.sigma)
    selem = disk(opts.kernel)
    closed = closing(edges, selem)
    filled = ndimage.binary_fill_holes(closed)

    w, h = img_gray.shape
    min_size = int(w * h * opts.minsize)
    max_size = int(w * h * opts.maxsize)
    label_objects, nb_labels = ndimage.label(filled)
    print nb_labels
    sizes = np.bincount(label_objects.ravel())
    print sizes
    logging.debug("Find objects with pixels between {0} and {1}"\
                    .format(min_size, max_size))
    mask_sizes = np.logical_and(sizes >= min_size, sizes <= max_size)
    cleaned = mask_sizes[label_objects]

    label_objects, nb_labels = ndimage.label(cleaned)
    print nb_labels

    ax1.set_title('Original picture')
    ax1.imshow(img)

    ax2.set_title(r'Edge detection (\sigma={0}, k={1})'.format(sigma, kernel))
    edges = color.gray2rgb(edges)
    cleaned = color.gray2rgb(cleaned)
    ax2.imshow(cleaned)

    ax3.set_title('Object detection')
    ax3.imshow(img)

    img_ravel = []
    for row in img:
        for r, g, b in row:
            img_ravel.append((r, g, b))
    label_ravel = label_objects.ravel()
    assert len(img_ravel) == len(label_ravel)
    data = []
    for i, region in enumerate(regionprops(label_objects)):
        if i >= opts.count:
            break
        i += 1
        pixels = [pix for pix, label in zip(img_ravel, label_ravel) \
                        if label == i]
        npixels, rgb_1q, rgb_m, rgb_3q = pixel_stats(pixels)
        data.append((i, npixels, rgb_m))
        # draw rectangle around segmented coins
        minr, minc, maxr, maxc = region.bbox
        rect = Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, ec='w', lw=1)
        ax3.add_patch(rect)
        mc, mr = (minc + maxc) / 2, (minr + maxr) / 2
        ax3.text(mc, mr, "\#{0}".format(i), color='w',
                    ha="center", va="center")

    yy = .7
    for i, npixels, rgbx in data:
        itag =  "\#{0}:".format(i)
        pixeltag = "{0} pixels".format(npixels)
        rgbx = [x / 255. for x in rgbx]
        hashtag = ",".join(str(x) for x in rgb2ints(rgbx))
        print >> sys.stderr, itag, pixeltag, hashtag
        ax4.text(.05, yy, itag, va="center")
        ax4.text(.2, yy, pixeltag, va="center")
        yy -= .04
        ax4.add_patch(Rectangle((.2, yy - .025), .25, .05, lw=0,
                      fc=rgb2hex(rgbx)))
        ax4.text(.5, yy, hashtag, va="center")
        yy -= .06
    normalize_axes(ax4)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
