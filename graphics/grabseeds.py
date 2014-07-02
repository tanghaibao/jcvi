#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import os.path as op
import sys
import logging
from math import sin, cos, pi

import numpy as np
from jcvi.graphics.base import plt, savefig, normalize_axes, Rectangle

from wand.image import Image
from scipy.ndimage import binary_fill_holes, label
from skimage import filter, color
from skimage.measure import regionprops
from skimage.morphology import disk, closing
from jcvi.utils.counter import Counter
from jcvi.utils.webcolors import rgb_to_hex, closest_color
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, need_update


np.seterr(all="ignore")


def main():

    actions = (
        ('seeds', 'extract seed color from images'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def p_round(n, precision=5):
    precision = int(precision)
    return int(round(n / float(precision))) * precision


def pixel_stats(img):
    img = [(p_round(r), p_round(g), p_round(b)) for r, g, b in img]
    c = Counter(img)
    imgx, count = c.most_common(1)[0]
    return imgx


def load_image(pngfile, resize=1000, format="jpeg", rotate=0):
    resizefile = pngfile.rsplit(".", 1)[0] + ".resize.jpg"
    if need_update(pngfile, resizefile):
        img = Image(filename=pngfile)
        if rotate:
            img.rotate(rotate)
        if resize:
            w, h = img.size
            nw, nh = resize, resize * h / w
            img.resize(nw, nh)
            logging.debug("Image resized from ({0} x {1}px) to ({2}px x {3}px)".\
                            format(w, h, nw, nh))
        img.format = format
        img.save(filename=resizefile)

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
    g2.add_option("--minsize", default=.001, type="float",
                help="Min ratio of object to image")
    g2.add_option("--maxsize", default=.5, type="float",
                help="Max ratio of object to image")
    g2.add_option("--count", default=10, type="int",
                help="Report max number of objects")
    p.add_option_group(g2)

    g3 = OptionGroup(p, "De-noise")
    valid_filters = ("canny", "roberts", "sobel")
    g3.add_option("--filter", default="canny", choices=valid_filters,
                help="Edge detection algorithm")
    g3.add_option("--sigma", default=1, type="int",
                help="Canny edge detection sigma, higher for noisy image")
    g3.add_option("--kernel", default=2, type="int",
                help="Edge closure, higher for noisy image")
    p.add_option_group(g3)
    opts, args, iopts = p.set_image_options(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    pf = op.basename(pngfile).split(".")[0]
    sigma, kernel = opts.sigma, opts.kernel
    ff = opts.filter

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

    # Edge detection
    img_gray = color.rgb2gray(img)
    logging.debug("Running {0} edge detection ...".format(ff))
    if ff == "canny":
        edges = filter.canny(img_gray, sigma=opts.sigma)
    elif ff == "roberts":
        edges = filter.roberts(img_gray)
    elif ff == "sobel":
        edges = filter.sobel(img_gray)
    selem = disk(opts.kernel)
    closed = closing(edges, selem)
    filled = binary_fill_holes(closed)

    # Object size filtering
    w, h = img_gray.shape
    min_size = int(w * h * opts.minsize)
    max_size = int(w * h * opts.maxsize)
    label_objects, nb_labels = label(filled)
    sizes = np.bincount(label_objects.ravel())
    logging.debug("Find objects with pixels between {0} and {1}"\
                    .format(min_size, max_size))
    mask_sizes = np.logical_and(sizes >= min_size, sizes <= max_size)
    cleaned = mask_sizes[label_objects]

    label_objects, nb_labels = label(cleaned)
    logging.debug("A total of {0} objects identified.".format(nb_labels))

    # Plotting
    ax1.set_title('Original picture')
    ax1.imshow(img)

    ax2.set_title('Edge detection\n({0}, $\sigma$={1}, $k$={2})'.\
                    format(ff, sigma, kernel))
    edges = color.gray2rgb(edges)
    cleaned = color.gray2rgb(cleaned)
    ax2.imshow(cleaned)

    ax3.set_title('Object detection')
    ax3.imshow(img)

    data = []
    # Calculate region properties
    rp = regionprops(label_objects)
    for i, props in enumerate(rp):
        i += 1
        if i >= opts.count:
            break

        y0, x0 = props.centroid
        orientation = props.orientation
        major, minor = props.major_axis_length, props.minor_axis_length
        major_dx = cos(orientation) * major / 2
        major_dy = sin(orientation) * major / 2
        minor_dx = sin(orientation) * minor / 2
        minor_dy = cos(orientation) * minor / 2
        ax2.plot((x0 - major_dx, x0 + major_dx),
                 (y0 + major_dy, y0 - major_dy), 'r-')
        ax2.plot((x0 - minor_dx, x0 + minor_dx),
                 (y0 - minor_dy, y0 + minor_dy), 'r-')

        npixels = int(props.area)
        # Sample the center of the blob for color
        d = int(round(minor / 2 * .7))
        square = img[(y0 - d):(y0 + d), (x0 - d):(x0 + d)]
        pixels = []
        for row in square:
            pixels.extend(row)
        logging.debug("Seed #{0}: {1} pixels ({2} sampled)".\
                        format(i, npixels, len(pixels)))

        rgbx = pixel_stats(pixels)
        data.append((i, npixels, rgbx, major, minor))
        minr, minc, maxr, maxc = props.bbox
        rect = Rectangle((minc, minr), maxc - minc, maxr - minr,
                                  fill=False, ec='w', lw=1)
        ax3.add_patch(rect)
        mc, mr = (minc + maxc) / 2, (minr + maxr) / 2
        ax3.text(mc, mr, "{0}".format(i), color='w',
                    ha="center", va="center", size=8)

    for ax in (ax2, ax3):
        ax.set_xlim(0, h)
        ax.set_ylim(w, 0)

    filename = op.basename(pngfile).replace('_', '\_')
    ax4.text(.1, .8, "File: {0}".format(filename), color='g')
    # Output identified seed stats
    yy = .7
    for i, npixels, rgbx, major, minor in data:
        itag =  "{0}:".format(i)
        ellipse_size = major * minor * pi / 4
        pixeltag = "length={0} width={1}".format(int(major), int(minor))
        sizetag = "size={0} ellipse={1}".format(npixels, int(ellipse_size))
        hashtag = ",".join(str(x) for x in rgbx)
        hashtag = "{0} {1}".format(hashtag, closest_color(rgbx))
        print pf, itag, pixeltag, sizetag, hashtag
        ax4.text(.05, yy, itag, va="center")
        ax4.text(.2, yy, pixeltag, va="center")
        yy -= .04
        ax4.add_patch(Rectangle((.2, yy - .025), .15, .05, lw=0,
                      fc=rgb_to_hex(rgbx)))
        ax4.text(.4, yy, hashtag, va="center")
        yy -= .06

    normalize_axes(ax4)

    for ax in (ax1, ax2, ax3):
        xticklabels = [int(x) for x in ax.get_xticks()]
        yticklabels = [int(x) for x in ax.get_yticks()]
        ax.set_xticklabels(xticklabels, family='Helvetica')
        ax.set_yticklabels(yticklabels, family='Helvetica')

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
