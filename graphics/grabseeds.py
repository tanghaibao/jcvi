#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import os.path as op
import sys
import logging
from math import sin, cos

import numpy as np
from jcvi.graphics.base import plt, savefig, normalize_axes, Rectangle

from Image import open as iopen
from wand.image import Image
from pytesseract import image_to_string
from scipy.ndimage import binary_fill_holes, label, distance_transform_edt
from skimage import filter, color
from skimage.feature import peak_local_max
from skimage.measure import regionprops
from skimage.morphology import disk, closing, watershed
from jcvi.utils.counter import Counter
from jcvi.utils.webcolors import rgb_to_hex, closest_color
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher


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


def slice(s, m):
    assert ':' in s
    ra, rb = s.split(':')
    ra = 0 if ra == '' else int(ra)
    rb = m if rb == '' else int(rb)
    return ra, rb


def convert_image(pngfile, resize=1000, format="jpeg", rotate=0,
                  rows=':', cols=':', labelrows=None, labelcols=None):
    pf = pngfile.rsplit(".", 1)[0]
    resizefile = pf + ".resize.jpg"
    mainfile = pf + ".main.jpg"
    labelfile = pf + ".label.jpg"
    img = Image(filename=pngfile)
    exif = dict((k, v) for k, v in img.metadata.items() if k.startswith('exif:'))

    # Rotation, slicing and cropping of main image
    if rotate:
        img.rotate(rotate)
    if resize:
        w, h = img.size
        nw, nh = resize, resize * h / w
        img.resize(nw, nh)
        logging.debug("Image resized from {0}px:{1}px to {2}px:{3}px".\
                        format(w, h, nw, nh))
    img.format = format
    img.save(filename=resizefile)

    rimg = img.clone()
    if rows != ':' or cols != ':':
        w, h = img.size
        ra, rb = slice(rows, h)
        ca, cb = slice(cols, w)
        # left, top, right, bottom
        img.crop(ca, ra, cb, rb)
        logging.debug("Crop image to {0}:{1} {2}:{3}".format(ra, rb, ca, cb))
        img.format = format
        img.save(filename=mainfile)
    else:
        mainfile = resizefile

    # Extract text labels from image
    if labelrows or labelcols:
        if labelrows and not labelcols:
            labelcols = ':'
        if labelcols and not labelrows:
            labelrows = ':'
        ra, rb = slice(labelrows, h)
        ca, cb = slice(labelcols, w)
        rimg.crop(ca, ra, cb, rb)
        logging.debug("Extract label from {0}:{1} {2}:{3}".format(ra, rb, ca, cb))
        rimg.format = format
        rimg.save(filename=labelfile)
    else:
        labelfile = None

    return resizefile, mainfile, labelfile, exif


def load_image(resizefile):
    img = plt.imread(resizefile)
    img = np.flipud(img)
    h, w, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(resizefile, w, h))
    return img


def seeds(args):
    """
    %prog seeds [pngfile|jpgfile]

    Extract seed color from [pngfile|jpgfile]. Use --rows and --cols to crop image.
    """
    p = OptionParser(seeds.__doc__)
    g1 = OptionGroup(p, "Image manipulation")
    g1.add_option("--rotate", default=0, type="int",
                help="Rotate degrees clockwise")
    g1.add_option("--rows", default=':',
                help="Crop rows e.g. `:800` from first 800 rows")
    g1.add_option("--cols", default=':',
                help="Crop cols e.g. `-800:` from last 800 cols")
    g1.add_option("--labelrows",
                help="Label rows e.g. `:800` from first 800 rows")
    g1.add_option("--labelcols",
                help="Label cols e.g. `-800: from last 800 rows")
    p.add_option_group(g1)

    g2 = OptionGroup(p, "Object recognition")
    g2.add_option("--minsize", default=.001, type="float",
                help="Min ratio of object to image")
    g2.add_option("--maxsize", default=.5, type="float",
                help="Max ratio of object to image")
    g2.add_option("--count", default=100, type="int",
                help="Report max number of objects")
    g2.add_option("--watershed", default=False, action="store_true",
                help="Run watershed to segment touching objects")
    p.add_option_group(g2)

    g3 = OptionGroup(p, "De-noise")
    valid_filters = ("canny", "roberts", "sobel")
    g3.add_option("--filter", default="canny", choices=valid_filters,
                help="Edge detection algorithm")
    g3.add_option("--sigma", default=1, type="int",
                help="Canny edge detection sigma, higher for noisy image")
    g3.add_option("--kernel", default=2, type="int",
                help="Edge closure, higher if the object edges are dull")
    p.add_option_group(g3)

    g4 = OptionGroup(p, "Output")
    g4.add_option("--edges", default=False, action="store_true",
                help="Visualize edges in middle PDF panel")
    p.add_option_group(g4)
    opts, args, iopts = p.set_image_options(args, figsize='12x6')

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    pf = op.basename(pngfile).split(".")[0]
    sigma, kernel = opts.sigma, opts.kernel
    rows, cols = opts.rows, opts.cols
    labelrows, labelcols = opts.labelrows, opts.labelcols
    ff = opts.filter

    resizefile, mainfile, labelfile, exif = convert_image(pngfile, rotate=opts.rotate,
                                    rows=rows, cols=cols,
                                    labelrows=labelrows, labelcols=labelcols)

    oimg = load_image(resizefile)
    img = load_image(mainfile)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1,
                                             figsize=(iopts.w, iopts.h))

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
    label_objects, nb_labels = label(filled)

    # Watershed algorithm
    if opts.watershed:
        distance = distance_transform_edt(filled)
        local_maxi = peak_local_max(distance, indices=False, footprint=np.ones((3, 3)),
                                    labels=filled)
        markers, nmarkers = label(local_maxi)
        logging.debug("Identified {0} watershed markers".format(nmarkers))
        label_objects = watershed(-distance, markers, mask=filled)

    # Object size filtering
    sizes = np.bincount(label_objects.ravel())
    w, h = img_gray.shape
    min_size = int(w * h * opts.minsize)
    max_size = int(w * h * opts.maxsize)
    logging.debug("Find objects with pixels between {0} and {1}"\
                    .format(min_size, max_size))
    mask_sizes = np.logical_and(sizes >= min_size, sizes <= max_size)
    cleaned = mask_sizes[label_objects]

    olabels, nb_labels = label(cleaned)
    logging.debug("A total of {0} objects identified.".format(nb_labels))

    # Plotting
    ax1.set_title('Original picture')
    ax1.imshow(oimg)

    ax2.set_title('Edge detection\n({0}, $\sigma$={1}, $k$={2})'.\
                    format(ff, sigma, kernel))
    edges = color.gray2rgb(edges)
    cleaned = color.gray2rgb(cleaned)
    ax2_img = edges if opts.edges else cleaned
    ax2.imshow(ax2_img)

    ax3.set_title('Object detection')
    ax3.imshow(img)

    data = []
    # Calculate region properties
    rp = regionprops(olabels)
    for i, props in enumerate(rp):
        i += 1
        if i > opts.count:
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
                    ha="center", va="center", size=6)

    for ax in (ax2, ax3):
        ax.set_xlim(0, h)
        ax.set_ylim(w, 0)

    filename = op.basename(pngfile).replace('_', '\_')
    if labelfile:
        accession = image_to_string(iopen(labelfile))
        accession = " ".join(accession.split())  # normalize spaces
        accession = accession.replace('_', '\_')
    else:
        accession = pf

    ax4.text(.1, .82, "File: {0}".format(filename), color='g')
    ax4.text(.1, .76, "Label: {0}".format(accession), color='m')
    # Output identified seed stats
    yy = .7
    for i, npixels, rgbx, major, minor in data:
        pixeltag = "length={0} width={1} size={2}".\
                        format(int(major), int(minor), npixels)
        hashtag = ",".join(str(x) for x in rgbx)
        hashtag = "{0} {1}".format(hashtag, closest_color(rgbx))
        print >> sys.stderr, accession, "seed", "{0}:".format(i), \
                        pixeltag, hashtag
        if i > 7:
            continue
        ax4.text(.01, yy, str(i), va="center", bbox=dict(fc='none', ec='k'))
        ax4.text(.1, yy, pixeltag, va="center")
        yy -= .04
        ax4.add_patch(Rectangle((.1, yy - .025), .12, .05, lw=0,
                      fc=rgb_to_hex(rgbx)))
        ax4.text(.27, yy, hashtag, va="center")
        yy -= .06

    ax4.text(.1 , yy, "(A total of {0} objects displayed)".format(nb_labels),
             color="darkslategrey")
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
