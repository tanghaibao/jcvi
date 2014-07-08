#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""

import os.path as op
import sys
import string
import logging
from math import sin, cos

import numpy as np
from jcvi.graphics.base import plt, savefig, normalize_axes, Rectangle, latex

from Image import open as iopen
from wand.image import Image
from scipy.ndimage import binary_fill_holes, distance_transform_edt
from skimage.color import gray2rgb, rgb2gray
from skimage.filter import canny, roberts, sobel
from skimage.feature import peak_local_max
from skimage.measure import regionprops, label
from skimage.morphology import disk, closing, watershed
from skimage.segmentation import clear_border
from jcvi.utils.counter import Counter
from jcvi.utils.webcolors import rgb_to_hex, closest_color
from jcvi.formats.base import must_open
from jcvi.apps.tesseract import image_to_string
from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, \
                iglob, mkdir


np.seterr(all="ignore")


class Seed (object):

    def __init__(self, imagename, accession, seedno, rgb, props, exif):
        self.imagename = imagename
        self.accession = accession
        self.seedno = seedno
        self.area = int(round(props.area))
        self.length = int(round(props.major_axis_length))
        self.width = int(round(props.minor_axis_length))
        self.props = props
        self.rgb = rgb
        self.colorname = closest_color(rgb)
        self.datetime = exif.get('exif:DateTimeOriginal', "none")
        self.rgbtag = ','.join(str(x) for x in rgb)
        self.pixeltag = "length={0} width={1} area={2}".\
                        format(self.length, self.width, self.area)
        self.hashtag = " ".join((self.rgbtag, self.colorname))

    def __str__(self):
        return " ".join(str(x) for x in (self.accession, "seed",
                        "{0}:".format(self.seedno), self.pixeltag, self.hashtag))

    @classmethod
    def header(cls):
        fields = "ImageName DateTime Accession SeedNum Area Length Width" \
                 " ColorName RGB"
        return "\t".join(fields.split())

    @property
    def tsvline(self):
        return "\t".join(str(x) for x in (self.imagename, self.datetime,
                        self.accession, self.seedno,
                        self.area, self.length, self.width,
                        self.colorname, self.rgbtag))


def main():

    actions = (
        ('batchseeds', 'extract seed metrics for each image in a directory'),
        ('seeds', 'extract seed metrics from one image'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add_seeds_options(p, args):
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
    g3.add_option("--border", default=10, type="int",
                help="Remove image border of certain pixels")
    p.add_option_group(g3)

    g4 = OptionGroup(p, "Output")
    g4.add_option("--edges", default=False, action="store_true",
                help="Visualize edges in middle PDF panel")
    g4.add_option("--outdir", default=".",
                help="Store intermediate images and PDF in folder")
    p.add_option_group(g4)
    opts, args, iopts = p.set_image_options(args, figsize='12x6')

    return opts, args, iopts


def batchseeds(args):
    """
    %prog batchseeds folder

    Extract seed metrics for each image in a directory.
    """
    from jcvi.formats.pdf import cat
    from jcvi.formats.excel import fromcsv

    xargs = args[1:]
    p = OptionParser(batchseeds.__doc__)
    opts, args, iopts = add_seeds_options(p, args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    folder = folder.rstrip('/')
    assert op.isdir(folder)
    images = []
    for im in iglob(folder, "*.jpg", "*.JPG", "*.png"):
        if im.endswith(".resize.jpg") or \
           im.endswith(".main.jpg") or \
           im.endswith(".label.jpg"):
            continue
        images.append(im)

    outdir = folder + "-debug"
    outfile = folder + "-output.tsv"
    fw = must_open(outfile, 'w')
    print >> fw, Seed.header()
    nseeds = 0
    for im in images:
        objects = seeds([im, "--outdir={0}".format(outdir)] + xargs)
        for o in objects:
            print >> fw, o.tsvline
        nseeds += len(objects)
    fw.close()
    logging.debug("Processed {0} images.".format(len(images)))
    excelfile = fromcsv([outfile, "--rgb=8"])
    logging.debug("A total of {0} objects written to `{1}`.".\
                    format(nseeds, excelfile))

    pdfs = iglob(outdir, "*.pdf")
    outpdf = folder + "-output.pdf"
    cat(pdfs + ["--outfile={0}".format(outpdf)])

    logging.debug("Debugging information written to `{0}`.".format(outpdf))
    return outfile


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


def convert_image(pngfile, outdir=".",
                  resize=1000, format="jpeg", rotate=0,
                  rows=':', cols=':', labelrows=None, labelcols=None):
    pf = op.basename(pngfile).rsplit(".", 1)[0]
    resizefile = op.join(outdir, pf + ".resize.jpg")
    mainfile = op.join(outdir, pf + ".main.jpg")
    labelfile = op.join(outdir, pf + ".label.jpg")
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
        logging.debug("Crop image to {0}:{1} {2}:{3}".format(ra, rb, ca, cb))
        img.crop(ca, ra, cb, rb)
        img.format = format
        img.save(filename=mainfile)
    else:
        mainfile = resizefile

    # Extract text labels from image
    if labelrows or labelcols:
        w, h = rimg.size
        if labelrows and not labelcols:
            labelcols = ':'
        if labelcols and not labelrows:
            labelrows = ':'
        ra, rb = slice(labelrows, h)
        ca, cb = slice(labelcols, w)
        logging.debug("Extract label from {0}:{1} {2}:{3}".format(ra, rb, ca, cb))
        rimg.crop(ca, ra, cb, rb)
        rimg.format = format
        rimg.save(filename=labelfile)
    else:
        labelfile = None

    return resizefile, mainfile, labelfile, exif


def extract_label(labelfile):
    accession = image_to_string(iopen(labelfile))
    accession = " ".join(accession.split())  # normalize spaces
    accession = "".join(x for x in accession if x in string.printable)
    if not accession:
        accession = "none"
    return accession


def load_image(resizefile):
    img = plt.imread(resizefile)
    img = np.flipud(img)
    h, w, c = img.shape
    logging.debug("Image `{0}` loaded ({1}px x {2}px).".format(resizefile, w, h))
    return img


def seeds(args):
    """
    %prog seeds [pngfile|jpgfile]

    Extract seed metrics from [pngfile|jpgfile]. Use --rows and --cols to crop image.
    """
    p = OptionParser(seeds.__doc__)
    p.set_outfile()
    opts, args, iopts = add_seeds_options(p, args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pngfile, = args
    pf = op.basename(pngfile).rsplit(".", 1)[0]
    sigma, kernel = opts.sigma, opts.kernel
    rows, cols = opts.rows, opts.cols
    labelrows, labelcols = opts.labelrows, opts.labelcols
    ff = opts.filter
    outdir = opts.outdir
    if outdir != '.':
        mkdir(outdir)

    resizefile, mainfile, labelfile, exif = \
                      convert_image(pngfile, outdir=outdir,
                                    rotate=opts.rotate,
                                    rows=rows, cols=cols,
                                    labelrows=labelrows, labelcols=labelcols)

    oimg = load_image(resizefile)
    img = load_image(mainfile)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1,
                                             figsize=(iopts.w, iopts.h))

    # Edge detection
    img_gray = rgb2gray(img)
    logging.debug("Running {0} edge detection ...".format(ff))
    if ff == "canny":
        edges = canny(img_gray, sigma=opts.sigma)
    elif ff == "roberts":
        edges = roberts(img_gray)
    elif ff == "sobel":
        edges = sobel(img_gray)
    edges = clear_border(edges, buffer_size=opts.border)
    selem = disk(kernel)
    closed = closing(edges, selem) if kernel else edges
    filled = binary_fill_holes(closed)

    # Watershed algorithm
    if opts.watershed:
        distance = distance_transform_edt(filled)
        local_maxi = peak_local_max(distance, threshold_rel=.05, indices=False)
        coordinates = peak_local_max(distance, threshold_rel=.05)
        markers, nmarkers = label(local_maxi, return_num=True)
        logging.debug("Identified {0} watershed markers".format(nmarkers))
        labels = watershed(closed, markers, mask=filled)
    else:
        labels = label(filled)

    # Object size filtering
    w, h = img_gray.shape
    min_size = int(w * h * opts.minsize)
    max_size = int(w * h * opts.maxsize)
    logging.debug("Find objects with pixels between {0} and {1}"\
                    .format(min_size, max_size))

    # Plotting
    ax1.set_title('Original picture')
    ax1.imshow(oimg)

    ax2.set_title('Edge detection\n({0}, $\sigma$={1}, $k$={2})'.\
                    format(ff, sigma, kernel))
    closed = gray2rgb(closed)
    ax2_img = labels
    if opts.edges:
        ax2_img = closed
    elif opts.watershed:
        ax2.plot(coordinates[:, 1], coordinates[:, 0], 'g.')
    ax2.imshow(ax2_img)

    ax3.set_title('Object detection')
    ax3.imshow(img)

    filename = op.basename(pngfile)
    if labelfile:
        accession = extract_label(labelfile)
    else:
        accession = pf

    # Calculate region properties
    rp = regionprops(labels)
    rp = [x for x in rp if min_size <= x.area <= max_size]
    nb_labels = len(rp)
    logging.debug("A total of {0} objects identified.".format(nb_labels))
    objects = []
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
        d = min(int(round(minor / 2 * .7)), 100)
        square = img[(y0 - d):(y0 + d), (x0 - d):(x0 + d)]
        pixels = []
        for row in square:
            pixels.extend(row)
        logging.debug("Seed #{0}: {1} pixels ({2} sampled)".\
                        format(i, npixels, len(pixels)))

        rgb = pixel_stats(pixels)
        objects.append(Seed(filename, accession, i, rgb, props, exif))
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

    # Output identified seed stats
    ax4.text(.1, .92, "File: {0}".format(latex(filename)), color='g')
    ax4.text(.1, .86, "Label: {0}".format(latex(accession)), color='m')
    yy = .8
    fw = must_open(opts.outfile, "w")
    for o in objects:
        print >> fw, o.tsvline
        i = o.seedno
        if i > 7:
            continue
        ax4.text(.01, yy, str(i), va="center", bbox=dict(fc='none', ec='k'))
        ax4.text(.1, yy, o.pixeltag, va="center")
        yy -= .04
        ax4.add_patch(Rectangle((.1, yy - .025), .12, .05, lw=0,
                      fc=rgb_to_hex(o.rgb)))
        ax4.text(.27, yy, o.hashtag, va="center")
        yy -= .06
    ax4.text(.1 , yy, "(A total of {0} objects displayed)".format(nb_labels),
             color="darkslategrey")
    normalize_axes(ax4)

    for ax in (ax1, ax2, ax3):
        xticklabels = [int(x) for x in ax.get_xticks()]
        yticklabels = [int(x) for x in ax.get_yticks()]
        ax.set_xticklabels(xticklabels, family='Helvetica')
        ax.set_yticklabels(yticklabels, family='Helvetica')

    image_name = op.join(outdir, pf + "." + iopts.format)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    return objects


if __name__ == '__main__':
    main()
