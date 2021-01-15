#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""
import json
import os.path as op
import sys
import string
import logging
from math import sin, cos, pi

from collections import Counter

import numpy as np
from jcvi.graphics.base import (
    plt,
    savefig,
    normalize_axes,
    Rectangle,
    latex,
    load_image,
)

from PIL.Image import open as iopen
from pytesseract import image_to_string
from wand.image import Image
from webcolors import rgb_to_hex, normalize_integer_triplet
from scipy.ndimage import binary_fill_holes, distance_transform_edt
from scipy.optimize import fmin_bfgs as fmin
from skimage.color import gray2rgb, rgb2gray
from skimage.filters import roberts, sobel
from skimage.feature import canny, peak_local_max
from skimage.measure import regionprops, label
from skimage.morphology import disk, closing, watershed
from skimage.segmentation import clear_border

from jcvi.utils.webcolors import closest_color
from jcvi.formats.base import must_open
from jcvi.algorithms.formula import reject_outliers, get_kmeans
from jcvi.apps.base import (
    OptionParser,
    OptionGroup,
    ActionDispatcher,
    iglob,
    mkdir,
    datadir,
)


np.seterr(all="ignore")


class Seed(object):
    def __init__(self, imagename, accession, seedno, rgb, props, exif):
        self.imagename = imagename
        self.accession = accession
        self.seedno = seedno
        y, x = props.centroid
        self.x, self.y = int(round(x)), int(round(y))
        self.location = "{0}|{1}".format(self.x, self.y)
        self.area = int(round(props.area))
        self.length = int(round(props.major_axis_length))
        self.width = int(round(props.minor_axis_length))
        self.props = props
        self.circularity = 4 * pi * props.area / props.perimeter ** 2
        self.rgb = rgb
        self.colorname = closest_color(rgb)
        self.datetime = exif.get("exif:DateTimeOriginal", "none")
        self.rgbtag = triplet_to_rgb(rgb)
        self.pixeltag = "length={0} width={1} area={2}".format(
            self.length, self.width, self.area
        )
        self.hashtag = " ".join((self.rgbtag, self.colorname))
        self.calibrated = False

    def __str__(self):
        fields = [
            self.imagename,
            self.datetime,
            self.accession,
            self.seedno,
            self.location,
            self.area,
            "{0:.2f}".format(self.circularity),
            self.length,
            self.width,
            self.colorname,
            self.rgbtag,
        ]
        if self.calibrated:
            fields += [
                self.pixelcmratio,
                self.rgbtransform,
                self.correctedlength,
                self.correctedwidth,
                self.correctedcolorname,
                self.correctedrgb,
            ]
        return "\t".join(str(x) for x in fields)

    @classmethod
    def header(cls, calibrate=False):
        fields = (
            "ImageName DateTime Accession SeedNum Location "
            "Area Circularity Length(px) Width(px) ColorName RGB".split()
        )
        if calibrate:
            fields += (
                "PixelCMratio RGBtransform Length(cm)"
                " Width(cm) CorrectedColorName CorrectedRGB".split()
            )
        return "\t".join(fields)

    def calibrate(self, pixel_cm_ratio, tr):
        self.pixelcmratio = "{0:.2f}".format(pixel_cm_ratio)
        self.rgbtransform = ",".join(["{0:.2f}".format(x) for x in tr.flatten()])
        self.correctedlength = "{0:.2f}".format(self.length / pixel_cm_ratio)
        self.correctedwidth = "{0:.2f}".format(self.width / pixel_cm_ratio)
        correctedrgb = np.dot(tr, np.array(self.rgb))
        self.correctedrgb = triplet_to_rgb(correctedrgb)
        self.correctedcolorname = closest_color(correctedrgb)
        self.calibrated = True


def rgb_to_triplet(rgb):
    return tuple([int(x) for x in rgb.split(",")])


def triplet_to_rgb(triplet):
    triplet = normalize_integer_triplet(triplet)
    return ",".join(str(int(round(x))) for x in triplet)


def main():

    actions = (
        ("batchseeds", "extract seed metrics for each image in a directory"),
        ("seeds", "extract seed metrics from one image"),
        ("calibrate", "calibrate pixel-inch ratio and color adjustment"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def total_error(x, colormap):
    xs = np.reshape(x, (3, 3))
    error_squared = sum([np.linalg.norm(np.dot(xs, o) - e) ** 2 for o, e in colormap])
    return error_squared ** 0.5


def calibrate(args):
    """
    %prog calibrate calibrate.JPG boxsize

    Calibrate pixel-inch ratio and color adjustment.
    - `calibrate.JPG` is the photo containig a colorchecker
    - `boxsize` is the measured size for the boxes on printed colorchecker, in
      squared centimeter (cm2) units
    """
    xargs = args[2:]
    p = OptionParser(calibrate.__doc__)
    opts, args, iopts = add_seeds_options(p, args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    imagefile, boxsize = args
    boxsize = float(boxsize)

    # Read in color checker
    colorcheckerfile = op.join(datadir, "colorchecker.txt")
    colorchecker = []
    expected = 0
    for row in open(colorcheckerfile):
        boxes = row.split()
        colorchecker.append(boxes)
        expected += len(boxes)

    folder = op.split(imagefile)[0]
    objects = seeds([imagefile, "--outdir={0}".format(folder)] + xargs)
    nseeds = len(objects)
    logging.debug("Found {0} boxes (expected={1})".format(nseeds, expected))
    assert (
        expected - 4 <= nseeds <= expected + 4
    ), "Number of boxes drastically different from {0}".format(expected)

    # Calculate pixel-cm ratio
    boxes = [t.area for t in objects]
    reject = reject_outliers(boxes)
    retained_boxes = [b for r, b in zip(reject, boxes) if not r]
    mbox = np.median(retained_boxes)  # in pixels
    pixel_cm_ratio = (mbox / boxsize) ** 0.5
    logging.debug(
        "Median box size: {0} pixels. Measured box size: {1} cm2".format(mbox, boxsize)
    )
    logging.debug("Pixel-cm ratio: {0}".format(pixel_cm_ratio))

    xs = [t.x for t in objects]
    ys = [t.y for t in objects]
    xs = [float(itemx) for itemx in xs]
    ys = [float(itemy) for itemy in ys]
    idx_xs = get_kmeans(xs, 6)
    idx_ys = get_kmeans(ys, 4)
    for xi, yi, s in zip(idx_xs, idx_ys, objects):
        s.rank = (yi, xi)

    objects.sort(key=lambda x: x.rank)

    colormap = []
    for s in objects:
        x, y = s.rank
        observed, expected = s.rgb, rgb_to_triplet(colorchecker[x][y])
        colormap.append((np.array(observed), np.array(expected)))

    # Color transfer
    tr0 = np.eye(3).flatten()
    print("Initial distance:", total_error(tr0, colormap), file=sys.stderr)
    tr = fmin(total_error, tr0, args=(colormap,))
    tr.resize((3, 3))
    print("RGB linear transform:\n", tr, file=sys.stderr)
    calib = {"PixelCMratio": pixel_cm_ratio, "RGBtransform": tr.tolist()}

    jsonfile = op.join(folder, "calibrate.json")
    fw = must_open(jsonfile, "w")
    print(json.dumps(calib, indent=4), file=fw)
    fw.close()
    logging.debug("Calibration specs written to `{0}`.".format(jsonfile))

    return jsonfile


def add_seeds_options(p, args):
    g1 = OptionGroup(p, "Image manipulation")
    g1.add_option("--rotate", default=0, type="int", help="Rotate degrees clockwise")
    g1.add_option(
        "--rows", default=":", help="Crop rows e.g. `:800` from first 800 rows"
    )
    g1.add_option(
        "--cols", default=":", help="Crop cols e.g. `-800:` from last 800 cols"
    )
    g1.add_option("--labelrows", help="Label rows e.g. `:800` from first 800 rows")
    g1.add_option("--labelcols", help="Label cols e.g. `-800: from last 800 rows")
    valid_colors = ("red", "green", "blue", "purple", "yellow", "orange", "INVERSE")
    g1.add_option(
        "--changeBackground",
        default=0,
        choices=valid_colors,
        help="Changes background color",
    )
    p.add_option_group(g1)

    g2 = OptionGroup(p, "Object recognition")
    g2.add_option(
        "--minsize",
        default=0.05,
        type="float",
        help="Min percentage of object to image",
    )
    g2.add_option(
        "--maxsize", default=50, type="float", help="Max percentage of object to image"
    )
    g2.add_option(
        "--count", default=100, type="int", help="Report max number of objects"
    )
    g2.add_option(
        "--watershed",
        default=False,
        action="store_true",
        help="Run watershed to segment touching objects",
    )
    p.add_option_group(g2)

    g3 = OptionGroup(p, "De-noise")
    valid_filters = ("canny", "roberts", "sobel")
    g3.add_option(
        "--filter",
        default="canny",
        choices=valid_filters,
        help="Edge detection algorithm",
    )
    g3.add_option(
        "--sigma",
        default=1,
        type="int",
        help="Canny edge detection sigma, higher for noisy image",
    )
    g3.add_option(
        "--kernel",
        default=2,
        type="int",
        help="Edge closure, higher if the object edges are dull",
    )
    g3.add_option(
        "--border", default=5, type="int", help="Remove image border of certain pixels"
    )
    p.add_option_group(g3)

    g4 = OptionGroup(p, "Output")
    g4.add_option("--calibrate", help="JSON file to correct distance and color")
    g4.add_option(
        "--edges",
        default=False,
        action="store_true",
        help="Visualize edges in middle PDF panel",
    )
    g4.add_option(
        "--outdir", default=".", help="Store intermediate images and PDF in folder"
    )
    g4.add_option("--prefix", help="Output prefix")
    g4.add_option(
        "--noheader", default=False, action="store_true", help="Do not print header"
    )
    p.add_option_group(g4)
    opts, args, iopts = p.set_image_options(args, figsize="12x6", style="white")

    return opts, args, iopts


def batchseeds(args):
    """
    %prog batchseeds folder

    Extract seed metrics for each image in a directory.
    """
    from jcvi.formats.pdf import cat

    xargs = args[1:]
    p = OptionParser(batchseeds.__doc__)
    opts, args, iopts = add_seeds_options(p, args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (folder,) = args
    folder = folder.rstrip("/")
    outdir = folder + "-debug"
    outfile = folder + "-output.tsv"
    assert op.isdir(folder)
    images = []
    jsonfile = opts.calibrate or op.join(folder, "calibrate.json")
    if not op.exists(jsonfile):
        jsonfile = None
    for im in iglob(folder, "*.jpg,*.JPG,*.png"):
        if im.endswith((".resize.jpg", ".main.jpg", ".label.jpg")):
            continue
        if op.basename(im).startswith("calibrate"):
            continue
        images.append(im)

    fw = must_open(outfile, "w")
    print(Seed.header(calibrate=jsonfile), file=fw)
    nseeds = 0
    for im in images:
        imargs = [im, "--noheader", "--outdir={0}".format(outdir)] + xargs
        if jsonfile:
            imargs += ["--calibrate={0}".format(jsonfile)]
        objects = seeds(imargs)
        for o in objects:
            print(o, file=fw)
        nseeds += len(objects)
    fw.close()
    logging.debug("Processed {0} images.".format(len(images)))
    logging.debug("A total of {0} objects written to `{1}`.".format(nseeds, outfile))

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
    assert ":" in s
    ra, rb = s.split(":")
    ra = 0 if ra == "" else int(ra)
    rb = m if rb == "" else int(rb)
    return ra, rb


def convert_background(pngfile, new_background):
    """Replace the background color with the specified background color, default is blue"""
    if new_background:
        _name, _ext = op.splitext(op.basename(pngfile))
        _name += "_bgxform"
        newfile = op.join(op.dirname(pngfile), _name + _ext)

        img = iopen(pngfile)
        pixels = list(img.getdata())
        h, w = img.size

        # Get Standard Deviation of RGB
        rgbArray = []
        for x in range(255):
            rgbArray.append(x)
        stdRGB = np.std(rgbArray) * 0.8

        # Get average color
        obcolor = [None, None, None]
        pixel_values = []
        for t in range(3):
            pixel_color = img.getdata(band=t)
            for pixel in pixel_color:
                if pixel > (stdRGB):
                    pixel_values.append(pixel)
            obcolor[t] = sum(pixel_values) / len(pixel_values)

        # Get background color using average color and standard deviation
        for t in range(3):
            pixel_color = img.getdata(band=t)
            seed_pixel_values = []
            for i in pixel_color:
                if (i > (obcolor[t] - stdRGB)) and (i < (obcolor[t] + stdRGB)):
                    seed_pixel_values.append(i)
            obcolor[t] = sum(seed_pixel_values) / len(seed_pixel_values)
        # Selection of colors based on option parser
        nbcolor = [None, None, None]
        if new_background == "INVERSE":
            nbcolor = [None, None, None]
            for t in range(3):
                nbcolor[t] = 255 - obcolor[t]
        elif new_background == "red":
            nbcolor = [255, 0, 0]

        elif new_background == "green":
            nbcolor = [0, 255, 0]

        elif new_background == "blue":
            nbcolor = [0, 0, 255]

        elif new_background == "yellow":
            nbcolor = [255, 255, 0]

        elif new_background == "purple":
            nbcolor = [255, 0, 255]

        elif new_background == "orange":
            nbcolor = [255, 165, 0]

        # Change Background Color
        obcolor = tuple(obcolor)
        nbcolor = tuple(nbcolor)
        for i in range(len(pixels)):
            r, g, b = pixels[i]
            if r >= (obcolor[0] - stdRGB) and r <= (obcolor[0] + stdRGB):
                if g >= (obcolor[1] - stdRGB) and g <= (obcolor[1] + stdRGB):
                    if b >= (obcolor[2] - stdRGB) and b <= (obcolor[2] + stdRGB):
                        pixels[i] = nbcolor
        img.putdata(pixels)
        img.save(newfile, "PNG")
        return newfile
    return pngfile


def convert_image(
    pngfile,
    pf,
    outdir=".",
    resize=1000,
    format="jpeg",
    rotate=0,
    rows=":",
    cols=":",
    labelrows=None,
    labelcols=None,
):
    resizefile = op.join(outdir, pf + ".resize.jpg")
    mainfile = op.join(outdir, pf + ".main.jpg")
    labelfile = op.join(outdir, pf + ".label.jpg")
    img = Image(filename=pngfile)
    exif = dict((k, v) for k, v in img.metadata.items() if k.startswith("exif:"))

    # Rotation, slicing and cropping of main image
    if rotate:
        img.rotate(rotate)
    if resize:
        w, h = img.size
        if min(w, h) > resize:
            if w < h:
                nw, nh = resize, resize * h / w
            else:
                nw, nh = resize * w / h, resize
            img.resize(nw, nh)
            logging.debug(
                "Image `{0}` resized from {1}px:{2}px to {3}px:{4}px".format(
                    pngfile, w, h, nw, nh
                )
            )
    img.format = format
    img.save(filename=resizefile)

    rimg = img.clone()
    if rows != ":" or cols != ":":
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
            labelcols = ":"
        if labelcols and not labelrows:
            labelrows = ":"
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

    (pngfile,) = args
    pf = opts.prefix or op.basename(pngfile).rsplit(".", 1)[0]
    sigma, kernel = opts.sigma, opts.kernel
    rows, cols = opts.rows, opts.cols
    labelrows, labelcols = opts.labelrows, opts.labelcols
    ff = opts.filter
    calib = opts.calibrate
    outdir = opts.outdir
    if outdir != ".":
        mkdir(outdir)
    if calib:
        calib = json.load(must_open(calib))
        pixel_cm_ratio, tr = calib["PixelCMratio"], calib["RGBtransform"]
        tr = np.array(tr)
    nbcolor = opts.changeBackground
    pngfile = convert_background(pngfile, nbcolor)
    resizefile, mainfile, labelfile, exif = convert_image(
        pngfile,
        pf,
        outdir=outdir,
        rotate=opts.rotate,
        rows=rows,
        cols=cols,
        labelrows=labelrows,
        labelcols=labelcols,
    )
    oimg = load_image(resizefile)
    img = load_image(mainfile)

    fig, (ax1, ax2, ax3, ax4) = plt.subplots(
        ncols=4, nrows=1, figsize=(iopts.w, iopts.h)
    )
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
        local_maxi = peak_local_max(distance, threshold_rel=0.05, indices=False)
        coordinates = peak_local_max(distance, threshold_rel=0.05)
        markers, nmarkers = label(local_maxi, return_num=True)
        logging.debug("Identified {0} watershed markers".format(nmarkers))
        labels = watershed(closed, markers, mask=filled)
    else:
        labels = label(filled)

    # Object size filtering
    w, h = img_gray.shape
    canvas_size = w * h
    min_size = int(round(canvas_size * opts.minsize / 100))
    max_size = int(round(canvas_size * opts.maxsize / 100))
    logging.debug(
        "Find objects with pixels between {0} ({1}%) and {2} ({3}%)".format(
            min_size, opts.minsize, max_size, opts.maxsize
        )
    )

    # Plotting
    ax1.set_title("Original picture")
    ax1.imshow(oimg)

    params = "{0}, $\sigma$={1}, $k$={2}".format(ff, sigma, kernel)
    if opts.watershed:
        params += ", watershed"
    ax2.set_title("Edge detection\n({0})".format(params))
    closed = gray2rgb(closed)
    ax2_img = labels
    if opts.edges:
        ax2_img = closed
    elif opts.watershed:
        ax2.plot(coordinates[:, 1], coordinates[:, 0], "g.")
    ax2.imshow(ax2_img, cmap=iopts.cmap)

    ax3.set_title("Object detection")
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
        ax2.plot((x0 - major_dx, x0 + major_dx), (y0 + major_dy, y0 - major_dy), "r-")
        ax2.plot((x0 - minor_dx, x0 + minor_dx), (y0 - minor_dy, y0 + minor_dy), "r-")

        npixels = int(props.area)
        # Sample the center of the blob for color
        d = min(int(round(minor / 2 * 0.35)) + 1, 50)
        x0d, y0d = int(round(x0)), int(round(y0))
        square = img[(y0d - d) : (y0d + d), (x0d - d) : (x0d + d)]
        pixels = []
        for row in square:
            pixels.extend(row)
        logging.debug(
            "Seed #{0}: {1} pixels ({2} sampled) - {3:.2f}%".format(
                i, npixels, len(pixels), 100.0 * npixels / canvas_size
            )
        )

        rgb = pixel_stats(pixels)
        objects.append(Seed(filename, accession, i, rgb, props, exif))
        minr, minc, maxr, maxc = props.bbox
        rect = Rectangle(
            (minc, minr), maxc - minc, maxr - minr, fill=False, ec="w", lw=1
        )
        ax3.add_patch(rect)
        mc, mr = (minc + maxc) / 2, (minr + maxr) / 2
        ax3.text(mc, mr, "{0}".format(i), color="w", ha="center", va="center", size=6)

    for ax in (ax2, ax3):
        ax.set_xlim(0, h)
        ax.set_ylim(w, 0)

    # Output identified seed stats
    ax4.text(0.1, 0.92, "File: {0}".format(latex(filename)), color="g")
    ax4.text(0.1, 0.86, "Label: {0}".format(latex(accession)), color="m")
    yy = 0.8
    fw = must_open(opts.outfile, "w")
    if not opts.noheader:
        print(Seed.header(calibrate=calib), file=fw)
    for o in objects:
        if calib:
            o.calibrate(pixel_cm_ratio, tr)
        print(o, file=fw)
        i = o.seedno
        if i > 7:
            continue
        ax4.text(0.01, yy, str(i), va="center", bbox=dict(fc="none", ec="k"))
        ax4.text(0.1, yy, o.pixeltag, va="center")
        yy -= 0.04
        ax4.add_patch(
            Rectangle((0.1, yy - 0.025), 0.12, 0.05, lw=0, fc=rgb_to_hex(o.rgb))
        )
        ax4.text(0.27, yy, o.hashtag, va="center")
        yy -= 0.06
    ax4.text(
        0.1,
        yy,
        "(A total of {0} objects displayed)".format(nb_labels),
        color="darkslategray",
    )
    normalize_axes(ax4)

    for ax in (ax1, ax2, ax3):
        xticklabels = [int(x) for x in ax.get_xticks()]
        yticklabels = [int(x) for x in ax.get_yticks()]
        ax.set_xticklabels(xticklabels, family="Helvetica", size=8)
        ax.set_yticklabels(yticklabels, family="Helvetica", size=8)

    image_name = op.join(outdir, pf + "." + iopts.format)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    return objects


if __name__ == "__main__":
    main()
