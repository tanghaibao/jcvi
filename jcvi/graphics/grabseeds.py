#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Image processing pipelines for phenotyping projects.
"""
import json
import os.path as op
import string
import sys

from collections import Counter
from math import cos, pi, sin
from typing import Any, List, Optional, Tuple

import numpy as np

from PIL.Image import open as iopen
from pyefd import elliptic_fourier_descriptors
from pytesseract import image_to_string
from scipy.ndimage import binary_fill_holes, distance_transform_edt
from scipy.optimize import fmin_bfgs as fmin
from skimage.color import gray2rgb, rgb2gray
from skimage.feature import canny, peak_local_max
from skimage.filters import roberts, sobel
from skimage.measure import find_contours, regionprops, label
from skimage.morphology import disk, closing
from skimage.segmentation import clear_border, watershed
from wand.image import Image
from webcolors import rgb_to_hex, normalize_integer_triplet

from ..algorithms.formula import get_kmeans, reject_outliers
from ..apps.base import (
    ActionDispatcher,
    OptionParser,
    datadir,
    logger,
    iglob,
    mkdir,
)
from ..formats.base import must_open
from ..formats.pdf import cat
from ..utils.webcolors import closest_color

from .base import (
    Rectangle,
    latex,
    load_image,
    normalize_axes,
    plt,
    savefig,
    set_helvetica_axis,
)


np.seterr(all="ignore")

RGBTuple = Tuple[int, int, int]


class Seed(object):
    """
    Seed object with metrics.
    """

    def __init__(
        self,
        imagename: str,
        accession: str,
        seedno: int,
        rgb: RGBTuple,
        props: Any,
        efds: np.ndarray,
        exif: dict,
    ):
        self.imagename = imagename
        self.accession = accession
        self.seedno = seedno
        y, x = props.centroid
        self.x, self.y = int(round(x)), int(round(y))
        self.location = f"{self.x}|{self.y}"
        self.area = int(round(props.area))
        self.length = int(round(props.major_axis_length))
        self.width = int(round(props.minor_axis_length))
        self.props = props
        self.efds = efds
        self.circularity = 4 * pi * props.area / props.perimeter**2
        self.rgb = rgb
        self.colorname = closest_color(rgb)
        self.datetime = exif.get("exif:DateTimeOriginal", "none")
        self.rgbtag = triplet_to_rgb(rgb)
        self.pixeltag = f"length={self.length} width={self.width} area={self.area}"
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
            f"{self.circularity:.2f}",
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
        fields += [",".join(f"{x:.3f}" for x in self.efds)]
        return "\t".join(str(x) for x in fields)

    @classmethod
    def header(cls, calibrated: bool = False) -> str:
        """
        Return header line for the TSV file.
        """
        fields = (
            "ImageName DateTime Accession SeedNum Location "
            "Area Circularity Length(px) Width(px) ColorName RGB".split()
        )
        if calibrated:
            fields += (
                "PixelCMratio RGBtransform Length(cm)"
                " Width(cm) CorrectedColorName CorrectedRGB".split()
            )
        fields += ["EllipticFourierDescriptors"]
        return "\t".join(fields)

    def calibrate(self, pixel_cm_ratio: float, tr: np.ndarray):
        """
        Calibrate pixel-inch ratio and color adjustment.
        """
        self.pixelcmratio = f"{pixel_cm_ratio:.2f}"
        self.rgbtransform = ",".join([f"{x:.2f}" for x in tr.flatten()])
        self.correctedlength = f"{self.length / pixel_cm_ratio:.2f}"
        self.correctedwidth = f"{self.width / pixel_cm_ratio:.2f}"
        correctedrgb = np.dot(tr, np.array(self.rgb))
        self.correctedrgb = triplet_to_rgb(correctedrgb)
        self.correctedcolorname = closest_color(correctedrgb)
        self.calibrated = True


def rgb_to_triplet(rgb: str) -> RGBTuple:
    """
    Convert RGB string to triplet.
    """
    return tuple([int(x) for x in rgb.split(",")][:3])


def triplet_to_rgb(triplet: RGBTuple) -> str:
    """
    Convert triplet to RGB string.
    """
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


def total_error(x: np.ndarray, colormap: Tuple[Tuple[np.ndarray, np.ndarray]]) -> float:
    """
    Calculate total error between observed and expected colors.
    """
    xs = np.reshape(x, (3, 3))
    error_squared = sum(np.linalg.norm(np.dot(xs, o) - e) ** 2 for o, e in colormap)
    return error_squared**0.5


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
    _, args, _ = add_seeds_options(p, args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    imagefile, boxsize = args
    boxsize = float(boxsize)

    # Read in color checker
    colorcheckerfile = op.join(datadir, "colorchecker.txt")
    colorchecker = []
    expected = 0
    with open(colorcheckerfile, encoding="utf-8") as file:
        for row in file:
            boxes = row.split()
            colorchecker.append(boxes)
            expected += len(boxes)

    folder = op.split(imagefile)[0]
    objects = seeds([imagefile, f"--outdir={folder}"] + xargs)
    nseeds = len(objects)
    logger.debug("Found %d boxes (expected=%d)", nseeds, expected)
    assert (
        expected - 4 <= nseeds <= expected + 4
    ), f"Number of boxes drastically different from {expected}"

    # Calculate pixel-cm ratio
    boxes = [t.area for t in objects]
    reject = reject_outliers(boxes)
    retained_boxes = [b for r, b in zip(reject, boxes) if not r]
    mbox = np.median(retained_boxes)  # in pixels
    pixel_cm_ratio = (mbox / boxsize) ** 0.5
    logger.debug("Median box size: %d pixels. Measured box size: %d cm2", mbox, boxsize)
    logger.debug("Pixel-cm ratio: %.2f", pixel_cm_ratio)

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
    logger.debug("Calibration specs written to `%s`.", jsonfile)

    return jsonfile


def add_seeds_options(p, args):
    """
    Add options to the OptionParser for seeds() and batchseeds() functions.
    """
    g1 = p.add_argument_group("Image manipulation")
    g1.add_argument("--rotate", default=0, type=int, help="Rotate degrees clockwise")
    g1.add_argument(
        "--rows", default=":", help="Crop rows e.g. `:800` from first 800 rows"
    )
    g1.add_argument(
        "--cols", default=":", help="Crop cols e.g. `-800:` from last 800 cols"
    )
    g1.add_argument("--labelrows", help="Label rows e.g. `:800` from first 800 rows")
    g1.add_argument("--labelcols", help="Label cols e.g. `-800: from last 800 rows")
    valid_colors = ("red", "green", "blue", "purple", "yellow", "orange", "INVERSE")
    g1.add_argument(
        "--changeBackground",
        default=0,
        choices=valid_colors,
        help="Changes background color",
    )

    g2 = p.add_argument_group("Object recognition")
    g2.add_argument(
        "--minsize",
        default=0.05,
        type=float,
        help="Min percentage of object to image",
    )
    g2.add_argument(
        "--maxsize", default=50, type=float, help="Max percentage of object to image"
    )
    g2.add_argument(
        "--count", default=100, type=int, help="Report max number of objects"
    )
    g2.add_argument(
        "--watershed",
        default=False,
        action="store_true",
        help="Run watershed to segment touching objects",
    )

    g3 = p.add_argument_group("De-noise")
    valid_filters = ("canny", "roberts", "sobel")
    g3.add_argument(
        "--filter",
        default="canny",
        choices=valid_filters,
        help="Edge detection algorithm",
    )
    g3.add_argument(
        "--sigma",
        default=1,
        type=int,
        help="Canny edge detection sigma, higher for noisy image",
    )
    g3.add_argument(
        "--kernel",
        default=2,
        type=int,
        help="Edge closure, higher if the object edges are dull",
    )
    g3.add_argument(
        "--border", default=5, type=int, help="Remove image border of certain pixels"
    )

    g4 = p.add_argument_group("Output")
    g4.add_argument("--calibrate", help="JSON file to correct distance and color")
    g4.add_argument(
        "--edges",
        default=False,
        action="store_true",
        help="Visualize edges in middle PDF panel",
    )
    g4.add_argument(
        "--outdir", default=".", help="Store intermediate images and PDF in folder"
    )
    g4.add_argument("--prefix", help="Output prefix")
    g4.add_argument(
        "--noheader", default=False, action="store_true", help="Do not print header"
    )
    opts, args, iopts = p.set_image_options(args, figsize="12x6", style="white")

    return opts, args, iopts


def batchseeds(args):
    """
    %prog batchseeds folder

    Extract seed metrics for each image in a directory.
    """
    xargs = args[1:]
    p = OptionParser(batchseeds.__doc__)
    opts, args, _ = add_seeds_options(p, args)

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
    print(Seed.header(calibrated=bool(jsonfile)), file=fw)
    nseeds = 0
    for im in images:
        imargs = [im, "--noheader", f"--outdir={outdir}"] + xargs
        if jsonfile:
            imargs += [f"--calibrate={jsonfile}"]
        objects = seeds(imargs)
        for o in objects:
            print(o, file=fw)
        nseeds += len(objects)
    fw.close()
    logger.debug("Processed %d images.", len(images))
    logger.debug("A total of %d objects written to `%s`.", nseeds, outfile)

    pdfs = iglob(outdir, "*.pdf")
    outpdf = folder + "-output.pdf"
    cat(pdfs + [f"--outfile={outpdf}"])

    logger.debug("Debugging information written to `%s`.", outpdf)
    return outfile


def p_round(n: int, precision: int = 5) -> int:
    """
    Round to the nearest precision.
    """
    precision = int(precision)
    return int(round(n / float(precision))) * precision


def pixel_stats(img: List[RGBTuple]) -> RGBTuple:
    """
    Get the most common pixel color.
    """
    img = [(p_round(r), p_round(g), p_round(b)) for r, g, b in img]
    c = Counter(img)
    imgx, _ = c.most_common(1)[0]
    return imgx


def slice_to_ints(s: str, m: int) -> Tuple[int, int]:
    """
    Parse slice string.
    """
    assert ":" in s
    ra, rb = s.split(":")
    ra = 0 if ra == "" else int(ra)
    rb = m if rb == "" else int(rb)
    return ra, rb


def convert_background(pngfile: str, new_background: str):
    """
    Replace the background color with the specified background color, default is
    blue.
    """
    if new_background:
        _name, _ext = op.splitext(op.basename(pngfile))
        _name += "_bgxform"
        newfile = op.join(op.dirname(pngfile), _name + _ext)

        img = iopen(pngfile)
        pixels = list(img.getdata())

        # Get Standard Deviation of RGB
        rgb_array = []
        for x in range(255):
            rgb_array.append(x)
        std_rgb = np.std(rgb_array) * 0.8

        # Get average color
        obcolor = [0, 0, 0]
        pixel_values = []
        for t in range(3):
            pixel_color = img.getdata(band=t)
            for pixel in pixel_color:
                if pixel > std_rgb:
                    pixel_values.append(pixel)
            obcolor[t] = sum(pixel_values) // len(pixel_values)

        # Get background color using average color and standard deviation
        for t in range(3):
            pixel_color = img.getdata(band=t)
            seed_pixel_values = []
            for i in pixel_color:
                if obcolor[t] - std_rgb < i < obcolor[t] + std_rgb:
                    seed_pixel_values.append(i)
            obcolor[t] = sum(seed_pixel_values) // len(seed_pixel_values)
        # Selection of colors based on option parser
        nbcolor = [0, 0, 0]
        if new_background == "INVERSE":
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
        for idx, pixel in enumerate(pixels):
            if all(o - std_rgb <= p <= o + std_rgb for o, p in zip(obcolor, pixel)):
                pixels[idx] = nbcolor
        img.putdata(pixels)
        img.save(newfile, "PNG")
        return newfile
    return pngfile


def convert_image(
    pngfile: str,
    pf: str,
    outdir: str = ".",
    resize: int = 1000,
    img_format: str = "jpeg",
    rotate: int = 0,
    rows: str = ":",
    cols: str = ":",
    labelrows: Optional[str] = None,
    labelcols: Optional[str] = None,
) -> Tuple[str, str, Optional[str], dict]:
    """
    Convert image to JPEG format and resize it.
    """
    resizefile = op.join(outdir, pf + ".resize.jpg")
    mainfile = op.join(outdir, pf + ".main.jpg")
    labelfile = op.join(outdir, pf + ".label.jpg")
    img = Image(filename=pngfile)
    exif = dict((k, img.metadata[k]) for k in img.metadata if k.startswith("exif:"))

    # Rotation, slicing and cropping of main image
    if rotate:
        img.rotate(rotate)
    if resize:
        w, h = img.size
        if min(w, h) > resize:
            if w < h:
                nw, nh = resize, resize * h // w
            else:
                nw, nh = resize * w // h, resize
            img.resize(nw, nh)
            logger.debug(
                "Image `%s` resized from %dpx:%dpx to %dpx:%dpx", pngfile, w, h, nw, nh
            )
    img.format = img_format
    img.save(filename=resizefile)

    rimg = img.clone()
    if rows != ":" or cols != ":":
        w, h = img.size
        ra, rb = slice_to_ints(rows, h)
        ca, cb = slice_to_ints(cols, w)
        # left, top, right, bottom
        logger.debug("Crop image to %d:%d %d:%d", ra, rb, ca, cb)
        img.crop(ca, ra, cb, rb)
        img.format = img_format
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
        ra, rb = slice_to_ints(labelrows, h)
        ca, cb = slice_to_ints(labelcols, w)
        logger.debug("Extract label from %d:%d %d:%d", ra, rb, ca, cb)
        rimg.crop(ca, ra, cb, rb)
        rimg.format = img_format
        rimg.save(filename=labelfile)
    else:
        labelfile = None

    return resizefile, mainfile, labelfile, exif


def extract_label(labelfile: str) -> str:
    """
    Extract accession number from label image.
    """
    accession = image_to_string(iopen(labelfile))
    accession = " ".join(accession.split())  # normalize spaces
    accession = "".join(x for x in accession if x in string.printable)
    if not accession:
        accession = "none"
    return accession


def efd_feature(contour: np.ndarray) -> np.ndarray:
    """
    To use EFD as features, one can write a small wrapper function.

    Based on: https://pyefd.readthedocs.io/en/latest
    """
    coeffs = elliptic_fourier_descriptors(contour, normalize=True)
    # skip the first three coefficients, which are always 1, 0, 0
    return coeffs.flatten()[3:]


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

    _, (ax1, ax2, ax3, ax4) = plt.subplots(ncols=4, nrows=1, figsize=(iopts.w, iopts.h))
    # Edge detection
    img_gray = rgb2gray(img)
    logger.debug("Running %s edge detection ...", ff)
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
        logger.debug("Identified %d watershed markers", nmarkers)
        labels = watershed(closed, markers, mask=filled)
    else:
        labels = label(filled)

    # Object size filtering
    w, h = img_gray.shape
    canvas_size = w * h
    min_size = int(round(canvas_size * opts.minsize / 100))
    max_size = int(round(canvas_size * opts.maxsize / 100))
    logger.debug(
        "Find objects with pixels between %d (%d%%) and %d (%d%%)",
        min_size,
        opts.minsize,
        max_size,
        opts.maxsize,
    )

    # Plotting
    ax1.set_title("Original picture")
    ax1.imshow(oimg)

    params = rf"{ff}, $\sigma$={sigma}, $k$={kernel}"
    if opts.watershed:
        params += ", watershed"
    ax2.set_title(f"Edge detection\n({params})")
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
    logger.debug("A total of %d objects identified.", nb_labels)
    objects = []
    for i, props in enumerate(rp):
        i += 1
        if i > opts.count:
            break

        contour = find_contours(labels == props.label, 0.5)[0]
        efds = efd_feature(contour)
        y0, x0 = props.centroid
        orientation = props.orientation
        major, minor = props.major_axis_length, props.minor_axis_length
        major_dx = sin(orientation) * major / 2
        major_dy = cos(orientation) * major / 2
        minor_dx = cos(orientation) * minor / 2
        minor_dy = -sin(orientation) * minor / 2
        ax2.plot((x0 - major_dx, x0 + major_dx), (y0 - major_dy, y0 + major_dy), "r-")
        ax2.plot((x0 - minor_dx, x0 + minor_dx), (y0 - minor_dy, y0 + minor_dy), "r-")
        ax2.plot(contour[:, 1], contour[:, 0], "y-")

        npixels = int(props.area)
        # Sample the center of the blob for color
        d = min(int(round(minor / 2 * 0.35)) + 1, 50)
        x0d, y0d = int(round(x0)), int(round(y0))
        square = img[(y0d - d) : (y0d + d), (x0d - d) : (x0d + d)]
        pixels = []
        for row in square:
            pixels.extend(row)
        logger.debug(
            "Seed #%d: %d pixels (%d sampled) - %.2f%%",
            i,
            npixels,
            len(pixels),
            100.0 * npixels / canvas_size,
        )

        rgb = pixel_stats(pixels)
        objects.append(Seed(filename, accession, i, rgb, props, efds, exif))
        minr, minc, maxr, maxc = props.bbox
        rect = Rectangle(
            (minc, minr), maxc - minc, maxr - minr, fill=False, ec="w", lw=1
        )
        ax3.add_patch(rect)
        mc, mr = (minc + maxc) // 2, (minr + maxr) // 2
        ax3.text(mc, mr, f"{i}", color="w", ha="center", va="center", size=6)

    for ax in (ax2, ax3):
        ax.set_xlim(0, h)
        ax.set_ylim(w, 0)

    # Output identified seed stats
    ax4.text(0.1, 0.92, f"File: {latex(filename)}", color="g")
    ax4.text(0.1, 0.86, f"Label: {latex(accession)}", color="m")
    yy = 0.8
    fw = must_open(opts.outfile, "w")
    if not opts.noheader:
        print(Seed.header(calibrated=calib), file=fw)
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
        f"(A total of {nb_labels} objects displayed)",
        color="darkslategray",
    )
    normalize_axes(ax4)

    for ax in (ax1, ax2, ax3):
        set_helvetica_axis(ax)

    image_name = op.join(outdir, pf + "." + iopts.format)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    return objects


if __name__ == "__main__":
    main()
