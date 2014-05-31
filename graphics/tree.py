#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import os.path as op
import logging

from ete2 import Tree

from jcvi.formats.sizes import Sizes
from jcvi.formats.base import DictFile
from jcvi.graphics.base import plt, savefig
from jcvi.graphics.glyph import ExonGlyph, get_setups
from jcvi.apps.base import OptionParser, glob


def truncate_name(name, rule=None):
    """
    shorten taxa names for tree display

    Options of rule. This only affects tree display.
    - headn (eg. head3 truncates first 3 chars)
    - oheadn (eg. ohead3 retains only the first 3 chars)
    - tailn (eg. tail3 truncates last 3 chars)
    - otailn (eg. otail3 retains only the last 3 chars)
    n = 1 ~ 99
    """
    import re

    if rule is None:
        return name

    k = re.search("(?<=^head)[0-9]{1,2}$", rule)
    if k:
        k = k.group(0)
        tname = name[int(k):]
    else:
        k = re.search("(?<=^ohead)[0-9]{1,2}$", rule)
        if k:
            k = k.group(0)
            tname = name[:int(k)]
        else:
            k = re.search("(?<=^tail)[0-9]{1,2}$", rule)
            if k:
                k = k.group(0)
                tname = name[:-int(k)]
            else:
                k = re.search("(?<=^otail)[0-9]{1,2}$", rule)
                if k:
                    k = k.group(0)
                    tname = name[-int(k):]
                else:
                    print >>sys.stderr, truncate_name.__doc__
                    raise ValueError('Wrong rule for truncation!')
    return tname


def decode_name(name, barcodemap):
    """
    rename seq/taxon name, typically for a tree display,
    according to a barcode map given in a dictionary

    By definition barcodes should be distinctive.
    """
    for barcode in barcodemap:
        if barcode in name:
            return barcodemap[barcode]

    return name


def draw_tree(ax, tx, rmargin=.3, leafcolor="k", supportcolor="k",
              outgroup=None, reroot=True, gffdir=None, sizes=None,
              trunc_name=None, SH=None, scutoff=0, barcodefile=None,
              leafcolorfile=None, leaffont=12):
    """
    main function for drawing phylogenetic tree
    """

    t = Tree(tx)
    if reroot:
        if outgroup:
            R = t.get_common_ancestor(*outgroup)
        else:
            # Calculate the midpoint node
            R = t.get_midpoint_outgroup()

        if R != t:
            t.set_outgroup(R)

    farthest, max_dist = t.get_farthest_leaf()

    margin = .05
    xstart = margin
    ystart = 1 - margin
    canvas = 1 - rmargin - 2 * margin
    tip = .005
    # scale the tree
    scale = canvas / max_dist

    num_leaves = len(t.get_leaf_names())
    yinterval = canvas / (num_leaves + 1)

    # get exons structures, if any
    structures = {}
    if gffdir:
        gffiles = glob("{0}/*.gff*".format(gffdir))
        setups, ratio = get_setups(gffiles, canvas=rmargin / 2, noUTR=True)
        structures = dict((a, (b, c)) for a, b, c in setups)

    if sizes:
        sizes = Sizes(sizes).mapping

    if barcodefile:
        barcodemap = DictFile(barcodefile, delimiter="\t")

    if leafcolorfile:
        leafcolors = DictFile(leafcolorfile, delimiter="\t")

    coords = {}
    i = 0
    for n in t.traverse("postorder"):
        dist = n.get_distance(t)
        xx = xstart + scale * dist

        if n.is_leaf():
            yy = ystart - i * yinterval
            i += 1

            if trunc_name:
                name = truncate_name(n.name, rule=trunc_name)
            else:
                name = n.name

            if barcodefile:
                name = decode_name(name, barcodemap)

            sname = name.replace("_", "-")

            try:
                lc = leafcolors[n.name]
            except Exception:
                lc = leafcolor
            else:
                # if color is given as "R,G,B"
                if "," in lc:
                    lc = map(float, lc.split(","))

            ax.text(xx + tip, yy, sname, va="center",
                    fontstyle="italic", size=leaffont, color=lc)

            gname = n.name.split("_")[0]
            if gname in structures:
                mrnabed, cdsbeds = structures[gname]
                ExonGlyph(ax, 1 - rmargin / 2, yy, mrnabed, cdsbeds,
                          align="right", ratio=ratio)
            if sizes and gname in sizes:
                size = sizes[gname]
                size = size / 3 - 1  # base pair converted to amino acid
                size = "{0}aa".format(size)
                ax.text(1 - rmargin / 2 + tip, yy, size, size=leaffont)

        else:
            children = [coords[x] for x in n.get_children()]
            children_x, children_y = zip(*children)
            min_y, max_y = min(children_y), max(children_y)
            # plot the vertical bar
            ax.plot((xx, xx), (min_y, max_y), "k-")
            # plot the horizontal bar
            for cx, cy in children:
                ax.plot((xx, cx), (cy, cy), "k-")
            yy = sum(children_y) * 1. / len(children_y)
            support = n.support
            if support > 1:
                support = support / 100.
            if not n.is_root():
                if support > scutoff / 100.:
                    ax.text(xx, yy+.005, "{0:d}".format(int(abs(support * 100))),
                        ha="right", size=leaffont, color=supportcolor)

        coords[n] = (xx, yy)

    # scale bar
    br = .1
    x1 = xstart + .1
    x2 = x1 + br * scale
    yy = ystart - i * yinterval
    ax.plot([x1, x1], [yy - tip, yy + tip], "k-")
    ax.plot([x2, x2], [yy - tip, yy + tip], "k-")
    ax.plot([x1, x2], [yy, yy], "k-")
    ax.text((x1 + x2) / 2, yy - tip, "{0:g}".format(br),
            va="top", ha="center", size=leaffont)

    if SH is not None:
        xs = x1
        ys = (margin + yy) / 2.
        ax.text(xs, ys, "SH test against ref tree: {0}"\
                .format(SH), ha="left", size=leaffont, color="g")

    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()


def read_trees(tree):
    from urlparse import parse_qs
    from jcvi.formats.base import read_block

    trees = []

    fp = open(tree)
    for header, tx in read_block(fp, "#"):
        header = parse_qs(header[1:])
        label = header["label"][0].strip("\"")
        outgroup = header["outgroup"]
        trees.append((label, outgroup, "".join(tx)))

    return trees


def main(args):
    """
    %prog newicktree

    Plot Newick formatted tree. The gene structure can be plotted along if
    --gffdir is given. The gff file needs to be `genename.gff`. If --sizes is
    on, also show the number of amino acids.

    With --barcode a mapping file can be provided to convert seq names to
    eg. species names, useful in unified tree display. This file should have
    distinctive barcodes in column1 and new names in column2, tab delimited.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--outgroup", help="Outgroup for rerooting the tree. " + \
                 "Use comma to separate multiple taxa.")
    p.add_option("--noreroot", default=False, action="store_true", \
                 help="Don't reroot the input tree [default: %default]")
    p.add_option("--rmargin", default=.3, type="float",
                 help="Set blank rmargin to the right [default: %default]")
    p.add_option("--gffdir", default=None,
                 help="The directory that contain GFF files [default: %default]")
    p.add_option("--sizes", default=None,
                 help="The FASTA file or the sizes file [default: %default]")
    p.add_option("--SH", default=None, type="string",
                 help="SH test p-value [default: %default]")
    p.add_option("--scutoff", default=0, type="int",
                 help="cutoff for displaying node support, 0-100 [default: %default]")
    p.add_option("--barcode", default=None,
                 help="path to seq names barcode mapping file: " \
                 "barcode<tab>new_name [default: %default]")
    p.add_option("--leafcolor", default="k",
                 help="Font color for the OTUs, or path to a file " \
                 "containing color mappings: leafname<tab>color [default: %default]")
    p.add_option("--leaffont", default=12,
                 help="Font size for the OTUs")

    opts, args, iopts = p.set_image_options(args, figsize="8x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    outgroup = None
    reroot = not opts.noreroot
    if opts.outgroup:
        outgroup = opts.outgroup.split(",")

    if datafile == "demo":
        tx = """(((Os02g0681100:0.1151,Sb04g031800:0.11220)1.0:0.0537,
        (Os04g0578800:0.04318,Sb06g026210:0.04798)-1.0:0.08870)1.0:0.06985,
        ((Os03g0124100:0.08845,Sb01g048930:0.09055)1.0:0.05332,
        (Os10g0534700:0.06592,Sb01g030630:0.04824)-1.0:0.07886):0.09389);"""
    else:
        logging.debug("Load tree file `{0}`.".format(datafile))
        tx = open(datafile).read()

    pf = datafile.rsplit(".", 1)[0]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    if op.isfile(opts.leafcolor):
        leafcolor = "k"
        leafcolorfile = opts.leafcolor
    else:
        leafcolor = opts.leafcolor
        leafcolorfile = None

    draw_tree(root, tx, rmargin=opts.rmargin, leafcolor=leafcolor, \
              outgroup=outgroup, reroot=reroot, gffdir=opts.gffdir, \
              sizes=opts.sizes, SH=opts.SH, scutoff=opts.scutoff, \
              barcodefile=opts.barcode, leafcolorfile=leafcolorfile,
              leaffont=opts.leaffont)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main(sys.argv[1:])
