#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import re
import logging

from glob import glob
from optparse import OptionParser

from ete2 import Tree
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, _, set_image_options, savefig
from jcvi.graphics.glyph import ExonGlyph, get_setups
from jcvi.apps.base import debug
debug()


def draw_tree(ax, tx, rmargin=.3, leafcolor="k", supportcolor="k",
              outgroup=None, reroot=True, gffdir=None, sizes=None,
              trunc_name=None):
    """
    main function for drawing phylogenetic tree

    Options for truncating seq names are below. This only affects display.
    - headn (eg. head3 means truncating first 3 chars)
    - oheadn (eg. ohead3 means only retain first 3 chars)
    - tailn (eg. tail3 means truncating last 3 chars)
    - otailn (eg. otail3 means only retain last 3 chars)
    n = 1 ~ 99
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

    margin = .05 if len(farthest.name)<10 else .1
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

    coords = {}
    i = 0
    for n in t.traverse("postorder"):
        dist = n.get_distance(t)
        xx = xstart + scale * dist

        if n.is_leaf():
            yy = ystart - i * yinterval
            i += 1

            if trunc_name:
                k = re.search("(?<=^head)[0-9]{1,2}$", trunc_name).group(0)
                if k:
                    name = n.name[int(k):]
                else:
                    k = re.search("(?<=^ohead)[0-9]{1,2}$", \
                        trunc_name).group(0)
                    if k:
                        name = n.name[:int(k)]
                    else:
                        k = re.search("(?<=^tail)[0-9]{1,2}$", \
                            trunc_name).group(0)
                        if k:
                            name = n.name[:-int(k)]
                        else:
                            k = re.search("(?<=^otail)[0-9]{1,2}$", \
                                trunc_name).group(0)
                            if k:
                                name = n.name[-int(k):]
                            else:
                                raise ValueError('Wrong trunc_name option.')
            else:
                name = n.name
            ax.text(xx + tip, yy, name.replace("_","-"), va="center",
                    fontstyle="italic", size=8, color=leafcolor)

            gname = n.name.split("_")[0]
            if gname in structures:
                mrnabed, cdsbeds = structures[gname]
                ExonGlyph(ax, 1 - rmargin / 2, yy, mrnabed, cdsbeds,
                          align="right", ratio=ratio)
            if sizes and gname in sizes:
                size = sizes[gname]
                size = size / 3 - 1  # base pair converted to amino acid
                size = _("{0}aa".format(size))
                ax.text(1 - rmargin / 2 + tip, yy, size)

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
                support = support/100.
            if not n.is_root():
                ax.text(xx, yy, _("{0:d}".format(int(abs(support * 100)))),
                        ha="right", size=8, color=supportcolor)

        coords[n] = (xx, yy)

    # scale bar
    br = .1
    x1 = xstart + .1
    x2 = x1 + br * scale
    yy = ystart - i * yinterval
    ax.plot([x1, x1], [yy - tip, yy + tip], "k-")
    ax.plot([x2, x2], [yy - tip, yy + tip], "k-")
    ax.plot([x1, x2], [yy, yy], "k-")
    ax.text((x1 + x2) / 2, yy - tip, _("{0:g}".format(br)),
            va="top", ha="center", size=10)

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


def main():
    """
    %prog newicktree

    Plot Newick formatted tree. The gene structure can be plotted along if
    --gffdir is given. The gff file needs to be `genename.gff`. If --sizes is
    on, also show the number of amino acids.
    """
    p = OptionParser(main.__doc__)
    p.add_option("--outgroup", help="Outgroup for rerooting the tree. " + \
                 "Use comma to separate multiple taxa.")
    p.add_option("--reroot", action="store_true", \
                 help="Reroot the input tree [default: %default]")
    p.add_option("--rmargin", default=.3, type="float",
                 help="Set blank rmargin to the right [default: %default]")
    p.add_option("--gffdir", default=None,
                 help="The directory that contain GFF files [default: %default]")
    p.add_option("--sizes", default=None,
                 help="The FASTA file or the sizes file [default: %default]")
    p.add_option("--leafcolor", default="k",
                 help="The font color for the OTUs [default: %default]")

    opts, args, iopts = set_image_options(p, figsize="8x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    outgroup = None
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

    draw_tree(root, tx, rmargin=opts.rmargin, leafcolor=opts.leafcolor, \
              outgroup=outgroup, reroot=opts.reroot, gffdir=opts.gffdir, \
              sizes=opts.sizes)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
