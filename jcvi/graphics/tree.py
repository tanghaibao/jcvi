#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import sys
import os.path as op
import logging

from ete3 import Tree

from jcvi.formats.sizes import Sizes
from jcvi.formats.base import DictFile, LineFile
from jcvi.graphics.base import Rectangle, plt, savefig, markup, normalize_axes
from jcvi.graphics.glyph import ExonGlyph, get_setups
from jcvi.apps.base import OptionParser, glob


class LeafInfoLine:
    def __init__(self, row, delimiter=","):
        args = [x.strip() for x in row.split(delimiter)]
        self.name = args[0]
        self.color = args[1]
        self.new_name = None
        if len(args) > 2:
            self.new_name = args[2]


class LeafInfoFile(LineFile):
    def __init__(self, filename, delimiter=","):
        super(LeafInfoFile, self).__init__(filename)
        self.cache = {}
        with open(filename) as fp:
            for row in fp:
                if row[0] == "#":
                    continue
                line = LeafInfoLine(row, delimiter=delimiter)
                self.cache[line.name] = line


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
        tname = name[int(k) :]
    else:
        k = re.search("(?<=^ohead)[0-9]{1,2}$", rule)
        if k:
            k = k.group(0)
            tname = name[: int(k)]
        else:
            k = re.search("(?<=^tail)[0-9]{1,2}$", rule)
            if k:
                k = k.group(0)
                tname = name[: -int(k)]
            else:
                k = re.search("(?<=^otail)[0-9]{1,2}$", rule)
                if k:
                    k = k.group(0)
                    tname = name[-int(k) :]
                else:
                    print(truncate_name.__doc__, file=sys.stderr)
                    raise ValueError("Wrong rule for truncation!")
    return tname


def draw_tree(
    ax,
    tx,
    margin=0.08,
    rmargin=0.2,
    tip=0.01,
    treecolor="k",
    supportcolor="k",
    outgroup=None,
    dashedoutgroup=False,
    reroot=True,
    gffdir=None,
    sizes=None,
    trunc_name=None,
    SH=None,
    scutoff=0,
    leafcolor="k",
    leaffont=12,
    leafinfo=None,
):
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

        if R is not t:
            t.set_outgroup(R)

        # By default, the distance to outgroup and non-outgroup is the same
        # we re-adjust the distances so that the outgroups will appear
        # farthest from everything else
        if dashedoutgroup:
            a, b = t.children
            # Avoid even split
            total = a.dist + b.dist
            newR = t.get_common_ancestor(*outgroup)
            a.dist = 0.9 * total
            b.dist = total - a.dist

    farthest, max_dist = t.get_farthest_leaf()

    xstart = margin
    ystart = 2 * margin
    canvas = 1 - 2 * margin
    # scale the tree
    scale = (canvas - rmargin) / max_dist

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
            yy = ystart + i * yinterval
            i += 1

            if trunc_name:
                name = truncate_name(n.name, rule=trunc_name)
            else:
                name = n.name

            if leafinfo and n.name in leafinfo:
                line = leafinfo[n.name]
                lc = line.color
                sname = line.new_name
            else:
                lc = leafcolor
                sname = None
            lc = lc or "k"
            sname = sname or name.replace("_", "-")
            # if color is given as "R,G,B"
            if "," in lc:
                lc = [float(x) for x in lc.split(",")]

            ax.text(
                xx + tip,
                yy,
                markup(sname),
                va="center",
                fontstyle="italic",
                size=leaffont,
                color=lc,
            )

            gname = n.name.split("_")[0]
            if gname in structures:
                mrnabed, cdsbeds = structures[gname]
                ExonGlyph(
                    ax,
                    1 - rmargin / 2,
                    yy,
                    mrnabed,
                    cdsbeds,
                    align="right",
                    ratio=ratio,
                )
            if sizes and gname in sizes:
                size = sizes[gname]
                size = size / 3 - 1  # base pair converted to amino acid
                size = "{0}aa".format(size)
                ax.text(1 - rmargin / 2 + tip, yy, size, size=leaffont)

        else:
            linestyle = "--" if (dashedoutgroup and n is t) else "-"
            children = [coords[x] for x in n.get_children()]
            children_x, children_y = zip(*children)
            min_y, max_y = min(children_y), max(children_y)
            # plot the vertical bar
            ax.plot((xx, xx), (min_y, max_y), linestyle, color=treecolor)
            # plot the horizontal bar
            for cx, cy in children:
                ax.plot((xx, cx), (cy, cy), linestyle, color=treecolor)
            yy = sum(children_y) * 1.0 / len(children_y)
            support = n.support
            if support > 1:
                support = support / 100.0
            if not n.is_root() and supportcolor:
                if support > scutoff / 100.0:
                    ax.text(
                        xx,
                        yy + 0.005,
                        "{0:d}".format(int(abs(support * 100))),
                        ha="right",
                        size=leaffont,
                        color=supportcolor,
                    )

        coords[n] = (xx, yy)

    # scale bar
    br = 0.1
    x1 = xstart + 0.1
    x2 = x1 + br * scale
    yy = margin
    ax.plot([x1, x1], [yy - tip, yy + tip], "-", color=treecolor)
    ax.plot([x2, x2], [yy - tip, yy + tip], "-", color=treecolor)
    ax.plot([x1, x2], [yy, yy], "-", color=treecolor)
    ax.text(
        (x1 + x2) / 2,
        yy - tip,
        "{0:g}".format(br),
        va="top",
        ha="center",
        size=leaffont,
        color=treecolor,
    )

    if SH is not None:
        xs = x1
        ys = (margin + yy) / 2.0
        ax.text(
            xs,
            ys,
            "SH test against ref tree: {0}".format(SH),
            ha="left",
            size=leaffont,
            color="g",
        )

    normalize_axes(ax)


def read_trees(tree):
    from six.moves.urllib.parse import parse_qs
    from jcvi.formats.base import read_block

    trees = []

    fp = open(tree)
    for header, tx in read_block(fp, "#"):
        header = parse_qs(header[1:])
        label = header["label"][0].strip('"')
        outgroup = header["outgroup"]
        color, = header.get("color", ["k"])
        trees.append((label, outgroup, color, "".join(tx)))

    return trees


def draw_geoscale(ax, minx=0, maxx=175):
    """
    Draw geological epoch on million year ago (mya) scale.
    """
    a, b = 0.1, 0.6  # Correspond to 200mya and 0mya

    def cv(x):
        return b - (x - b) / (maxx - minx) * (b - a)

    ax.plot((a, b), (0.5, 0.5), "k-")
    tick = 0.015
    for mya in xrange(maxx - 25, 0, -25):
        p = cv(mya)
        ax.plot((p, p), (0.5, 0.5 - tick), "k-")
        ax.text(p, 0.5 - 2.5 * tick, str(mya), ha="center", va="center")
    ax.text(
        (a + b) / 2,
        0.5 - 5 * tick,
        "Time before present (million years)",
        ha="center",
        va="center",
    )

    # Source:
    # http://www.weston.org/schools/ms/biologyweb/evolution/handouts/GSAchron09.jpg
    Geo = (
        ("Neogene", 2.6, 23.0, "#fee400"),
        ("Paleogene", 23.0, 65.5, "#ff9a65"),
        ("Cretaceous", 65.5, 145.5, "#80ff40"),
        ("Jurassic", 145.5, 201.6, "#33fff3"),
    )
    h = 0.05
    for era, start, end, color in Geo:
        start, end = cv(start), cv(end)
        end = max(a, end)
        p = Rectangle(
            (end, 0.5 + tick / 2), abs(start - end), h, lw=1, ec="w", fc=color
        )
        ax.text(
            (start + end) / 2,
            0.5 + (tick + h) / 2,
            era,
            ha="center",
            va="center",
            size=9,
        )
        ax.add_patch(p)


def main(args):
    """
    %prog newicktree

    Plot Newick formatted tree. The gene structure can be plotted along if
    --gffdir is given. The gff file needs to be `genename.gff`. If --sizes is
    on, also show the number of amino acids.
    """
    p = OptionParser(main.__doc__)
    p.add_option(
        "--outgroup",
        help="Outgroup for rerooting the tree. "
        + "Use comma to separate multiple taxa.",
    )
    p.add_option(
        "--noreroot",
        default=False,
        action="store_true",
        help="Don't reroot the input tree",
    )
    p.add_option(
        "--rmargin", default=0.2, type="float", help="Set blank rmargin to the right"
    )
    p.add_option("--gffdir", default=None, help="The directory that contain GFF files")
    p.add_option("--sizes", default=None, help="The FASTA file or the sizes file")
    p.add_option("--SH", default=None, type="string", help="SH test p-value")
    p.add_option(
        "--geoscale", default=False, action="store_true", help="Plot geological scale"
    )

    group = p.add_option_group("Node style")
    group.add_option("--leafcolor", default="k", help="Font color for the OTUs")
    group.add_option("--leaffont", default=12, help="Font size for the OTUs")
    group.add_option(
        "--leafinfo", help="CSV file specifying the leaves: name,color,new_name"
    )
    group.add_option(
        "--scutoff",
        default=0,
        type="int",
        help="cutoff for displaying node support, 0-100",
    )
    group.add_option(
        "--no_support",
        dest="support",
        default=True,
        action="store_false",
        help="Do not print node support values",
    )

    group = p.add_option_group("Edge style")
    group.add_option(
        "--dashedoutgroup",
        default=False,
        action="store_true",
        help="Gray out the edges connecting outgroup and non-outgroup",
    )

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
        logging.debug("Load tree file `{0}`".format(datafile))
        tx = open(datafile).read()

    pf = datafile.rsplit(".", 1)[0]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    supportcolor = "k" if opts.support else None

    if opts.geoscale:
        draw_geoscale(root)

    else:
        leafinfo = LeafInfoFile(opts.leafinfo).cache if opts.leafinfo else None

        draw_tree(
            root,
            tx,
            rmargin=opts.rmargin,
            supportcolor=supportcolor,
            outgroup=outgroup,
            dashedoutgroup=opts.dashedoutgroup,
            reroot=reroot,
            gffdir=opts.gffdir,
            sizes=opts.sizes,
            SH=opts.SH,
            scutoff=opts.scutoff,
            leafcolor=opts.leafcolor,
            leaffont=opts.leaffont,
            leafinfo=leafinfo,
        )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main(sys.argv[1:])
