#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import logging
import sys

from collections import defaultdict
from itertools import groupby

from ete3 import Tree

from jcvi.apps.base import OptionParser, glob
from jcvi.formats.base import LineFile
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import (
    FancyBboxPatch,
    Rectangle,
    linear_shade,
    markup,
    normalize_axes,
    plt,
    savefig,
    set3_n,
)
from jcvi.graphics.glyph import ExonGlyph, TextCircle, get_setups


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


class WGDInfoLine:
    def __init__(self, row, delimiter=",", defaultcolor="#7fc97f"):
        args = [x.strip() for x in row.split(delimiter)]
        self.node_name = args[0]
        self.divergence = float(args[1]) / 100
        self.name = args[2]
        self.color = args[3] or defaultcolor
        self.style = args[4]


class WGDInfoFile(LineFile):
    def __init__(self, filename, delimiter=","):
        super(WGDInfoFile, self).__init__(filename)
        self.cache = defaultdict(list)
        with open(filename) as fp:
            for row in fp:
                if row[0] == "#":
                    continue
                line = WGDInfoLine(row, delimiter=delimiter)
                self.cache[line.node_name].append(line)


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


def draw_wgd_xy(ax, xx, yy, wgdline):
    """Draw WGD at (xx, yy) position

    Args:
        ax (axis): Matplotlib axes
        xx (float): x position
        yy (float): y position
        wgdline (WGDInfo): WGDInfoLines that contains the styling information
    """
    TextCircle(
        ax,
        xx,
        yy,
        wgdline.name,
        fc=wgdline.color,
        radius=0.0225,
        color="k",
        fontweight="bold",
    )


def draw_wgd(ax, y, rescale, name, wgdcache):
    """Draw WGD given a name and the WGDInfo cache.

    Args:
        ax (matplotlib.axes): matplotlib axes
        y (float): y position
        rescale (function): Rescale function to generate x position
        name (str): Name of the line (usually the taxon/internal name)
        wgdcache (Dict): Dictionary containing WGDInfoLines
    """
    if not wgdcache or name not in wgdcache:
        return
    for line in wgdcache[name]:
        draw_wgd_xy(ax, rescale(line.divergence), y, line)


def draw_tree(
    ax,
    t,
    hpd=None,
    margin=0.1,
    rmargin=0.2,
    ymargin=0.1,
    tip=0.01,
    treecolor="k",
    supportcolor="k",
    internal=True,
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
    wgdinfo=None,
    geoscale=False,
    groups=[],
):
    """
    main function for drawing phylogenetic tree
    """

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
    print("max_dist = {}".format(max_dist), file=sys.stderr)

    xstart = margin
    ystart = 2 * ymargin
    # scale the tree
    scale = (1 - margin - rmargin) / max_dist

    def rescale(dist):
        return xstart + scale * dist

    def rescale_divergence(divergence):
        return rescale(max_dist - divergence)

    num_leaves = len(t.get_leaf_names())
    yinterval = (1 - ystart) / num_leaves
    ytop = ystart + (num_leaves - 0.5) * yinterval

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
    color_groups = []  # Used to plot groups to the right of the tree
    for n in t.traverse("postorder"):
        dist = n.get_distance(t)
        xx = rescale(dist)

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
            color_groups.append((lc, yy, xx))

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
            # plot HPD if exists
            if hpd and n.name in hpd:
                a, b = hpd[n.name]
                ax.plot(
                    (rescale_divergence(a), rescale_divergence(b)),
                    (yy, yy),
                    "-",
                    color="darkslategray",
                    alpha=0.4,
                    lw=2,
                )
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
            if internal and n.name:
                TextCircle(ax, xx, yy, n.name, size=9)
            else:  # Just a dot
                TextCircle(ax, xx, yy, None, radius=0.002)

        coords[n] = (xx, yy)
        # WGD info
        draw_wgd(ax, yy, rescale_divergence, n.name, wgdinfo)

    # scale bar
    if geoscale:
        draw_geoscale(
            ax, ytop, margin=margin, rmargin=rmargin, yy=ymargin, max_dist=max_dist
        )
    else:
        br = 0.1
        x1 = xstart + 0.1
        x2 = x1 + br * scale
        yy = ymargin
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

    # Groupings on the right, often to used to show groups such as phylogenetic
    # clades
    if groups:
        color_groups.sort()
        group_extents = []
        for color, group in groupby(color_groups, key=lambda x: x[0]):
            group = list(group)
            _, min_yy, xx = min(group)
            _, max_yy, xx = max(group)
            group_extents.append((min_yy, max_yy, xx, color))
        group_extents.sort(reverse=True)

        for group_name, (min_yy, max_yy, xx, color) in zip(groups, group_extents):
            group_color = linear_shade(color, fraction=0.85)
            ax.add_patch(
                FancyBboxPatch(
                    (xx, min_yy - yinterval / 2),
                    rmargin - 0.01,
                    max_yy - min_yy + yinterval,
                    boxstyle="round,pad=-0.002,rounding_size=0.005",
                    fc=group_color,
                    ec=group_color,
                )
            )
            # Add the group label
            horizontal = (max_yy - min_yy) < 0.2
            mid_yy = (min_yy + max_yy) / 2
            label_rightend = 0.98
            if horizontal:
                ax.text(
                    label_rightend,
                    mid_yy,
                    markup(group_name),
                    color="darkslategray",
                    ha="right",
                    va="center",
                )
            else:
                ax.text(
                    label_rightend,
                    mid_yy,
                    markup(group_name),
                    color="darkslategray",
                    ha="right",
                    va="center",
                    rotation=-90,
                )

    if SH is not None:
        xs = x1
        ys = (ymargin + yy) / 2.0
        ax.text(
            xs,
            ys,
            "SH test against ref tree: {0}".format(SH),
            ha="left",
            size=leaffont,
            color="g",
        )


def read_trees(tree):
    from urllib.parse import parse_qs
    from jcvi.formats.base import read_block

    trees = []

    fp = open(tree)
    for header, tx in read_block(fp, "#"):
        header = parse_qs(header[1:])
        label = header["label"][0].strip('"')
        outgroup = header["outgroup"]
        (color,) = header.get("color", ["k"])
        trees.append((label, outgroup, color, "".join(tx)))

    return trees


def draw_geoscale(
    ax, ytop, margin=0.1, rmargin=0.2, yy=0.1, max_dist=3.0, contrast_epochs=True
):
    """
    Draw geological epoch on million year ago (mya) scale.
    max_dist = 3.0 => max is 300 mya
    """
    import math

    a, b = margin, 1 - rmargin  # Correspond to 300mya and 0mya
    minx, maxx = 0, int(max_dist * 100)

    def cv(x):
        return b - (x - b) / (maxx - minx) * (b - a)

    ax.plot((a, b), (yy, yy), "k-")
    tick = 0.0125
    scale_start = int(math.ceil(maxx / 25) * 25)
    for mya in range(scale_start - 25, 0, -25):
        p = cv(mya)
        ax.plot((p, p), (yy, yy - tick), "k-")
        ax.text(p, yy - 2.5 * tick, str(mya), ha="center", va="center")

    ax.text(
        (a + b) / 2,
        yy - 5 * tick,
        "Time before present (million years)",
        ha="center",
        va="center",
    )

    # Source:
    # https://en.wikipedia.org/wiki/Geological_period
    Geo = (
        ("Neogene", 2.588, 23.03),
        ("Paleogene", 23.03, 66.0),
        ("Cretaceous", 66.0, 145.5),
        ("Jurassic", 145.5, 201.3),
        ("Triassic", 201.3, 252.17),
        ("Permian", 252.17, 298.9),
        ("Carboniferous", 298.9, 358.9),
    )
    h = 0.05
    for (era, start, end), color in zip(Geo, set3_n(len(Geo))):
        if maxx - start < 10:  # not visible enough
            continue
        start, end = cv(start), cv(end)
        end = max(a, end)
        p = Rectangle((end, yy + tick / 2), abs(start - end), h, lw=1, ec="w", fc=color)
        ax.text(
            (start + end) / 2,
            yy + (tick + h) / 2,
            era,
            ha="center",
            va="center",
            size=8,
        )
        ax.add_patch(p)

    # We highlight recent epochs for better visualization, we just highlight
    # Neogene and Cretaceous as these are more relevant for most phylogeny
    if contrast_epochs:
        for era, start, end in Geo:
            if not era in ("Neogene", "Cretaceous"):
                continue

            # Make a beige patch
            start, end = cv(start), cv(end)
            ax.add_patch(
                Rectangle(
                    (end, yy + tick + h),
                    abs(start - end),
                    ytop - yy - tick - h,
                    fc="beige",
                    ec="beige",
                )
            )


def parse_tree(infile):
    """Parse newick formatted tree file and returns a tuple consisted of a
    Tree object, and a HPD dictionary if 95%HPD is found in the newick string,
    otherwise None

    Args:
        infile (str): Path to the tree file
    """
    import re

    with open(infile) as fp:
        treedata = fp.read()
    hpd_re = re.compile(r"( \[&95%HPD=[^[]*\])")

    def repl(match):
        repl.count += 1
        name = "N{}".format(repl.count)
        lb, ub = re.findall(r"HPD=\{(.*), (.*)\}", match.group(0))[0]
        repl.hpd[name] = (float(lb), float(ub))
        return name

    repl.count = 0
    repl.hpd = {}

    treedata, changed = re.subn(hpd_re, repl, treedata)
    if repl.hpd:
        print(repl.hpd, file=sys.stderr)

    return (Tree(treedata, format=1), repl.hpd) if changed else (Tree(treedata), None)


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
    group.add_option(
        "--no_internal",
        dest="internal",
        default=True,
        action="store_false",
        help="Do not show internal nodes",
    )

    group = p.add_option_group("Edge style")
    group.add_option(
        "--dashedoutgroup",
        default=False,
        action="store_true",
        help="Gray out the edges connecting outgroup and non-outgroup",
    )

    group = p.add_option_group("Additional annotations")
    group.add_option(
        "--geoscale",
        default=False,
        action="store_true",
        help="Plot geological scale",
    )
    group.add_option(
        "--wgdinfo", help="CSV specifying the position and style of WGD events"
    )
    group.add_option(
        "--groups",
        help="Group names from top to bottom, to the right of the tree. "
        "Each distinct color in --leafinfo is considered part of the same group. "
        "Separate the names with comma, such as 'eudicots,,monocots,'. "
        "Empty names will be ignored for that specific group. ",
    )

    opts, args, iopts = p.set_image_options(args, figsize="10x7")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (datafile,) = args
    outgroup = None
    reroot = not opts.noreroot
    if opts.outgroup:
        outgroup = opts.outgroup.split(",")

    hpd = None
    if datafile == "demo":
        t = Tree(
            """(((Os02g0681100:0.1151,Sb04g031800:0.11220)1.0:0.0537,
        (Os04g0578800:0.04318,Sb06g026210:0.04798)-1.0:0.08870)1.0:0.06985,
        ((Os03g0124100:0.08845,Sb01g048930:0.09055)1.0:0.05332,
        (Os10g0534700:0.06592,Sb01g030630:0.04824)-1.0:0.07886):0.09389);"""
        )
    else:
        logging.debug("Load tree file `{0}`".format(datafile))
        t, hpd = parse_tree(datafile)

    pf = datafile.rsplit(".", 1)[0]

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    supportcolor = "k" if opts.support else None
    margin, rmargin = 0.1, opts.rmargin  # Left and right margin
    leafinfo = LeafInfoFile(opts.leafinfo).cache if opts.leafinfo else None
    wgdinfo = WGDInfoFile(opts.wgdinfo).cache if opts.wgdinfo else None

    draw_tree(
        root,
        t,
        hpd=hpd,
        margin=margin,
        rmargin=rmargin,
        ymargin=margin,
        supportcolor=supportcolor,
        internal=opts.internal,
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
        wgdinfo=wgdinfo,
        geoscale=opts.geoscale,
        groups=opts.groups.split(",") if opts.groups else [],
    )

    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main(sys.argv[1:])
