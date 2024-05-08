#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog mcscan.txt all.bed layout.csv

Illustrate MCscan multiple collinearity alignments. Use layout.csv to indicate
the positions of tracks. For example:

#x, y, rotation, ha, va, color, ratio
0.5, 0.6, 0, left, center, g
0.25, 0.7, 45, center, center, m

With the row ordering corresponding to the column ordering in the MCscan output.

For "ha" (horizontal alignment), accepted values are: left|right|leftalign|rightalign|center|""
For "va" (vertical alignment), accepted values are: top|bottom|center|""(empty)
"""

import sys

from typing import List, Optional

import numpy as np

from matplotlib import transforms
from matplotlib.path import Path

from ..apps.base import OptionParser, logger
from ..compara.synteny import BlockFile
from ..formats.base import DictFile
from ..formats.bed import Bed
from ..utils.cbook import human_size
from ..utils.validator import validate_in_choices, validate_in_range

from .base import (
    AbstractLayout,
    PathPatch,
    markup,
    plt,
    savefig,
)
from .glyph import (
    BasePalette,
    Glyph,
    OrientationPalette,
    OrthoGroupPalette,
    RoundLabel,
)
from .tree import draw_tree, read_trees


HorizontalAlignments = ("left", "right", "leftalign", "rightalign", "center", "")
VerticalAlignments = ("top", "bottom", "center", "")
CANVAS_SIZE = 0.65


class LayoutLine(object):
    """
    Parse a line in the layout file. The line is in the following format:

    *0.5, 0.6, 0, left, center, g, 1, chr1
    """

    def __init__(self, row, delimiter=","):
        self.hidden = row[0] == "*"
        if self.hidden:
            row = row[1:]
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.x = float(args[0])
        validate_in_range(self.x, 0, 1, "XPosition(x) column")
        self.y = float(args[1])
        validate_in_range(self.y, 0, 1, "YPosition(y) column")
        self.rotation = int(args[2])
        self.ha = args[3]
        validate_in_choices(
            self.ha, HorizontalAlignments, "HorizontaAlignment(ha) column"
        )
        self.va = args[4]
        validate_in_choices(self.va, VerticalAlignments, "VerticalAlignment(va) column")
        self.color = args[5]
        self.ratio = 1
        if len(args) > 6:
            self.ratio = float(args[6])
        if len(args) > 7:
            self.label = args[7].strip()
        else:
            self.label = None
        if len(args) > 8:
            self.label_fontsize = float(args[8])
        else:
            self.label_fontsize = 10


class Layout(AbstractLayout):
    """
    Parse the layout file.
    """

    def __init__(self, filename, delimiter=",", seed: Optional[int] = None):
        super(Layout, self).__init__(filename)
        fp = open(filename, encoding="utf-8")
        self.edges = []
        for row in fp:
            if row[0] == "#":
                continue
            if row[0] == "e":
                args = row.rstrip().split(delimiter)
                args = [x.strip() for x in args]
                a, b = args[1:3]
                if len(args) >= 4 and args[3]:
                    blockcolor = args[3]
                else:
                    blockcolor = None
                if len(args) >= 5 and args[4]:
                    samearc = args[4]
                else:
                    samearc = None
                a, b = int(a), int(b)
                assert args[0] == "e"
                self.edges.append((a, b, blockcolor, samearc))
            else:
                self.append(LayoutLine(row, delimiter=delimiter))

        self.assign_colors(seed=seed)


class Shade(object):
    """
    Draw a shade between two tracks.
    """

    Styles = ("curve", "line")

    def __init__(
        self,
        ax,
        a,
        b,
        ymid_pad: float = 0.0,
        highlight=False,
        style="curve",
        ec="k",
        fc="k",
        alpha=0.2,
        lw=1,
        zorder=1,
    ):
        """Create syntenic wedges between tracks.

        Args:
            ax: matplotlib Axes
            a (tuple of floats): ((start_x, start_y), (end_x, end_y))
            b (tuple of floats): ((start_x, start_y), (end_x, end_y))
            ymid_pad (float): Adjustment to y-mid position of Bezier controls, curve style only
            highlight (bool, optional): Plot this shade if color is specified. Defaults to False.
            style (str, optional): Style. Defaults to "curve", must be one of
            ("curve", "line")
            ec (str, optional): Edge color. Defaults to "k".
            fc (str, optional): Face color. Defaults to "k".
            alpha (float, optional): Transparency. Defaults to 0.2.
            lw (int, optional): Line width. Defaults to 1.
            zorder (int, optional): Z-order. Defaults to 1.
        """
        fc = fc or "gainsboro"  # Default block color is grayish
        assert style in self.Styles, f"style must be one of {self.Styles}"
        a1, a2 = a
        b1, b2 = b
        ax1, ay1 = a1
        ax2, ay2 = a2
        bx1, by1 = b1
        bx2, by2 = b2
        if ax1 is None or ax2 is None or bx1 is None or bx2 is None:
            logger.warning("Shade: None found in coordinates, skipping")
            return
        M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
        if style == "curve":
            ymid1 = (ay1 + by1) / 2 + ymid_pad
            ymid2 = (ay2 + by2) / 2 + ymid_pad
            pathdata = [
                (M, a1),
                (C4, (ax1, ymid1)),
                (C4, (bx1, ymid1)),
                (C4, b1),
                (L, b2),
                (C4, (bx2, ymid2)),
                (C4, (ax2, ymid2)),
                (C4, a2),
                (CP, a1),
            ]
        else:
            pathdata = [(M, a1), (L, b1), (L, b2), (L, a2), (CP, a1)]
        codes, verts = zip(*pathdata)
        path = Path(verts, codes)
        if highlight:
            ec = fc = highlight

        pp = PathPatch(path, ec=ec, fc=fc, alpha=alpha, lw=lw, zorder=zorder)
        ax.add_patch(pp)


class Region(object):
    """
    Draw a region of synteny.
    """

    def __init__(
        self,
        ax,
        ext,
        layout,
        bed,
        scale,
        switch=None,
        chr_label=True,
        loc_label=True,
        gene_labels: Optional[set] = None,
        genelabelsize=0,
        genelabelrotation=25,
        pad=0.05,
        vpad=0.015,
        extra_features=None,
        glyphstyle="box",
        glyphcolor: BasePalette = OrientationPalette(),
    ):
        x, y = layout.x, layout.y
        ratio = layout.ratio
        scale /= ratio
        self.y = y
        lr = layout.rotation
        tr = transforms.Affine2D().rotate_deg_around(x, y, lr) + ax.transAxes
        inv = ax.transAxes.inverted()

        start, end, si, ei, chrom, orientation, span = ext
        flank = span / scale / 2
        xstart, xend = x - flank, x + flank
        self.xstart, self.xend = xstart, xend

        cv = lambda t: xstart + abs(t - startbp) / scale
        hidden = layout.hidden

        # Chromosome
        if not hidden:
            ax.plot((xstart, xend), (y, y), color="gray", transform=tr, lw=2, zorder=1)

        self.genes = genes = bed[si : ei + 1]
        startbp, endbp = start.start, end.end
        if orientation == "-":
            startbp, endbp = endbp, startbp

        if switch:
            chrom = switch.get(chrom, chrom)
        if layout.label:
            chrom = layout.label

        label = "-".join(
            (
                human_size(startbp, target="Mb", precision=2)[:-2],
                human_size(endbp, target="Mb", precision=2),
            )
        )

        height = 0.012
        self.gg = {}
        # Genes
        for g in genes:
            gstart, gend = g.start, g.end
            strand = g.strand
            if strand == "-":
                gstart, gend = gend, gstart
            if orientation == "-":
                strand = "+" if strand == "-" else "-"

            x1, x2, a, b = self.get_coordinates(gstart, gend, y, cv, tr, inv)
            gene_name = g.accn
            self.gg[gene_name] = (a, b)

            color, zorder = (
                glyphcolor.get_color_and_zorder(strand)
                if isinstance(glyphcolor, OrientationPalette)
                else glyphcolor.get_color_and_zorder(gene_name)
            )

            if hidden:
                continue
            gp = Glyph(
                ax,
                x1,
                x2,
                y,
                height,
                gradient=False,
                fc=color,
                style=glyphstyle,
                zorder=zorder,
            )
            gp.set_transform(tr)
            if genelabelsize and (not gene_labels or gene_name in gene_labels):
                if genelabelrotation == 0:
                    text_x = x1 if x1 > x2 else x2
                    text_y = y
                else:
                    text_x = (x1 + x2) / 2
                    text_y = y + height / 2 + genelabelsize * vpad / 3
                ax.text(
                    text_x,
                    text_y,
                    markup(gene_name),
                    size=genelabelsize,
                    rotation=genelabelrotation,
                    ha="left",
                    va="center",
                    color="lightslategray",
                )

        # Extra features (like repeats)
        if extra_features:
            for g in extra_features:
                gstart, gend = g.start, g.end
                x1, x2, a, b = self.get_coordinates(gstart, gend, y, cv, tr, inv)
                gp = Glyph(
                    ax,
                    x1,
                    x2,
                    y,
                    height * 3 / 4,
                    gradient=False,
                    fc="#ff7f00",
                    style=glyphstyle,
                    zorder=2,
                )
                gp.set_transform(tr)

        ha, va = layout.ha, layout.va

        hpad = 0.02
        if ha == "left":
            xx = xstart - hpad
            ha = "right"
        elif ha == "leftalign":
            xx = 0.5 - CANVAS_SIZE / 2 - hpad
            ha = "right"
        elif ha == "right":
            xx = xend + hpad
            ha = "left"
        elif ha == "rightalign":
            xx = 0.5 + CANVAS_SIZE / 2 + hpad
            ha = "left"
        else:
            xx = x
            ha = "center"

        # Tentative solution to labels stick into glyph
        magic = 40.0
        cc = abs(lr) / magic if abs(lr) > magic else 1
        if va == "top":
            yy = y + cc * pad
        elif va == "bottom":
            yy = y - cc * pad
        else:
            yy = y

        l = np.array((xx, yy))
        trans_angle = ax.transAxes.transform_angles(np.array((lr,)), l.reshape((1, 2)))[
            0
        ]
        lx, ly = l
        if not hidden:
            bbox = dict(boxstyle="round", fc="w", ec="w", alpha=0.5)
            kwargs = dict(
                ha=ha, va="center", rotation=trans_angle, bbox=bbox, zorder=10
            )

            chr_label = markup(chrom) if chr_label else None
            loc_label = label if loc_label else None
            if chr_label:
                if loc_label:
                    ax.text(
                        lx,
                        ly + vpad,
                        chr_label,
                        size=layout.label_fontsize,
                        color=layout.color,
                        **kwargs,
                    )
                    ax.text(
                        lx,
                        ly - vpad,
                        loc_label,
                        color="lightslategrey",
                        size=layout.label_fontsize,
                        **kwargs,
                    )
                else:
                    ax.text(lx, ly, chr_label, color=layout.color, **kwargs)

    def get_coordinates(self, gstart, gend, y, cv, tr, inv):
        """
        Get coordinates of a gene.
        """
        x1, x2 = cv(gstart), cv(gend)
        a, b = tr.transform((x1, y)), tr.transform((x2, y))
        a, b = inv.transform(a), inv.transform(b)
        return x1, x2, a, b


def ymid_offset(samearc: Optional[str], pad: float = 0.05):
    """
    Adjustment to ymid, this is useful to adjust the appearance of the Bezier
    curves between the tracks.
    """
    if samearc == "above":
        return 2 * pad
    if samearc == "above2":
        return 4 * pad
    if samearc == "below":
        return -2 * pad
    if samearc == "below2":
        return -4 * pad
    return 0


class Synteny(object):
    """
    Draw the synteny plot.
    """

    def __init__(
        self,
        fig,
        root,
        datafile,
        bedfile,
        layoutfile,
        switch=None,
        tree=None,
        extra_features=None,
        chr_label: bool = True,
        loc_label: bool = True,
        gene_labels: Optional[set] = None,
        genelabelsize: int = 0,
        genelabelrotation: int = 25,
        pad: float = 0.05,
        vpad: float = 0.015,
        scalebar: bool = False,
        shadestyle: str = "curve",
        glyphstyle: str = "arrow",
        glyphcolor: str = "orientation",
        seed: Optional[int] = None,
        prune_features=True,
    ):
        _, h = fig.get_figwidth(), fig.get_figheight()
        bed = Bed(bedfile)
        order = bed.order
        bf = BlockFile(datafile)
        self.layout = lo = Layout(layoutfile, seed=seed)
        switch = DictFile(switch, delimiter="\t") if switch else None
        if extra_features:
            extra_features = Bed(extra_features)

        exts = []
        extras = []
        for i in range(bf.ncols):
            ext = bf.get_extent(i, order)
            exts.append(ext)
            if extra_features:
                start, end, _, _, chrom, _, span = ext
                start, end = start.start, end.end  # start, end coordinates
                ef = list(extra_features.extract(chrom, start, end))

                # Pruning removes minor features with < 0.1% of the region
                if prune_features:
                    ef_pruned = [x for x in ef if x.span >= span / 1000]
                    logger.info(
                        "Extracted %d features (%d after pruning)",
                        len(ef),
                        len(ef_pruned),
                    )
                    extras.append(ef_pruned)
                else:
                    logger.info("Extracted %d features", len(ef))
                    extras.append(ef)

        maxspan = max(exts, key=lambda x: x[-1])[-1]
        scale = maxspan / CANVAS_SIZE

        self.gg = gg = {}
        self.rr = []
        ymids = []
        glyphcolor = (
            OrientationPalette()
            if glyphcolor == "orientation"
            else OrthoGroupPalette(bf.grouper())
        )
        for i in range(bf.ncols):
            ext = exts[i]
            ef = extras[i] if extras else None
            r = Region(
                root,
                ext,
                lo[i],
                bed,
                scale,
                switch,
                gene_labels=gene_labels,
                genelabelsize=genelabelsize,
                genelabelrotation=genelabelrotation,
                chr_label=chr_label,
                loc_label=loc_label,
                vpad=vpad,
                extra_features=ef,
                glyphstyle=glyphstyle,
                glyphcolor=glyphcolor,
            )
            self.rr.append(r)
            # Use tid and accn to store gene positions
            gg.update(dict(((i, k), v) for k, v in r.gg.items()))
            ymids.append(r.y)

        for i, j, blockcolor, samearc in lo.edges:
            ymid_pad = ymid_offset(samearc, pad)
            for ga, gb, h in bf.iter_pairs(i, j):
                a, b = gg[(i, ga)], gg[(j, gb)]
                Shade(
                    root, a, b, ymid_pad, fc=blockcolor, lw=0, alpha=1, style=shadestyle
                )

            for ga, gb, h in bf.iter_pairs(i, j, highlight=True):
                a, b = gg[(i, ga)], gg[(j, gb)]
                Shade(
                    root,
                    a,
                    b,
                    ymid_pad,
                    alpha=1,
                    highlight=h,
                    zorder=2,
                    style=shadestyle,
                )

        if scalebar:
            logger.info("Build scalebar (scale=%.3f)", scale)
            # Find the best length of the scalebar
            ar = [1, 2, 5]
            candidates = (
                [1000 * x for x in ar]
                + [10000 * x for x in ar]
                + [100000 * x for x in ar]
            )
            # Find the one that's close to an optimal canvas size
            dists = [(abs(x / scale - 0.12), x) for x in candidates]
            dist, candidate = min(dists)
            dist = candidate / scale
            x, y, yp = 0.22, 0.92, 0.005
            a, b = x - dist / 2, x + dist / 2
            lsg = "lightslategrey"
            root.plot([a, a], [y - yp, y + yp], "-", lw=2, color=lsg)
            root.plot([b, b], [y - yp, y + yp], "-", lw=2, color=lsg)
            root.plot([a, b], [y, y], "-", lw=2, color=lsg)
            root.text(
                x,
                y + 0.02,
                human_size(candidate, precision=0),
                ha="center",
                va="center",
            )

        if tree:
            trees = read_trees(tree)
            ntrees = len(trees)
            logger.debug("A total of %d trees imported.", ntrees)
            xiv = 1.0 / ntrees
            yiv = 0.3
            xstart = 0
            ystart = min(ymids) - 0.4
            for i in range(ntrees):
                ax = fig.add_axes([xstart, ystart, xiv, yiv])
                label, outgroup, color, tx = trees[i]
                draw_tree(
                    ax,
                    tx,
                    outgroup=outgroup,
                    rmargin=0.4,
                    leaffont=11,
                    treecolor=color,
                    supportcolor=color,
                    leafcolor=color,
                )
                xstart += xiv
                RoundLabel(ax, 0.5, 0.3, label, fill=True, fc="lavender", color=color)


def draw_gene_legend(
    ax,
    x1: float,
    x2: float,
    ytop: float,
    d: float = 0.04,
    text: bool = False,
    repeat: bool = False,
    glyphstyle: str = "box",
):
    """
    Draw a legend for gene glyphs.
    """
    forward, backward = OrientationPalette.forward, OrientationPalette.backward
    ax.plot([x1, x1 + d], [ytop, ytop], ":", color=forward, lw=2)
    ax.plot([x1 + d], [ytop], ">", color=forward, mec=forward)
    ax.plot([x2, x2 + d], [ytop, ytop], ":", color=backward, lw=2)
    ax.plot([x2], [ytop], "<", color=backward, mec="g")
    if text:
        ax.text(x1 + d / 2, ytop + d / 2, "gene (+)", ha="center")
        ax.text(x2 + d / 2, ytop + d / 2, "gene (-)", ha="center")
    if repeat:
        xr = (x1 + x2 + d) / 2
        Glyph(
            ax,
            xr - d / 2,
            xr + d / 2,
            ytop,
            0.012 * 3 / 4,
            gradient=False,
            fc="#ff7f00",
            style=glyphstyle,
            zorder=2,
        )
        ax.text(xr, ytop + d / 2, "repeat", ha="center")


def main(args: List[str]):
    p = OptionParser(__doc__)
    p.add_argument("--switch", help="Rename the seqid with two-column file")
    p.add_argument("--tree", help="Display trees on the bottom of the figure")
    p.add_argument("--extra", help="Extra features in BED format")
    p.add_argument(
        "--genelabels",
        help='Show only these gene labels, separated by comma. Example: "At1g12340,At5g54690"',
    )
    p.add_argument(
        "--genelabelsize",
        default=0,
        type=int,
        help="Show gene labels at this font size, useful for debugging. "
        + "However, plot may appear visually crowded. "
        + "Reasonably good values are 2 to 6 [Default: disabled]",
    )
    p.add_argument(
        "--genelabelrotation",
        default=25,
        type=int,
        help="Rotate gene labels at this angle (anti-clockwise), useful for debugging.",
    )
    p.add_argument(
        "--scalebar",
        default=False,
        action="store_true",
        help="Add scale bar to the plot",
    )
    p.add_argument(
        "--glyphstyle",
        default="box",
        choices=Glyph.Styles,
        help="Style of feature glyphs",
    )
    p.add_argument(
        "--glyphcolor",
        default="orientation",
        choices=Glyph.Palette,
        help="Glyph coloring based on",
    )
    p.add_argument(
        "--shadestyle",
        default="curve",
        choices=Shade.Styles,
        help="Style of syntenic wedges",
    )
    p.add_argument(
        "--outputprefix",
        default="",
        help="Prefix for the output file",
    )
    p.add_argument(
        "--noprune",
        default=False,
        action="store_true",
        help="If set, do not exclude small features from annotation track (<1%% of region)",
    )
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 3:
        sys.exit(not p.print_help())

    datafile, bedfile, layoutfile = args
    switch = opts.switch
    tree = opts.tree
    gene_labels = None if not opts.genelabels else set(opts.genelabels.split(","))
    prune_features = not opts.noprune

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))
    Synteny(
        fig,
        root,
        datafile,
        bedfile,
        layoutfile,
        switch=switch,
        tree=tree,
        extra_features=opts.extra,
        gene_labels=gene_labels,
        genelabelsize=opts.genelabelsize,
        genelabelrotation=opts.genelabelrotation,
        scalebar=opts.scalebar,
        shadestyle=opts.shadestyle,
        glyphstyle=opts.glyphstyle,
        glyphcolor=opts.glyphcolor,
        seed=iopts.seed,
        prune_features=prune_features,
    )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    outputprefix = opts.outputprefix
    if outputprefix:
        pf = outputprefix
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)

    return image_name


if __name__ == "__main__":
    main(sys.argv[1:])
