#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
%prog alignment.blocks alignment.bed layout.csv

Illustrate shared identity between aligned chromosomes. Use layout.csv to indicate
the positions of tracks. For example:

# x, y, rotation, ha, va, color, ratio, label, chrmName, rStart, rEnd, chrmMax
0.5, 0.6, 0, top, center, g, 1, Spp 1 Chrm 8, Chromosome_08, 2000000, 2108098, 2108098
0.5, 0.4, 0, right, center, m, 1, Spp 2 Chrm 1.4 , Scaffold_1.4, 50000, 200000, 2787645
0.5, 0.2, 45, bottom, center, b, 1, Spp3 Chr 5, SppX_Chromosome_05, 60000, 460000, 2046703
# edges
e, 0, 1
e, 1, 2

With the row ordering corresponding to the column ordering in the alignment blockfile.
"""

import numpy as np
import sys

# from jcvi.formats.base import DictFile
from jcvi.apps.base import OptionParser, logger
from jcvi.compara.synteny import BlockFile
from jcvi.formats.bed import Bed
from jcvi.graphics.base import (
    AbstractLayout,
    markup,
    mpl,
    Path,
    PathPatch,
    plt,
    savefig,
)
from jcvi.graphics.glyph import Glyph, RoundLabel  # , GeneGlyph
from jcvi.utils.cbook import human_size


# Default colours for ribbons with different orientations
forward, backward = "#1f77b4", "#2ca02c"


class LayoutLine(object):
    # 0-11: x,y,rotation,ha,va,color,ratio,label,chrmName,rStart,rEnd,chrmMax
    def __init__(self, row, delimiter=","):
        self.hidden = row[0] == "*"
        if self.hidden:
            row = row[1:]
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.x = float(args[0])
        self.y = float(args[1])
        self.rotation = int(args[2])
        self.ha = args[3]
        self.va = args[4]
        self.color = args[5]
        self.ratio = 1
        if len(args) > 6:
            self.ratio = float(args[6])
        if len(args) > 7:
            self.label = args[7].strip()
        else:
            self.label = None
        if len(args) > 8:
            self.chrmName = str(args[8])
        if len(args) > 9:
            self.rStart = int(args[9])
            self.rEnd = int(args[10])
        else:
            self.rStart = None
            self.rEnd = None
        if len(args) > 10:
            self.chrmMax = int(args[10])
        else:
            self.chrmMax = None


class Layout(AbstractLayout):
    def __init__(self, filename, delimiter=","):
        super(Layout, self).__init__(filename)
        fp = open(filename)
        self.edges = []
        for row in fp:
            # Skip blank lines
            if row.strip():
                # Skip header rows
                if row[0] == "#":
                    continue
                # Import edges
                if row[0] == "e":
                    args = row.rstrip().split(delimiter)
                    args = [x.strip() for x in args]
                    # From index 1 and 2
                    a, b = args[1:3]
                    a, b = int(a), int(b)
                    assert args[0] == "e"
                    self.edges.append((a, b))
                else:
                    # Parse all other lines as sequence tracks
                    self.append(LayoutLine(row, delimiter=delimiter))

        # Track number of layout lines added
        self.lines = len(self)

        if 2 <= len(self) <= 8:
            self.assign_colors()


class Shade(object):
    Styles = ("curve", "line")

    def __init__(
        self,
        ax,
        a,
        b,
        ymid,
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
            # Note: a & b customised from synteny for ribbon module. Description is from synteny version.
            #a (tuple of floats): ((start_x, start_y), (end_x, end_y))
            #b (tuple of floats): ((start_x, start_y), (end_x, end_y))
            ymid (float): y-mid position for curve style only
            highlight (bool, optional): Plot this shade if color is specified. Defaults to False.
            style (str, optional): Style. Defaults to "curve", must be one of
            ("curve", "line")
            ec (str, optional): Edge color. Defaults to "k".
            fc (str, optional): Face color. Defaults to "k".
            alpha (float, optional): Transparency. Defaults to 0.2.
            lw (int, optional): Line width. Defaults to 1.
            zorder (int, optional): Z-order. Defaults to 1.
        """
        assert style in Shade.Styles, "style must be one of {}".format(self.Styles)

        # a1, a2 = a
        # b1, b2 = b
        a1, a2, a_start, a_end, a_strand, a_score = a
        b1, b2, b_start, b_end, b_strand, b_score = b
        ax1, ay1 = a1
        ax2, ay2 = a2
        bx1, by1 = b1
        bx2, by2 = b2
        M, C4, L, CP = Path.MOVETO, Path.CURVE4, Path.LINETO, Path.CLOSEPOLY
        if style == "curve":
            pathdata = [
                (M, a1),
                (C4, (ax1, ymid)),
                (C4, (bx1, ymid)),
                (C4, b1),
                (L, b2),
                (C4, (bx2, ymid)),
                (C4, (ax2, ymid)),
                (C4, a2),
                (CP, a1),
            ]
        else:
            pathdata = [(M, a1), (L, b1), (L, b2), (L, a2), (CP, a1)]
        codes, verts = list(zip(*pathdata))
        path = Path(verts, codes)
        if highlight:
            ec = fc = highlight

        pp = PathPatch(path, ec=ec, fc=fc, alpha=alpha, lw=lw, zorder=zorder)
        ax.add_patch(pp)


class Region(object):
    def __init__(
        self,
        ax,
        ext,
        layout,
        bed,
        scale,
        orientation=None,
        parent_ori=None,
        chr_label=True,
        loc_label=True,
        pad=0.05,
        vpad=0.015,
        features=None,
        plotRibbonBlocks=False,
        annotcolor="g",
    ):
        x, y = layout.x, layout.y
        ratio = layout.ratio
        scale /= ratio
        self.y = y
        lr = layout.rotation
        tr = mpl.transforms.Affine2D().rotate_deg_around(x, y, lr) + ax.transAxes
        inv = ax.transAxes.inverted()

        # Import standard track extent (by max feature range)
        start, end, si, ei, chrm, trackOri, span = ext
        startbp, endbp = (
            start.start,
            end.end,
        )  # Get start/end of first/last features from blockfile
        # Default all input tracks to Fwd orientation
        trackOri = "+"

        # Override track extent if custom range provided
        if layout.chrmName and layout.rStart and layout.rEnd:
            startbp, endbp, chrm = (layout.rStart, layout.rEnd, layout.chrmName)
            span = endbp - (startbp + 1)

        # Set flank size
        flank = span / scale / 2

        # Adjust x-coords to flank size
        xstart, xend = x - flank, x + flank
        self.xstart, self.xend = xstart, xend

        def cv(t):
            return xstart + abs(t - startbp) / scale

        hidden = layout.hidden

        # Plot Chromosome Bar
        if not hidden:
            ax.plot(
                (xstart, xend), (y, y), color="dimgray", transform=tr, lw=4, zorder=2
            )

        # Compose Chromosome label
        if orientation == "R" or trackOri == "-":
            label_startbp, label_endbp = endbp, startbp
        else:
            label_startbp, label_endbp = startbp, endbp

        if layout.label:
            chrm = layout.label

        label = "-".join(
            (
                human_size(label_startbp, target="Mb", precision=2)[:-2],
                human_size(label_endbp, target="Mb", precision=2),
            )
        )

        # Plot annotation tracks
        height = 0.012
        self.gg = {}

        # Process and plot ribbon intervals as annotation features
        self.ribbonBlocks = ribbonBlocks = bed[si : ei + 1]
        for g in ribbonBlocks:
            gstart, gend, outofrange = self.clip2range(g.start, g.end, startbp, endbp)
            if outofrange:
                continue
            # Set feature orientation - Note: Not really 'strand' just orientation relative to Parent where Parent is always '+'
            strand = g.strand
            if orientation == "R":
                gstart, gend = self.flip(gstart, gend, startbp, endbp)
                assert gstart <= gend
                # Flip feature strand (relative to parent) if self track has been flipped
                strand = "+" if strand == "-" else "-"
            # Invert gene start/end positions if feature on '-' strand
            if strand == "-":
                gstart, gend = gend, gstart
                assert gstart >= gend
            x1, x2, a, b = self.get_coordinates(gstart, gend, y, cv, tr, inv)
            # Note: gg is later used to plot ribbons
            # TODO: Create separate dict for (gstart, gend, strand, g.score) so that we don't have to tweak downstream functions that only want (a,b).
            self.gg[g.accn] = (a, b, gstart, gend, strand, g.score)
            # Filp feature orientation if Parent track has been flipped
            if parent_ori == "R":
                strand = "+" if strand == "-" else "-"
            color = forward if strand == "+" else backward
            if not hidden and plotRibbonBlocks:
                gp = Glyph(ax, x1, x2, y, height, gradient=False, fc=color, zorder=3)
                gp.set_transform(tr)

        # Set default feature colour
        feat_col = annotcolor.strip()

        # Add annotation track offset
        offset = 0.005

        # Set default annotation feature height
        feat_height = height * 0.3

        # Plot feature track (genes)
        if features:
            for g in features:
                # Bedline attributes =  seqid,start,end,accn,score,strand, extra = [color,offset,height]
                gstart, gend, outofrange = self.clip2range(
                    g.start, g.end, startbp, endbp
                )
                if outofrange:
                    continue
                if orientation == "R":
                    gstart, gend = self.flip(gstart, gend, startbp, endbp)
                x1, x2, a, b = self.get_coordinates(gstart, gend, y, cv, tr, inv)
                # Set custom annotation color / y-offset / height
                annot_col, annot_offset, annot_height = self.annotFormat(
                    g.extra, offset, feat_col, feat_height
                )
                # Note: Use GeneGlyph instead of Glyph for terminal CDS exon or for genes.
                # Need to store feature direction, and either flag terminal exons OR unique parent feature names (to determine last exon).
                # GeneGlyph(self, ax, x1, x2, y, height, gradient=True, tip=.0025, color="k", shadow=False)
                gp = Glyph(
                    ax,
                    x1,
                    x2,
                    y + annot_offset,
                    annot_height,
                    gradient=False,
                    fc=annot_col,
                    zorder=4,
                )
                gp.set_transform(tr)

        # Position and apply chromosome labels
        ha, va = layout.ha, layout.va

        hpad = 0.02
        if ha == "left":
            xx = xstart - hpad
            ha = "right"
        elif ha == "right":
            xx = xend + hpad
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

            # Format Chrm label for LaTeX
            chr_label = markup(chrm) if chr_label else None
            # Add seq coords to label
            loc_label = label if loc_label else None
            if chr_label:
                if loc_label:
                    ax.text(lx, ly + vpad, chr_label, color=layout.color, **kwargs)
                    ax.text(
                        lx,
                        ly - vpad,
                        loc_label,
                        color="lightslategrey",
                        size=10,
                        **kwargs,
                    )
                else:
                    ax.text(lx, ly, chr_label, color=layout.color, **kwargs)

    def annotFormat(self, extra, offset, feat_col, feat_height):
        annot_col = feat_col
        annot_offset = offset
        annot_height = feat_height
        # extra [colour,offset multiplier, height multiplier]
        if extra:
            if len(extra) >= 1 and str(extra[0]) != ".":
                annot_col = str(extra[0])
            if len(extra) >= 2 and str(extra[1]) != ".":
                annot_offset = offset * float(extra[1])
            if len(extra) >= 3 and str(extra[2]) != "." and float(extra[2]) > 0:
                annot_height = feat_height * float(extra[2])
        return (annot_col, annot_offset, annot_height)

    def get_coordinates(self, gstart, gend, y, cv, tr, inv):
        x1, x2 = cv(gstart), cv(gend)
        a, b = tr.transform((x1, y)), tr.transform((x2, y))
        a, b = inv.transform(a), inv.transform(b)
        return x1, x2, a, b

    def clip2range(self, f_start, f_end, start, end):
        # Return unchanged if annotation within plot range
        if f_start >= start and f_end <= end:
            return f_start, f_end, False
        # Clip if plot range fully within annotation
        if f_start <= start and f_end >= end:
            return start, end, False
        # If annotation starts before plot range
        if f_start <= start:
            # And ends before plot range, exclude annotation
            if f_end <= start:
                return None, None, True
            # If end is within plot range, clip to range start
            elif f_end >= start:
                return start, f_end, False
        # If annotation ends after plot range
        if f_end >= end:
            # And begins after range end, exclude annotation
            if f_start >= end:
                return None, None, True
            # If annotation start is within plot range, clip to range end
            elif f_start <= end:
                return f_start, end, False

    def flip(self, gstart, gend, startbp, endbp):
        upstream = gstart - startbp
        downstream = endbp - gend
        newEnd = endbp - upstream
        newStart = startbp + downstream
        return newStart, newEnd


class Synteny(object):
    def __init__(
        self,
        fig,
        root,
        blocks,
        bedfile,
        layoutfile,
        orientation=None,
        tree=None,
        features=None,
        chr_label=True,
        loc_label=True,
        pad=0.05,
        vpad=0.015,
        scalebar=False,
        paintbyscore=False,
        paintbystrand=False,
        prune_features=True,
        plotRibbonBlocks=False,
        annotcolor="g",
        scaleAlpha=False,
        noline=False,
        shadestyle="curve",
    ):
        # Set figure dimensions
        w, h = fig.get_figwidth(), fig.get_figheight()
        # Import primary annotations from bedfile
        bed = Bed(bedfile)
        # Order bed features as dict {'ppa004748m': (42353, scaffold_5    1666610 1670773 ppa004748m  0   +),..}
        # genename: (orderedPosition, seqname, start, end, genename, score, strand)
        order = bed.order
        # Import synteny links (pairs of source > target IDs for features in BED file)
        bf = BlockFile(blocks, defaultcolor="hide")
        # Import layout config
        self.layout = lo = Layout(layoutfile)

        # Check correct number of orientation values
        if orientation:
            if lo.lines != len(orientation):
                logger.error("Incorrect number of orientation instructions")
                sys.exit(0)

        # Import features track
        if features:
            features = Bed(features)

        # Init list for track extent, features, and extra features
        exts = []
        feats = []

        # Collect span for any tracks with custom range coords
        customSpans = []

        # For each track (by index)
        for i in range(bf.ncols):
            # Set range coord containing all annotations
            ext = bf.get_extent(i, order)
            # Add to ext index
            exts.append(ext)
            # ext has format:
            # (  chr1  10730   27033   GSVIVT01012261001   0   +,
            #   chr1 558160  562819  GSVIVT01012208001   0   +,
            #   0,
            #   49,
            #   'chr1',
            #   '+',
            #   552088)
            # = first feature, last feature, bf idx ff, bf idx lf, chr name,
            #   track orientation, span (range from ff to lf)

            # Override block feature extents if custom range provided in layout
            if lo[i].chrmName and lo[i].rStart and lo[i].rEnd:
                start, end, chrm = (lo[i].rStart, lo[i].rEnd, lo[i].chrmName)
                span = end - (start + 1)
                customSpans.append(span)
            else:
                start, end, si, ei, chrm, orientation, span = ext
                # Get start/end of first/last features from blockfile
                start, end = start.start, end.end
                customSpans.append(span)

            if features:
                # Unpack coords from 'features' bed object
                # Filter for features within range to be plotted
                fe = list(features.extract(chrm, start, end))
                # Pruning removes minor features with < 0.1% of the region
                if prune_features:
                    fe_pruned = [x for x in fe if x.span >= span / 1000]
                    logger.info(
                        "Extracted {0} features "
                        "({1} after pruning)".format(len(fe), len(fe_pruned)),
                        file=sys.stderr,
                    )
                    feats.append(fe_pruned)
                else:
                    fe_all = [x for x in fe]
                    feats.append(fe_all)

        # Find largest coord range for any track
        # maxspan = max(exts, key=lambda x: x[-1])[-1]
        maxspan = max(customSpans)
        scale = maxspan / 0.65

        self.gg = gg = {}
        self.rr = []
        ymids = []

        # Plot annotations
        for i in range(bf.ncols):
            ext = exts[i]
            fe = feats[i] if feats else None
            ori = orientation[i] if orientation else None
            parent_ori = orientation[i - 1] if orientation and i >= 1 else None
            r = Region(
                root,
                ext,
                lo[i],
                bed,
                scale,
                orientation=ori,
                parent_ori=parent_ori,
                chr_label=chr_label,
                loc_label=loc_label,
                vpad=vpad,
                features=fe,
                plotRibbonBlocks=plotRibbonBlocks,
                annotcolor=annotcolor,
            )
            self.rr.append(r)
            # Use tid and accn to store gene positions
            gg.update(dict(((i, k), v) for k, v in list(r.gg.items())))
            ymids.append(r.y)

        # Plot ribbons
        for i, j in lo.edges:
            # geneA, geneB, highlight colour
            for ga, gb, h in bf.iter_pairs(i, j):
                # Coords in gg should be correctly trimmed from Region processing
                # Need to skip annotations in blockfile that are outside plotted region in layout file
                if (i, ga) not in list(gg.keys()) or (j, gb) not in list(gg.keys()):
                    continue
                # Skip ribbons that will later be plotted with a hightlight color
                if h:
                    continue
                a, b = gg[(i, ga)], gg[(j, gb)]
                ymid = (ymids[i] + ymids[j]) / 2
                # Draw synteny ribbons
                ribbonColor = "gainsboro"
                baseAlpha = 0.8
                if paintbystrand:
                    ribbonColor = self.inversionCheck(a, b)
                if paintbyscore:
                    baseAlpha = self.scoreCheck(a, scaleAlpha)
                if noline:
                    lw = 0.0
                else:
                    lw = 0.1
                Shade(
                    root,
                    a,
                    b,
                    ymid,
                    fc=ribbonColor,
                    lw=lw,
                    alpha=baseAlpha,
                    style=shadestyle,
                )
                # Add second low alpha ribbon track with thin outlines. Helps visualise very thin bands.
                # Shade(root, a, b, ymid, fc=ribbonColor, lw=0.1, alpha=0.125, style=shadestyle)

            # Paint ribbons for which a highlight colour was set
            for ga, gb, h in bf.iter_pairs(i, j, highlight=True):
                if (i, ga) not in list(gg.keys()) or (j, gb) not in list(gg.keys()):
                    continue
                if h == "hide":
                    continue
                a, b = gg[(i, ga)], gg[(j, gb)]
                ymid = (ymids[i] + ymids[j]) / 2
                baseAlpha = 1
                if paintbyscore:
                    baseAlpha = self.scoreCheck(a, scaleAlpha)
                Shade(
                    root,
                    a,
                    b,
                    ymid,
                    alpha=baseAlpha,
                    highlight=h,
                    zorder=1,
                    style=shadestyle,
                )

        if scalebar:
            logger.info("Build scalebar (scale={})".format(scale), file=sys.stderr)
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
            x, y, yp = 0.2, 0.96, 0.005
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
            from jcvi.graphics.tree import draw_tree, read_trees

            trees = read_trees(tree)
            ntrees = len(trees)
            logger.debug("A total of {0} trees imported.".format(ntrees))
            xiv = 1.0 / ntrees
            yiv = 0.3
            xstart = 0
            ystart = min(ymids) - 0.4
            for i in range(ntrees):
                ax = fig.add_axes([xstart, ystart, xiv, yiv])
                label, outgroup, tx = trees[i]
                draw_tree(ax, tx, outgroup=outgroup, rmargin=0.4, leaffont=11)
                xstart += xiv
                RoundLabel(ax, 0.5, 0.3, label, fill=True, fc="lavender", color="r")

    def inversionCheck(self, a, b):
        a1, a2, a_start, a_end, a_strand, a_score = a
        b1, b2, b_start, b_end, b_strand, b_score = b
        aDir = "F" if a_start < a_end else "R"
        bDir = "F" if b_start < b_end else "R"
        if aDir == bDir:
            return forward
        else:
            return backward

    def scoreCheck(self, a, scaleAlpha):
        a1, a2, a_start, a_end, a_strand, a_score = a
        assert 0 <= float(a_score) <= 100
        if scaleAlpha:
            OldValue = float(a_score)
            OldMin = 50
            OldMax = 100
            NewMax = 90
            NewMin = 10
            newScore = (
                ((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)
            ) + NewMin
            if newScore < NewMin:
                newScore = NewMin
            return float(newScore) / 100
        else:
            return float(a_score) / 100


def draw_gene_legend(ax, x1, x2, ytop, d=0.04, text=False, repeat=False):
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
            zorder=2,
        )
        ax.text(xr, ytop + d / 2, "repeat", ha="center")


def set_strand_colors(colorCodes):
    global forward
    global backward
    forward, backward = colorCodes.strip().split(",")


def main():
    # Get cmd line args
    p = OptionParser(__doc__)
    p.add_option("--outfile", default=None, help="Prefix for output graphic.")
    p.add_option(
        "--annotations",
        help="Feature annotations in BED format. \n \
        [1] seqid \n \
        [2] start \n \
        [3] end \n \
        [4] accn \n \
        [5] score \n \
        [6] strand \n \
        [7] Custom color. i.e. '.' = use default color, else any pyplot compatible color code: 'g', 'green, '#02ab2e', etc. \n \
        [8] Vertical offset multiplier. i.e. 0 = plot on chrm line, -1 = plot below chrm, 1 = plot above chrm \n \
        [9] Feature height multiplier. i.e. 1 = default, 2 = double height ",
    )
    p.add_option(
        "--noprune",
        default=True,
        action="store_false",
        help="If set, do not exclude small features from annotation track. ",
    )
    p.add_option(
        "--reorient",
        default=None,
        help="Comma delimited string of 'F' or 'R' characters. Must be same number and order as tracks in layout file. \n i.e. F,F,R will flip the third track. Default: All Forward.",
    )
    p.add_option(
        "--tree", help="Display trees on the bottom of the figure [default: %default]"
    )
    p.add_option(
        "--scalebar",
        default=False,
        action="store_true",
        help="Add scale bar to the plot",
    )
    p.add_option(
        "--paintbyscore",
        default=False,
        action="store_true",
        help="Set ribbon transparancy using score column from ribbon bedfile.",
    )
    p.add_option(
        "--paintbystrand",
        default=False,
        action="store_true",
        help="Set ribbon colour by alignment orientation. Red=Inverted, Blue=Same",
    )
    p.add_option(
        "--plotHits",
        default=False,
        action="store_true",
        help="If set, plot ribbon features from blockfile also as annotations.",
    )
    p.add_option(
        "--strandcolors",
        default=None,
        help="Comma delimited string of color codes for forward or inverted orientation ribbons. Used by paintbystrand. Default: 'b,g' for forward,inverted.",
    )
    p.add_option(
        "--annotcolor",
        default="g",
        help="Comma delimited string of color code for 'annotation' feature tracks. Default: 'g'.",
    )
    p.add_option(
        "--scaleAlpha",
        default=False,
        action="store_true",
        help="If set, ribbon alpha values will be rescaled from 0.5-1 to 0.15-0.95 Note: Assumes 50% min identity in alignments.",
    )
    p.add_option(
        "--transparent",
        default=False,
        action="store_true",
        help="If set, save image with transparent background.",
    )
    p.add_option(
        "--noline",
        default=False,
        action="store_true",
        help="If set, do not draw outline on ribbons.",
    )
    p.add_option(
        "--shadestyle",
        default="curve",
        choices=Shade.Styles,
        help="Style of syntenic wedges",
    )
    # Unpack options
    # opts = {'style': 'darkgrid', 'format': 'pdf', 'scalebar': False, 'extra': 'grape_peach_cacao.bed', 'tree': None, 'diverge': 'PiYG', 'cmap': 'jet', 'figsize': '8x7', 'font': 'Helvetica', 'dpi': 300}
    # args = positional args, data files
    # iopts = (2400px x 2100px)
    opts, args, iopts = p.set_image_options(figsize="20x20")

    # Check for data files
    if len(args) != 3:
        logger.error("Requires 3 data file args.")
        sys.exit(not p.print_help())

    # Unpack data file paths
    blocks, bedfile, layoutfile = args

    # Set strand colors as global vars
    if opts.strandcolors:
        set_strand_colors(opts.strandcolors)

    # Get custom track orientations
    if opts.reorient:
        flip = [x for x in opts.reorient.strip().split(",") if x in ("F", "R")]
    else:
        flip = None

    # Create base image
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Plot all the things!
    Synteny(
        fig,
        root,
        blocks,
        bedfile,
        layoutfile,
        orientation=flip,
        tree=opts.tree,
        features=opts.annotations,
        scalebar=opts.scalebar,
        paintbyscore=opts.paintbyscore,
        paintbystrand=opts.paintbystrand,
        prune_features=opts.noprune,
        plotRibbonBlocks=opts.plotHits,
        annotcolor=opts.annotcolor,
        scaleAlpha=opts.scaleAlpha,
        noline=opts.noline,
        shadestyle=opts.shadestyle,
    )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    if opts.outfile:
        image_name = str(opts.outfile) + "." + iopts.format
    else:
        pf = blocks.rsplit(".", 1)[0]
        image_name = pf + "." + iopts.format

    savefig(image_name, dpi=iopts.dpi, iopts=iopts, transparent=opts.transparent)


if __name__ == "__main__":
    main()
