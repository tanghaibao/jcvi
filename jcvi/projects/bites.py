#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the Brapa bites paper

Tang et al. (2012) Altered Patterns of Fractionation and Exon Deletions in
Brassica rapa Support a Two-Step Model of Paleohexaploidy. Genetics.
<http://www.genetics.org/content/190/4/1563.short>
"""
from more_itertools import pairwise

from jcvi.graphics.base import plt, Rectangle, Polygon, CirclePolygon, savefig
from jcvi.graphics.glyph import RoundLabel, arrowprops, TextCircle
from jcvi.graphics.chromosome import Chromosome
from jcvi.apps.base import OptionParser, ActionDispatcher, fname


def main():

    actions = (
        ("excision", "show intra-chromosomal recombination"),
        ("bites", "show the bites calling pipeline"),
        ("scenario", "show step-wise genome merger events in brapa"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def excision(args):
    """
    %prog excision

    Illustrate the mechanism of illegitimate recombination.
    """
    p = OptionParser(__doc__)
    opts, args = p.parse_args(args)

    fig = plt.figure(1, (5, 5))
    root = fig.add_axes([0, 0, 1, 1])

    plt.plot((0.2, 0.8), (0.6, 0.6), "r-", lw=3)
    plt.plot((0.4, 0.6), (0.6, 0.6), "b>-", mfc="g", mec="w", ms=12, lw=3)
    plt.plot((0.3, 0.7), (0.5, 0.5), "r-", lw=3)
    plt.plot((0.5,), (0.5,), "b>-", mfc="g", mec="w", ms=12, lw=3)

    # Circle excision
    plt.plot((0.5,), (0.45,), "b>-", mfc="g", mec="w", ms=12, lw=3)
    circle = CirclePolygon((0.5, 0.4), 0.05, fill=False, lw=3, ec="b")
    root.add_patch(circle)

    arrow_dist = 0.07
    ar_xpos, ar_ypos = 0.5, 0.52
    root.annotate(
        " ", (ar_xpos, ar_ypos), (ar_xpos, ar_ypos + arrow_dist), arrowprops=arrowprops
    )

    RoundLabel(root, 0.2, 0.64, "Gene")
    RoundLabel(root, 0.3, 0.54, "Excision")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


def bites(args):
    """
    %prog bites

    Illustrate the pipeline for automated bite discovery.
    """

    p = OptionParser(__doc__)
    opts, args = p.parse_args()

    fig = plt.figure(1, (6, 6))
    root = fig.add_axes([0, 0, 1, 1])

    # HSP pairs
    hsps = (
        ((50, 150), (60, 180)),
        ((190, 250), (160, 235)),
        ((300, 360), (270, 330)),
        ((430, 470), (450, 490)),
        ((570, 620), (493, 543)),
        ((540, 555), (370, 385)),  # non-collinear hsps
    )

    titlepos = (0.9, 0.65, 0.4)
    titles = ("Compare orthologous region", "Find collinear HSPs", "Scan paired gaps")
    ytip = 0.01
    mrange = 650.0
    m = lambda x: x / mrange * 0.7 + 0.1
    for i, (ya, title) in enumerate(zip(titlepos, titles)):
        yb = ya - 0.1
        plt.plot((0.1, 0.8), (ya, ya), "-", color="gray", lw=2, zorder=1)
        plt.plot((0.1, 0.8), (yb, yb), "-", color="gray", lw=2, zorder=1)
        RoundLabel(root, 0.5, ya + 4 * ytip, title)
        root.text(0.9, ya, "A. thaliana", ha="center", va="center")
        root.text(0.9, yb, "B. rapa", ha="center", va="center")
        myhsps = hsps
        if i >= 1:
            myhsps = hsps[:-1]
        for (a, b), (c, d) in myhsps:
            a, b, c, d = [m(x) for x in (a, b, c, d)]
            r1 = Rectangle((a, ya - ytip), b - a, 2 * ytip, fc="r", lw=0, zorder=2)
            r2 = Rectangle((c, yb - ytip), d - c, 2 * ytip, fc="r", lw=0, zorder=2)
            r3 = Rectangle((a, ya - ytip), b - a, 2 * ytip, fill=False, zorder=3)
            r4 = Rectangle((c, yb - ytip), d - c, 2 * ytip, fill=False, zorder=3)
            r5 = Polygon(
                ((a, ya - ytip), (c, yb + ytip), (d, yb + ytip), (b, ya - ytip)),
                fc="r",
                alpha=0.2,
            )
            rr = (r1, r2, r3, r4, r5)
            if i == 2:
                rr = rr[:-1]
            for r in rr:
                root.add_patch(r)

    # Gap pairs
    hspa, hspb = zip(*myhsps)
    gapa, gapb = [], []
    for (a, b), (c, d) in pairwise(hspa):
        gapa.append((b + 1, c - 1))
    for (a, b), (c, d) in pairwise(hspb):
        gapb.append((b + 1, c - 1))
    gaps = zip(gapa, gapb)
    tpos = titlepos[-1]

    yy = tpos - 0.05
    for i, ((a, b), (c, d)) in enumerate(gaps):
        i += 1
        a, b, c, d = [m(x) for x in (a, b, c, d)]
        xx = (a + b + c + d) / 4
        TextCircle(root, xx, yy, str(i))

    # Bites
    ystart = 0.24
    ytip = 0.05
    bites = (
        ("Bite(40=>-15)", True),
        ("Bite(50=>35)", False),
        ("Bite(70=>120)", False),
        ("Bite(100=>3)", True),
    )
    for i, (bite, selected) in enumerate(bites):
        xx = 0.15 if (i % 2 == 0) else 0.55
        yy = ystart - i / 2 * ytip
        i += 1
        TextCircle(root, xx, yy, str(i))
        color = "k" if selected else "gray"
        root.text(xx + ytip, yy, bite, size=10, color=color, va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


def scenario(args):
    """
    %prog scenario

    Illustration of the two-step genome merger process for B. rapa companion paper.
    """
    p = OptionParser(__doc__)
    p.parse_args()

    fig = plt.figure(1, (5, 5))
    root = fig.add_axes([0, 0, 1, 1])

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    # Layout format: (x, y, label, (chr lengths))
    anc = (0.5, 0.9, "Ancestor", (1,))
    s1 = (0.2, 0.6, "Genome I", (1,))
    s2 = (0.5, 0.6, "Genome II", (1,))
    s3 = (0.8, 0.6, "Genome III", (1,))
    tetra = (0.35, 0.4, "Tetraploid I / II", (0.5, 0.9))
    hexa = (0.5, 0.1, "Hexaploid I / II / III", (0.36, 0.46, 0.9))
    labels = (anc, s1, s2, s3, tetra, hexa)
    connections = (
        (anc, s1),
        (anc, s2),
        (anc, s3),
        (s1, tetra),
        (s2, tetra),
        (tetra, hexa),
        (s3, hexa),
    )

    xinterval = 0.02
    yratio = 0.05
    for xx, yy, label, chrl in labels:
        # RoundLabel(root, xx, yy, label)
        root.text(xx, yy, label, ha="center", va="center")
        offset = len(label) * 0.012
        for i, c in enumerate(chrl):
            ya = yy + yratio * c
            yb = yy - yratio * c
            Chromosome(root, xx - offset + i * xinterval, ya, yb, width=0.01)

    # Comments
    comments = ((0.15, 0.33, "II dominant"), (0.25, 0.03, "III dominant"))

    for xx, yy, c in comments:
        root.text(xx, yy, c, size=9, ha="center", va="center")

    # Branches
    tip = 0.04
    for a, b in connections:
        xa, ya, la, chra = a
        xb, yb, lb, chrb = b
        plt.plot((xa, xb), (ya - tip, yb + 2 * tip), "k-", lw=2, alpha=0.5)

    figname = fname() + ".pdf"
    savefig(figname, dpi=300)


if __name__ == "__main__":
    main()
