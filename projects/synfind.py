#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""

from copy import deepcopy

from jcvi.graphics.base import FancyArrow, plt, savefig, panel_labels, markup
from jcvi.graphics.glyph import CartoonRegion, RoundRect
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('cartoon', 'generate cartoon illustration of SynFind'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def plot_diagram(ax, x, y, A, B, tag, label):
    ax.text(x, y + .14, "{0}: {1}".format(tag, label), ha="center")
    strip = tag != 'G'
    A.draw(ax, x, y + .06, gene_len=.02, strip=strip)
    B.draw(ax, x, y, gene_len=.02, strip=strip)


def cartoon(args):
    """
    %prog synteny.py

    Generate cartoon illustration of SynFind.
    """
    p = OptionParser(cartoon.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Panel A
    A = CartoonRegion(41)
    A.draw(root, .35, .85, strip=False, color=False)
    x1, x2 = A.x1, A.x2
    lsg = "lightslategray"
    pad = .01
    xc, yc = .35, .88
    arrowlen = x2 - xc - pad
    arrowprops = dict(length_includes_head=True, width=.01, fc=lsg, lw=0,
                      head_length=arrowlen * .15, head_width=.03)
    p = FancyArrow(xc - pad, yc, -arrowlen, 0, shape="left", **arrowprops)
    root.add_patch(p)
    p = FancyArrow(xc + pad, yc, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    yt = yc + 4 * pad
    root.text((x1 + xc) / 2, yt, "20 genes upstream", ha="center")
    root.text((x2 + xc) / 2, yt, "20 genes downstream", ha="center")
    root.plot((xc,), (yc,), "o", mfc='w', mec=lsg, mew=2, lw=2, color=lsg)
    root.text(xc, yt, "Query gene", ha="center")

    # Panel B
    A.draw(root, .35, .7, strip=False)

    RoundRect(root, (.07, .49), .56, .14, fc='y', alpha=.2)
    a = deepcopy(A)
    a.evolve(mode='S', target=10)
    a.draw(root, .35, .6)
    b = deepcopy(A)
    b.evolve(mode='F', target=8)
    b.draw(root, .35, .56)
    c = deepcopy(A)
    c.evolve(mode='G', target=6)
    c.draw(root, .35, .52)

    for x in (a, b, c):
        root.text(.64, x.y, "Score={0}".format(x.nonwhites), va="center")

    # Panel C
    A.truncate_between_flankers()
    a.truncate_between_flankers()
    b.truncate_between_flankers()
    c.truncate_between_flankers(target=6)

    plot_diagram(root, .14, .2, A, a, "S", "syntenic")
    plot_diagram(root, .37, .2, A, b, "F", "missing, with both flankers")
    plot_diagram(root, .6, .2, A, c, "G", "missing, with one flanker")

    labels = ((.04, .95, 'A'), (.04, .75, 'B'), (.04, .4, 'C'))
    panel_labels(root, labels)

    # Descriptions
    xt = .85
    desc = ("Extract neighborhood",
            "of *window* size",
            "Generate anchors",
            "Group anchors within *window* apart",
            "Find regions above *score* cutoff",
            "Identify flankers",
            "Annotate syntelog class"
            )
    for yt, t in zip((.88, .84, .64, .6, .56, .3, .26), desc):
        root.text(xt, yt, markup(t), ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
