#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""

from copy import deepcopy

from jcvi.graphics.base import FancyArrow, plt, savefig, panel_labels
from jcvi.graphics.glyph import CartoonRegion
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('cartoon', 'generate cartoon illustration of SynFind'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def cartoon(args):
    """
    %prog synteny.py

    Generate cartoon illustration of SynFind.
    """
    from jcvi.projects.misc import plot_diagram

    p = OptionParser(cartoon.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x6")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Panel A
    A = CartoonRegion(41)
    A.draw(root, .35, .85)
    x1, x2 = A.x1, A.x2
    lsg = "lightslategray"
    p = FancyArrow(.35, .88, x1 - .35, 0, fc=lsg, lw=0, shape="left",
                   length_includes_head=True, width=.01,
                   head_length=abs(x1 - .35) * .15, head_width=.03)
    root.add_patch(p)
    p = FancyArrow(.35, .88, x2 - .35, 0, fc=lsg, lw=0, shape="right",
                   length_includes_head=True, width=.01,
                   head_length=abs(x2 - .35) * .15, head_width=.03)
    root.add_patch(p)

    # Panel B
    A.draw(root, .35, .7)

    a = deepcopy(A)
    a.evolve(target=10)
    a.draw(root, .35, .6, strip=True)
    b = deepcopy(A)
    b.evolve(target=8)
    b.draw(root, .35, .56, strip=True)
    c = deepcopy(A)
    c.evolve(target=6)
    c.draw(root, .35, .52, strip=True)

    # Panel C
    plot_diagram(root, .2, .2, "S", "syntenic", gradient=False)
    plot_diagram(root, .4, .2, "F", "missing, with both flankers",
                 gradient=False)
    plot_diagram(root, .6, .2, "G", "missing, with one flanker", gradient=False)

    labels = ((.04, .95, 'A'), (.04, .75, 'B'), (.04, .4, 'C'))
    panel_labels(root, labels)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
