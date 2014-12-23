#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run quality control (QC) on gene annotation. MAKER output was used during
testing. Several aspects of annotation QC are implemented in this script.

- Trim UTRs. MAKER sometimes predict UTRs that extend into other genes.
- Remove overlapping models.
"""

import sys

from jcvi.formats.gff import Gff, make_index
from jcvi.formats.base import must_open
from jcvi.utils.range import range_minmax
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('trimUTR', 'remove UTRs in the annotation set'),
        ('uniq', 'remove overlapping gene models'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_cds_minmax(g, cid, level=2):
    cds = [x for x in g.children(cid, level) if x.featuretype == "CDS"]
    cdsranges = [(x.start, x.end) for x in cds]
    return range_minmax(cdsranges)


def trim(c, start, end):
    cstart, cend = c.start, c.end
    # Trim coordinates for feature c based on overlap to start and end
    c.start, c.end = max(cstart, start), min(cend, end)
    if c.start != cstart or c.end != cend:
        print >> sys.stderr, c.id, \
                "[{0}, {1}] => [{2}, {3}]".format(cstart, cend, c.start, c.end)
    else:
        print >> sys.stderr, c.id, "no change"


def trimUTR(args):
    """
    %prog trimUTR gffile

    Remove UTRs in the annotation set.
    """
    p = OptionParser(trimUTR.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    g = make_index(gffile)
    gff = Gff(gffile)
    mRNA_register = {}
    fw = must_open(opts.outfile, "w")
    for c in gff:
        cid, ctype = c.accn, c.type
        if ctype == "gene":
            start, end = get_cds_minmax(g, cid)
            trim(c, start, end)
        elif ctype == "mRNA":
            start, end = get_cds_minmax(g, cid, level=1)
            trim(c, start, end)
            mRNA_register[cid] = (start, end)
        elif ctype != "CDS":
            start, end = mRNA_register[c.parent]
            trim(c, start, end)
        if c.start > c.end:
            print >> sys.stderr, cid, \
                    "destroyed [{0} > {1}]".format(c.start, c.end)
        else:
            print >> fw, c


if __name__ == '__main__':
    main()
