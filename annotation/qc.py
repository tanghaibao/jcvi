#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run quality control (QC) on gene annotation. MAKER output was used during
testing. Several aspects of annotation QC are implemented in this script.

- Trim UTRs. MAKER sometimes predict UTRs that extend into other genes.
- Remove overlapping models.
"""

import sys

from jcvi.formats.gff import Gff, get_piles, make_index, import_feats, \
            populate_children
from jcvi.formats.base import must_open
from jcvi.formats.sizes import Sizes
from jcvi.utils.range import Range, range_minmax, range_chain
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('trimUTR', 'remove UTRs in the annotation set'),
        ('uniq', 'remove overlapping gene models'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def uniq(args):
    """
    %prog uniq gffile cdsfasta

    Remove overlapping gene models. Similar to formats.gff.uniq(), overlapping
    'piles' are processed, one by one.

    Here, we use a different algorithm, that retains the best non-overlapping
    subset witin each pile, rather than single best model. Scoring function is
    also different, rather than based on score or span, we optimize for the
    subset that show the best combined score. Score is defined by:

    score = (1 - AED) * length
    """
    p = OptionParser(uniq.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gffile, cdsfasta = args
    gff = Gff(gffile)
    sizes = Sizes(cdsfasta).mapping
    gene_register = {}
    for g in gff:
        if g.type != "mRNA":
            continue
        aed = float(g.attributes["_AED"][0])
        gene_register[g.parent] = (1 - aed) * sizes[g.accn]

    allgenes = import_feats(gffile)
    g = get_piles(allgenes)

    bestids = set()
    for group in g:
        ranges = [Range(x.seqid, x.start, x.end, \
                    gene_register[x.accn], x.accn) for x in group]
        selected_chain, score = range_chain(ranges)
        bestids |= set(x.id for x in selected_chain)

    removed = set(x.accn for x in allgenes) - bestids
    fw = open("removed.ids", "w")
    print >> fw, "\n".join(sorted(removed))
    fw.close()
    populate_children(opts.outfile, bestids, gffile, "gene")


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
