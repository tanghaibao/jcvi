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
        ('nmd', 'identify transcript variant candidates for nonsense-mediated decay')
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


def trim(c, start, end, trim5=False, trim3=False, both=True):
    cstart, cend = c.start, c.end
    # Trim coordinates for feature c based on overlap to start and end
    if ((trim5 or both) and c.strand == '+') \
            or ((trim3 or both) and c.strand == '-'):
        c.start = max(cstart, start)
    if ((trim3 or both) and c.strand == '+') \
            or ((trim5 or both) and c.strand == '-'):
        c.end = min(cend, end)

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
    from jcvi.formats.base import SetFile

    p = OptionParser(trimUTR.__doc__)
    p.add_option("--trim5", default=None, type="str", \
        help="File containing gene list for 5' UTR trimming")
    p.add_option("--trim3", default=None, type="str", \
        help="File containing gene list for 3' UTR trimming")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    g = make_index(gffile)
    gff = Gff(gffile)

    trim_both = False if (opts.trim5 or opts.trim3) else True
    trim5 = SetFile(opts.trim5) if opts.trim5 else set()
    trim3 = SetFile(opts.trim3) if opts.trim3 else set()

    mRNA_register = {}
    fw = must_open(opts.outfile, "w")
    for c in gff:
        cid, ctype = c.accn, c.type
        t5, t3 = False, False
        if ctype == "gene":
            t5 = True if cid in trim5 else False
            t3 = True if cid in trim3 else False
            start, end = get_cds_minmax(g, cid)
            trim(c, start, end, trim5=t5, trim3=t3, both=trim_both)
        elif ctype == "mRNA":
            if any(id in trim5 for id in (cid, c.parent)):
                t5 = True
                trim5.add(cid)
            if any(id in trim3 for id in (cid, c.parent)):
                t3 = True
                trim3.add(cid)
            start, end = get_cds_minmax(g, cid, level=1)
            trim(c, start, end, trim5=t5, trim3=t3, both=trim_both)
            mRNA_register[cid] = (start, end)
        elif ctype != "CDS":
            t5 = True if c.parent in trim5 else False
            t3 = True if c.parent in trim3 else False
            start, end = mRNA_register[c.parent]
            trim(c, start, end, trim5=t5, trim3=t3, both=trim_both)
        if c.start > c.end:
            print >> sys.stderr, cid, \
                    "destroyed [{0} > {1}]".format(c.start, c.end)
        else:
            print >> fw, c


def nmd(args):
    """
    %prog nmd gffile

    Identify transcript variants which might be candidates for nonsense
    mediated decay (NMD)

    A transcript is considered to be a candidate for NMD when the CDS stop
    codon is located more than 50nt upstream of terminal splice site donor

    References:
    http://www.nature.com/horizon/rna/highlights/figures/s2_spec1_f3.html
    http://www.biomedcentral.com/1741-7007/7/23/figure/F1
    """
    import __builtin__
    from jcvi.utils.cbook import enumerate_reversed

    p = OptionParser(nmd.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gffile, = args
    gff = make_index(gffile)

    fw = must_open(opts.outfile, "w")
    for gene in gff.features_of_type('gene', order_by=('seqid', 'start')):
        _enumerate = __builtin__.enumerate if gene.strand == "-" else enumerate_reversed
        for mrna in gff.children(gene, featuretype='mRNA', order_by=('start')):
            tracker = dict()
            tracker['exon'] = list(gff.children(mrna, featuretype='exon', order_by=('start')))
            tracker['cds'] = [None] * len(tracker['exon'])

            tcds_pos = None
            for i, exon in _enumerate(tracker['exon']):
                for cds in gff.region(region=exon, featuretype='CDS', completely_within=True):
                    if mrna.id in cds['Parent']:
                        tracker['cds'][i] = cds
                        tcds_pos = i
                        break
                if tcds_pos: break

            NMD, distance = False, 0
            if (mrna.strand == "+" and tcds_pos + 1 < len(tracker['exon'])) \
                or (mrna.strand == "-" and tcds_pos - 1 >= 0):
                tcds = tracker['cds'][tcds_pos]
                texon = tracker['exon'][tcds_pos]

                PTC = tcds.end if mrna.strand == '+' else tcds.start
                TDSS = texon.end if mrna.strand == '+' else texon.start
                distance = abs(TDSS - PTC)
                NMD = True if distance > 50 else False

            print >> fw, "\t".join(str(x) for x in (gene.id, mrna.id, \
                gff.children_bp(mrna, child_featuretype='CDS'), distance, NMD))

    fw.close()


if __name__ == '__main__':
    main()
