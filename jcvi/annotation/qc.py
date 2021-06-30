#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run quality control (QC) on gene annotation. MAKER output was used during
testing. Several aspects of annotation QC are implemented in this script.

- Trim UTRs. MAKER sometimes predict UTRs that extend into other genes.
- Remove overlapping models.
"""
import sys

from jcvi.formats.gff import (
    Gff,
    get_piles,
    make_index,
    import_feats,
    populate_children,
    to_range,
)
from jcvi.formats.base import must_open
from jcvi.formats.sizes import Sizes
from jcvi.utils.range import range_minmax, range_chain, range_overlap
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ("trimUTR", "remove UTRs in the annotation set"),
        ("uniq", "remove overlapping gene models"),
        ("nmd", "identify transcript variant candidates for nonsense-mediated decay"),
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
        ranges = [to_range(x, score=gene_register[x.accn], id=x.accn) for x in group]
        selected_chain, score = range_chain(ranges)
        bestids |= set(x.id for x in selected_chain)

    removed = set(x.accn for x in allgenes) - bestids
    fw = open("removed.ids", "w")
    print("\n".join(sorted(removed)), file=fw)
    fw.close()
    populate_children(opts.outfile, bestids, gffile, "gene")


def get_cds_minmax(g, cid, level=2):
    cds = [x for x in g.children(cid, level) if x.featuretype == "CDS"]
    cdsranges = [(x.start, x.end) for x in cds]
    return range_minmax(cdsranges)


def trim(c, start, end, trim5=False, trim3=False, both=True):
    cstart, cend = c.start, c.end
    # Trim coordinates for feature c based on overlap to start and end
    if ((trim5 or both) and c.strand == "+") or ((trim3 or both) and c.strand == "-"):
        c.start = max(cstart, start)
    if ((trim3 or both) and c.strand == "+") or ((trim5 or both) and c.strand == "-"):
        c.end = min(cend, end)

    if c.start != cstart or c.end != cend:
        print(
            c.id,
            "trimmed [{0}, {1}] => [{2}, {3}]".format(cstart, cend, c.start, c.end),
            file=sys.stderr,
        )
    else:
        print(c.id, "no change", file=sys.stderr)


def reinstate(c, rc, trim5=False, trim3=False, both=True):
    cstart, cend = c.start, c.end
    # reinstate coordinates for feature `c` based on reference feature `refc`
    if ((trim5 or both) and c.strand == "+") or ((trim3 or both) and c.strand == "-"):
        c.start = rc.start
    if ((trim3 or both) and c.strand == "+") or ((trim5 or both) and c.strand == "-"):
        c.end = rc.end

    if c.start != cstart or c.end != cend:
        print(
            c.id,
            "reinstated [{0}, {1}] => [{2}, {3}]".format(cstart, cend, c.start, c.end),
            file=sys.stderr,
        )
    else:
        print(c.id, "no change", file=sys.stderr)


def cmp_children(cid, gff, refgff, cftype="CDS"):
    start, end = get_cds_minmax(gff, cid, level=1)
    rstart, rend = get_cds_minmax(refgff, cid, level=1)
    return (
        ((start == rstart) and (end == rend))
        and (
            len(list(gff.children(cid, featuretype=cftype)))
            == len(list(refgff.children(cid, featuretype=cftype)))
        )
        and (
            gff.children_bp(cid, child_featuretype=cftype)
            == refgff.children_bp(cid, child_featuretype=cftype)
        )
    )


def fprint(c, fw):
    if c.start > c.end:
        print(c.id, "destroyed [{0} > {1}]".format(c.start, c.end), file=sys.stderr)
    else:
        print(c, file=fw)


def trimUTR(args):
    """
    %prog trimUTR gffile

    Remove UTRs in the annotation set.

    If reference GFF3 is provided, reinstate UTRs from reference
    transcripts after trimming.

    Note: After running trimUTR, it is advised to also run
    `python -m jcvi.formats.gff fixboundaries` on the resultant GFF3
    to adjust the boundaries of all parent 'gene' features
    """
    import gffutils
    from jcvi.formats.base import SetFile

    p = OptionParser(trimUTR.__doc__)
    p.add_option(
        "--trim5",
        default=None,
        type="str",
        help="File containing gene list for 5' UTR trimming",
    )
    p.add_option(
        "--trim3",
        default=None,
        type="str",
        help="File containing gene list for 3' UTR trimming",
    )
    p.add_option(
        "--trimrange",
        default=None,
        type="str",
        help="File containing gene list for UTR trim back"
        + "based on suggested (start, stop) coordinate range",
    )
    p.add_option(
        "--refgff",
        default=None,
        type="str",
        help="Reference GFF3 used as fallback to replace UTRs",
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    gff = make_index(gffile)

    trim_both = False if (opts.trim5 or opts.trim3) else True
    trim5 = SetFile(opts.trim5) if opts.trim5 else set()
    trim3 = SetFile(opts.trim3) if opts.trim3 else set()
    trimrange = dict()
    if opts.trimrange:
        trf = must_open(opts.trimrange)
        for tr in trf:
            assert (
                len(tr.split("\t")) == 3
            ), "Must specify (start, stop) coordinate range"
            id, start, stop = tr.split("\t")
            trimrange[id] = (int(start), int(stop))
        trf.close()

    refgff = make_index(opts.refgff) if opts.refgff else None

    fw = must_open(opts.outfile, "w")
    for feat in gff.iter_by_parent_childs(
        featuretype="gene", order_by=("seqid", "start"), level=1
    ):
        for c in feat:
            cid, ctype, cparent = (
                c.id,
                c.featuretype,
                c.attributes.get("Parent", [None])[0],
            )
            t5, t3 = False, False
            if ctype == "gene":
                t5 = True if cid in trim5 else False
                t3 = True if cid in trim3 else False
                start, end = get_cds_minmax(gff, cid)
                trim(c, start, end, trim5=t5, trim3=t3, both=trim_both)
                fprint(c, fw)
            elif ctype == "mRNA":
                utr_types, extras = [], set()
                if any(id in trim5 for id in (cid, cparent)):
                    t5 = True
                    trim5.add(cid)
                if any(id in trim3 for id in (cid, cparent)):
                    t3 = True
                    trim3.add(cid)
                refc = None
                if refgff:
                    try:
                        refc = refgff[cid]
                        refctype = refc.featuretype
                        refptype = refgff[refc.attributes["Parent"][0]].featuretype
                        if refctype == "mRNA" and refptype == "gene":
                            if cmp_children(cid, gff, refgff, cftype="CDS"):
                                reinstate(c, refc, trim5=t5, trim3=t3, both=trim_both)
                                if t5:
                                    utr_types.append("five_prime_UTR")
                                if t3:
                                    utr_types.append("three_prime_UTR")
                                for utr_type in utr_types:
                                    for utr in refgff.children(
                                        refc, featuretype=utr_type
                                    ):
                                        extras.add(utr)
                                        for exon in refgff.region(
                                            region=utr, featuretype="exon"
                                        ):
                                            if exon.attributes["Parent"][0] == cid:
                                                extras.add(exon)
                        else:
                            refc = None
                    except gffutils.exceptions.FeatureNotFoundError:
                        pass
                start, end = get_cds_minmax(gff, cid, level=1)
                if cid in trimrange:
                    start, end = range_minmax([trimrange[cid], (start, end)])
                if not refc:
                    trim(c, start, end, trim5=t5, trim3=t3, both=trim_both)
                fprint(c, fw)
                for cc in gff.children(cid, order_by="start"):
                    _ctype = cc.featuretype
                    if _ctype not in utr_types:
                        if _ctype != "CDS":
                            if _ctype == "exon":
                                eskip = [
                                    range_overlap(to_range(cc), to_range(x))
                                    for x in extras
                                    if x.featuretype == "exon"
                                ]
                                if any(eskip):
                                    continue
                            trim(cc, start, end, trim5=t5, trim3=t3, both=trim_both)
                            fprint(cc, fw)
                        else:
                            fprint(cc, fw)
                for x in extras:
                    fprint(x, fw)
    fw.close()


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
    from jcvi.utils.cbook import enumerate_reversed

    p = OptionParser(nmd.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (gffile,) = args
    gff = make_index(gffile)

    fw = must_open(opts.outfile, "w")
    for gene in gff.features_of_type("gene", order_by=("seqid", "start")):
        _enumerate = enumerate if gene.strand == "-" else enumerate_reversed
        for mrna in gff.children(gene, featuretype="mRNA", order_by="start"):
            tracker = dict()
            tracker["exon"] = list(
                gff.children(mrna, featuretype="exon", order_by="start")
            )
            tracker["cds"] = [None] * len(tracker["exon"])

            tcds_pos = None
            for i, exon in _enumerate(tracker["exon"]):
                for cds in gff.region(
                    region=exon, featuretype="CDS", completely_within=True
                ):
                    if mrna.id in cds["Parent"]:
                        tracker["cds"][i] = cds
                        tcds_pos = i
                        break
                if tcds_pos:
                    break

            NMD, distance = False, 0
            if (mrna.strand == "+" and tcds_pos + 1 < len(tracker["exon"])) or (
                mrna.strand == "-" and tcds_pos - 1 >= 0
            ):
                tcds = tracker["cds"][tcds_pos]
                texon = tracker["exon"][tcds_pos]

                PTC = tcds.end if mrna.strand == "+" else tcds.start
                TDSS = texon.end if mrna.strand == "+" else texon.start
                distance = abs(TDSS - PTC)
                NMD = True if distance > 50 else False

            print(
                "\t".join(
                    str(x)
                    for x in (
                        gene.id,
                        mrna.id,
                        gff.children_bp(mrna, child_featuretype="CDS"),
                        distance,
                        NMD,
                    )
                ),
                file=fw,
            )

    fw.close()


if __name__ == "__main__":
    main()
