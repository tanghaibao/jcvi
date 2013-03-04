#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Catalog gene losses, and bites within genes.
"""

import sys

from optparse import OptionParser
from itertools import groupby

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed
from jcvi.utils.range import range_minmax, range_overlap
from jcvi.utils.cbook import gene_name
from jcvi.algorithms.synteny import add_beds, check_beds
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('loss', 'extract likely gene loss candidates'),
        ('genestatus', 'tag genes based on translation from GMAP models'),
        ('summary', 'provide summary of fractionation'),
        ('gaps', 'check gene locations against gaps'),
        # Specific study (requires specific datasets)
        ('napus', 'extract napus gene loss vs diploid ancestors'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def gaps(args):
    """
    %prog gaps idsfile fractionationfile gapsbed

    Check gene locations against gaps. `idsfile` contains a list of IDs to query
    into `fractionationfile` in order to get expected locations.
    """
    from jcvi.formats.base import DictFile
    from jcvi.apps.base import popen
    from jcvi.utils.cbook import percentage

    p = OptionParser(gaps.__doc__)
    p.add_option("--bdist", default=0, type="int",
                 help="Base pair distance [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    idsfile, frfile, gapsbed = args
    bdist = opts.bdist
    d =  DictFile(frfile, keypos=1, valuepos=2)
    bedfile = idsfile + ".bed"
    fw = open(bedfile, "w")
    fp = open(idsfile)
    total = 0
    for row in fp:
        id = row.strip()
        hit = d[id]
        tag, pos = get_tag(hit, None)
        seqid, start, end = pos
        start -= bdist
        end += bdist
        print >> fw, "\t".join(str(x) for x in (seqid, start - 1, end, id))
        total += 1
    fw.close()

    cmd = "intersectBed -a {0} -b {1} -v | wc -l".format(bedfile, gapsbed)
    not_in_gaps = popen(cmd).read()
    not_in_gaps = int(not_in_gaps)
    in_gaps = total - not_in_gaps
    print >> sys.stderr, "Ids in gaps: {1}".\
            format(total, percentage(in_gaps, total))


def genestatus(args):
    """
    %prog genestatus diploid.gff3.exon.ids

    Tag genes based on translation from GMAP models, using fasta.translate()
    --ids.
    """
    from itertools import groupby
    p = OptionParser(genestatus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    idsfile, = args
    fp = open(idsfile)
    data = []
    for row in fp:
        mRNA, label = row.split()
        labelatoms = label.split(",")
        if label == "complete" or label == "contain_ns,complete":
            tag = "complete"
        if "cannot_translate" in labelatoms:
            tag = "pseudogene"
        elif "five_prime_missing" in labelatoms or \
             "three_prime_missing" in labelatoms:
            tag = "partial"
        data.append((mRNA, label, tag))

    key = lambda x: x[0].split(".")[0]
    for gene, cc in groupby(data, key=key):
        cc = list(cc)
        tags = [x[-1] for x in cc]
        if "complete" in tags:
            tag = "complete"
        elif "partial" in tags:
            tag = "partial"
        else:
            tag = "pseudogene"
        print "\t".join((gene, tag))


def summary(args):
    """
    %prog summary diploid.napus.fractionation gmap.status

    Provide summary of fractionation. `fractionation` file is generated with
    loss(). `gmap.status` is generated with genestatus().
    """
    from jcvi.formats.base import DictFile
    from jcvi.utils.cbook import percentage

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    frfile, statusfile = args
    status = DictFile(statusfile)
    fp = open(frfile)
    inside, outside = 0, 0
    syntenic, non_syntenic = 0, 0
    s, ns, nf = 0, 0, 0
    complete, pseudogene, partial, gmap_fail = 0, 0, 0, 0
    partial_deletions = []
    for row in fp:
        seqid, gene, tag = row.split()
        if tag == '.':
            outside += 1
        else:
            inside += 1
            if tag[0] == '[':
                non_syntenic += 1
                if tag.startswith("[S]"):
                    s += 1
                    gstatus = status.get(gene, None)
                    if gstatus == 'complete':
                        complete += 1
                    elif gstatus == 'pseudogene':
                        pseudogene += 1
                    elif gstatus == 'partial':
                        partial_deletions.append(gene)
                        partial += 1
                    else:
                        gmap_fail += 1
                elif tag.startswith("[NS]"):
                    ns += 1
                elif tag.startswith("[NF]"):
                    nf += 1
            else:
                syntenic += 1

    m = "{0} inside synteny blocks\n".format(inside)
    m += "{0} outside synteny blocks\n".format(outside)
    m += "{0} has syntenic gene\n".format(syntenic)
    m += "{0} lack syntenic gene\n".format(non_syntenic)
    m += "{0} has sequence match in syntenic location\n".format(s)
    m += "{0} has sequence match in non-syntenic location\n".format(ns)
    m += "{0} has no sequence match\n".format(nf)
    m += "{0} syntenic sequence - complete model\n".format(percentage(complete, s))
    m += "{0} syntenic sequence - partial model\n".format(percentage(partial, s))
    m += "{0} syntenic sequence - pseudogene\n".format(percentage(pseudogene, s))
    m += "{0} syntenic sequence - gmap fail\n".format(percentage(gmap_fail, s))
    print >> sys.stderr, m

    fw = open("partial_deletions", "w")
    print >> fw, "\n".join(partial_deletions)
    fw.close()


def get_tag(name, order):
    if name[0] == '[':
        tag, tname = name[1:].split(']')
        seqid, se = tname.split(":")
        start, end = se.split("-")
        start, end = int(start), int(end)
    else:
        tag = None
        xi, x = order[name]
        seqid, start, end = x.seqid, x.start, x.end
    return tag, (seqid, start, end)


def napus(args):
    """
    %prog napus napus.bed brapa.boleracea.i1.blocks diploid.napus.fractionation

    Extract napus gene loss vs diploid ancestors. We are looking specifically
    for anything that has the pattern:

        BR - BO    or     BR - BO
        |                       |
        AN                     CN

    Step 1: extract BR - BO syntenic pairs
    Step 2: get diploid gene retention patterns from BR or BO as query
    Step 3: look for if AN or CN is NS(non-syntenic) or NF(not found) and
    specifically with NS, the NS location is actually the homeologous site.
    """
    p = OptionParser(napus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    napusbed, brbo, dpnp = args
    retention = {}
    fp = open(dpnp)
    for row in fp:
        seqid, query, hit = row.split()
        retention[query] = hit

    order = Bed(napusbed).order

    fp = open(brbo)
    for row in fp:
        br, bo = row.split()
        if '.' in (br, bo):
            continue
        an, cn = retention[br], retention[bo]
        row = "\t".join((br, bo, an, cn))
        if '.' in (an, cn):
            #print row
            continue

        # label loss candidates
        antag, anrange = get_tag(an, order)
        cntag, cnrange = get_tag(cn, order)

        if range_overlap(anrange, cnrange):
            if (antag, cntag) == ("NS", None):
                row = row + "\tAN LOST|{0}".format(br)
            if (antag, cntag) == (None, "NS"):
                row = row + "\tCN LOST|{0}".format(bo)

        print row


def region_str(region):
    return "{0}:{1}-{2}".format(*region)


def loss(args):
    """
    %prog loss a.b.i1.blocks a.b-genomic.blast

    Extract likely gene loss candidates between genome a and b.
    """
    p = OptionParser(loss.__doc__)
    p.add_option("--bed", default=False, action="store_true",
                 help="Genomic BLAST is in bed format [default: %default]")
    p.add_option("--gdist", default=20, type="int",
                 help="Gene distance [default: %default]")
    p.add_option("--bdist", default=20000, type="int",
                 help="Base pair distance [default: %default]")
    add_beds(p)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blocksfile, genomicblast = args
    gdist, bdist = opts.gdist, opts.bdist
    qbed, sbed, qorder, sorder, is_self = check_beds(blocksfile, p, opts)
    blocks = []
    fp = open(blocksfile)
    genetrack = {}
    proxytrack = {}
    for row in fp:
        a, b = row.split()
        genetrack[a] = b
        blocks.append((a, b))

    data = []
    for key, rows in groupby(blocks, key=lambda x: x[-1]):
        rows = list(rows)
        data.append((key, rows))

    imax = len(data) - 1
    for i, (key, rows) in enumerate(data):
        if i == 0 or i == imax:
            continue
        if key != '.':
            continue

        before, br = data[i - 1]
        after, ar = data[i + 1]
        bi, bx = sorder[before]
        ai, ax = sorder[after]
        dist = abs(bi - ai)
        if bx.seqid != ax.seqid or dist > gdist:
            continue

        start, end = range_minmax(((bx.start, bx.end), (ax.start, ax.end)))
        proxy = (bx.seqid, start - bdist, end + bdist)
        for a, b in rows:
            proxytrack[a] = proxy

    tags = {}
    if opts.bed:
        bed = Bed(genomicblast, sorted=False)
        key = lambda x: gene_name(x.accn.rsplit(".", 1)[0])
        for query, bb in groupby(bed, key=key):
            bb = list(bb)
            if query not in proxytrack:
                continue

            proxy = proxytrack[query]
            tag = "NS"
            best_b = bb[0]
            for b in bb:
                hsp = (b.seqid, b.start, b.end)
                if range_overlap(proxy, hsp):
                    tag = "S"
                    best_b = b
                    break

            hsp = (best_b.seqid, best_b.start, best_b.end)
            proxytrack[query] = hsp
            tags[query] = tag

    else:
        blast = Blast(genomicblast)
        for query, bb in blast.iter_hits():
            bb = list(bb)
            query = gene_name(query)
            if query not in proxytrack:
                continue

            proxy = proxytrack[query]
            tag = "NS"
            best_b = bb[0]
            for b in bb:
                hsp = (b.subject, b.sstart, b.sstop)
                if range_overlap(proxy, hsp):
                    tag = "S"
                    best_b = b
                    break

            hsp = (best_b.subject, best_b.sstart, best_b.sstop)
            proxytrack[query] = hsp
            tags[query] = tag

    for b in qbed:
        accn = b.accn
        target_region = genetrack[accn]
        if accn in proxytrack:
            target_region = region_str(proxytrack[accn])
            if accn in tags:
                ptag = "[{0}]".format(tags[accn])
            else:
                ptag = "[NF]"
            target_region = ptag + target_region

        print "\t".join((b.seqid, accn, target_region))


if __name__ == '__main__':
    main()
