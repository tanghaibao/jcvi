#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Catalog gene losses, and bites within genes.
"""
import sys
import logging

from itertools import groupby

from jcvi.formats.blast import Blast
from jcvi.formats.bed import Bed
from jcvi.utils.range import range_minmax, range_overlap, range_distance
from jcvi.utils.cbook import gene_name
from jcvi.utils.grouper import Grouper
from jcvi.compara.synteny import check_beds
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (
        # Identify true gene loss
        ("loss", "extract likely gene loss candidates"),
        ("validate", "confirm synteny loss against CDS bed overlaps"),
        ("summary", "provide summary of fractionation"),
        ("gaps", "check gene locations against gaps"),
        # Gene specific status
        ("gffselect", "dump gff for the missing genes"),
        ("genestatus", "tag genes based on translation from GMAP models"),
        # Specific study for napus (requires specific datasets)
        ("napus", "extract gene loss vs diploid ancestors (napus)"),
        ("merge", "merge protein quartets table with registry (napus)"),
        ("segment", "merge adjacent gene loss into segmental loss (napus)"),
        ("offdiag", "find gene pairs that are off diagonal"),
        ("diff", "calculate diff of size of syntenic regions"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def offdiag(args):
    """
    %prog offdiag diploid.napus.1x1.lifted.anchors

    Find gene pairs that are off diagnoal. "Off diagonal" are the pairs that are
    not on the orthologous chromosomes. For example, napus chrA01 and brapa A01.
    """
    p = OptionParser(offdiag.__doc__)
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorsfile,) = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts)

    fp = open(anchorsfile)
    pf = "-".join(anchorsfile.split(".")[:2])
    header = "Block-id|Napus|Diploid|Napus-chr|Diploid-chr|RBH?".split("|")
    print("\t".join(header))
    i = -1
    for row in fp:
        if row[0] == "#":
            i += 1
            continue
        q, s, score = row.split()
        rbh = "no" if score[-1] == "L" else "yes"
        qi, qq = qorder[q]
        si, ss = sorder[s]
        oqseqid = qseqid = qq.seqid
        osseqid = sseqid = ss.seqid
        sseqid = sseqid.split("_")[0][-3:]
        if qseqid[0] == "A":
            qseqid = qseqid[-3:]  # A09 => A09
        elif qseqid[0] == "C":
            qseqid = "C0" + qseqid[-1]  # C9 => C09
        else:
            continue
        if qseqid == sseqid or sseqid[-2:] == "nn":
            continue
        block_id = pf + "-block-{0}".format(i)
        print("\t".join((block_id, q, s, oqseqid, osseqid, rbh)))


def diff(args):
    """
    %prog diff simplefile

    Calculate difference of pairwise syntenic regions.
    """
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(diff.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (simplefile,) = args
    fp = open(simplefile)
    data = [x.split() for x in fp]
    spans = []
    for block_id, ab in groupby(data[1:], key=lambda x: x[0]):
        a, b = list(ab)
        aspan, bspan = a[4], b[4]
        aspan, bspan = int(aspan), int(bspan)
        spans.append((aspan, bspan))
    aspans, bspans = zip(*spans)
    dspans = [b - a for a, b, in spans]
    s = SummaryStats(dspans)
    print("For a total of {0} blocks:".format(len(dspans)), file=sys.stderr)
    print("Sum of A: {0}".format(sum(aspans)), file=sys.stderr)
    print("Sum of B: {0}".format(sum(bspans)), file=sys.stderr)
    print("Sum of Delta: {0} ({1})".format(sum(dspans), s), file=sys.stderr)


def estimate_size(accns, bed, order, conservative=True):
    """
    Estimate the bp length for the deletion tracks, indicated by the gene accns.
    True different levels of estimates vary on conservativeness.
    """
    accns = [order[x] for x in accns]
    ii, bb = zip(*accns)
    mini, maxi = min(ii), max(ii)
    if not conservative:  # extend one gene
        mini -= 1
        maxi += 1
    minb = bed[mini]
    maxb = bed[maxi]
    assert minb.seqid == maxb.seqid
    distmode = "ss" if conservative else "ee"
    ra = (minb.seqid, minb.start, minb.end, "+")
    rb = (maxb.seqid, maxb.start, maxb.end, "+")

    dist, orientation = range_distance(ra, rb, distmode=distmode)
    assert dist != -1
    return dist


def segment(args):
    """
    %prog segment loss.ids bedfile

    Merge adjacent gene loss into segmental loss.

    Then based on the segmental loss, estimate amount of DNA loss in base pairs.
    Two estimates can be given:
    - conservative: just within the start and end of a single gene
    - aggressive: extend the deletion track to the next gene

    The real deletion size is within these estimates.
    """
    from jcvi.formats.base import SetFile

    p = OptionParser(segment.__doc__)
    p.add_option(
        "--chain",
        default=1,
        type="int",
        help="Allow next N genes to be chained",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    idsfile, bedfile = args
    bed = Bed(bedfile)
    order = bed.order
    ids = SetFile(idsfile)
    losses = Grouper()
    skip = opts.chain
    for i, a in enumerate(bed):
        a = a.accn
        for j in range(i + 1, i + 1 + skip):
            if j >= len(bed):
                break
            b = bed[j].accn
            if a in ids:
                losses.join(a, a)
            if a in ids and b in ids:
                losses.join(a, b)

    losses = list(losses)
    singletons = [x for x in losses if len(x) == 1]
    segments = [x for x in losses if len(x) > 1]
    ns, nm, nt = len(singletons), len(segments), len(losses)
    assert ns + nm == nt

    # Summary for all segments
    for x in sorted(singletons) + sorted(segments):
        print(
            "\t".join(
                str(x)
                for x in ("|".join(sorted(x)), len(x), estimate_size(x, bed, order))
            )
        )

    # Find longest segment stretch
    if segments:
        mx, maxsegment = max([(len(x), x) for x in segments])
        print("Longest stretch: run of {0} genes".format(mx), file=sys.stderr)
        print("  {0}".format("|".join(sorted(maxsegment))), file=sys.stderr)
        seg_asize = sum(estimate_size(x, bed, order) for x in segments)
        seg_bsize = sum(
            estimate_size(x, bed, order, conservative=False) for x in segments
        )
    else:
        seg_asize = seg_bsize = 0

    sing_asize = sum(estimate_size(x, bed, order) for x in singletons)
    sing_bsize = sum(
        estimate_size(x, bed, order, conservative=False) for x in singletons
    )
    total_asize = sing_asize + seg_asize
    total_bsize = sing_bsize + seg_bsize
    print(
        "Singleton ({0}): {1} - {2} bp".format(ns, sing_asize, sing_bsize),
        file=sys.stderr,
    )
    print(
        "Segment ({0}): {1} - {2} bp".format(nm, seg_asize, seg_bsize), file=sys.stderr
    )
    print(
        "Total ({0}): {1} - {2} bp".format(nt, total_asize, total_bsize),
        file=sys.stderr,
    )
    print(
        "Average ({0}): {1} bp".format(nt, (total_asize + total_bsize) / 2),
        file=sys.stderr,
    )


def merge(args):
    """
    %prog merge protein-quartets registry LOST

    Merge protein quartets table with dna quartets registry. This is specific
    to the napus project.
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(merge.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    quartets, registry, lost = args
    qq = DictFile(registry, keypos=1, valuepos=3)
    lost = DictFile(lost, keypos=1, valuepos=0, delimiter="|")
    qq.update(lost)
    fp = open(quartets)
    cases = {
        "AN,CN": 4,
        "BO,AN,CN": 8,
        "BO,CN": 2,
        "BR,AN": 1,
        "BR,AN,CN": 6,
        "BR,BO": 3,
        "BR,BO,AN": 5,
        "BR,BO,AN,CN": 9,
        "BR,BO,CN": 7,
    }
    ip = {
        "syntenic_model": "Syntenic_model_excluded_by_OMG",
        "complete": "Predictable",
        "partial": "Truncated",
        "pseudogene": "Pseudogene",
        "random": "Match_random",
        "real_ns": "Transposed",
        "gmap_fail": "GMAP_fail",
        "AN LOST": "AN_LOST",
        "CN LOST": "CN_LOST",
        "BR LOST": "BR_LOST",
        "BO LOST": "BO_LOST",
        "outside": "Outside_synteny_blocks",
        "[NF]": "Not_found",
    }
    for row in fp:
        atoms = row.strip().split("\t")
        genes = atoms[:4]
        tag = atoms[4]
        a, b, c, d = [qq.get(x, ".").rsplit("-", 1)[-1] for x in genes]
        qqs = [c, d, a, b]
        for i, q in enumerate(qqs):
            if atoms[i] != ".":
                qqs[i] = "syntenic_model"
        # Make comment
        comment = "Case{0}".format(cases[tag])
        dots = sum([1 for x in genes if x == "."])
        if dots == 1:
            idx = genes.index(".")
            status = qqs[idx]
            status = ip[status]
            comment += "-" + status
        print(row.strip() + "\t" + "\t".join(qqs + [comment]))


def gffselect(args):
    """
    %prog gffselect gmaplocation.bed expectedlocation.bed translated.ids tag

    Try to match up the expected location and gmap locations for particular
    genes. translated.ids was generated by fasta.translate --ids. tag must be
    one of "complete|pseudogene|partial".
    """
    from jcvi.formats.bed import intersectBed_wao

    p = OptionParser(gffselect.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    gmapped, expected, idsfile, tag = args
    data = get_tags(idsfile)
    completeness = dict((a.replace("mrna", "path"), c) for (a, b, c) in data)

    seen = set()
    idsfile = expected.rsplit(".", 1)[0] + ".ids"
    fw = open(idsfile, "w")
    cnt = 0
    for a, b in intersectBed_wao(expected, gmapped):
        if b is None:
            continue
        aname, bbname = a.accn, b.accn
        bname = bbname.split(".")[0]
        if completeness[bbname] != tag:
            continue
        if aname == bname:
            if bname in seen:
                continue
            seen.add(bname)
            print(bbname, file=fw)
            cnt += 1
    fw.close()

    logging.debug("Total {0} records written to `{1}`.".format(cnt, idsfile))


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
    p.add_option("--bdist", default=0, type="int", help="Base pair distance")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    idsfile, frfile, gapsbed = args
    bdist = opts.bdist
    d = DictFile(frfile, keypos=1, valuepos=2)
    bedfile = idsfile + ".bed"
    fw = open(bedfile, "w")
    fp = open(idsfile)
    total = 0
    for row in fp:
        id = row.strip()
        hit = d[id]
        tag, pos = get_tag(hit, None)
        seqid, start, end = pos
        start, end = max(start - bdist, 1), end + bdist
        print("\t".join(str(x) for x in (seqid, start - 1, end, id)), file=fw)
        total += 1
    fw.close()

    cmd = "intersectBed -a {0} -b {1} -v | wc -l".format(bedfile, gapsbed)
    not_in_gaps = popen(cmd).read()
    not_in_gaps = int(not_in_gaps)
    in_gaps = total - not_in_gaps
    print("Ids in gaps: {1}".format(total, percentage(in_gaps, total)), file=sys.stderr)


def get_tags(idsfile):
    fp = open(idsfile)
    data = []
    for row in fp:
        mRNA, label = row.split()
        labelatoms = label.split(",")
        if label == "complete" or label == "contain_ns,complete":
            tag = "complete"
        if "cannot_translate" in labelatoms:
            tag = "pseudogene"
        elif "five_prime_missing" in labelatoms or "three_prime_missing" in labelatoms:
            tag = "partial"
        data.append((mRNA, label, tag))
    return data


def genestatus(args):
    """
    %prog genestatus diploid.gff3.exon.ids

    Tag genes based on translation from GMAP models, using fasta.translate()
    --ids.
    """
    p = OptionParser(genestatus.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (idsfile,) = args
    data = get_tags(idsfile)
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
        print("\t".join((gene, tag)))


def summary(args):
    """
    %prog summary diploid.napus.fractionation gmap.status

    Provide summary of fractionation. `fractionation` file is generated with
    loss(). `gmap.status` is generated with genestatus().
    """
    from jcvi.formats.base import DictFile
    from jcvi.utils.cbook import percentage, Registry

    p = OptionParser(summary.__doc__)
    p.add_option("--extra", help="Cross with extra tsv file")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    frfile, statusfile = args
    status = DictFile(statusfile)
    fp = open(frfile)
    registry = Registry()  # keeps all the tags for any given gene
    for row in fp:
        seqid, gene, tag = row.split()
        if tag == ".":
            registry[gene].append("outside")
        else:
            registry[gene].append("inside")
            if tag[0] == "[":
                registry[gene].append("no_syntenic_model")
                if tag.startswith("[S]"):
                    registry[gene].append("[S]")
                    gstatus = status.get(gene, None)
                    if gstatus == "complete":
                        registry[gene].append("complete")
                    elif gstatus == "pseudogene":
                        registry[gene].append("pseudogene")
                    elif gstatus == "partial":
                        registry[gene].append("partial")
                    else:
                        registry[gene].append("gmap_fail")
                elif tag.startswith("[NS]"):
                    registry[gene].append("[NS]")
                    if "random" in tag or "Scaffold" in tag:
                        registry[gene].append("random")
                    else:
                        registry[gene].append("real_ns")
                elif tag.startswith("[NF]"):
                    registry[gene].append("[NF]")
            else:
                registry[gene].append("syntenic_model")

    inside = registry.count("inside")
    outside = registry.count("outside")
    syntenic = registry.count("syntenic_model")
    non_syntenic = registry.count("no_syntenic_model")
    s = registry.count("[S]")
    ns = registry.count("[NS]")
    nf = registry.count("[NF]")
    complete = registry.count("complete")
    pseudogene = registry.count("pseudogene")
    partial = registry.count("partial")
    gmap_fail = registry.count("gmap_fail")
    random = registry.count("random")
    real_ns = registry.count("real_ns")

    complete_models = registry.get_tag("complete")
    pseudogenes = registry.get_tag("pseudogene")
    partial_deletions = registry.get_tag("partial")

    m = "{0} inside synteny blocks\n".format(inside)
    m += "{0} outside synteny blocks\n".format(outside)
    m += "{0} has syntenic gene\n".format(syntenic)
    m += "{0} lack syntenic gene\n".format(non_syntenic)
    m += "{0} has sequence match in syntenic location\n".format(s)
    m += "{0} has sequence match in non-syntenic location\n".format(ns)
    m += "{0} has sequence match in un-ordered scaffolds\n".format(random)
    m += "{0} has sequence match in real non-syntenic location\n".format(real_ns)
    m += "{0} has no sequence match\n".format(nf)
    m += "{0} syntenic sequence - complete model\n".format(percentage(complete, s))
    m += "{0} syntenic sequence - partial model\n".format(percentage(partial, s))
    m += "{0} syntenic sequence - pseudogene\n".format(percentage(pseudogene, s))
    m += "{0} syntenic sequence - gmap fail\n".format(percentage(gmap_fail, s))
    print(m, file=sys.stderr)

    aa = ["complete_models", "partial_deletions", "pseudogenes"]
    bb = [complete_models, partial_deletions, pseudogenes]
    for a, b in zip(aa, bb):
        fw = open(a, "w")
        print("\n".join(b), file=fw)
        fw.close()

    extra = opts.extra
    if extra:
        registry.update_from(extra)

    fp.seek(0)
    fw = open("registry", "w")
    for row in fp:
        seqid, gene, tag = row.split()
        ts = registry[gene]
        print("\t".join((seqid, gene, tag, "-".join(ts))), file=fw)
    fw.close()

    logging.debug("Registry written.")


def get_tag(name, order):
    if name[0] == "[":
        tag, tname = name[1:].split("]")
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
    Step 4: categorize gene losses into singleton, or segmental (defined as
    consecutive losses with a maximum skip of 1
    """
    from jcvi.utils.cbook import SummaryStats

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

    quartetsfile = "quartets"
    fp = open(brbo)
    fw = open(quartetsfile, "w")
    AL = "AN LOST"
    CL = "CN LOST"
    for row in fp:
        br, bo = row.split()
        if "." in (br, bo):
            continue
        an, cn = retention[br], retention[bo]
        row = "\t".join((br, bo, an, cn))
        if "." in (an, cn):
            # print row
            continue

        # label loss candidates
        antag, anrange = get_tag(an, order)
        cntag, cnrange = get_tag(cn, order)

        if range_overlap(anrange, cnrange):
            if (antag, cntag) == ("NS", None):
                row = row + "\t{0}|{1}".format(AL, br)
            if (antag, cntag) == (None, "NS"):
                row = row + "\t{0}|{1}".format(CL, bo)

        print(row, file=fw)
    fw.close()

    logging.debug("Quartets and gene losses written to `{0}`.".format(quartetsfile))

    # Parse the quartets file to extract singletons vs.segmental losses
    fp = open(quartetsfile)
    fw = open(quartetsfile + ".summary", "w")
    data = [x.rstrip().split("\t") for x in fp]
    skip = 1  # max distance between losses

    g = Grouper()
    losses = [(len(x) == 5) for x in data]
    for i, d in enumerate(losses):
        if not d:
            continue
        g.join(i, i)
        itag = data[i][-1].split("|")[0]
        for j in range(i + 1, i + skip + 1):
            jtag = data[j][-1].split("|")[0]
            if j < len(losses) and losses[j] and itag == jtag:
                g.join(i, j)

    losses = list(g)
    singletons = [x for x in losses if len(x) == 1]
    segments = [x for x in losses if len(x) > 1]
    ns, nm = len(singletons), len(segments)
    assert len(losses) == ns + nm

    grab_tag = lambda pool, tag: [
        x for x in pool if all(data[z][-1].startswith(tag) for z in x)
    ]

    an_loss_singletons = grab_tag(singletons, AL)
    cn_loss_singletons = grab_tag(singletons, CL)
    als, cls = len(an_loss_singletons), len(cn_loss_singletons)

    an_loss_segments = grab_tag(segments, AL)
    cn_loss_segments = grab_tag(segments, CL)
    alm, clm = len(an_loss_segments), len(cn_loss_segments)
    mixed = len(segments) - alm - clm
    assert mixed == 0

    logging.debug("Singletons: {0} (AN LOSS: {1}, CN LOSS: {2})".format(ns, als, cls))
    logging.debug("Segments: {0} (AN LOSS: {1}, CN LOSS: {2})".format(nm, alm, clm))
    print(SummaryStats([len(x) for x in losses]), file=sys.stderr)

    for x in singletons + segments:
        print("### LENGTH =", len(x), file=fw)
        for i in x:
            print("\t".join(data[i]), file=fw)
    fw.close()


def region_str(region):
    return "{0}:{1}-{2}".format(*region)


def loss(args):
    """
    %prog loss a.b.i1.blocks [a.b-genomic.blast]

    Extract likely gene loss candidates between genome a and b.
    """
    p = OptionParser(loss.__doc__)
    p.add_option(
        "--bed",
        default=False,
        action="store_true",
        help="Genomic BLAST is in bed format",
    )
    p.add_option("--gdist", default=20, type="int", help="Gene distance")
    p.add_option(
        "--bdist",
        default=20000,
        type="int",
        help="Base pair distance",
    )
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    blocksfile = args[0]
    emptyblast = len(args) == 1
    if emptyblast:
        genomicblast = "empty.blast"
        sh("touch {0}".format(genomicblast))
    else:
        genomicblast = args[1]

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
        if key != ".":
            continue

        before, br = data[i - 1]
        after, ar = data[i + 1]
        bi, bx = sorder[before]
        ai, ax = sorder[after]
        dist = abs(bi - ai)
        if bx.seqid != ax.seqid or dist > gdist:
            continue

        start, end = range_minmax(((bx.start, bx.end), (ax.start, ax.end)))
        start, end = max(start - bdist, 1), end + bdist
        proxy = (bx.seqid, start, end)
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

        print("\t".join((b.seqid, accn, target_region)))

    if emptyblast:
        sh("rm -f {0}".format(genomicblast))


def validate(args):
    """
    %prog validate diploid.napus.fractionation cds.bed

    Check whether [S] intervals overlap with CDS.
    """
    from jcvi.formats.bed import intersectBed_wao

    p = OptionParser(validate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fractionation, cdsbed = args
    fp = open(fractionation)

    sbed = "S.bed"
    fw = open(sbed, "w")
    for row in fp:
        a, b, c = row.split()
        if not c.startswith("[S]"):
            continue

        tag, (seqid, start, end) = get_tag(c, None)
        print("\t".join(str(x) for x in (seqid, start - 1, end, b)), file=fw)

    fw.close()

    pairs = {}
    for a, b in intersectBed_wao(sbed, cdsbed):
        if b is None:
            continue
        pairs[a.accn] = b.accn

    validated = fractionation + ".validated"
    fw = open(validated, "w")
    fp.seek(0)
    fixed = 0
    for row in fp:
        a, b, c = row.split()
        if b in pairs:
            assert c.startswith("[S]")
            c = pairs[b]
            fixed += 1

        print("\t".join((a, b, c)), file=fw)

    logging.debug("Fixed {0} [S] cases in `{1}`.".format(fixed, validated))
    fw.close()


if __name__ == "__main__":
    main()
