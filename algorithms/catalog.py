#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging
import string

from optparse import OptionParser

from jcvi.formats.blast import BlastLine
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed
from jcvi.formats.base import must_open
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.algorithms.synteny import AnchorFile, add_beds, check_beds
from jcvi.apps.base import debug, set_outfile, set_stripnames, \
        ActionDispatcher, need_update, sh
debug()


def main():

    actions = (
        ('tandem', 'identify tandem gene groups within certain distance'),
        ('ortholog', 'run a combined synteny and RBH pipeline to call orthologs'),
        ('group', 'cluster the anchors into ortho-groups'),
        ('omgprepare', 'prepare to run Sankoff OMG algorithm'),
        ('omg', 'run Sankoff OMG algorithm to get orthologs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def group(args):
    """
    %prog group anchorfiles

    Group the anchors into ortho-groups. Can input multiple anchor files.
    """
    p = OptionParser(group.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    anchorfiles = args
    groups = Grouper()
    for anchorfile in anchorfiles:
        ac = AnchorFile(anchorfile)
        for a, b, idx in ac.iter_pairs():
            groups.join(a, b)

    ngroups = len(groups)
    nmembers = sum(len(x) for x in groups)
    logging.debug("Created {0} groups with {1} members.".\
                  format(ngroups, nmembers))

    outfile = opts.outfile
    fw = must_open(outfile, "w")
    for g in groups:
        print >> fw, ",".join(sorted(g))
    fw.close()

    return outfile


def omg(args):
    """
    %prog omg weightsfile

    Run Sankoff's OMG algorithm to get orthologs. Download OMG code at:
    <http://137.122.149.195/IsbraSoftware/OMGMec.html>
    """
    from glob import glob
    from collections import defaultdict

    p = OptionParser(omg.__doc__)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    weightsfiles = args
    groupfile = group(weightsfiles + ["--outfile=groups"])

    weights = defaultdict(list)
    for wf in weightsfiles:
        fp = open(wf)
        for row in fp:
            a, b, c = row.split()
            weights[a].append((a, b, c))

    infofiles = glob("*.info")
    info = {}
    for infof in infofiles:
        fp = open(infof)
        for row in fp:
            a = row.split()[0]
            info[a] = row.rstrip()

    fp = open(groupfile)
    for row in fp:
        genes = row.rstrip().split(",")
        ngenes = len(genes)
        if ngenes < 10:
            continue

        header = "a group of genes  :length ={0}".format(ngenes)
        print header
        for gene in genes:
            gene_pairs = weights[gene]
            for a, b, c in gene_pairs:
                if b not in genes:
                    continue

                print "weight", c
                print info[a]
                print info[b]
                print


def geneinfo(bed, order, genomeidx, ploidy):
    bedfile = bed.filename
    p = bedfile.split(".")[0]
    idx = genomeidx[p]
    pd = ploidy[p]
    infofile = p + ".info"

    if not need_update(bedfile, infofile):
        return infofile

    fwinfo = open(infofile, "w")

    for s in bed:
        chr = "".join(x for x in s.seqid if x in string.digits)
        chr = int(chr)
        print >> fwinfo, "\t".join(str(x) for x in \
                    (s.accn, chr, s.start, s.end, s.strand, idx, pd))
    fwinfo.close()

    logging.debug("Update info file `{0}`.".format(infofile))

    return infofile


def omgprepare(args):
    """
    %prog omgprepare anchorsfile blastfile ploidy

    Prepare to run Sankoff's OMG algorithm to get orthologs.
    """
    from jcvi.formats.blast import cscore
    from jcvi.formats.base import DictFile

    p = OptionParser(omgprepare.__doc__)
    set_stripnames(p)
    add_beds(p)

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    anchorfile, blastfile, ploidy = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    fp = open(ploidy)
    genomeidx = dict((x.split()[0], i) for i, x in enumerate(fp))
    fp.close()

    ploidy = DictFile(ploidy)
    qp = qbed.filename.split(".")[0]
    sp = sbed.filename.split(".")[0]

    geneinfo(qbed, qorder, genomeidx, ploidy)
    geneinfo(sbed, sorder, genomeidx, ploidy)

    pf = blastfile.rsplit(".", 1)[0]
    cscorefile = pf + ".cscore"
    #cscore([blastfile, "-o", cscorefile, "--cutoff=.5"])
    ac = AnchorFile(anchorfile)
    pairs = set((a, b) for a, b, i in ac.iter_pairs())
    logging.debug("Imported {0} pairs from `{1}`.".format(len(pairs), anchorfile))

    weightsfile = pf + ".weights"
    fp = open(cscorefile)
    fw = open(weightsfile, "w")
    npairs = 0
    for row in fp:
        a, b, c = row.split()
        c = float(c)
        if (a, b) not in pairs and c < .9:
            continue

        c = int(c * 100)
        qi, q = qorder[a]
        si, s = sorder[b]

        print >> fw, "\t".join((a, b, str(c)))
        npairs += 1
    fw.close()

    logging.debug("Write {0} pairs to `{1}`.".format(npairs, weightsfile))


def make_ortholog(blocksfile, rbhfile, orthofile):
    from jcvi.formats.base import DictFile

    # Generate mapping both ways
    adict = DictFile(rbhfile)
    bdict = DictFile(rbhfile, keypos=1, valuepos=0)
    adict.update(bdict)

    fp = open(blocksfile)
    fw = open(orthofile, "w")
    nrecruited = 0
    for row in fp:
        a, b = row.split()
        if b == '.':
            if a in adict:
                b = adict[a]
                nrecruited += 1
                b += "'"
        print >> fw, "\t".join((a, b))

    logging.debug("Recruited {0} pairs from RBH.".format(nrecruited))
    fp.close()
    fw.close()


def ortholog(args):
    """
    %prog ortholog a.cds a.bed b.cds b.bed

    Run a sensitive pipeline to find orthologs between two species a and b.
    The pipeline runs LAST and 1-to-1 quota synteny blocks as the backbone of
    such predictions. Extra orthologs will be recruited from reciprocal best
    match (RBH).
    """
    from jcvi.apps.last import main as last_main
    from jcvi.apps.blastfilter import main as blastfilter_main
    from jcvi.algorithms.quota import main as quota_main
    from jcvi.algorithms.synteny import scan, screen, mcscan
    from jcvi.formats.blast import cscore

    p = OptionParser(ortholog.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    afasta, abed, bfasta, bbed = args
    aprefix = afasta.split(".")[0]
    bprefix = bfasta.split(".")[0]
    pprefix = ".".join((aprefix, bprefix))
    qprefix = ".".join((bprefix, aprefix))
    last = pprefix + ".last"
    if need_update((afasta, bfasta), last):
        last_main([bfasta, afasta, "-o", last, "-a", "16"])

    filtered_last = last + ".filtered"
    bstring = ["--qbed=" + abed, "--sbed=" + bbed]
    if need_update(last, filtered_last):
        blastfilter_main([last, "--cscore=.7",
                          "--tandem_Nmax=10"] + bstring)

    anchors = pprefix + ".anchors"
    lifted_anchors = pprefix + ".lifted.anchors"
    if need_update(filtered_last, lifted_anchors):
        scan([filtered_last, anchors, "--dist=20",
              "--liftover=" + last] + bstring)

    blockids = pprefix + ".1x1.ids"
    if need_update(lifted_anchors, blockids):
        quota_main([lifted_anchors, "--quota=1:1",
                   "-o", blockids] + bstring)

    ooanchors = lifted_anchors + ".1x1"
    if need_update(blockids, ooanchors):
        screen([lifted_anchors, ooanchors, "--ids=" + blockids,
                "--minspan=0"] + bstring)

    pblocks = pprefix + ".1x1.blocks"
    qblocks = qprefix + ".1x1.blocks"
    if need_update(ooanchors, [pblocks, qblocks]):
        mcscan([abed, ooanchors, "--iter=1", "-o", pblocks])
        mcscan([bbed, ooanchors, "--iter=1", "-o", qblocks])

    rbh = pprefix + ".rbh"
    if need_update(last, rbh):
        cscore([last, "-o", rbh])

    portho = pprefix + ".ortholog"
    qortho = qprefix + ".ortholog"
    if need_update([pblocks, qblocks, rbh], [portho, qortho]):
        make_ortholog(pblocks, rbh, portho)
        make_ortholog(qblocks, rbh, qortho)


def tandem_main(blast_file, cds_file, bed_file, N=3, P=50, is_self=True, \
    strip_name=".", ofile=sys.stderr):

    # get the sizes for the CDS first
    f = Fasta(cds_file)
    sizes = dict(f.itersizes())

    # retrieve the locations
    bed = Bed(bed_file)
    order = bed.order

    if is_self:
        # filter the blast file
        g = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            qi, q = order[query]
            si, s = order[subject]

            if q.seqid == s.seqid and abs(qi - si) <= N:
                g.join(query, subject)

    else:
        homologs = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len)*P/100.:
                continue

            query = gene_name(b.query, strip_name)
            subject = gene_name(b.subject, strip_name)
            homologs.join(query, subject)

        g = Grouper()
        for i, atom in enumerate(bed):
            for x in range(1, N+1):
                if all([i-x >= 0, bed[i-x].seqid == atom.seqid, \
                    homologs.joined(bed[i-x].accn, atom.accn)]):
                    g.join(bed[i-x].accn, atom.accn)

    # dump the grouper
    fw = must_open(ofile, "w")
    ngenes, nfamilies = 0, 0
    families = []
    for group in sorted(g):
        if len(group) >= 2:
            print >>fw, ",".join(sorted(group))
            ngenes += len(group)
            nfamilies += 1
            families.append(sorted(group))

    longest_family = max(families, key=lambda x: len(x))

    # generate reports
    print >>sys.stderr, "Proximal paralogues (dist=%d):" % N
    print >>sys.stderr, "Total %d genes in %d families" % (ngenes, nfamilies)
    print >>sys.stderr, "Longest families (%d): %s" % (len(longest_family),
        ",".join(longest_family))

    return families


def tandem(args):
    """
    %prog tandem blast_file cds_file bed_file [options]

    Find tandem gene clusters that are separated by N genes, based on filtered
    blast_file by enforcing alignments between any two genes at least 50%
    (or user specified value) of either gene.

    pep_file can also be used in same manner.
    """
    p = OptionParser(tandem.__doc__)
    p.add_option("--tandem_Nmax", dest="tandem_Nmax", type="int", default=3,
               help="merge tandem genes within distance [default: %default]")
    p.add_option("--percent_overlap", type="int", default=50,
               help="tandem genes have >=x% aligned sequence, x=0-100 \
               [default: %default]")
    p.add_option("--not_self", default=False, action="store_true",
                 help="provided is not self blast file [default: %default]")
    p.add_option("--strip_gene_name", dest="sep", type="string", default=".",
               help="strip alternative splicing. Use None for no stripping. \
               [default: %default]")
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    blast_file, cds_file, bed_file = args
    N = opts.tandem_Nmax
    P = opts.percent_overlap
    is_self = not opts.not_self
    sep = opts.sep
    ofile = opts.outfile

    tandem_main(blast_file, cds_file, bed_file, N=N, P=P, is_self=is_self, \
        strip_name=sep, ofile=ofile)


if __name__ == '__main__':
    main()
