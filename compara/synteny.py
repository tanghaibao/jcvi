#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import os.path as op
import sys
import logging

import numpy as np
from collections import defaultdict

from jcvi.algorithms.lis import heaviest_increasing_subsequence as his
from jcvi.formats.bed import Bed
from jcvi.formats.blast import BlastLine
from jcvi.formats.base import BaseFile, SetFile, read_block, must_open
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name, human_size
from jcvi.utils.range import Range, range_chain
from jcvi.apps.base import OptionParser, ActionDispatcher


class AnchorFile (BaseFile):

    def __init__(self, filename, minsize=0):
        super(AnchorFile, self).__init__(filename)
        self.blocks = list(self.iter_blocks(minsize=minsize))

    def iter_blocks(self, minsize=0):
        fp = open(self.filename)
        for header, lines in read_block(fp, "#"):
            lines = [x.split() for x in lines]
            if len(lines) >= minsize:
                yield lines

    def iter_pairs(self, minsize=0):
        block_id = -1
        for rows in self.iter_blocks(minsize=minsize):
            block_id += 1
            for row in rows:
                a, b = row[:2]
                yield a, b, block_id

    def print_to_file(self, filename="stdout", accepted=None):
        fw = must_open(filename, "w")
        blocks = self.blocks
        nremoved = 0
        ncorrected = 0
        for block in blocks:
            print >> fw, "###"
            for line in block:
                a, b, score = line
                pair = (a, b)
                if accepted:
                    if pair not in accepted:
                        nremoved += 1
                        continue
                    av = accepted[pair]
                    if score != av and score != av + 'L':
                        score = av
                        ncorrected += 1
                print >> fw, "\t".join((a, b, score))
        fw.close()

        logging.debug("Removed {0} existing anchors.".format(nremoved))
        logging.debug("Corrected scores for {0} anchors.".format(ncorrected))
        logging.debug("Anchors written to `{0}`.".format(filename))

    def blast(self, blastfile=None, outfile=None):
        """
        convert anchor file to 12 col blast file
        """
        from jcvi.formats.blast import BlastSlow, BlastLineByConversion

        if not outfile:
            outfile = self.filename + ".blast"

        if blastfile is not None:
            blasts = BlastSlow(blastfile).to_dict()
        else:
            blasts = None

        fw = must_open(outfile, "w", checkexists=True)
        nlines = 0
        for a, b, id in self.iter_pairs():
            if (a, b) in blasts:
                bline = blasts[(a, b)]
            elif (b, a) in blasts:
                bline = blasts[(b, a)]
            else:
                line = "\t".join((a, b))
                bline = BlastLineByConversion(line, mode="110000000000")

            print >> fw, bline
            nlines += 1
        fw.close()

        logging.debug("A total of {0} BLAST lines written to `{1}`."\
                        .format(nlines, outfile))

        return outfile


class BlockFile (BaseFile):

    def __init__(self, filename, header=False):
        super(BlockFile, self).__init__(filename)
        fp = must_open(filename)
        hd = fp.next().rstrip().split("\t")
        ncols = len(hd)
        if header:
            self.header = hd
        else:
            fp.seek(0)
            self.header = range(ncols)

        data = []
        highlight = []
        for row in fp:
            hl = ("*" in row)
            # r* highlights the block in red color
            if hl:
                hl, row = row.split("*", 1)
                hl = hl or "r"
            atoms = row.rstrip().split("\t")
            atoms = [x.strip() for x in atoms]
            atoms = ["." if x == "" else x for x in atoms]
            if len(atoms) > ncols:
                atoms = atoms[:ncols]
            elif len(atoms) < ncols:
                atoms = atoms + ["."] * (ncols - len(atoms))
            data.append(atoms)
            highlight.append(hl)

        self.data = data
        self.highlight = highlight
        self.columns = zip(*data)
        self.ncols = ncols

    def get_extent(self, i, order, debug=True):
        col = self.columns[i]
        ocol = [order[x] for x in col if x in order]
        orientation = '+' if ocol[0][0] <= ocol[-1][0] else '-'
        si, start = min(ocol)
        ei, end = max(ocol)
        same_chr = (start.seqid == end.seqid)
        chr = start.seqid if same_chr else None
        ngenes = ei - si + 1
        if debug:
            r = "{0}:{1}-{2}".format(chr, start.start, end.end)
            print >> sys.stderr, "Column {0}: {1} - {2} ({3})".\
                    format(i, start.accn, end.accn, r)
            print >> sys.stderr, "  {0} .. {1} ({2}) features .. {3}".\
                    format(chr, ngenes, len(ocol), orientation)

        span = abs(start.start - end.end)

        return start, end, si, ei, chr, orientation, span

    def iter_pairs(self, i, j, highlight=False):
        for h, d in zip(self.highlight, self.data):
            if highlight and not h:
                continue

            a, b = d[i], d[j]
            if "." in (a, b) or "" in (a, b):
                continue

            yield a, b, h

    def iter_all_pairs(self):
        ncols = self.ncols
        for i in xrange(ncols):
            for j in xrange(i + 1, ncols):
                for a, b, h in self.iter_pairs(i, j):
                    yield a, b, h

    def iter_gene_col(self):
        for hd, col in zip(self.header, self.columns):
            for g in col:
                if g not in (".", ""):
                    yield g, hd

    def query_gene(self, gene, color=None, invert=False):
        """
        Used in mcscanq() for query
        """
        qi = self.columns[0].index(gene)
        ndata = len(self.data)
        for col in self.columns[1:]:
            upstream_dist = downstream_dist = 1000
            # search upstream
            for i in xrange(qi - 1, -1, -1):
                if col[i] not in (".", ""):
                    upstream = col[i]
                    upstream_dist = qi - i
                    break
            # search downstream
            for i in xrange(qi, ndata):
                if col[i] not in (".", ""):
                    downstream = col[i]
                    downstream_dist = i - qi
                    break
            closest = upstream if upstream_dist < downstream_dist \
                    else downstream
            # output in .simple format
            if invert:
                line = "\t".join(str(x) for x in \
                        (closest, closest, gene, gene, 0, "+"))
            else:
                line = "\t".join(str(x) for x in \
                        (gene, gene, closest, closest, 0, "+"))
            if color is not None:
                line = color + "*" + line
            yield line


class SimpleFile (object):

    def __init__(self, simplefile, defaultcolor='#fb8072', order=None):
        # Sometimes the simplefile has query and subject wrong
        fp = open(simplefile)
        self.blocks = []
        for row in fp:
            if row[:2] == "##" or row.startswith("StartGeneA"):
                continue
            hl = ("*" in row)
            if hl:
                hl, row = row.split("*", 1)
                hl = hl or defaultcolor
            a, b, c, d, score, orientation = row.split()
            if order and a not in order:
                a, b, c, d = c, d, a, b
            if orientation == '-':
                c, d = d, c
            score = int(score)
            self.blocks.append((a, b, c, d, score, orientation, hl))


def _score(cluster):
    """
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y = zip(*cluster)[:2]
    return min(len(set(x)), len(set(y)))


def group_hits(blasts):
    # grouping the hits based on chromosome pair
    all_hits = defaultdict(list)
    for b in blasts:
        all_hits[(b.qseqid, b.sseqid)].append((b.qi, b.si, b.score))

    return all_hits


def read_blast(blast_file, qorder, sorder, is_self=False, ostrip=True):
    """
    read the blast and convert name into coordinates
    """
    fp = open(blast_file)
    filtered_blast = []
    seen = set()
    for row in fp:
        b = BlastLine(row)
        query, subject = b.query, b.subject
        if query == subject:
            continue
        if ostrip:
            query, subject = gene_name(query), gene_name(subject)
        if query not in qorder or subject not in sorder:
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]

        if is_self and qi > si:
            # remove redundant a<->b to one side when doing self-self BLAST
            query, subject = subject, query
            qi, si = si, qi
            q, s = s, q

        key = query, subject
        if key in seen:
            continue
        seen.add(key)

        b.qseqid, b.sseqid = q.seqid, s.seqid
        b.qi, b.si = qi, si
        b.query, b.subject = query, subject

        filtered_blast.append(b)

    logging.debug("A total of {0} BLAST imported from `{1}`.".\
                  format(len(filtered_blast), blast_file))

    return filtered_blast


def read_anchors(ac, qorder, sorder, minsize=0):
    """
    anchors file are just (geneA, geneB) pairs (with possible deflines)
    """
    all_anchors = defaultdict(list)
    nanchors = 0
    anchor_to_block = {}

    for a, b, idx in ac.iter_pairs(minsize=minsize):
        if a not in qorder or b not in sorder:
            continue
        qi, q = qorder[a]
        si, s = sorder[b]
        pair = (qi, si)

        all_anchors[(q.seqid, s.seqid)].append(pair)
        anchor_to_block[pair] = idx
        nanchors += 1

    logging.debug("A total of {0} anchors imported.".format(nanchors))
    assert nanchors == len(anchor_to_block)

    return all_anchors, anchor_to_block


def synteny_scan(points, xdist, ydist, N):
    """
    This is the core single linkage algorithm which behaves in O(n):
    iterate through the pairs, foreach pair we look back on the
    adjacent pairs to find links
    """
    clusters = Grouper()
    n = len(points)
    points.sort()
    for i in xrange(n):
        for j in xrange(i - 1, -1, -1):
            # x-axis distance
            del_x = points[i][0] - points[j][0]
            if del_x > xdist:
                break
            # y-axis distance
            del_y = points[i][1] - points[j][1]
            if abs(del_y) > ydist:
                continue
            # otherwise join
            clusters.join(points[i], points[j])

    # select clusters that are at least >=N
    clusters = [sorted(cluster) for cluster in list(clusters) \
            if _score(cluster) >= N]

    return clusters


def batch_scan(points, xdist=20, ydist=20, N=5):
    """
    runs synteny_scan() per chromosome pair
    """
    chr_pair_points = group_hits(points)

    clusters = []
    for chr_pair in sorted(chr_pair_points.keys()):
        points = chr_pair_points[chr_pair]
        clusters.extend(synteny_scan(points, xdist, ydist, N))

    return clusters


def synteny_liftover(points, anchors, dist):
    """
    This is to get the nearest anchors for all the points (useful for the
    `liftover` operation below).
    """
    from scipy.spatial import cKDTree

    points = np.array(points, dtype=int)
    ppoints = points[:, :2] if points.shape[1] > 2 else points
    tree = cKDTree(anchors, leafsize=16)
    #print tree.data
    dists, idxs = tree.query(ppoints, p=1, distance_upper_bound=dist)

    for point, dist, idx in zip(points, dists, idxs):
        if idx == tree.n:  # nearest is out of range
            continue
        if dist == 0:  # already in anchors
            continue

        yield point, tuple(anchors[idx])


def check_beds(hintfile, p, opts):

    wd, hintfile = op.split(hintfile)
    if not (opts.qbed and opts.sbed):
        try:
            q, s = hintfile.split(".", 2)[:2]
            opts.qbed = op.join(wd, q + ".bed")
            opts.sbed = op.join(wd, s + ".bed")
            logging.debug("Assuming --qbed={0} --sbed={1}".\
                         format(opts.qbed, opts.sbed))
        except:
            print >> sys.stderr, "Options --qbed and --sbed are required"
            sys.exit(not p.print_help())

    qbed_file, sbed_file = opts.qbed, opts.sbed
    # is this a self-self blast?
    is_self = (qbed_file == sbed_file)
    if is_self:
        logging.debug("Looks like self-self comparison.")

    qbed = Bed(opts.qbed)
    sbed = Bed(opts.sbed)
    qorder = qbed.order
    sorder = sbed.order

    return qbed, sbed, qorder, sorder, is_self


def add_options(p, args, dist=10):
    """
    scan and liftover has similar interfaces, so share common options
    returns opts, files
    """
    p.set_beds()
    p.add_option("--dist", default=dist, type="int",
            help="Extent of flanking regions to search [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blast_file, anchor_file = args

    return blast_file, anchor_file, opts.dist, opts


def main():

    actions = (
        ('scan', 'get anchor list using single-linkage algorithm'),
        ('summary', 'provide statistics for pairwise blocks'),
        ('liftover', 'given anchor list, pull adjacent pairs from blast file'),
        ('mcscan', 'stack synteny blocks on a reference bed'),
        ('mcscanq', 'query multiple synteny blocks'),
        ('screen', 'extract subset of blocks from anchorfile'),
        ('simple', 'convert anchorfile to simple block descriptions'),
        ('stats', 'provide statistics for mscan blocks'),
        ('depth', 'calculate the depths in the two genomes in comparison'),
        ('breakpoint', 'identify breakpoints where collinearity ends'),
        ('matrix', 'make oxford grid based on anchors file'),
        ('coge', 'convert CoGe file to anchors file'),
        ('spa', 'convert chr ordering from SPA to simple lists'),
        ('layout', 'compute layout based on .simple file'),
        ('rebuild', 'rebuild anchors file from prebuilt blocks file'),
        ('fromaligns', 'convert aligns file to anchors file'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def colinear_evaluate_weights(tour, data):
    tour = dict((s, i) for i, s in enumerate(tour))
    data = [(tour[x], score) for x, score in data if x in tour]
    return his(data)[-1],


def layout(args):
    """
    %prog layout query.subject.simple query.seqids subject.seqids

    Compute optimal seqids order in a second genome, based on seqids on one
    genome, given the pairwise blocks in .simple format.
    """
    from jcvi.algorithms.ec import GA_setup, GA_run

    p = OptionParser(layout.__doc__)
    p.set_beds()
    p.set_cpus(cpus=32)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    simplefile, qseqids, sseqids = args
    qbed, sbed, qorder, sorder, is_self = check_beds(simplefile, p, opts)

    qseqids = qseqids.strip().split(",")
    sseqids = sseqids.strip().split(",")
    qseqids_ii = dict((s, i) for i, s in enumerate(qseqids))
    sseqids_ii = dict((s, i) for i, s in enumerate(sseqids))

    blocks = SimpleFile(simplefile).blocks
    scores = defaultdict(int)
    for a, b, c, d, score, orientation, hl in blocks:
        qi, q = qorder[a]
        si, s = sorder[c]
        qseqid, sseqid = q.seqid, s.seqid
        if sseqid not in sseqids:
            continue
        scores[sseqids_ii[sseqid], qseqid] += score

    data = []
    for (a, b), score in sorted(scores.items()):
        if b not in qseqids_ii:
            continue
        data.append((qseqids_ii[b], score))

    tour = range(len(qseqids))
    toolbox = GA_setup(tour)
    toolbox.register("evaluate", colinear_evaluate_weights, data=data)
    tour, fitness = GA_run(toolbox, ngen=100, npop=100, cpus=opts.cpus)
    tour = [qseqids[x] for x in tour]

    print ",".join(tour)


def fromaligns(args):
    """
    %prog fromaligns out.aligns

    Convert aligns file (old MCscan output) to anchors file.
    """
    p = OptionParser(fromaligns.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    alignsfile, = args
    fp = must_open(alignsfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        if row.startswith("## Alignment"):
            print >> fw, "###"
            continue
        if row[0] == '#' or not row.strip():
            continue
        atoms = row.split(':')[-1].split()
        print >> fw, "\t".join(atoms[:2])
    fw.close()


def mcscanq(args):
    """
    %prog mcscanq query.ids blocksfile

    Query multiple synteny blocks to get the closest alignment feature. Mostly
    used for 'highlighting' the lines in the synteny plot, drawn by
    graphics.karyotype and graphics.synteny.
    """
    p = OptionParser(mcscanq.__doc__)
    p.add_option("--color", help="Add color highlight, used in plotting")
    p.add_option("--invert", default=False, action="store_true",
                 help="Invert query and subject [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    qids, blocksfile = args
    b = BlockFile(blocksfile)
    fp = open(qids)
    for gene in fp:
        gene = gene.strip()
        for line in b.query_gene(gene, color=opts.color, invert=opts.invert):
            print line


def spa(args):
    """
    %prog spa spafiles

    Convert chromosome ordering from SPA to simple lists. First column is the
    reference order.
    """
    from jcvi.algorithms.graph import merge_paths
    from jcvi.utils.cbook import uniqify

    p = OptionParser(spa.__doc__)
    p.add_option("--unmapped", default=False, action="store_true",
                 help="Include unmapped scaffolds in the list [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    spafiles = args
    paths = []
    mappings = []
    missings = []
    for spafile in spafiles:
        fp = open(spafile)
        path = []
        mapping = []
        missing = []
        for row in fp:
            if row[0] == '#' or not row.strip():
                continue

            atoms = row.rstrip().split('\t')
            if len(atoms) == 2:
                a, c2 = atoms
                assert a == "unmapped"
                missing.append(c2)
                continue

            c1, c2, orientation = atoms
            path.append(c1)
            mapping.append(c2)

        paths.append(uniqify(path))
        mappings.append(mapping)
        missings.append(missing)

    ref = merge_paths(paths)
    print "ref", len(ref), ",".join(ref)
    for spafile, mapping, missing in zip(spafiles, mappings, missings):
        mapping = [x for x in mapping if "random" not in x]
        mapping = uniqify(mapping)
        if len(mapping) < 50 and opts.unmapped:
            mapping = uniqify(mapping + missing)

        print spafile, len(mapping), ",".join(mapping)


def rebuild(args):
    """
    %prog rebuild blocksfile blastfile

    Rebuild anchors file from pre-built blocks file.
    """
    p = OptionParser(rebuild.__doc__)
    p.add_option("--header", default=False, action="store_true",
                 help="First line is header [default: %default]")
    p.add_option("--write_blast", default=False, action="store_true",
                 help="Get blast records of rebuilt anchors [default: %default]")
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blocksfile, blastfile = args
    bk = BlockFile(blocksfile, header=opts.header)
    fw = open("pairs", "w")
    for a, b, h in bk.iter_all_pairs():
        print >> fw, "\t".join((a, b))
    fw.close()

    if opts.write_blast:
        AnchorFile("pairs").blast(blastfile, "pairs.blast")

    fw = open("tracks", "w")
    for g, col in bk.iter_gene_col():
        print >> fw, "\t".join(str(x) for x in (g, col))
    fw.close()


def coge(args):
    """
    %prog coge cogefile

    Convert CoGe file to anchors file.
    """
    p = OptionParser(coge.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    cogefile, = args
    fp = must_open(cogefile)
    cogefile = cogefile.replace(".gz", "")
    ksfile = cogefile + ".ks"
    anchorsfile = cogefile + ".anchors"
    fw_ks = must_open(ksfile, "w")
    fw_ac = must_open(anchorsfile, "w")

    tag = "###"
    print >> fw_ks, tag
    for header, lines in read_block(fp, tag):
        print >> fw_ac, tag
        lines = list(lines)
        for line in lines:
            if line[0] == '#':
                continue
            ks, ka, achr, a, astart, astop, bchr, \
                    b, bstart, bstop, ev, ss = line.split()
            a = a.split("||")[3]
            b = b.split("||")[3]
            print >> fw_ac, "\t".join((a, b, ev))
            print >> fw_ks, ",".join((";".join((a, b)), ks, ka, ks, ka))

    fw_ks.close()
    fw_ac.close()


def matrix(args):
    """
    %prog matrix all.bed anchorfile matrixfile

    Make oxford grid based on anchors file.
    """

    p = OptionParser(matrix.__doc__)
    p.add_option("--seqids", help="File with seqids [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, anchorfile, matrixfile = args
    ac = AnchorFile(anchorfile)
    seqidsfile = opts.seqids
    if seqidsfile:
        seqids = SetFile(seqidsfile, delimiter=',')

    order = Bed(bedfile).order
    blocks = ac.blocks
    m = defaultdict(int)
    fw = open(matrixfile, "w")
    aseqids = set()
    bseqids = set()
    for block in blocks:
        a, b, scores = zip(*block)
        ai, af = order[a[0]]
        bi, bf = order[b[0]]
        aseqid = af.seqid
        bseqid = bf.seqid
        if seqidsfile:
            if (aseqid not in seqids) or (bseqid not in seqids):
                continue
        m[(aseqid, bseqid)] += len(block)
        aseqids.add(aseqid)
        bseqids.add(bseqid)

    aseqids = list(aseqids)
    bseqids = list(bseqids)
    print >> fw, "\t".join(["o"] + bseqids)
    for aseqid in aseqids:
        print >> fw, "\t".join([aseqid] + \
                    [str(m[(aseqid, x)]) for x in bseqids])


def get_boundary_bases(start, end, order):

    from jcvi.utils.range import range_minmax

    (i, s), (j, e) = order[start], order[end]
    seqid = s.seqid
    assert seqid == e.seqid

    startbase, endbase = range_minmax([(s.start, s.end), (e.start, e.end)])

    return seqid, startbase, endbase


def simple(args):
    """
    %prog simple anchorfile --qbed=qbedfile --sbed=sbedfile [options]

    Write the block ends for each block in the anchorfile.
    GeneA1    GeneA2    GeneB1    GeneB2   +/-      score

    Optional additional columns:
    orderA1   orderA2   orderB1   orderB2  sizeA    sizeB   size    block_id

    With base coordinates (--coords):
    block_id  seqidA    startA    endA     bpSpanA  GeneA1   GeneA2  geneSpanA
    block_id  seqidB    startB    endB     bpSpanB  GeneB1   GeneB2  geneSpanB
    """
    p = OptionParser(simple.__doc__)
    p.add_option("--rich", default=False, action="store_true", \
                help="Output additional columns [default: %default]")
    p.add_option("--coords", default=False, action="store_true",
                help="Output columns with base coordinates [default: %default]")
    p.add_option("--noheader", default=False, action="store_true",
                help="Don't output header [default: %default]")
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    additional = opts.rich
    coords = opts.coords
    header = not opts.noheader

    ac = AnchorFile(anchorfile)
    simplefile = anchorfile.rsplit(".", 1)[0] + ".simple"

    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    pf = "-".join(anchorfile.split(".", 2)[:2])
    blocks = ac.blocks

    if coords:
        h = "Block|Chr|Start|End|Span|StartGene|EndGene|GeneSpan|Orientation"
    else:
        h = "StartGeneA|EndGeneA|StartGeneB|EndGeneB|Orientation|Score"
        if additional:
            h += "|StartOrderA|EndOrderA|StartOrderB|EndOrderB|"\
                  "SizeA|SizeB|Size|Block"

    fws = open(simplefile, "w")
    if header:
        print >> fws, "\t".join(h.split("|"))

    atotalbase = btotalbase = 0
    for i, block in enumerate(blocks):

        a, b, scores = zip(*block)
        a = [qorder[x] for x in a]
        b = [sorder[x] for x in b]
        ia, oa = zip(*a)
        ib, ob = zip(*b)

        astarti, aendi = min(ia), max(ia)
        bstarti, bendi = min(ib), max(ib)
        astart, aend = min(a)[1].accn, max(a)[1].accn
        bstart, bend = min(b)[1].accn, max(b)[1].accn

        sizeA = len(set(ia))
        sizeB = len(set(ib))
        size = len(block)

        slope, intercept = np.polyfit(ia, ib, 1)
        orientation = "+" if slope >= 0 else '-'
        aspan = aendi - astarti + 1
        bspan = bendi - bstarti + 1
        score = int((aspan * bspan) ** .5)
        score = str(score)
        block_id = pf + "-block-{0}".format(i)

        if coords:

            aseqid, astartbase, aendbase = \
                    get_boundary_bases(astart, aend, qorder)
            bseqid, bstartbase, bendbase = \
                    get_boundary_bases(bstart, bend, sorder)
            abase = aendbase - astartbase + 1
            bbase = bendbase - bstartbase + 1
            atotalbase += abase
            btotalbase += bbase

            # Write dual lines
            aargs = [block_id, aseqid, astartbase, aendbase,
                     abase, astart, aend, aspan, "+"]
            bargs = [block_id, bseqid, bstartbase, bendbase,
                     bbase, bstart, bend, bspan, orientation]

            for args in (aargs, bargs):
                print >> fws, "\t".join(str(x) for x in args)
            continue

        args = [astart, aend, bstart, bend, score, orientation]
        if additional:
            args += [astarti, aendi, bstarti, bendi,
                     sizeA, sizeB, size, block_id]
        print >> fws, "\t".join(str(x) for x in args)

    fws.close()
    logging.debug("A total of {0} blocks written to `{1}`.".format(i + 1, simplefile))

    if coords:
        print >> sys.stderr, "Total block span in {0}: {1}".format(qbed.filename, \
                        human_size(atotalbase, precision=0))
        print >> sys.stderr, "Total block span in {0}: {1}".format(sbed.filename, \
                        human_size(btotalbase, precision=0))
        print >> sys.stderr, "Ratio: {0:.1f}x".format(\
                        max(atotalbase, btotalbase) * 1. / min(atotalbase, btotalbase))


def screen(args):
    """
    %prog screen anchorfile newanchorfile --qbed=qbedfile --sbed=sbedfile [options]

    Extract subset of blocks from anchorfile. Provide several options:

    1. Option --ids: a file with IDs, 0-based, comma separated, all in one line.
    2. Option --seqids: only allow seqids in this file.
    3. Option --seqpairs: only allow seqpairs in this file, one per line, e.g. "Chr01,Chr05".
    4. Option --minspan: remove blocks with less span than this.
    5. Option --minsize: remove blocks with less number of anchors than this.
    """
    p = OptionParser(screen.__doc__)
    p.set_beds()

    p.add_option("--ids", help="File with block IDs (0-based) [default: %default]")
    p.add_option("--seqids", help="File with seqids [default: %default]")
    p.add_option("--seqpairs", help="File with seqpairs [default: %default]")
    p.add_option("--nointra", action="store_true",
                 help="Remove intra-chromosomal blocks [default: %default]")
    p.add_option("--minspan", default=0, type="int",
                 help="Only blocks with span >= [default: %default]")
    p.add_option("--minsize", default=0, type="int",
                 help="Only blocks with anchors >= [default: %default]")
    p.add_option("--simple", action="store_true",
                 help="Write simple anchorfile with block ends [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    anchorfile, newanchorfile = args
    ac = AnchorFile(anchorfile)
    idsfile = opts.ids
    seqidsfile = opts.seqids
    seqpairsfile = opts.seqpairs
    minspan = opts.minspan
    minsize = opts.minsize
    osimple = opts.simple
    nointra = opts.nointra
    ids, seqids, seqpairs = None, None, None

    if idsfile:
        ids = SetFile(idsfile, delimiter=',')
        ids = set(int(x) for x in ids)
    if seqidsfile:
        seqids = SetFile(seqidsfile, delimiter=',')
    if seqpairsfile:
        fp = open(seqpairsfile)
        seqpairs = set()
        for row in fp:
            a, b = row.strip().split(",")
            seqpairs.add((a, b))
            seqpairs.add((b, a))

    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    blocks = ac.blocks
    selected = 0
    fw = open(newanchorfile, "w")

    for i, block in enumerate(blocks):
        if ids and i not in ids:
            continue

        a, b, scores = zip(*block)
        a = [qorder[x] for x in a]
        b = [sorder[x] for x in b]
        ia, oa = zip(*a)
        ib, ob = zip(*b)
        aspan = max(ia) - min(ia) + 1
        bspan = max(ib) - min(ib) + 1
        aseqid = oa[0].seqid
        bseqid = ob[0].seqid

        if seqids:
            if (aseqid not in seqids) or (bseqid not in seqids):
                continue

        if seqpairs:
            if (aseqid, bseqid) not in seqpairs:
                continue

        if nointra and aseqid == bseqid:
            continue

        if minsize:
            if len(block) < minsize:
                continue

        if minspan:
            if aspan < minspan or bspan < minspan:
                continue

        selected += 1
        print >> fw, "###"
        for line in block:
            print >> fw, "\t".join(line)

    fw.close()

    if osimple:
        simple([newanchorfile, "--noheader", \
                "--qbed=" + qbed.filename, "--sbed=" + sbed.filename])

    logging.debug("Before: {0} blocks, After: {1} blocks".\
                  format(len(blocks), selected))


def summary(args):
    """
    %prog summary anchorfile

    Provide statistics for pairwise blocks.
    """
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(summary.__doc__)
    p.add_option("--prefix", help="Generate per block stats [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    ac = AnchorFile(anchorfile)
    clusters = ac.blocks

    nclusters = len(clusters)
    nanchors = [len(c) for c in clusters]
    nranchors = [_score(c) for c in clusters]  # non-redundant anchors
    print >> sys.stderr, "A total of {0} (NR:{1}) anchors found in {2} clusters.".\
                  format(sum(nanchors), sum(nranchors), nclusters)
    print >> sys.stderr, "Stats:", SummaryStats(nanchors)
    print >> sys.stderr, "NR stats:", SummaryStats(nranchors)

    prefix = opts.prefix
    if prefix:
        pad = len(str(nclusters))
        for i, c in enumerate(clusters):
            block_id = "{0}{1:0{2}d}".format(prefix, i + 1, pad)
            print "\t".join((block_id, str(len(c))))


def stats(args):
    """
    %prog stats blocksfile

    Provide statistics for MCscan-style blocks. The count of homologs in each
    pivot gene is recorded.
    """
    from jcvi.utils.cbook import percentage

    p = OptionParser(stats.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    blocksfile, = args
    fp = open(blocksfile)
    counts = defaultdict(int)
    total = orthologous = 0
    for row in fp:
        atoms = row.rstrip().split("\t")
        hits = [x for x in atoms[1:] if x != '.']
        counts[len(hits)] += 1
        total += 1
        if atoms[1] != '.':
            orthologous += 1

    print >> sys.stderr, "Total lines: {0}".format(total)
    for i, n in sorted(counts.items()):
        print >> sys.stderr, "Count {0}: {1}".format(i, percentage(n, total))

    print >> sys.stderr

    matches = sum(n for i, n in counts.items() if i != 0)
    print >> sys.stderr, "Total lines with matches: {0}".\
                format(percentage(matches, total))
    for i, n in sorted(counts.items()):
        if i == 0:
            continue

        print >> sys.stderr, "Count {0}: {1}".format(i, percentage(n, matches))

    print >> sys.stderr
    print >> sys.stderr, "Orthologous matches: {0}".\
                format(percentage(orthologous, matches))


def get_best_pair(qs, ss, ts):
    pairs = {}
    for q, s, t in zip(qs, ss, ts):
        t = long(t)
        if q not in pairs or pairs[q][1] < t:
            pairs[q] = (s, t)

    # Discard score
    spairs = dict((q, s) for q, (s, t) in pairs.items())
    return spairs


def get_range(q, s, t, i, order, block_pairs, clip=10):
    pairs = get_best_pair(q, s, t)
    score = len(pairs)
    block_pairs[i].update(pairs)

    q = [order[x][0] for x in q]
    q.sort()
    qmin = q[0]
    qmax = q[-1]
    if qmax - qmin >= 2 * clip:
        qmin += clip / 2
        qmax -= clip / 2

    return Range("0", qmin, qmax, score=score, id=i)


def mcscan(args):
    """
    %prog mcscan bedfile anchorfile [options]

    Stack synteny blocks on a reference bed, MCSCAN style. The first column in
    the output is the reference order, given in the bedfile. Then each column
    next to it are separate 'tracks'.

    If --mergetandem=tandem_file is specified, tandem_file should have each
    tandem cluster as one line, tab separated.
    """
    p = OptionParser(mcscan.__doc__)
    p.add_option("--iter", default=100, type="int",
                 help="Max number of chains to output [default: %default]")
    p.add_option("--ascii", default=False, action="store_true",
                 help="Output symbols rather than gene names [default: %default]")
    p.add_option("--Nm", default=10, type="int",
                 help="Clip block ends to allow slight overlaps [default: %default]")
    p.add_option("--trackids", action="store_true",
                 help="Track block IDs in separate file [default: %default]")
    p.add_option("--mergetandem", default=None,
                 help="merge tandems genes in output acoording to PATH-TO-TANDEM_FILE, "\
                 "cannot be used with --ascii")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, anchorfile = args
    ascii = opts.ascii
    clip = opts.Nm
    trackids = opts.trackids
    ofile = opts.outfile
    mergetandem = opts.mergetandem
    bed = Bed(bedfile)
    order = bed.order

    if trackids:
        olog = ofile + ".tracks"
        fwlog = must_open(olog, "w")

    if mergetandem:
        assert not ascii
        tandems = {}
        for row in file(mergetandem):
            row = row.split()
            s = ";".join(row)
            for atom in row:
                tandems[atom] = s

    ac = AnchorFile(anchorfile)
    ranges = []
    block_pairs = defaultdict(dict)
    blocks = ac.blocks
    for i, ib in enumerate(blocks):
        q, s, t = zip(*ib)
        if q[0] not in order:
            q, s = s, q

        r = get_range(q, s, t, i, order, block_pairs, clip=clip)
        ranges.append(r)

        assert q[0] in order
        if s[0] not in order:
            continue

        # is_self comparison
        q, s = s, q
        r = get_range(q, s, t, i, order, block_pairs, clip=clip)
        ranges.append(r)

    fw = must_open(ofile, "w")

    tracks = []
    print >> sys.stderr, "Chain started: {0} blocks".format(len(ranges))
    iteration = 0
    while ranges:
        if iteration >= opts.iter:
            break

        selected, score = range_chain(ranges)
        tracks.append(selected)
        selected = set(x.id for x in selected)
        if trackids:
            print >> fwlog, ",".join(str(x) for x in sorted(selected))

        ranges = [x for x in ranges if x.id not in selected]
        msg = "Chain {0}: score={1}".format(iteration, score)
        if ranges:
            msg += " {0} blocks remained..".format(len(ranges))
        else:
            msg += " done!"

        print >> sys.stderr, msg
        iteration += 1

    mbed = []
    for b in bed:
        id = b.accn
        atoms = []
        for track in tracks:
            track_ids = [x.id for x in track]
            for tid in track_ids:
                pairs = block_pairs[tid]
                anchor = pairs.get(id, ".")
                if anchor != ".":
                    break
            if ascii and anchor != ".":
                anchor = "x"
            atoms.append(anchor)
        mbed.append((id, atoms))

    for id, atoms in mbed:
        sep = "" if ascii else "\t"
        if mergetandem:
            for i, atom in enumerate(atoms):
                atoms[i] = tandems.get(atom, atom)
        print >> fw, "\t".join((id, sep.join(atoms)))

    logging.debug("MCscan blocks written to `{0}`.".format(ofile))
    if trackids:
        logging.debug("Block IDs written to `{0}`.".format(olog))


def write_details(fw, details, bed):
    """
    Write per gene depth to file
    """
    for a, b, depth in details:
        for i in xrange(a, b):
            gi = bed[i].accn
            print >> fw, "\t".join((gi, str(depth)))


def depth(args):
    """
    %prog depth anchorfile --qbed qbedfile --sbed sbedfile

    Calculate the depths in the two genomes in comparison, given in --qbed and
    --sbed. The synteny blocks will be layered on the genomes, and the
    multiplicity will be summarized to stderr.
    """
    from jcvi.utils.range import range_depth

    p = OptionParser(depth.__doc__)
    p.add_option("--depthfile",
                 help="Generate file with gene and depth [default: %default]")
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    depthfile = opts.depthfile
    ac = AnchorFile(anchorfile)
    qranges = []
    sranges = []
    blocks = ac.blocks
    for ib in blocks:
        q, s, t = zip(*ib)
        q = [qorder[x] for x in q]
        s = [sorder[x] for x in s]
        qrange = (min(q)[0], max(q)[0])
        srange = (min(s)[0], max(s)[0])
        qranges.append(qrange)
        sranges.append(srange)
        if is_self:
            qranges.append(srange)

    qgenome = op.basename(qbed.filename).split(".")[0]
    sgenome = op.basename(sbed.filename).split(".")[0]
    print >> sys.stderr, "Genome {0} depths:".format(qgenome)
    ds, details = range_depth(qranges, len(qbed))
    if depthfile:
        fw = open(depthfile, "w")
        write_details(fw, details, qbed)

    if is_self:
        return

    print >> sys.stderr, "Genome {0} depths:".format(sgenome)
    ds, details = range_depth(sranges, len(sbed))
    if depthfile:
        write_details(fw, details, sbed)
        fw.close()
        logging.debug("Depth written to `{0}`.".format(depthfile))


def get_blocks(scaffold, bs, order, xdist=20, ydist=20, N=6):
    points = []
    for b in bs:
        accn = b.accn.rsplit(".", 1)[0]
        if accn not in order:
            continue
        x, xx = order[accn]
        y = (b.start + b.end) / 2
        points.append((x, y))

    #print scaffold, points
    blocks = synteny_scan(points, xdist, ydist, N)
    return blocks


def breakpoint(args):
    """
    %prog breakpoint blastfile bedfile [options]

    Identify breakpoints where collinearity ends. `blastfile` contains mapping
    from markers (query) to scaffolds (subject). `bedfile` contains marker
    locations in the related species.
    """
    from jcvi.formats.blast import bed
    from jcvi.utils.range import range_interleave

    p = OptionParser(breakpoint.__doc__)
    p.add_option("--xdist", type="int", default=20,
                 help="xdist (in related genome) cutoff [default: %default]")
    p.add_option("--ydist", type="int", default=200000,
                 help="ydist (in current genome) cutoff [default: %default]")
    p.add_option("-n", type="int", default=5,
                 help="number of markers in a block [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blastfile, bedfile = args
    order = Bed(bedfile).order
    blastbedfile = bed([blastfile])
    bbed = Bed(blastbedfile)
    for scaffold, bs in bbed.sub_beds():
        blocks = get_blocks(scaffold, bs, order,
                            xdist=opts.xdist, ydist=opts.ydist, N=opts.n)
        sblocks = []
        for block in blocks:
            xx, yy = zip(*block)
            sblocks.append((scaffold, min(yy), max(yy)))
        iblocks = range_interleave(sblocks)
        for ib in iblocks:
            ch, start, end = ib
            print "{0}\t{1}\t{2}".format(ch, start - 1, end)


def scan(args):
    """
    %prog scan blastfile anchor_file [options]

    pull out syntenic anchors from blastfile based on single-linkage algorithm
    """
    p = OptionParser(scan.__doc__)
    p.add_option("-n", "--min_size", dest="n", type="int", default=4,
            help="minimum number of anchors in a cluster [default: %default]")
    p.add_option("--liftover",
            help="Scan BLAST file to find extra anchors [default: %default]")
    p.set_stripnames()

    blast_file, anchor_file, dist, opts = add_options(p, args, dist=20)
    qbed, sbed, qorder, sorder, is_self = check_beds(blast_file, p, opts)

    filtered_blast = read_blast(blast_file, qorder, sorder, \
                                is_self=is_self, ostrip=False)

    fw = open(anchor_file, "w")
    logging.debug("Chaining distance = {0}".format(dist))

    clusters = batch_scan(filtered_blast, xdist=dist, ydist=dist, N=opts.n)
    for cluster in clusters:
        print >>fw, "###"
        for qi, si, score in cluster:
            query, subject = qbed[qi].accn, sbed[si].accn
            print >>fw, "\t".join((query, subject, str(int(score))))

    fw.close()
    summary([anchor_file])

    lo = opts.liftover
    if not lo:
        return anchor_file

    bedopts = ["--qbed=" + opts.qbed, "--sbed=" + opts.sbed]
    ostrip = [] if opts.strip_names else ["--no_strip_names"]
    newanchorfile = liftover([lo, anchor_file] + bedopts + ostrip)
    return newanchorfile


def liftover(args):
    """
    %prog liftover blastfile anchorfile [options]

    Typical use for this program is given a list of anchors (syntennic
    genes), choose from the blastfile the pairs that are close to the anchors.

    Anchorfile has the following format, each row defines a pair.

        geneA geneB
        geneC geneD
    """
    p = OptionParser(liftover.__doc__)
    p.set_stripnames()

    blast_file, anchor_file, dist, opts = add_options(p, args)
    qbed, sbed, qorder, sorder, is_self = check_beds(blast_file, p, opts)

    filtered_blast = read_blast(blast_file, qorder, sorder,
                            is_self=is_self, ostrip=opts.strip_names)
    blast_to_score = dict(((b.qi, b.si), int(b.score)) for b in filtered_blast)
    accepted = dict(((b.query, b.subject), str(int(b.score))) \
                     for b in filtered_blast)

    ac = AnchorFile(anchor_file)
    all_hits = group_hits(filtered_blast)
    all_anchors, anchor_to_block = read_anchors(ac, qorder, sorder)

    # select hits that are close to the anchor list
    lifted = 0
    for chr_pair in sorted(all_anchors.keys()):
        hits = np.array(all_hits[chr_pair])
        anchors = np.array(all_anchors[chr_pair])

        if not len(hits):
            continue

        for point, nearest in synteny_liftover(hits, anchors, dist):
            qi, si = point[:2]
            block_id = anchor_to_block[nearest]
            query, subject = qbed[qi].accn, sbed[si].accn
            score = blast_to_score[(qi, si)]

            ac.blocks[block_id].append((query, subject, str(score) + "L"))
            lifted += 1

    logging.debug("{0} new pairs found.".format(lifted))
    newanchorfile = anchor_file.rsplit(".", 1)[0] + ".lifted.anchors"
    ac.print_to_file(filename=newanchorfile, accepted=accepted)
    summary([newanchorfile])

    return newanchorfile


if __name__ == '__main__':
    main()
