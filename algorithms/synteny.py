#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys
import logging

import numpy as np
from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.blast import BlastLine
from jcvi.formats.base import BaseFile, read_block
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.utils.range import Range, range_chain
from jcvi.formats.base import must_open
from jcvi.apps.base import ActionDispatcher, debug
debug()


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

    def iter_pairs(self):
        fp = open(self.filename)
        block_id = -1
        for row in fp:
            if row[0] == '#':
                block_id += 1
                continue
            a, b = row.split()[:2]
            yield a, b, block_id

    def print_to_file(self, filename="stdout", accepted=None):
        fw = must_open(filename, "w")
        blocks = self.blocks
        nremoved = 0
        for block in blocks:
            print >> fw, "###"
            for line in block:
                a, b = line[:2]
                if accepted and (a, b) not in accepted:
                    nremoved += 1
                    continue
                print >> fw, "\t".join(str(x) for x in line)
        fw.close()

        logging.debug("Removed {0} existing anchors.".format(nremoved))
        logging.debug("Anchors written to `{0}`.".format(filename))


class BlockFile (BaseFile):

    def __init__(self, filename):
        super(BlockFile, self).__init__(filename)
        fp = must_open(filename)
        data = []
        for row in fp:
            atoms = row.rstrip().split("\t")
            data.append(atoms)

        self.data = data
        self.columns = zip(*data)
        self.ncols = len(self.columns)

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
            print >> sys.stderr, "Column {0}: {1} - {2}".\
                    format(i, start.accn, end.accn)
            print >> sys.stderr, "  {0} .. {1} ({2}) features .. {3}".\
                    format(chr, ngenes, len(ocol), orientation)

        span = abs(start.start - end.end)

        return start, end, si, ei, chr, orientation, span

    def iter_pairs(self, i, j):
        for d in self.data:
            a, b = d[i], d[j]
            if "." in (a, b):
                continue

            yield a, b


def _score(cluster):
    """
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y, scores = zip(*cluster)
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

    logging.debug("A total of {0} BLAST matches imported.".format(len(filtered_blast)))

    return filtered_blast


def read_anchors(ac, qorder, sorder):
    """
    anchors file are just (geneA, geneB) pairs (with possible deflines)
    """
    all_anchors = defaultdict(list)
    nanchors = 0
    anchor_to_block = {}
    block_extent = defaultdict(list)

    for a, b, idx in ac.iter_pairs():
        if a not in qorder or b not in sorder:
            continue
        qi, q = qorder[a]
        si, s = sorder[b]
        pair = (qi, si)
        block_extent[idx].append(pair)

        all_anchors[(q.seqid, s.seqid)].append(pair)
        anchor_to_block[pair] = idx
        nanchors += 1

    logging.debug("A total of {0} anchors imported.".format(nanchors))
    assert nanchors == len(anchor_to_block)

    for idx, dd in block_extent.items():
        xs, ys = zip(*dd)
        block_extent[idx] = (min(xs), max(xs), min(ys), max(ys))

    return all_anchors, anchor_to_block, block_extent


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


def add_beds(p):

    p.add_option("--qbed", help="Path to qbed (required)")
    p.add_option("--sbed", help="Path to sbed (required)")


def check_beds(p, opts):

    if not (opts.qbed and opts.sbed):
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


def add_options(p, args):
    """
    scan and liftover has similar interfaces, so share common options
    returns opts, files
    """
    add_beds(p)
    p.add_option("--dist", default=10, type="int",
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
        ('liftover', 'given anchor list, pull adjancent pairs from blast file'),
        ('mcscan', 'stack synteny blocks on a reference bed'),
        ('stats', 'provide statistics for mscan blocks'),
        ('depth', 'calculate the depths in the two genomes in comparison'),
        ('group', 'cluster the anchors into ortho-groups'),
        ('breakpoint', 'identify breakpoints where collinearity ends'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary anchorfile

    Provide statistics for pairwise blocks.
    """
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
    ac = AnchorFile(anchorfile)
    clusters = ac.blocks

    nclusters = len(clusters)
    nanchors = [len(c) for c in clusters]
    print >> sys.stderr, "A total of {0} anchors found in {1} clusters.".\
                  format(sum(nanchors), nclusters)
    print >> sys.stderr, SummaryStats(nanchors)


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
    %prog mcscan bedfile anchorfile

    Stack synteny blocks on a reference bed, MCSCAN style. The first column in
    the output is the reference order, given in the bedfile. Then each column
    next to it are separate 'tracks'.
    """
    p = OptionParser(mcscan.__doc__)
    p.add_option("--iter", default=100, type="int",
                 help="Max number of chains to output [default: %default]")
    p.add_option("--ascii", default=False, action="store_true",
                 help="Output symbols rather than gene names [default: %default]")
    p.add_option("--Nm", default=10, type="int",
                 help="Clip block ends to allow slight overlaps [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, anchorfile = args
    ascii = opts.ascii
    clip = opts.Nm
    bed = Bed(bedfile)
    order = bed.order

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

    tracks = []
    print >> sys.stderr, "Chain started: {0} blocks".format(len(ranges))
    iteration = 0
    while ranges:
        if iteration >= opts.iter:
            break

        selected, score = range_chain(ranges)
        tracks.append(selected)
        selected = set(x.id for x in selected)
        ranges = [x for x in ranges if x.id not in selected]
        msg = "Chain {0}: score={1}".format(iteration, score)
        if ranges:
            msg += " {0} blocks remained..".format(len(ranges))
        else:
            msg += " done!"

        print >> sys.stderr, msg
        iteration += 1

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

        sep = "" if ascii else "\t"
        print "\t".join((id, sep.join(atoms)))


def group(args):
    """
    %prog group anchorfiles

    Group the anchors into ortho-groups. Can input multiple anchor files.
    """
    p = OptionParser(group.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    anchorfiles = args
    groups = Grouper()
    for anchorfile in anchorfiles:
        ac = AnchorFile(anchorfile)
        for a, b in ac.iter_pairs():
            groups.join(a, b)

    ngroups = len(groups)
    nmembers = sum(len(x) for x in groups)
    logging.debug("Created {0} groups with {1} members.".\
                  format(ngroups, nmembers))

    for g in groups:
        print ",".join(sorted(g))


def depth(args):
    """
    %prog depth anchorfile --qbed qbedfile --sbed sbedfile

    Calculate the depths in the two genomes in comparison, given in --qbed and
    --sbed. The synteny blocks will be layered on the genomes, and the
    multiplicity will be summarized to stderr.
    """
    from jcvi.utils.range import range_depth

    p = OptionParser(depth.__doc__)
    add_beds(p)

    opts, args = p.parse_args(args)
    qbed, sbed, qorder, sorder, is_self = check_beds(p, opts)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorfile, = args
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

    qgenome = qbed.filename.split(".")[0]
    sgenome = sbed.filename.split(".")[0]
    print >> sys.stderr, "Genome {0} depths:".format(qgenome)
    range_depth(qranges, len(qbed))
    if is_self:
        return

    print >> sys.stderr, "Genome {0} depths:".format(sgenome)
    range_depth(sranges, len(sbed))


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
    %prog breakpoint blastfile bedfile

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
    key = lambda x: x[1]
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
    p.add_option("-n", type="int", default=5,
            help="minimum number of anchors in a cluster [default: %default]")
    p.add_option("--liftover",
            help="Scan BLAST file to find extra anchors [default: %default]")

    blast_file, anchor_file, dist, opts = add_options(p, args)
    qbed, sbed, qorder, sorder, is_self = check_beds(p, opts)

    filtered_blast = read_blast(blast_file, qorder, sorder, is_self=is_self)

    fw = open(anchor_file, "w")
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
    newanchorfile = liftover([lo, anchor_file] + bedopts)
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

    p.add_option("--no_strip_names", dest="strip_names",
            action="store_false", default=True,
            help="do not strip alternative splicing "
            "(e.g. At5g06540.1 -> At5g06540)")

    blast_file, anchor_file, dist, opts = add_options(p, args)
    qbed, sbed, qorder, sorder, is_self = check_beds(p, opts)

    filtered_blast = read_blast(blast_file, qorder, sorder,
                            is_self=is_self, ostrip=opts.strip_names)
    blast_to_score = dict(((b.qi, b.si), int(b.score)) for b in filtered_blast)
    accepted = set((b.query, b.subject) for b in filtered_blast)

    ac = AnchorFile(anchor_file)
    all_hits = group_hits(filtered_blast)
    all_anchors, anchor_to_block, block_extent = read_anchors(ac, qorder, sorder)

    # select hits that are close to the anchor list
    fw = sys.stdout
    lifted = 0
    for chr_pair in sorted(all_anchors.keys()):
        hits = np.array(all_hits[chr_pair])
        anchors = np.array(all_anchors[chr_pair])

        #logging.debug("%s: %d" % (chr_pair, len(anchors)))
        if not len(hits):
            continue

        for point, nearest in synteny_liftover(hits, anchors, dist):
            qi, si = point[:2]
            block_id = anchor_to_block[nearest]
            query, subject = qbed[qi].accn, sbed[si].accn
            score = blast_to_score[(qi, si)]

            xmin, xmax, ymin, ymax = block_extent[block_id]
            # Lifted pairs cannot be outside the bounding box
            if qi < xmin or qi > xmax or si < ymin or si > ymax:
                continue

            ac.blocks[block_id].append((query, subject, str(score) + "L"))
            lifted += 1

    logging.debug("{0} new pairs found.".format(lifted))
    newanchorfile = anchor_file.rsplit(".", 1)[0] + ".lifted.anchors"
    ac.print_to_file(filename=newanchorfile, accepted=accepted)
    summary([newanchorfile])

    return newanchorfile


if __name__ == '__main__':
    main()
