#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""Syntenty inference in comparative genomics
"""

import os.path as op
import sys
import logging

import numpy as np
from collections import defaultdict

try:
    from collections.abc import Iterable
except ImportError:
    from collections import Iterable

from jcvi.algorithms.lis import heaviest_increasing_subsequence as his
from jcvi.formats.bed import Bed, BedLine
from jcvi.formats.blast import Blast
from jcvi.formats.base import BaseFile, SetFile, read_block, must_open
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name, human_size
from jcvi.utils.range import Range, range_chain
from jcvi.apps.base import OptionParser, ActionDispatcher


class AnchorFile(BaseFile):
    def __init__(self, filename, minsize=0):
        super(AnchorFile, self).__init__(filename)
        self.blocks = list(self.iter_blocks(minsize=minsize))

    def iter_blocks(self, minsize=0):
        fp = open(self.filename)
        for _, lines in read_block(fp, "#"):
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

    def make_ranges(self, order, clip=10):
        """Prepare anchors information into a set of ranges for chaining"""
        ranges = []
        block_pairs = defaultdict(dict)
        blocks = self.blocks
        for i, ib in enumerate(blocks):
            q, s, t = zip(*ib)
            if q[0] not in order:
                q, s = s, q

            r = make_range(q, s, t, i, order, block_pairs, clip=clip)
            ranges.append(r)

            assert q[0] in order
            if s[0] not in order:
                continue

            # is_self comparison
            q, s = s, q
            r = make_range(q, s, t, i, order, block_pairs, clip=clip)
            ranges.append(r)
        return ranges, block_pairs

    def print_to_file(self, filename="stdout", accepted=None):
        fw = must_open(filename, "w")
        blocks = self.blocks
        nremoved = 0
        ncorrected = 0
        for block in blocks:
            print("###", file=fw)
            for line in block:
                a, b, score = line
                pair = (a, b)
                if accepted:
                    if pair not in accepted:
                        nremoved += 1
                        continue
                    av = accepted[pair]
                    if score != av and score != av + "L":
                        score = av
                        ncorrected += 1
                print("\t".join((a, b, score)), file=fw)
        fw.close()

        logging.debug("Removed %d existing anchors.", nremoved)
        logging.debug("Corrected scores for %d anchors.", ncorrected)
        logging.debug("Anchors written to `%s`.", filename)

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
        for a, b, _ in self.iter_pairs():
            if (a, b) in blasts:
                bline = blasts[(a, b)]
            elif (b, a) in blasts:
                bline = blasts[(b, a)]
            else:
                line = "\t".join((a, b))
                bline = BlastLineByConversion(line, mode="110000000000")

            print(bline, file=fw)
            nlines += 1
        fw.close()

        logging.debug("A total of %d BLAST lines written to `%s`.", nlines, outfile)

        return outfile

    @property
    def is_empty(self):
        blocks = self.blocks
        return not blocks or not blocks[0]


class BlockFile(BaseFile):
    """Parse .blocks file which is the mcscan output with multiple columns as 'tracks'"""

    def __init__(self, filename, defaultcolor="#fb8072", header=False):
        super(BlockFile, self).__init__(filename)
        fp = must_open(filename)
        hd = next(fp).rstrip().split("\t")
        ncols = len(hd)
        if header:
            self.header = hd
        else:
            fp.seek(0)
            self.header = range(ncols)

        data = []
        highlight = []
        for row in fp:
            hl = "*" in row
            # r* highlights the block in red color
            if hl:
                hl, row = row.split("*", 1)
                hl = hl or defaultcolor
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
        self.columns = list(zip(*data))
        self.ncols = ncols

    def get_extent(self, i, order, debug=True):
        # Some blocks file, such as ones manually edited, will have garbled
        # order, which prompts the hack below
        acol = [order[x][0] for x in self.columns[0] if x in order]
        bcol = [order[x][0] for x in self.columns[i] if x in order]
        elen = min(len(acol), len(bcol))
        ia, ib = acol[:elen], bcol[:elen]
        orientation = get_orientation(ia, ib)

        ocol = [order[x] for x in self.columns[i] if x in order]
        # orientation = '+' if ocol[0][0] <= ocol[-1][0] else '-'
        si, start = min(ocol)
        ei, end = max(ocol)
        same_chr = start.seqid == end.seqid
        chr = start.seqid if same_chr else None
        ngenes = ei - si + 1
        if debug:
            r = "{0}:{1}-{2}".format(chr, start.start, end.end)
            print(
                "Column {0}: {1} - {2} ({3})".format(i, start.accn, end.accn, r),
                file=sys.stderr,
            )
            print(
                "  {0} .. {1} ({2}) features .. {3}".format(
                    chr, ngenes, len(ocol), orientation
                ),
                file=sys.stderr,
            )

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
        for i in range(ncols):
            for j in range(i + 1, ncols):
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
            for i in range(qi - 1, -1, -1):
                if col[i] not in (".", ""):
                    upstream = col[i]
                    upstream_dist = qi - i
                    break
            # search downstream
            for i in range(qi, ndata):
                if col[i] not in (".", ""):
                    downstream = col[i]
                    downstream_dist = i - qi
                    break
            closest = upstream if upstream_dist < downstream_dist else downstream
            # output in .simple format
            if invert:
                line = "\t".join(str(x) for x in (closest, closest, gene, gene, 0, "+"))
            else:
                line = "\t".join(str(x) for x in (gene, gene, closest, closest, 0, "+"))
            if color is not None:
                line = color + "*" + line
            yield line

    def grouper(self) -> Grouper:
        """Build orthogroup based on the gene matches."""
        grouper = Grouper()
        for row in self.data:
            if "." not in row:
                grouper.join(*row)
        logging.debug("A total of %d orthogroups formed", len(grouper))
        return grouper


class SimpleFile(object):
    def __init__(self, simplefile, defaultcolor="#fb8072", order=None):
        # Sometimes the simplefile has query and subject wrong
        fp = open(simplefile)
        self.blocks = []
        check = False
        for row in fp:
            if row[:2] == "##" or row.startswith("StartGeneA"):
                continue
            hl = "*" in row
            if hl:
                hl, row = row.split("*", 1)
                hl = hl or defaultcolor
            a, b, c, d, score, orientation = row.split()
            if order and a not in order:
                if c not in order:
                    check = True
                    print(
                        """{} {} {} {} can not found in bed files.""".format(
                            a, b, c, d
                        ),
                        file=sys.stderr,
                    )
                else:
                    a, b, c, d = c, d, a, b
            if orientation == "-":
                c, d = d, c
            score = int(score)
            self.blocks.append((a, b, c, d, score, orientation, hl))
        if check:
            print(
                "Error: some genes in blocks can't be found, please rerun after making sure that bed file agree with simple file.",
                file=sys.stderr,
            )
            exit(1)


def _score(cluster):
    """
    score of the cluster, in this case, is the number of non-repetitive matches
    """
    x, y = list(zip(*cluster))[:2]
    return min(len(set(x)), len(set(y)))


def get_orientation(ia, ib):
    """Infer the orientation of a pairwise block.

    Args:
        ia (List[int]): List a
        ib (List[int]): List b

    Returns:
        str: plus (+) or minus (-)
    """
    if len(ia) != len(ib) or len(ia) < 2:
        return "+"  # Just return a default orientation

    slope, _ = np.polyfit(ia, ib, 1)
    return "+" if slope >= 0 else "-"


def group_hits(blasts):
    if not blasts:
        return {"": []}

    # Already in the form of (qi, si, score)
    if isinstance(blasts[0], Iterable) and len(blasts[0]) == 3:
        return {"": blasts}

    # grouping the hits based on chromosome pair
    all_hits = defaultdict(list)
    for b in blasts:
        all_hits[(b.qseqid, b.sseqid)].append((b.qi, b.si, b.score))

    return all_hits


def read_blast(blast_file, qorder, sorder, is_self=False, ostrip=True):
    """Read the blast and convert name into coordinates"""
    filtered_blast = []
    seen = set()
    bl = Blast(blast_file)
    for b in bl:
        query, subject = b.query, b.subject
        if is_self and query == subject:
            continue
        if ostrip:
            query, subject = gene_name(query), gene_name(subject)
        if query not in qorder or subject not in sorder:
            continue

        qi, q = qorder[query]
        si, s = sorder[subject]

        if is_self:
            # remove redundant a<->b to one side when doing self-self BLAST
            if qi > si:
                query, subject = subject, query
                qi, si = si, qi
                q, s = s, q
            # Too close to diagonal! possible tandem repeats
            if q.seqid == s.seqid and si - qi < 40:
                continue

        key = query, subject
        if key in seen:
            continue
        seen.add(key)

        b.qseqid, b.sseqid = q.seqid, s.seqid
        b.qi, b.si = qi, si
        b.query, b.subject = query, subject

        filtered_blast.append(b)

    logging.debug(
        "A total of %d BLAST imported from `%s`.", len(filtered_blast), blast_file
    )

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


def synteny_scan(points, xdist, ydist, N, is_self=False, intrabound=300):
    """
    This is the core single linkage algorithm which behaves in O(n):
    iterate through the pairs, foreach pair we look back on the
    adjacent pairs to find links
    """
    clusters = Grouper()
    n = len(points)
    points.sort()
    for i in range(n):
        for j in range(i - 1, -1, -1):
            # x-axis distance
            del_x = points[i][0] - points[j][0]
            if del_x > xdist:
                break
            # y-axis distance
            del_y = points[i][1] - points[j][1]
            if abs(del_y) > ydist:
                continue
            # In self-comparison, ignore the anchors that are too close to the diagonal
            if is_self:
                intradist = min(
                    abs(points[i][0] - points[i][1]), abs(points[j][0] - points[j][1])
                )
                if intradist < intrabound:
                    continue
            # otherwise join
            clusters.join(points[i], points[j])

    # select clusters that are at least >=N
    clusters = [sorted(cluster) for cluster in list(clusters) if _score(cluster) >= N]

    return clusters


def batch_scan(points, xdist=20, ydist=20, N=5, is_self=False, intrabound=300):
    """
    runs synteny_scan() per chromosome pair
    """
    chr_pair_points = group_hits(points)

    clusters = []
    for chr_pair in sorted(chr_pair_points.keys()):
        points = chr_pair_points[chr_pair]
        clusters.extend(
            synteny_scan(
                points, xdist, ydist, N, is_self=is_self, intrabound=intrabound
            )
        )

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
    dists, idxs = tree.query(ppoints, p=1, distance_upper_bound=dist)

    for point, dist, idx in zip(points, dists, idxs):
        if idx == tree.n:  # nearest is out of range
            continue
        if dist == 0:  # already in anchors
            continue

        yield point, tuple(anchors[idx])


def get_bed_filenames(hintfile, p, opts):
    wd, hintfile = op.split(hintfile)
    if not (opts.qbed and opts.sbed):
        try:
            q, s = hintfile.split(".", 2)[:2]
            opts.qbed = op.join(wd, q + ".bed")
            opts.sbed = op.join(wd, s + ".bed")
            logging.debug("Assuming --qbed={0} --sbed={1}".format(opts.qbed, opts.sbed))
        except:
            print("Options --qbed and --sbed are required", file=sys.stderr)
            sys.exit(not p.print_help())

    return opts.qbed, opts.sbed


def check_beds(hintfile, p, opts, sorted=True):
    qbed_file, sbed_file = get_bed_filenames(hintfile, p, opts)
    # is this a self-self blast?
    is_self = qbed_file == sbed_file
    if is_self:
        logging.debug("Looks like self-self comparison.")

    qbed = Bed(opts.qbed, sorted=sorted)
    sbed = Bed(opts.sbed, sorted=sorted)
    qorder = qbed.order
    sorder = sbed.order

    return qbed, sbed, qorder, sorder, is_self


def add_options(p, args, dist=10):
    """
    scan and liftover has similar interfaces, so share common options
    returns opts, files
    """
    p.set_beds()
    p.add_option(
        "--dist", default=dist, type="int", help="Extent of flanking regions to search"
    )

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blast_file, anchor_file = args

    return blast_file, anchor_file, opts.dist, opts


def main():

    actions = (
        ("scan", "get anchor list using single-linkage algorithm"),
        ("summary", "provide statistics for pairwise blocks"),
        ("liftover", "given anchor list, pull adjacent pairs from blast file"),
        # Multiple synteny blocks inference
        ("mcscan", "stack synteny blocks on a reference bed"),
        ("mcscanq", "query multiple synteny blocks"),
        # Assemble multiple synteny blocks
        ("query", "collect matching region based on the query region"),
        ("assemble", "build blocks from regions defined by start and end"),
        # Filter synteny blocks
        ("screen", "extract subset of blocks from anchorfile"),
        ("simple", "convert anchorfile to simple block descriptions"),
        ("stats", "provide statistics for mscan blocks"),
        ("depth", "calculate the depths in the two genomes in comparison"),
        ("breakpoint", "identify breakpoints where collinearity ends"),
        ("matrix", "make oxford grid based on anchors file"),
        ("coge", "convert CoGe file to anchors file"),
        ("spa", "convert chr ordering from SPA to simple lists"),
        ("layout", "compute layout based on .simple file"),
        ("rebuild", "rebuild anchors file from prebuilt blocks file"),
        # Formatting
        ("fromaligns", "convert aligns file to anchors file"),
        ("toaligns", "convert anchors file to aligns file"),
    )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_region_size(region, bed, order):
    """Get a summary of a syntenic region, how many anchors it has and
    how many genes it spans.

    Args:
        region (List[str]): List of gene ids
        order (Dict[str, BedLine]): Bed order to retrieve the positions

    Returns:
        Tuple of three strs and two ints, start / end gene / seqid of the
        region and total anchor counts and the span (number of genes)
    """
    ris = [order[x] for x in region]
    min_ri, min_r = min(ris)
    max_ri, max_r = max(ris)
    anchor_count = len(region)
    span = max_ri - min_ri + 1
    min_seqid = min_r.seqid
    max_seqid = max_r.seqid
    assert min_seqid == max_seqid, "SeqId do not match, region invalid"
    return min_r.accn, max_r.accn, min_seqid, span, anchor_count


def query(args):
    """
    %prog query anchorsfile startGeneId endGeneId

    Collect matching region based on query region as given by startGeneId to
    endGeneId. This can be considered a local version of mcscan(). The bedfile
    must contain the range from startGeneId to endGeneId.

    Typical pipeline is to extract a set of pairwise syntenic regions to the
    selected region of interest and then assemble them into .blocks file for
    plotting purposes.
    """
    p = OptionParser(query.__doc__)
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    anchorsfile, start_gene_id, end_gene_id = args
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts)

    # Guess which is qbed, which is sbed
    if start_gene_id in sorder:  # flip query and subject
        qbed, sbed = sbed, qbed
        qorder, sorder = sorder, qorder

    ac = AnchorFile(anchorsfile)
    blocks = ac.blocks
    si, s = qorder[start_gene_id]
    ei, e = qorder[end_gene_id]
    target_region = qbed[si : ei + 1]
    target_genes = set(x.accn for x in target_region)

    # Go through all the blocks and pick out all matching regions
    regions = []
    for block in blocks:
        matching_region = set()
        for a, b, score in block:
            if not (a in target_genes or b in target_genes):
                continue
            if a in target_genes:
                matching_region.add(b)
            else:
                matching_region.add(a)
        if len(matching_region) < 2:
            continue
        # Print a summary of the matching region
        regions.append(get_region_size(matching_region, sbed, sorder))

    for min_accn, max_accn, seqid, span, anchor_count in sorted(
        regions, key=lambda x: (-x[-1], -x[-2])  # Sort by (anchor_count, span) DESC
    ):
        print(
            "{} {} ({}): span {}, anchors {}".format(
                min_accn, max_accn, seqid, span, anchor_count
            )
        )


def assemble(args):
    """
    %prog assemble regionsfile all.bed all.cds

    Assemble blocks file based on regions file. Regions file may look like:

    amborella evm_27.model.AmTr_v1.0_scaffold00004.87 evm_27.model.AmTr_v1.0_scaffold00004.204
    apostasia Ash010455 Ash010479 (fragScaff_scaffold_5)
    apostasia Ash018328 Ash018367 (original_scaffold_2912)
    apostasia Ash007533 Ash007562 (fragScaff_scaffold_132)
    apostasia Ash002281 Ash002299 (fragScaff_scaffold_86)

    Where each line lists a region, starting with the species name (species.bed
    must be present in the current directory). Followed by start and end gene.
    Contents after the 3rd field (end gene) are ignored. Using the example
    above, the final .blocks file will contain 5 columns, one column for each line.
    """
    import shutil
    from tempfile import mkdtemp, mkstemp

    from jcvi.apps.align import last
    from jcvi.formats.fasta import some

    p = OptionParser(assemble.__doc__)
    p.add_option(
        "--no_strip_names",
        default=False,
        action="store_true",
        help="Do not strip alternative splicing (e.g. At5g06540.1 -> At5g06540)",
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    strip_names = not opts.no_strip_names
    regionsfile, bedfile, cdsfile = args
    species_beds = {}
    column_genes = []
    pivot = None
    with open(regionsfile) as fp:
        for row in fp:
            species, start, end = row.split()[:3]
            if pivot is None:
                pivot = species
            if species not in species_beds:
                species_beds[species] = Bed(species + ".bed")
            bed = species_beds[species]
            order = bed.order
            si, s = order[start]
            ei, e = order[end]
            genes = set(x.accn for x in bed[si : ei + 1])
            column_genes.append(genes)

    # Write gene ids
    workdir = mkdtemp()
    fd, idsfile = mkstemp(dir=workdir)
    with open(idsfile, "w") as fw:
        for genes in column_genes:
            print(" ".join(genes), file=fw)

        logging.debug("Gene ids written to `{}`".format(idsfile))

    # Extract FASTA
    fd, fastafile = mkstemp(dir=workdir)
    some_args = [cdsfile, idsfile, fastafile]
    if not strip_names:
        some_args += ["--no_strip_names"]
    some(some_args)

    # Perform self-comparison and collect all pairs
    last_output = last([fastafile, fastafile, "--outdir", workdir])
    blast = Blast(last_output)
    pairs = set()
    for b in blast:
        query, subject = b.query, b.subject
        if strip_names:
            query, subject = gene_name(query), gene_name(subject)
        pairs.add((query, subject))
    logging.debug("Extracted {} gene pairs from `{}`".format(len(pairs), last_output))

    # Sort the pairs into columns
    N = len(column_genes)
    all_slots = []
    for i in range(N):
        for j in range(i + 1, N):
            genes_i = column_genes[i]
            genes_j = column_genes[j]
            for a, b in pairs:
                if not (a in genes_i and b in genes_j):
                    continue
                slots = ["."] * N
                slots[i] = a
                slots[j] = b
                all_slots.append(slots)

    # Compress the pairwise results and merge when possible
    # TODO: This is currently not optimized and inefficient
    def is_compatible(slots1, slots2):
        # At least intersects for one gene
        assert len(slots1) == len(slots2)
        flag = False
        for a, b in zip(slots1, slots2):
            if "." in (a, b):
                continue
            if a == b:
                flag = True
            else:
                return False
        return flag

    def merge(slots, processed):
        for i, a in enumerate(slots):
            if processed[i] == "." and a != ".":
                processed[i] = a

    processed_slots = []
    all_slots.sort()
    for slots in all_slots:
        merged = False
        for processed in processed_slots:
            if is_compatible(slots, processed):
                merge(slots, processed)  # Merge into that line
                merged = True
                break
        if not merged:  # New information
            processed_slots.append(slots)

    logging.debug(
        "Before compression: {}, After compression: {}".format(
            len(all_slots), len(processed_slots)
        )
    )

    pivot_order = species_beds[pivot].order
    pivot_max = len(species_beds[pivot])
    pivot_sort_key = lambda x: pivot_order[x[0]][0] if x[0] != "." else pivot_max
    processed_slots.sort(key=pivot_sort_key)

    with must_open(opts.outfile, "w") as fw:
        for slots in processed_slots:
            print("\t".join(slots), file=fw)

    # Cleanup
    shutil.rmtree(workdir)


def colinear_evaluate_weights(tour, data):
    tour = dict((s, i) for i, s in enumerate(tour))
    data = [(tour[x], score) for x, score in data if x in tour]
    return (his(data)[-1],)


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

    print(",".join(tour))


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

    (alignsfile,) = args
    fp = must_open(alignsfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        if row.startswith("## Alignment"):
            print("###", file=fw)
            continue
        if row[0] == "#" or not row.strip():
            continue
        atoms = row.split(":")[-1].split()
        print("\t".join(atoms[:2]), file=fw)
    fw.close()


def toaligns(args):
    """
    %prog fromaligns input.anchors

    Convert anchors file to tab-separated aligns file, adding the first column
    with the Block ID.
    """
    p = OptionParser(toaligns.__doc__)
    p.add_option("--prefix", default="b", help="Prefix to the block id")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (anchorfile,) = args
    ac = AnchorFile(anchorfile)
    logging.debug("A total of {} blocks imported".format(len(ac.blocks)))
    max_block_id_len = len(str(len(ac.blocks) - 1))
    header = "\t".join(("#Block ID", "Gene 1", "Gene 2"))

    with must_open(opts.outfile, "w") as fw:
        print(header, file=fw)
        for a, b, block_id in ac.iter_pairs():
            block_id = "{}{:0{}d}".format(opts.prefix, block_id, max_block_id_len)
            print("\t".join((block_id, a, b)), file=fw)


def mcscanq(args):
    """
    %prog mcscanq query.ids blocksfile

    Query multiple synteny blocks to get the closest alignment feature. Mostly
    used for 'highlighting' the lines in the synteny plot, drawn by
    graphics.karyotype and graphics.synteny.
    """
    p = OptionParser(mcscanq.__doc__)
    p.add_option("--color", help="Add color highlight, used in plotting")
    p.add_option(
        "--invert", default=False, action="store_true", help="Invert query and subject"
    )
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    qids, blocksfile = args
    b = BlockFile(blocksfile)
    fp = open(qids)
    for gene in fp:
        gene = gene.strip()
        for line in b.query_gene(gene, color=opts.color, invert=opts.invert):
            print(line)


def spa(args):
    """
    %prog spa spafiles

    Convert chromosome ordering from SPA to simple lists. First column is the
    reference order.
    """
    from jcvi.algorithms.graph import merge_paths
    from jcvi.utils.cbook import uniqify

    p = OptionParser(spa.__doc__)
    p.add_option(
        "--unmapped",
        default=False,
        action="store_true",
        help="Include unmapped scaffolds in the list",
    )
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
            if row[0] == "#" or not row.strip():
                continue

            atoms = row.rstrip().split("\t")
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
    print("ref", len(ref), ",".join(ref))
    for spafile, mapping, missing in zip(spafiles, mappings, missings):
        mapping = [x for x in mapping if "random" not in x]
        mapping = uniqify(mapping)
        if len(mapping) < 50 and opts.unmapped:
            mapping = uniqify(mapping + missing)

        print(spafile, len(mapping), ",".join(mapping))


def rebuild(args):
    """
    %prog rebuild blocksfile blastfile

    Rebuild anchors file from pre-built blocks file.
    """
    p = OptionParser(rebuild.__doc__)
    p.add_option(
        "--header", default=False, action="store_true", help="First line is header"
    )
    p.add_option(
        "--write_blast",
        default=False,
        action="store_true",
        help="Get blast records of rebuilt anchors",
    )
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blocksfile, blastfile = args
    bk = BlockFile(blocksfile, header=opts.header)
    fw = open("pairs", "w")
    for a, b, h in bk.iter_all_pairs():
        print("\t".join((a, b)), file=fw)
    fw.close()

    if opts.write_blast:
        AnchorFile("pairs").blast(blastfile, "pairs.blast")

    fw = open("tracks", "w")
    for g, col in bk.iter_gene_col():
        print("\t".join(str(x) for x in (g, col)), file=fw)
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

    (cogefile,) = args
    fp = must_open(cogefile)
    cogefile = cogefile.replace(".gz", "")
    ksfile = cogefile + ".ks"
    anchorsfile = cogefile + ".anchors"
    fw_ks = must_open(ksfile, "w")
    fw_ac = must_open(anchorsfile, "w")

    tag = "###"
    print(tag, file=fw_ks)
    for header, lines in read_block(fp, tag):
        print(tag, file=fw_ac)
        lines = list(lines)
        for line in lines:
            if line[0] == "#":
                continue
            (
                ks,
                ka,
                achr,
                a,
                astart,
                astop,
                bchr,
                b,
                bstart,
                bstop,
                ev,
                ss,
            ) = line.split()
            a = a.split("||")[3]
            b = b.split("||")[3]
            print("\t".join((a, b, ev)), file=fw_ac)
            print(",".join((";".join((a, b)), ks, ka, ks, ka)), file=fw_ks)

    fw_ks.close()
    fw_ac.close()


def matrix(args):
    """
    %prog matrix all.bed anchorfile matrixfile

    Make oxford grid based on anchors file.
    """

    p = OptionParser(matrix.__doc__)
    p.add_option("--seqids", help="File with seqids")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    bedfile, anchorfile, matrixfile = args
    ac = AnchorFile(anchorfile)
    seqidsfile = opts.seqids
    if seqidsfile:
        seqids = SetFile(seqidsfile, delimiter=",")

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
    print("\t".join(["o"] + bseqids), file=fw)
    for aseqid in aseqids:
        print("\t".join([aseqid] + [str(m[(aseqid, x)]) for x in bseqids]), file=fw)


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
    p.add_option(
        "--rich", default=False, action="store_true", help="Output additional columns"
    )
    p.add_option(
        "--coords",
        default=False,
        action="store_true",
        help="Output columns with base coordinates",
    )
    p.add_option(
        "--bed",
        default=False,
        action="store_true",
        help="Generate BED file for the blocks",
    )
    p.add_option(
        "--noheader", default=False, action="store_true", help="Don't output header"
    )
    p.set_beds()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorfile,) = args
    additional = opts.rich
    coords = opts.coords
    header = not opts.noheader
    bed = opts.bed
    if bed:
        coords = True
        bbed = Bed()

    ac = AnchorFile(anchorfile)
    simplefile = anchorfile.rsplit(".", 1)[0] + ".simple"

    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)
    pf = "-".join(anchorfile.split(".", 2)[:2])
    if ac.is_empty:
        logging.error("No blocks found in `%s`. Aborting ..", anchorfile)
        return

    if coords:
        h = "Block|Chr|Start|End|Span|StartGene|EndGene|GeneSpan|Orientation"
    else:
        h = "StartGeneA|EndGeneA|StartGeneB|EndGeneB|Orientation|Score"
        if additional:
            h += "|StartOrderA|EndOrderA|StartOrderB|EndOrderB|SizeA|SizeB|Size|Block"

    fws = open(simplefile, "w")
    if header:
        print("\t".join(h.split("|")), file=fws)

    blocks = ac.blocks
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

        orientation = get_orientation(ia, ib)
        aspan = aendi - astarti + 1
        bspan = bendi - bstarti + 1
        score = int((aspan * bspan) ** 0.5)
        score = str(score)
        block_id = pf + "-block-{0}".format(i)

        if coords:

            aseqid, astartbase, aendbase = get_boundary_bases(astart, aend, qorder)
            bseqid, bstartbase, bendbase = get_boundary_bases(bstart, bend, sorder)
            abase = aendbase - astartbase + 1
            bbase = bendbase - bstartbase + 1
            atotalbase += abase
            btotalbase += bbase

            # Write dual lines
            aargs = [
                block_id,
                aseqid,
                astartbase,
                aendbase,
                abase,
                astart,
                aend,
                aspan,
                "+",
            ]
            bargs = [
                block_id,
                bseqid,
                bstartbase,
                bendbase,
                bbase,
                bstart,
                bend,
                bspan,
                orientation,
            ]

            if bed:
                bbed.append(
                    BedLine(
                        "\t".join(
                            str(x)
                            for x in (
                                bseqid,
                                bstartbase - 1,
                                bendbase,
                                "{}:{}-{}".format(aseqid, astartbase, aendbase),
                                size,
                                orientation,
                            )
                        )
                    )
                )

            for args in (aargs, bargs):
                print("\t".join(str(x) for x in args), file=fws)
            continue

        args = [astart, aend, bstart, bend, score, orientation]
        if additional:
            args += [astarti, aendi, bstarti, bendi, sizeA, sizeB, size, block_id]
        print("\t".join(str(x) for x in args), file=fws)

    fws.close()
    logging.debug("A total of {0} blocks written to `{1}`.".format(i + 1, simplefile))

    if coords:
        print(
            "Total block span in {0}: {1}".format(
                qbed.filename, human_size(atotalbase, precision=2)
            ),
            file=sys.stderr,
        )
        print(
            "Total block span in {0}: {1}".format(
                sbed.filename, human_size(btotalbase, precision=2)
            ),
            file=sys.stderr,
        )
        print(
            "Ratio: {0:.1f}x".format(
                max(atotalbase, btotalbase) * 1.0 / min(atotalbase, btotalbase)
            ),
            file=sys.stderr,
        )

    if bed:
        bedfile = simplefile + ".bed"
        bbed.print_to_file(filename=bedfile, sorted=True)
        logging.debug("Bed file written to `{}`".format(bedfile))


def screen(args):
    """
    %prog screen anchorfile newanchorfile --qbed=qbedfile --sbed=sbedfile [options]

    Extract subset of blocks from anchorfile. Provide several options:

    1. Option --ids: a file with IDs, 0-based, comma separated, all in one line.
    2. Option --seqids: only allow seqids in this file.
    3. Option --seqpairs: only allow seqpairs in this file, one per line, e.g. "Chr01,Chr05".
    4. Option --minspan: remove blocks with less span than this.
    5. Option --minsize: remove blocks with less number of anchors than this.
    6. Option --intrabound: remove blocks that are too close to the diagonal on
       self dot plot that are typically artifacts
    """
    from jcvi.utils.range import range_distance

    p = OptionParser(screen.__doc__)
    p.set_beds()
    p.add_option("--ids", help="File with block IDs (0-based)")
    p.add_option("--seqids", help="File with seqids")
    p.add_option("--seqpairs", help="File with seqpairs")
    p.add_option(
        "--intrabound",
        default=300,
        type="int",
        help="Lower bound of intra-chromosomal blocks (only for self comparison)",
    )
    p.add_option("--minspan", default=0, type="int", help="Only blocks with span >=")
    p.add_option("--minsize", default=0, type="int", help="Only blocks with anchors >=")
    p.add_option(
        "--simple", action="store_true", help="Write simple anchorfile with block ends"
    )
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
    intrabound = opts.intrabound
    ids, seqids, seqpairs = None, None, None

    if idsfile:
        ids = SetFile(idsfile, delimiter=",")
        ids = set(int(x) for x in ids)
    if seqidsfile:
        seqids = SetFile(seqidsfile, delimiter=",")
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
        min_ia, max_ia = min(ia), max(ia)
        min_ib, max_ib = min(ib), max(ib)
        aspan = max_ia - min_ia + 1
        bspan = max_ib - min_ib + 1
        aseqid = oa[0].seqid
        bseqid = ob[0].seqid

        if seqids:
            if (aseqid not in seqids) or (bseqid not in seqids):
                continue

        if seqpairs:
            if (aseqid, bseqid) not in seqpairs:
                continue

        same_chromosome = is_self and (aseqid == bseqid)

        if same_chromosome:
            dist, _ = range_distance(
                (aseqid, min_ia, max_ia, "?"), (bseqid, min_ib, max_ib, "?")
            )
            if dist < intrabound:
                continue

        if minsize:
            if len(block) < minsize:
                continue

        if minspan:
            if aspan < minspan or bspan < minspan:
                continue

        selected += 1
        print("###", file=fw)
        for line in block:
            print("\t".join(line), file=fw)

    fw.close()

    if osimple:
        simple(
            [
                newanchorfile,
                "--noheader",
                "--qbed=" + qbed.filename,
                "--sbed=" + sbed.filename,
            ]
        )

    logging.debug("Before: {0} blocks, After: {1} blocks".format(len(blocks), selected))


def summary(args):
    """
    %prog summary anchorfile

    Provide statistics for pairwise blocks.
    """
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(summary.__doc__)
    p.add_option("--prefix", help="Generate per block stats")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorfile,) = args
    ac = AnchorFile(anchorfile)
    clusters = ac.blocks
    if clusters == [[]]:
        logging.debug("A total of 0 anchor was found. Aborted.")
        raise ValueError("A total of 0 anchor was found. Aborted.")

    nclusters = len(clusters)
    nanchors = [len(c) for c in clusters]
    nranchors = [_score(c) for c in clusters]  # non-redundant anchors
    print(
        "A total of {0} (NR:{1}) anchors found in {2} clusters.".format(
            sum(nanchors), sum(nranchors), nclusters
        ),
        file=sys.stderr,
    )
    print("Stats:", SummaryStats(nanchors), file=sys.stderr)
    print("NR stats:", SummaryStats(nranchors), file=sys.stderr)

    prefix = opts.prefix
    if prefix:
        pad = len(str(nclusters))
        for i, c in enumerate(clusters):
            block_id = "{0}{1:0{2}d}".format(prefix, i + 1, pad)
            print("\t".join((block_id, str(len(c)))))


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

    (blocksfile,) = args
    fp = open(blocksfile)
    counts = defaultdict(int)
    total = orthologous = 0
    for row in fp:
        atoms = row.rstrip().split("\t")
        hits = [x for x in atoms[1:] if x != "."]
        counts[len(hits)] += 1
        total += 1
        if atoms[1] != ".":
            orthologous += 1

    print("Total lines: {0}".format(total), file=sys.stderr)
    for i, n in sorted(counts.items()):
        print("Count {0}: {1}".format(i, percentage(n, total)), file=sys.stderr)

    print(file=sys.stderr)

    matches = sum(n for i, n in counts.items() if i != 0)
    print(
        "Total lines with matches: {0}".format(percentage(matches, total)),
        file=sys.stderr,
    )
    for i, n in sorted(counts.items()):
        if i == 0:
            continue

        print("Count {0}: {1}".format(i, percentage(n, matches)), file=sys.stderr)

    print(file=sys.stderr)
    print(
        "Orthologous matches: {0}".format(percentage(orthologous, matches)),
        file=sys.stderr,
    )


def get_best_pair(qs, ss, ts):
    pairs = {}
    for q, s, t in zip(qs, ss, ts):
        t = int(t[:-1]) if t[-1] == "L" else int(t)
        if q not in pairs or pairs[q][1] < t:
            pairs[q] = (s, t)

    # Discard score
    spairs = dict((q, s) for q, (s, t) in pairs.items())
    return spairs


def make_range(q, s, t, i, order, block_pairs, clip=10):
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
    p.add_option(
        "--iter", default=100, type="int", help="Max number of chains to output"
    )
    p.add_option(
        "--ascii",
        default=False,
        action="store_true",
        help="Output symbols rather than gene names",
    )
    p.add_option(
        "--Nm", default=10, type="int", help="Clip block ends to allow slight overlaps"
    )
    p.add_option(
        "--trackids", action="store_true", help="Track block IDs in separate file"
    )
    p.add_option(
        "--mergetandem",
        default=None,
        help="merge tandems genes in output acoording to PATH-TO-TANDEM_FILE, "
        "cannot be used with --ascii",
    )
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
        for row in open(mergetandem):
            row = row.split()
            s = ";".join(row)
            for atom in row:
                tandems[atom] = s

    ac = AnchorFile(anchorfile)
    ranges, block_pairs = ac.make_ranges(order, clip=clip)

    fw = must_open(ofile, "w")

    tracks = []
    print("Chain started: {0} blocks".format(len(ranges)), file=sys.stderr)
    iteration = 0
    while ranges:
        if iteration >= opts.iter:
            break

        selected, score = range_chain(ranges)
        tracks.append(selected)
        selected = set(x.id for x in selected)
        if trackids:
            print(",".join(str(x) for x in sorted(selected)), file=fwlog)

        ranges = [x for x in ranges if x.id not in selected]
        msg = "Chain {0}: score={1}".format(iteration, score)
        if ranges:
            msg += " {0} blocks remained..".format(len(ranges))
        else:
            msg += " done!"

        print(msg, file=sys.stderr)
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
        print("\t".join((id, sep.join(atoms))), file=fw)

    logging.debug("MCscan blocks written to `{0}`.".format(ofile))
    if trackids:
        logging.debug("Block IDs written to `{0}`.".format(olog))


def write_details(fw, details, bed):
    """
    Write per gene depth to file
    """
    for a, b, depth in details:
        for i in range(a, b):
            gi = bed[i].accn
            print("\t".join((gi, str(depth))), file=fw)


def depth(args):
    """
    %prog depth anchorfile --qbed qbedfile --sbed sbedfile

    Calculate the depths in the two genomes in comparison, given in --qbed and
    --sbed. The synteny blocks will be layered on the genomes, and the
    multiplicity will be summarized to stderr.
    """
    from jcvi.utils.range import range_depth
    from jcvi.graphics.base import latex

    p = OptionParser(depth.__doc__)
    p.add_option("--depthfile", help="Generate file with gene and depth")
    p.add_option(
        "--histogram", default=False, action="store_true", help="Plot histograms in PDF"
    )
    p.add_option("--xmax", type="int", help="x-axis maximum to display in plot")
    p.add_option("--title", default=None, help="Title to display in plot")
    p.add_option("--quota", help="Force to use this quota, e.g. 1:1, 1:2 ...")
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (anchorfile,) = args
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
    qtag = "Genome {0} depths".format(qgenome)
    print("{}:".format(qtag), file=sys.stderr)
    dsq, details = range_depth(qranges, len(qbed))
    if depthfile:
        fw = open(depthfile, "w")
        write_details(fw, details, qbed)

    if is_self:
        return

    stag = "Genome {0} depths".format(sgenome)
    print("{}:".format(stag), file=sys.stderr)
    dss, details = range_depth(sranges, len(sbed))
    if depthfile:
        write_details(fw, details, sbed)
        fw.close()
        logging.debug("Depth written to `{0}`.".format(depthfile))

    if not opts.histogram:
        return

    from jcvi.graphics.base import plt, quickplot_ax, savefig, normalize_axes

    # Plot two histograms one for query genome, one for subject genome
    plt.figure(1, (6, 3))
    f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)

    xmax = opts.xmax or max(4, max(list(dsq.keys()) + list(dss.keys())))
    if opts.quota:
        speak, qpeak = opts.quota.split(":")
        qpeak, speak = int(qpeak), int(speak)
    else:
        qpeak = find_peak(dsq)
        speak = find_peak(dss)

    qtag = "# of {} blocks per {} gene".format(sgenome, qgenome)
    stag = "# of {} blocks per {} gene".format(qgenome, sgenome)
    quickplot_ax(
        ax1,
        dss,
        0,
        xmax,
        stag,
        ylabel="Percentage of genome",
        highlight=range(1, speak + 1),
    )
    quickplot_ax(ax2, dsq, 0, xmax, qtag, ylabel=None, highlight=range(1, qpeak + 1))

    title = opts.title or "{} vs {} syntenic depths\n{}:{} pattern".format(
        qgenome, sgenome, speak, qpeak
    )
    root = f.add_axes([0, 0, 1, 1])
    vs, pattern = latex(title).split("\n")
    root.text(0.5, 0.97, vs, ha="center", va="center", color="darkslategray")
    root.text(0.5, 0.925, pattern, ha="center", va="center", color="tomato", size=16)
    print(title, file=sys.stderr)

    normalize_axes(root)

    pf = anchorfile.rsplit(".", 1)[0] + ".depth"
    image_name = pf + ".pdf"
    savefig(image_name)


def find_peak(data, cutoff=0.9):
    """
    This will look for the point where cumulative cutoff is reached. For
    example:

    >>> find_peak({0: 27, 1: 71, 2: 1})
    1
    """
    total_length = sum(data.values())
    count_cutoff = cutoff * total_length
    cum_sum = 0
    for i, count in sorted(data.items()):
        cum_sum += count
        if cum_sum > count_cutoff:
            return i


def get_blocks(scaffold, bs, order, xdist=20, ydist=20, N=6):
    points = []
    for b in bs:
        accn = b.accn.rsplit(".", 1)[0]
        if accn not in order:
            continue
        x, xx = order[accn]
        y = (b.start + b.end) / 2
        points.append((x, y))

    # print scaffold, points
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
    p.add_option(
        "--xdist", type="int", default=20, help="xdist (in related genome) cutoff"
    )
    p.add_option(
        "--ydist", type="int", default=200000, help="ydist (in current genome) cutoff"
    )
    p.add_option("-n", type="int", default=5, help="number of markers in a block")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blastfile, bedfile = args
    order = Bed(bedfile).order
    blastbedfile = bed([blastfile])
    bbed = Bed(blastbedfile)
    for scaffold, bs in bbed.sub_beds():
        blocks = get_blocks(
            scaffold, bs, order, xdist=opts.xdist, ydist=opts.ydist, N=opts.n
        )
        sblocks = []
        for block in blocks:
            xx, yy = zip(*block)
            sblocks.append((scaffold, min(yy), max(yy)))
        iblocks = range_interleave(sblocks)
        for ib in iblocks:
            ch, start, end = ib
            print("{0}\t{1}\t{2}".format(ch, start - 1, end))


def scan(args):
    """
    %prog scan blastfile anchor_file [options]

    pull out syntenic anchors from blastfile based on single-linkage algorithm
    """
    p = OptionParser(scan.__doc__)
    p.add_option(
        "-n",
        "--min_size",
        dest="n",
        type="int",
        default=4,
        help="minimum number of anchors in a cluster",
    )
    p.add_option(
        "--intrabound",
        default=300,
        type="int",
        help="Lower bound of intra-chromosomal blocks (only for self comparison)",
    )
    p.add_option("--liftover", help="Scan BLAST file to find extra anchors")
    p.add_option(
        "--liftover_dist",
        type="int",
        help="Distance to extend from liftover. Defaults to half of --dist",
    )
    p.set_stripnames()

    blast_file, anchor_file, dist, opts = add_options(p, args, dist=20)
    qbed, sbed, qorder, sorder, is_self = check_beds(blast_file, p, opts)

    intrabound = opts.intrabound
    filtered_blast = read_blast(
        blast_file, qorder, sorder, is_self=is_self, ostrip=False
    )

    fw = open(anchor_file, "w")
    logging.debug("Chaining distance = {0}".format(dist))

    clusters = batch_scan(
        filtered_blast,
        xdist=dist,
        ydist=dist,
        N=opts.n,
        is_self=is_self,
        intrabound=intrabound,
    )
    for cluster in clusters:
        print("###", file=fw)
        for qi, si, score in cluster:
            query, subject = qbed[qi].accn, sbed[si].accn
            print("\t".join((query, subject, str(int(score)))), file=fw)

    fw.close()
    summary([anchor_file])

    lo = opts.liftover
    if not lo:
        return anchor_file

    dargs = ["--qbed=" + opts.qbed, "--sbed=" + opts.sbed]
    if not opts.strip_names:
        dargs += ["--no_strip_names"]
    liftover_dist = opts.liftover_dist or dist // 2
    dargs += ["--dist={}".format(liftover_dist)]
    newanchorfile = liftover([lo, anchor_file] + dargs)
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

    filtered_blast = read_blast(
        blast_file, qorder, sorder, is_self=is_self, ostrip=opts.strip_names
    )
    blast_to_score = dict(((b.qi, b.si), int(b.score)) for b in filtered_blast)
    accepted = dict(((b.query, b.subject), str(int(b.score))) for b in filtered_blast)

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

    logging.debug("{} new pairs found (dist={}).".format(lifted, dist))
    newanchorfile = anchor_file.rsplit(".", 1)[0] + ".lifted.anchors"
    ac.print_to_file(filename=newanchorfile, accepted=accepted)
    summary([newanchorfile])

    return newanchorfile


if __name__ == "__main__":
    main()
