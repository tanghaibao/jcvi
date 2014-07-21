#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Syntenic path assembly.
"""

import sys
import logging

from string import maketrans
from itertools import groupby, combinations

from jcvi.formats.blast import BlastSlow
from jcvi.formats.sizes import Sizes
from jcvi.utils.iter import pairwise
from jcvi.algorithms.graph import BiGraph, BiEdge
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ("bed", "convert ANCHORS file to BED format"),
        ('fromblast', 'Generate path from BLAST file'),
        ('happy', 'Make graph from happy mapping data'),
        ('partition', 'Make individual graphs partitioned by happy mapping'),
        ('merge', 'Merge multiple graphs together and visualize'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed anchorsfile

    Convert ANCHORS file to BED format.
    """
    from collections import defaultdict
    from jcvi.compara.synteny import AnchorFile, check_beds
    from jcvi.formats.bed import Bed, BedLine
    from jcvi.formats.base import get_number

    p = OptionParser(bed.__doc__)
    p.add_option("--switch", default=False, action="store_true",
                 help="Switch reference and aligned map elements")
    p.add_option("--scale", type="float",
                 help="Scale the aligned map distance by factor")
    p.set_beds()
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    anchorsfile, = args
    switch = opts.switch
    scale = opts.scale
    ac = AnchorFile(anchorsfile)
    pairs = defaultdict(list)
    for a, b, block_id in ac.iter_pairs():
        pairs[a].append(b)

    qbed, sbed, qorder, sorder, is_self = check_beds(anchorsfile, p, opts)
    bd = Bed()
    for q in qbed:
        qseqid, qstart, qend, qaccn = q.seqid, q.start, q.end, q.accn
        if qaccn not in pairs:
            continue
        for s in pairs[qaccn]:
            si, s = sorder[s]
            sseqid, sstart, send, saccn = s.seqid, s.start, s.end, s.accn
        if switch:
            qseqid, sseqid = sseqid, qseqid
            qstart, sstart = sstart, qstart
            qend, send = send, qend
            qaccn, saccn = saccn, qaccn
        if scale:
            sstart /= scale
        bedline = "\t".join(str(x) for x in (qseqid, qstart - 1, qend,
                            "{0}:{1}".format(get_number(sseqid), sstart)))
        bd.append(BedLine(bedline))

    bd.print_to_file(filename=opts.outfile, sorted=True)


def happy_nodes(row, prefix=None):
    row = row.translate(None, "[](){}+-")
    scfs = [x.strip() for x in row.split(":")]
    if prefix:
        scfs = [prefix + x for x in scfs]
    return scfs


def happy_edges(row, prefix=None):
    """
    Convert a row in HAPPY file and yield edges.
    """
    trans = maketrans("[](){}", "      ")
    row = row.strip().strip("+")
    row = row.translate(trans)
    scfs = [x.strip("+") for x in row.split(":")]
    for a, b in pairwise(scfs):
        oa = '<' if a.strip()[0] == '-' else '>'
        ob = '<' if b.strip()[0] == '-' else '>'

        is_uncertain = a[-1] == ' ' or b[0] == ' '

        a = a.strip().strip('-')
        b = b.strip().strip('-')

        if prefix:
            a = prefix + a
            b = prefix + b

        e = BiEdge(a, b, oa, ob)
        yield e, is_uncertain


def partition(args):
    """
    %prog partition happy.txt synteny.graph

    Select edges from another graph and merge it with the certain edges built
    from the HAPPY mapping data.
    """
    allowed_format = ("png", "ps")
    p = OptionParser(partition.__doc__)
    p.add_option("--prefix", help="Add prefix to the name [default: %default]")
    p.add_option("--namestart", default=0, type="int",
                 help="Use a shorter name, starting index [default: %default]")
    p.add_option("--format", default="png", choices=allowed_format,
            help="Generate image of format [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    happyfile, graphfile = args
    bg = BiGraph()
    bg.read(graphfile, color="red")
    prefix = opts.prefix
    fp = open(happyfile)
    for i, row in enumerate(fp):
        nns = happy_nodes(row, prefix=prefix)
        nodes = set(nns)
        edges = happy_edges(row, prefix=prefix)

        small_graph = BiGraph()
        for e, is_uncertain in edges:
            if is_uncertain:
                e.color = "gray"
            small_graph.add_edge(e)

        for (u, v), e in bg.edges.items():
            # Grab edge if both vertices are on the same line
            if u in nodes and v in nodes:
                uv = (str(u), str(v))
                if uv in small_graph.edges:
                    e = small_graph.edges[uv]
                    e.color = "blue"  # supported by both evidences
                else:
                    small_graph.add_edge(e)

        print >> sys.stderr, small_graph

        pngfile = "A{0:02d}.{1}".format(i + 1, opts.format)
        telomeres = (nns[0], nns[-1])
        small_graph.draw(pngfile, namestart=opts.namestart,
                         nodehighlight=telomeres, dpi=72)

    legend = ["Edge colors:"]
    legend.append("[BLUE] Experimental + Synteny")
    legend.append("[BLACK] Experimental certain")
    legend.append("[GRAY] Experimental uncertain")
    legend.append("[RED] Synteny only")
    legend.append("Rectangle nodes are telomeres.")
    print >> sys.stderr, "\n".join(legend)


def merge(args):
    """
    %prog merge graphs

    Merge multiple graphs together and visualize.
    """
    p = OptionParser(merge.__doc__)
    p.add_option("--colorlist", default="black,red,pink,blue,green",
                 help="The color palette [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    colorlist = opts.colorlist.split(",")
    assert len(colorlist) >= len(args), "Need more colors in --colorlist"

    g = BiGraph()
    for a, c in zip(args, colorlist):
        g.read(a, color=c)

    g.draw("merged.png")


def happy(args):
    """
    %prog happy happy.txt

    Make bi-directed graph from HAPPY mapping data. JCVI encodes uncertainties
    in the order of the contigs / scaffolds.

    : separates scaffolds
    + means telomere (though the telomere repeats may not show because the
    telomere-adjacent sequence is missing)
    - means that the scaffold is in reverse orientation to that shown in the 2003
    TIGR scaffolds.

    Ambiguities are represented as follows, using Paul Dear.s description:
    [ ] means undetermined orientation. error quite possible (70% confidence?)
    ( ) means uncertain orientation. small chance of error (90% confidence?)
    { } means uncertain order.

    Example:
    +-8254707:8254647:-8254690:{[8254694]:[8254713]:[8254531]:[8254797]}:8254802:8254788+
    """
    p = OptionParser(happy.__doc__)
    p.add_option("--prefix", help="Add prefix to the name [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    happyfile, = args

    certain = "certain.graph"
    uncertain = "uncertain.graph"
    fw1 = open(certain, "w")
    fw2 = open(uncertain, "w")

    fp = open(happyfile)
    for row in fp:
        for e, is_uncertain in happy_edges(row, prefix=opts.prefix):
            fw = fw2 if is_uncertain else fw1
            print >> fw, e

    logging.debug("Edges written to `{0}`".format(",".join((certain, uncertain))))


def fromblast(args):
    """
    %prog fromblast blastfile subject.fasta

    Generate path from BLAST file. If multiple subjects map to the same query,
    an edge is constructed between them (with the link provided by the query).

    The BLAST file MUST be filtered, chained, supermapped.
    """
    from jcvi.formats.blast import sort
    from jcvi.utils.range import range_distance

    p = OptionParser(fromblast.__doc__)
    p.add_option("--clique", default=False, action="store_true",
                 help="Populate clique instead of linear path [default: %default]")
    p.add_option("--maxdist", default=100000, type="int",
                 help="Create edge within certain distance [default: %default]")
    p.set_verbose(help="Print verbose reports to stdout")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    blastfile, subjectfasta = args
    clique = opts.clique
    maxdist = opts.maxdist
    sort([blastfile, "--query"])
    blast = BlastSlow(blastfile, sorted=True)
    g = BiGraph()
    for query, blines in groupby(blast, key=lambda x: x.query):
        blines = list(blines)
        iterator = combinations(blines, 2) if clique else pairwise(blines)
        for a, b in iterator:
            asub, bsub = a.subject, b.subject
            if asub == bsub:
                continue

            arange = (a.query, a.qstart, a.qstop, "+")
            brange = (b.query, b.qstart, b.qstop, "+")
            dist, oo = range_distance(arange, brange, distmode="ee")
            if dist > maxdist:
                continue

            atag = ">" if a.orientation == "+" else "<"
            btag = ">" if b.orientation == "+" else "<"
            g.add_edge(BiEdge(asub, bsub, atag, btag))

    graph_to_agp(g, blastfile, subjectfasta, verbose=opts.verbose)


def graph_to_agp(g, blastfile, subjectfasta, verbose=False):

    from jcvi.formats.agp import order_to_agp

    logging.debug(str(g))
    g.write("graph.txt")
    #g.draw("graph.pdf")

    paths = []
    for path in g.iter_paths():
        m, oo = g.path(path)
        if len(oo) == 1:  # Singleton path
            continue
        paths.append(oo)
        if verbose:
            print m
            print oo

    npaths = len(paths)
    ntigs = sum(len(x) for x in paths)
    logging.debug("Graph decomposed to {0} paths with {1} components.".\
                  format(npaths, ntigs))

    agpfile = blastfile + ".agp"
    sizes = Sizes(subjectfasta)
    fwagp = open(agpfile, "w")
    scaffolded = set()
    for i, oo in enumerate(paths):
        ctgorder = [(str(ctg), ("+" if strand else "-")) \
                     for ctg, strand in oo]
        scaffolded |= set(ctg for ctg, strand in ctgorder)
        object = "pmol_{0:04d}".format(i)
        order_to_agp(object, ctgorder, sizes.mapping, fwagp)

    # Get the singletons as well
    nsingletons = 0
    for ctg, size in sizes.iter_sizes():
        if ctg in scaffolded:
            continue

        ctgorder = [(ctg, "+")]
        object = ctg
        order_to_agp(object, ctgorder, sizes.mapping, fwagp)
        nsingletons += 1
    logging.debug("Written {0} unscaffolded singletons.".format(nsingletons))

    fwagp.close()
    logging.debug("AGP file written to `{0}`.".format(agpfile))


if __name__ == '__main__':
    main()
