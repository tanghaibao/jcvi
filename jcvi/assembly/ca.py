#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Prepare input files for Celera Assembler, dispatch based on file suffix::

*.fasta: convert-fasta-to-v2.pl
*.sff: sffToCA
*.fastq: fastqToCA
"""
import logging
import os.path as op
import sys
from pickle import dump, load

import networkx as nx
from collections import Counter
from random import choice
from Bio import SeqIO

from jcvi.formats.base import must_open
from jcvi.formats.fasta import Fasta, SeqRecord, filter, format, parse_fasta
from jcvi.formats.blast import Blast
from jcvi.utils.console import printf
from jcvi.utils.range import range_minmax
from jcvi.algorithms.graph import graph_stats, graph_local_neighborhood
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    sh,
    need_update,
    glob,
    get_abs_path,
    popen,
)


def main():

    actions = (
        ("tracedb", "convert trace archive files to frg file"),
        ("clr", "prepare vector clear range file based on BLAST to vectors"),
        ("fasta", "convert fasta to frg file"),
        ("sff", "convert 454 reads to frg file"),
        ("fastq", "convert Illumina reads to frg file"),
        ("shred", "shred contigs into pseudo-reads"),
        ("astat", "generate the coverage-rho scatter plot"),
        ("unitigs", "output uniquely extended unitigs based on best.edges"),
        ("merger", "merge reads into unitigs offline"),
        ("removecontains", "remove contained reads from gkpStore"),
        ("graph", "visualize best.edges"),
        ("prune", "prune overlap graph"),
        ("overlap", "visualize overlaps for a given fragment"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


frgTemplate = """{{FRG
act:A
acc:{fragID}
rnd:1
sta:G
lib:{libID}
pla:0
loc:0
src:
.
seq:
{seq}
.
qlt:
{qvs}
.
hps:
.
clr:{clr_beg},{clr_end}
}}"""

headerTemplate = """{{VER
ver:2
}}
{{LIB
act:A
acc:{libID}
ori:U
mea:0.000
std:0.000
src:
.
nft:17
fea:
forceBOGunitigger=1
isNotRandom=0
doNotTrustHomopolymerRuns=0
doTrim_initialNone=1
doTrim_initialMerBased=0
doTrim_initialFlowBased=0
doTrim_initialQualityBased=0
doRemoveDuplicateReads=1
doTrim_finalLargestCovered=1
doTrim_finalEvidenceBased=0
doTrim_finalBestEdge=0
doRemoveSpurReads=1
doRemoveChimericReads=1
doCheckForSubReads=0
doConsensusCorrection=0
forceShortReadFormat=0
constantInsertSize=0
.
}}"""


class OverlapLine(object):

    # See doc: http://wgs-assembler.sourceforge.net/wiki/index.php/OverlapStore
    def __init__(self, line):
        args = line.split()
        self.aid = int(args[0])
        self.bid = int(args[1])
        self.orientation = args[2]
        self.ahang = int(args[3])
        self.bhang = int(args[4])
        self.erate = float(args[5])
        self.erate_adj = float(args[6])


def add_graph_options(p):
    p.add_option("--maxerr", default=100, type="int", help="Maximum error rate")
    p.add_option(
        "--frgctg",
        default="../9-terminator/asm.posmap.frgctg",
        help="Annotate graph with contig membership",
    )


def prune(args):
    """
    %prog prune best.edges

    Prune overlap graph.
    """
    from collections import defaultdict

    p = OptionParser(prune.__doc__)
    add_graph_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bestedges,) = args
    G = read_graph(bestedges, maxerr=opts.maxerr)
    reads_to_ctgs = parse_ctgs(bestedges, opts.frgctg)
    edges = defaultdict(int)
    r = defaultdict(int)
    for a, b, d in G.edges_iter(data=True):
        ua, ub = reads_to_ctgs.get(a), reads_to_ctgs.get(b)
        nn = (ua, ub).count(None)
        if nn == 0:
            if ua == ub:
                r["Same tigs"] += 1
            else:
                r["Diff tigs"] += 1
                if ua > ub:
                    ua, ub = ub, ua
                edges[(ua, ub)] += 1
        elif nn == 1:
            r["One null"] += 1
        else:
            assert nn == 2
            r["Two nulls"] += 1

    U = nx.Graph()
    difftigs = "diff_tigs.txt"
    neighbors = defaultdict(list)
    fw = open(difftigs, "w")
    for (ua, ub), count in edges.items():
        print("\t".join((ua, ub, str(count))), file=fw)
        U.add_edge(ua, ub, weight=count)
        neighbors[ua].append((ub, count))
        neighbors[ub].append((ua, count))
    fw.close()

    print("[Unitig edge property]", file=sys.stderr)
    for k, v in r.items():
        print(": ".join((k, str(v))), file=sys.stderr)
    print("Total: {0}".format(sum(r.values())), file=sys.stderr)

    print("[Unitig degree distribution]", file=sys.stderr)
    degrees = U.degree()
    degree_counter = Counter(degrees.values())
    for degree, count in sorted(degree_counter.items()):
        print("{0}\t{1}".format(degree, count), file=sys.stderr)

    # To find associative contigs, one look for a contig that is connected and
    # only connected to another single contig - and do that recursively until no
    # more contigs can be found
    associative = {}
    for ua, ubs in neighbors.items():
        if len(ubs) == 1:  # Only one neighbor
            ub, count = ubs[0]
            if count >= 2:  # Bubble
                associative[ua] = (ub, count)
    print(
        "A total of {0} associative contigs found".format(len(associative)),
        file=sys.stderr,
    )

    # Keep only one for mutual associative
    for ua, ub in associative.items():
        if ub in associative and ua < ub:
            print(ua, "mutually associative with", ub, file=sys.stderr)
            del associative[ub]
    print(
        "A total of {0} associative contigs retained".format(len(associative)),
        file=sys.stderr,
    )

    assids = "associative.ids"
    fw = open(assids, "w")
    for ua, (ub, count) in sorted(associative.items(), key=lambda x: (x[1], x[0])):
        print("\t".join((ua, ub, str(count))), file=fw)
    fw.close()
    logging.debug("Associative contigs written to `{0}`".format(assids))


def removecontains(args):
    """
    %prog removecontains 4-unitigger/best.contains asm.gkpStore

    Remove contained reads from gkpStore. This will improve assembly contiguity
    without sacrificing accuracy, when using bogart unitigger.
    """
    p = OptionParser(removecontains.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    contains, gkpStore = args

    s = set()
    fp = open(contains)
    for row in fp:
        if row[0] == "#":
            continue
        iid = int(row.split()[0])
        s.add(iid)

    cmd = "gatekeeper -dumpfragments -lastfragiid {}".format(gkpStore)
    gkpmsg = popen(cmd).read()
    last_iid = int(gkpmsg.strip().split()[-1])

    ndeleted = 0
    editfile = "delete.edit"
    fw = open(editfile, "w")
    for iid in range(1, last_iid + 1):
        if iid in s:
            print("frg iid {0} isdeleted 1".format(iid), file=fw)
            ndeleted += 1

    fw.close()
    assert len(s) == ndeleted
    logging.debug("A total of {0} contained reads flagged as deleted.".format(ndeleted))
    print("Now you can run:", file=sys.stderr)
    print("$ gatekeeper --edit {0} {1}".format(editfile, gkpStore), file=sys.stderr)


def overlap(args):
    """
    %prog overlap best.contains iid

    Visualize overlaps for a given fragment. Must be run in 4-unitigger. All
    overlaps for iid were retrieved, excluding the ones matching best.contains.
    """
    p = OptionParser(overlap.__doc__)
    p.add_option("--maxerr", default=2, type="int", help="Maximum error rate")
    p.add_option("--canvas", default=100, type="int", help="Canvas size")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bestcontains, iid = args
    canvas = opts.canvas

    bestcontainscache = bestcontains + ".cache"
    if need_update(bestcontains, bestcontainscache):
        fp = open(bestcontains)
        fw = open(bestcontainscache, "w")
        exclude = set()
        for row in fp:
            if row[0] == "#":
                continue
            j = int(row.split()[0])
            exclude.add(j)
        dump(exclude, fw)
        fw.close()

    exclude = load(open(bestcontainscache))
    logging.debug("A total of {0} reads to exclude".format(len(exclude)))

    cmd = "overlapStore -d ../asm.ovlStore -b {0} -e {0}".format(iid)
    cmd += " -E {0}".format(opts.maxerr)
    frags = []
    for row in popen(cmd):
        r = OverlapLine(row)
        if r.bid in exclude:
            continue
        frags.append(r)

    # Also include to query fragment
    frags.append(OverlapLine("{0} {0} N 0 0 0 0".format(iid)))
    frags.sort(key=lambda x: x.ahang)

    # Determine size of the query fragment
    cmd = "gatekeeper -b {0} -e {0}".format(iid)
    cmd += " -tabular -dumpfragments ../asm.gkpStore"
    fp = popen(cmd)
    row = next(fp)
    size = int(next(fp).split()[-1])

    # Determine size of canvas
    xmin = min(x.ahang for x in frags)
    xmax = max(x.bhang for x in frags)
    xsize = -xmin + size + xmax
    ratio = xsize / canvas

    for f in frags:
        fsize = -f.ahang + size + f.bhang
        a = (f.ahang - xmin) / ratio
        b = fsize / ratio
        t = "-" * b
        if f.orientation == "N":
            t = t[:-1] + ">"
        else:
            t = "<" + t[1:]
        if f.ahang == 0 and f.bhang == 0:
            t = "[green]{}".format(t)
        c = canvas - a - b
        printf(
            "{}{}{}{} ({})".format(
                " " * a, t, " " * c, str(f.bid).rjust(10), f.erate_adj
            ),
        )


def parse_ctgs(bestedges, frgtoctg):
    cache = "frgtoctg.cache"
    if need_update(frgtoctg, cache):
        reads_to_ctgs = {}
        frgtodeg = frgtoctg.replace(".frgctg", ".frgdeg")
        iidtouid = frgtoctg.replace(".posmap.frgctg", ".iidtouid")
        fp = open(iidtouid)
        frgstore = {}
        for row in fp:
            tag, iid, uid = row.split()
            if tag == "FRG":
                frgstore[uid] = int(iid)

        for pf, f in zip(("ctg", "deg"), (frgtoctg, frgtodeg)):
            fp = open(f)
            logging.debug("Parse posmap file `{0}`".format(f))
            for row in fp:
                frg, ctg = row.split()[:2]
                frg = frgstore[frg]
                reads_to_ctgs[frg] = pf + ctg
            logging.debug("Loaded mapping: {0}".format(len(reads_to_ctgs)))

        fw = open(cache, "w")
        dump(reads_to_ctgs, fw)
        fw.close()
        logging.debug("Contig mapping written to `{0}`".format(cache))

    reads_to_ctgs = load(open(cache))
    logging.debug("Contig mapping loaded from `{0}`".format(cache))
    return reads_to_ctgs


def read_graph(bestedges, maxerr=100, directed=False):
    logging.debug("Max error = {0}%".format(maxerr))
    tag = "dir." if directed else ""
    bestgraph = bestedges.split(".")[0] + ".err{0}.{1}graph".format(maxerr, tag)
    if need_update(bestedges, bestgraph):
        G = {} if directed else nx.Graph()
        fp = open(bestedges)
        best_store = {}
        for row in fp:
            if row[0] == "#":
                continue
            id1, lib_id, best5, o5, best3, o3, j1, j2 = row.split()
            id1, best5, best3 = int(id1), int(best5), int(best3)
            j1, j2 = float(j1), float(j2)
            if j1 <= maxerr or j2 <= maxerr:
                if not directed:
                    G.add_node(id1)
                id1p5, id1p3 = "{0}-5'".format(id1), "{0}-3'".format(id1)
                best5o5 = "{0}-{1}".format(best5, o5)
                best3o3 = "{0}-{1}".format(best3, o3)
                best_store[id1p5] = best5o5
                best_store[id1p3] = best3o3
            if best5 and j1 <= maxerr:
                if directed:
                    G[id1p5] = best5o5
                else:
                    G.add_edge(best5, id1, weight=10)
            if best3 and j2 <= maxerr:
                if directed:
                    G[id1p3] = best3o3
                else:
                    G.add_edge(id1, best3, weight=10)

        # Annotate edge weight for mutual best link, note that edge weights are
        # (11) set close to 10, to minimize impact to layout (Yifan Hu's
        # multilevel)
        nmutuals = 0
        for k, v in best_store.items():
            if best_store.get(v) == k and k < v:
                k, v = int(k.split("-")[0]), int(v.split("-")[0])
                G[k][v]["weight"] = 11
                nmutuals += 1
        logging.debug("Mutual best edges: {0}".format(nmutuals))

        if directed:
            fw = open(bestgraph, "w")
            dump(G, fw)
            fw.close()
        else:
            nx.write_gpickle(G, bestgraph)
        logging.debug("Graph pickled to `{0}`".format(bestgraph))

        # Compute node degree histogram and save in (degree, counts) tab file
        degrees = G.degree()
        degree_counter = Counter(degrees.values())
        degreesfile = "degrees.txt"
        fw = open(degreesfile, "w")
        for degree, count in sorted(degree_counter.items()):
            print("{0}\t{1}".format(degree, count), file=fw)
        fw.close()
        logging.debug("Node degree distribution saved to `{0}`".format(degreesfile))

        # Save high degree (top 1%) nodes in save in (node, degree) tab file
        percentile = sorted(degrees.values(), reverse=True)[len(degrees) / 1000]
        logging.debug("Top 0.1% has degree of at least {0}".format(percentile))
        hubs = [(k, v) for k, v in degrees.items() if v >= percentile]
        hubs.sort(key=lambda x: x[1], reverse=True)  # degress descending
        hubsfile = "hubs.txt"
        fw = open(hubsfile, "w")
        for node, degree in hubs:
            print("{0}\t{1}".format(node, degree), file=fw)
        fw.close()
        logging.debug("Hubs saved to `{0}`".format(hubsfile))

    logging.debug("Read graph from `{0}`".format(bestgraph))
    if directed:
        G = load(open(bestgraph))
    else:
        G = nx.read_gpickle(bestgraph)
        graph_stats(G)
    return G


def merger(args):
    """
    %prog merger layout gkpStore contigs.fasta

    Merge reads into one contig.
    """
    p = OptionParser(merger.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    layout, gkpstore, contigs = args
    fp = open(layout)
    pf = "0"
    iidfile = pf + ".iids"
    for i, row in enumerate(fp):
        logging.debug("Read unitig {0}".format(i))
        fw = open(iidfile, "w")
        layout = row.split("|")
        print("\n".join(layout), file=fw)
        fw.close()
        cmd = "gatekeeper -iid {0}.iids -dumpfasta {0} {1}".format(pf, gkpstore)
        sh(cmd)

        fastafile = "{0}.fasta".format(pf)
        newfastafile = "{0}.new.fasta".format(pf)
        format(
            [
                fastafile,
                newfastafile,
                "--sequential=replace",
                "--sequentialoffset=1",
                "--nodesc",
            ]
        )
        fasta([newfastafile])

        sh("rm -rf {0}".format(pf))
        cmd = "runCA {0}.frg -p {0} -d {0} consensus=pbutgcns".format(pf)
        cmd += " unitigger=bogart doFragmentCorrection=0 doUnitigSplitting=0"
        sh(cmd)
        outdir = "{0}/9-terminator".format(pf)

        cmd = "cat {0}/{1}.ctg.fasta {0}/{1}.deg.fasta {0}/{1}.singleton.fasta".format(
            outdir, pf
        )
        sh(cmd, outfile=contigs, append=True)


def unitigs(args):
    """
    %prog unitigs best.edges

    Reads Celera Assembler's "best.edges" and extract all unitigs.
    """
    p = OptionParser(unitigs.__doc__)
    p.add_option("--maxerr", default=2, type="int", help="Maximum error rate")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bestedges,) = args
    G = read_graph(bestedges, maxerr=opts.maxerr, directed=True)
    H = nx.Graph()
    intconv = lambda x: int(x.split("-")[0])
    for k, v in G.iteritems():
        if k == G.get(v, None):
            H.add_edge(intconv(k), intconv(v))

    nunitigs = nreads = 0
    for h in nx.connected_component_subgraphs(H, copy=False):
        st = [x for x in h if h.degree(x) == 1]
        if len(st) != 2:
            continue
        src, target = st
        path = list(nx.all_simple_paths(h, src, target))
        assert len(path) == 1
        (path,) = path
        print("|".join(str(x) for x in path))
        nunitigs += 1
        nreads += len(path)
    logging.debug(
        "A total of {0} unitigs built from {1} reads.".format(nunitigs, nreads)
    )


def graph(args):
    """
    %prog graph best.edges

    Convert Celera Assembler's "best.edges" to a GEXF which can be used to
    feed into Gephi to check the topology of the best overlapping graph. Mutual
    best edges are represented as thicker edges.

    Reference:
    https://github.com/PacificBiosciences/Bioinformatics-Training/blob/master/scripts/CeleraToGephi.py
    """
    p = OptionParser(graph.__doc__)
    p.add_option(
        "--query",
        default=-1,
        type="int",
        help="Search from node, -1 to select random node, 0 to disable",
    )
    p.add_option("--contig", help="Search from contigs, use comma to separate")
    p.add_option(
        "--largest", default=0, type="int", help="Only show largest components"
    )
    p.add_option("--maxsize", default=500, type="int", help="Max graph size")
    p.add_option(
        "--nomutualbest",
        default=False,
        action="store_true",
        help="Do not plot mutual best edges as heavy",
    )
    add_graph_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bestedges,) = args
    query = opts.query
    contig = opts.contig
    largest = opts.largest
    frgctg = opts.frgctg
    edgeweight = not opts.nomutualbest
    G = read_graph(bestedges, maxerr=opts.maxerr)

    if largest:
        H = list(nx.connected_component_subgraphs(G))
        c = min(len(H), largest)
        logging.debug("{0} components found, {1} retained".format(len(H), c))
        G = nx.Graph()
        for x in H[:c]:
            G.add_edges_from(x.edges())

    if query:
        if query == -1:
            query = choice(G.nodes())
        reads_to_ctgs = parse_ctgs(bestedges, frgctg)
        if contig:
            contigs = set(contig.split(","))
            core = [k for k, v in reads_to_ctgs.items() if v in contigs]
        else:
            ctg = reads_to_ctgs.get(query)
            core = [k for k, v in reads_to_ctgs.items() if v == ctg]
            logging.debug(
                "Reads ({0}) extended from the same contig {1}".format(len(core), ctg)
            )

        # Extract a local neighborhood
        SG = nx.Graph()
        H = graph_local_neighborhood(G, query=core, maxsize=opts.maxsize)
        SG.add_edges_from(H.edges(data=edgeweight))
        G = SG

        seen = []
        for n, attrib in G.nodes_iter(data=True):
            contig = reads_to_ctgs.get(n, "na")
            attrib["label"] = contig
            seen.append(contig)
        c = Counter(seen)
        cc = ["{0}({1})".format(k, v) for k, v in c.most_common()]
        print("Contigs: {0}".format(" ".join(cc)), file=sys.stderr)

    gexf = "best"
    if query >= 0:
        gexf += ".{0}".format(query)
    gexf += ".gexf"
    nx.write_gexf(G, gexf)
    logging.debug(
        "Graph written to `{0}` (|V|={1}, |E|={2})".format(gexf, len(G), G.size())
    )


def astat(args):
    """
    %prog astat coverage.log

    Create coverage-rho scatter plot.
    """
    p = OptionParser(astat.__doc__)
    p.add_option("--cutoff", default=1000, type="int", help="Length cutoff")
    p.add_option("--genome", default="", help="Genome name")
    p.add_option(
        "--arrDist",
        default=False,
        action="store_true",
        help="Use arrDist instead",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (covfile,) = args
    cutoff = opts.cutoff
    genome = opts.genome
    plot_arrDist = opts.arrDist

    suffix = ".{0}".format(cutoff)
    small_covfile = covfile + suffix
    update_covfile = need_update(covfile, small_covfile)
    if update_covfile:
        fw = open(small_covfile, "w")
    else:
        logging.debug("Found `{0}`, will use this one".format(small_covfile))
        covfile = small_covfile

    fp = open(covfile)
    header = next(fp)
    if update_covfile:
        fw.write(header)

    data = []
    msg = "{0} tigs scanned ..."
    for row in fp:
        tigID, rho, covStat, arrDist = row.split()
        tigID = int(tigID)
        if tigID % 1000000 == 0:
            sys.stderr.write(msg.format(tigID) + "\r")

        rho, covStat, arrDist = [float(x) for x in (rho, covStat, arrDist)]
        if rho < cutoff:
            continue

        if update_covfile:
            fw.write(row)
        data.append((tigID, rho, covStat, arrDist))

    print(msg.format(tigID), file=sys.stderr)

    from jcvi.graphics.base import plt, savefig

    logging.debug("Plotting {0} data points.".format(len(data)))
    tigID, rho, covStat, arrDist = zip(*data)

    y = arrDist if plot_arrDist else covStat
    ytag = "arrDist" if plot_arrDist else "covStat"

    fig = plt.figure(1, (7, 7))
    ax = fig.add_axes([0.12, 0.1, 0.8, 0.8])
    ax.plot(rho, y, ".", color="lightslategrey")

    xtag = "rho"
    info = (genome, xtag, ytag)
    title = "{0} {1} vs. {2}".format(*info)
    ax.set_title(title)
    ax.set_xlabel(xtag)
    ax.set_ylabel(ytag)

    if plot_arrDist:
        ax.set_yscale("log")

    imagename = "{0}.png".format(".".join(info))
    savefig(imagename, dpi=150)


def emitFragment(fw, fragID, libID, shredded_seq, clr=None, qvchar="l", fasta=False):
    """
    Print out the shredded sequence.
    """
    if fasta:
        s = SeqRecord(shredded_seq, id=fragID, description="")
        SeqIO.write([s], fw, "fasta")
        return

    seq = str(shredded_seq)
    slen = len(seq)
    qvs = qvchar * slen  # shredded reads have default low qv

    if clr is None:
        clr_beg, clr_end = 0, slen
    else:
        clr_beg, clr_end = clr

    print(
        frgTemplate.format(
            fragID=fragID,
            libID=libID,
            seq=seq,
            qvs=qvs,
            clr_beg=clr_beg,
            clr_end=clr_end,
        ),
        file=fw,
    )


def shred(args):
    """
    %prog shred fastafile

    Similar to the method of `shredContig` in runCA script. The contigs are
    shredded into pseudo-reads with certain length and depth.
    """
    p = OptionParser(shred.__doc__)
    p.set_depth(depth=2)
    p.add_option(
        "--readlen",
        default=1000,
        type="int",
        help="Desired length of the reads",
    )
    p.add_option(
        "--minctglen",
        default=0,
        type="int",
        help="Ignore contig sequence less than",
    )
    p.add_option(
        "--shift",
        default=50,
        type="int",
        help="Overlap between reads must be at least",
    )
    p.add_option(
        "--fasta",
        default=False,
        action="store_true",
        help="Output shredded reads as FASTA sequences",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    libID = fastafile.split(".")[0]
    depth = opts.depth
    readlen = opts.readlen
    shift = opts.shift

    outfile = libID + ".depth{0}".format(depth)
    if opts.fasta:
        outfile += ".fasta"
    else:
        outfile += ".frg"
    f = Fasta(fastafile, lazy=True)

    fw = must_open(outfile, "w", checkexists=True)
    if not opts.fasta:
        print(headerTemplate.format(libID=libID), file=fw)

    """
    Taken from runCA:

                    |*********|
                    |###################|
    |--------------------------------------------------|
     ---------------1---------------
               ---------------2---------------
                         ---------------3---------------
    *** - center_increments
    ### - center_range_width
    """
    for ctgID, (name, rec) in enumerate(f.iteritems_ordered()):
        seq = rec.seq
        seqlen = len(seq)
        if seqlen < opts.minctglen:
            continue

        shredlen = min(seqlen - shift, readlen)
        numreads = max(seqlen * depth / shredlen, 1)
        center_range_width = seqlen - shredlen

        ranges = []
        if depth == 1:
            if seqlen < readlen:
                ranges.append((0, seqlen))
            else:
                for begin in range(0, seqlen, readlen - shift):
                    end = min(seqlen, begin + readlen)
                    ranges.append((begin, end))
        else:
            if numreads == 1:
                ranges.append((0, shredlen))
            else:
                prev_begin = -1
                center_increments = center_range_width * 1.0 / (numreads - 1)
                for i in range(numreads):
                    begin = center_increments * i
                    end = begin + shredlen
                    begin, end = int(begin), int(end)

                    if begin == prev_begin:
                        continue

                    ranges.append((begin, end))
                    prev_begin = begin

        for shredID, (begin, end) in enumerate(ranges):
            shredded_seq = seq[begin:end]
            fragID = "{0}.{1}.frag{2}.{3}-{4}".format(libID, ctgID, shredID, begin, end)
            emitFragment(fw, fragID, libID, shredded_seq, fasta=opts.fasta)

    fw.close()
    logging.debug("Shredded reads are written to `{0}`.".format(outfile))
    return outfile


def tracedb(args):
    """
    %prog tracedb <xml|lib|frg>

    Run `tracedb-to-frg.pl` within current folder.
    """
    p = OptionParser(tracedb.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (action,) = args
    assert action in ("xml", "lib", "frg")

    CMD = "tracedb-to-frg.pl"
    xmls = glob("xml*")

    if action == "xml":
        for xml in xmls:
            cmd = CMD + " -xml {0}".format(xml)
            sh(cmd, outfile="/dev/null", errfile="/dev/null", background=True)

    elif action == "lib":
        cmd = CMD + " -lib {0}".format(" ".join(xmls))
        sh(cmd)

    elif action == "frg":
        for xml in xmls:
            cmd = CMD + " -frg {0}".format(xml)
            sh(cmd, background=True)


def make_matepairs(fastafile):
    """
    Assumes the mates are adjacent sequence records
    """
    assert op.exists(fastafile)

    matefile = fastafile.rsplit(".", 1)[0] + ".mates"
    if op.exists(matefile):
        logging.debug("matepairs file `{0}` found".format(matefile))
    else:
        logging.debug("parsing matepairs from `{0}`".format(fastafile))
        matefw = open(matefile, "w")
        it = SeqIO.parse(fastafile, "fasta")
        for fwd, rev in zip(it, it):
            print("{0}\t{1}".format(fwd.id, rev.id), file=matefw)

        matefw.close()

    return matefile


get_mean_sv = lambda size: (size, size / 5)


def split_fastafile(fastafile, maxreadlen=32000):
    pf = fastafile.split(".")[0]
    smallfastafile = pf + "-small.fasta"
    bigfastafile = pf + "-big.fasta"
    shredfastafile = pf + "-big.depth1.fasta"

    if need_update(fastafile, (smallfastafile, shredfastafile)):
        filter([fastafile, str(maxreadlen), "--less", "-o", smallfastafile])
        filter([fastafile, str(maxreadlen), "-o", bigfastafile])
        shred(
            [
                "--depth=1",
                "--shift={0}".format(maxreadlen / 100),
                "--readlen={0}".format(maxreadlen),
                "--fasta",
                bigfastafile,
            ]
        )

    return smallfastafile, shredfastafile


def fasta(args):
    """
    %prog fasta fastafile

    Convert reads formatted as FASTA file, and convert to CA frg file. If .qual
    file is found, then use it, otherwise just make a fake qual file. Mates are
    assumed as adjacent sequence records (i.e. /1, /2, /1, /2 ...) unless a
    matefile is given.
    """
    from jcvi.formats.fasta import clean, make_qual

    p = OptionParser(fasta.__doc__)
    p.add_option(
        "--clean",
        default=False,
        action="store_true",
        help="Clean up irregular chars in seq",
    )
    p.add_option("--matefile", help="Matepairs file")
    p.add_option(
        "--maxreadlen", default=262143, type="int", help="Maximum read length allowed"
    )
    p.add_option(
        "--minreadlen", default=1000, type="int", help="Minimum read length allowed"
    )
    p.add_option(
        "--sequential",
        default=False,
        action="store_true",
        help="Overwrite read name (e.g. long Pacbio name)",
    )
    p.set_size()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    maxreadlen = opts.maxreadlen
    minreadlen = opts.minreadlen
    if maxreadlen > 0:
        split = False
        f = Fasta(fastafile, lazy=True)
        for id, size in f.itersizes_ordered():
            if size > maxreadlen:
                logging.debug(
                    "Sequence {0} (size={1}) longer than max read len {2}".format(
                        id, size, maxreadlen
                    )
                )
                split = True
                break

        if split:
            for f in split_fastafile(fastafile, maxreadlen=maxreadlen):
                fasta([f, "--maxreadlen=0"])
            return

    plate = op.basename(fastafile).split(".")[0]

    mated = opts.size != 0
    mean, sv = get_mean_sv(opts.size)

    if mated:
        libname = "Sanger{0}Kb-".format(opts.size / 1000) + plate
    else:
        libname = plate

    frgfile = libname + ".frg"

    if opts.clean:
        cleanfasta = fastafile.rsplit(".", 1)[0] + ".clean.fasta"
        if need_update(fastafile, cleanfasta):
            clean([fastafile, "--canonical", "-o", cleanfasta])
        fastafile = cleanfasta

    if mated:
        qualfile = make_qual(fastafile, score=21)
        if opts.matefile:
            matefile = opts.matefile
            assert op.exists(matefile)
        else:
            matefile = make_matepairs(fastafile)

        cmd = "convert-fasta-to-v2.pl"
        cmd += " -l {0} -s {1} -q {2} ".format(libname, fastafile, qualfile)
        if mated:
            cmd += "-mean {0} -stddev {1} -m {2} ".format(mean, sv, matefile)

        sh(cmd, outfile=frgfile)
        return

    fw = must_open(frgfile, "w")
    print(headerTemplate.format(libID=libname), file=fw)

    sequential = opts.sequential
    i = j = 0
    for fragID, seq in parse_fasta(fastafile):
        if len(seq) < minreadlen:
            j += 1
            continue
        i += 1
        if sequential:
            fragID = libname + str(100000000 + i)
        emitFragment(fw, fragID, libname, seq)
    fw.close()

    logging.debug(
        "A total of {0} fragments written to `{1}` ({2} discarded).".format(
            i, frgfile, j
        )
    )


def sff(args):
    """
    %prog sff sffiles

    Convert reads formatted as 454 SFF file, and convert to CA frg file.
    Turn --nodedup on if another deduplication mechanism is used (e.g.
    CD-HIT-454). See assembly.sff.deduplicate().
    """
    p = OptionParser(sff.__doc__)
    p.add_option(
        "--prefix", dest="prefix", default=None, help="Output frg filename prefix"
    )
    p.add_option(
        "--nodedup",
        default=False,
        action="store_true",
        help="Do not remove duplicates",
    )
    p.set_size()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())

    sffiles = args
    plates = [x.split(".")[0].split("_")[-1] for x in sffiles]

    mated = opts.size != 0
    mean, sv = get_mean_sv(opts.size)

    if len(plates) > 1:
        plate = plates[0][:-1] + "X"
    else:
        plate = "_".join(plates)

    if mated:
        libname = "Titan{0}Kb-".format(opts.size / 1000) + plate
    else:
        libname = "TitanFrags-" + plate

    if opts.prefix:
        libname = opts.prefix

    cmd = "sffToCA"
    cmd += " -libraryname {0} -output {0} ".format(libname)
    cmd += " -clear 454 -trim chop "
    if mated:
        cmd += " -linker titanium -insertsize {0} {1} ".format(mean, sv)
    if opts.nodedup:
        cmd += " -nodedup "

    cmd += " ".join(sffiles)

    sh(cmd)


def fastq(args):
    """
    %prog fastq fastqfile

    Convert reads formatted as FASTQ file, and convert to CA frg file.
    """
    from jcvi.formats.fastq import guessoffset

    p = OptionParser(fastq.__doc__)
    p.add_option(
        "--outtie",
        dest="outtie",
        default=False,
        action="store_true",
        help="Are these outie reads?",
    )
    p.set_phred()
    p.set_size()

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(p.print_help())

    fastqfiles = [get_abs_path(x) for x in args]
    size = opts.size
    outtie = opts.outtie
    if size > 1000 and (not outtie):
        logging.debug("[warn] long insert size {0} but not outtie".format(size))

    mated = size != 0
    libname = op.basename(args[0]).split(".")[0]
    libname = libname.replace("_1_sequence", "")

    frgfile = libname + ".frg"
    mean, sv = get_mean_sv(opts.size)

    cmd = "fastqToCA"
    cmd += " -libraryname {0} ".format(libname)
    fastqs = " ".join("-reads {0}".format(x) for x in fastqfiles)
    if mated:
        assert len(args) in (1, 2), "you need one or two fastq files for mated library"
        fastqs = "-mates {0}".format(",".join(fastqfiles))
        cmd += "-insertsize {0} {1} ".format(mean, sv)
    cmd += fastqs

    offset = int(opts.phred) if opts.phred else guessoffset([fastqfiles[0]])
    illumina = offset == 64
    if illumina:
        cmd += " -type illumina"
    if outtie:
        cmd += " -outtie"

    sh(cmd, outfile=frgfile)


def clr(args):
    """
    %prog blastfile fastafiles

    Calculate the vector clear range file based BLAST to the vectors.
    """
    p = OptionParser(clr.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    blastfile = args[0]
    fastafiles = args[1:]

    sizes = {}
    for fa in fastafiles:
        f = Fasta(fa)
        sizes.update(f.itersizes())

    b = Blast(blastfile)
    for query, hits in b.iter_hits():

        qsize = sizes[query]
        vectors = list((x.qstart, x.qstop) for x in hits)
        vmin, vmax = range_minmax(vectors)

        left_size = vmin - 1
        right_size = qsize - vmax

        if left_size > right_size:
            clr_start, clr_end = 0, vmin
        else:
            clr_start, clr_end = vmax, qsize

        print("\t".join(str(x) for x in (query, clr_start, clr_end)))
        del sizes[query]

    for q, size in sorted(sizes.items()):
        print("\t".join(str(x) for x in (q, 0, size)))


if __name__ == "__main__":
    main()
