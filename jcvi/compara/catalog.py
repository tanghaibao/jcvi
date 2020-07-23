#!/usr/bin/env python
# -*- coding: UTF-8 -*-

from __future__ import print_function
import os.path as op
import sys
import logging
import string

from collections import defaultdict
from itertools import product, combinations

from jcvi.formats.blast import BlastLine
from jcvi.formats.fasta import Fasta
from jcvi.formats.bed import Bed
from jcvi.formats.base import must_open, BaseFile
from jcvi.utils.grouper import Grouper
from jcvi.utils.cbook import gene_name
from jcvi.compara.synteny import AnchorFile, check_beds
from jcvi.apps.base import (
    OptionParser,
    OptionGroup,
    glob,
    ActionDispatcher,
    need_update,
    sh,
    mkdir,
)


class OMGFile(BaseFile):
    def __init__(self, filename):
        super(OMGFile, self).__init__(filename)
        fp = open(filename)
        inblock = False
        components = []
        component = []
        for row in fp:
            if inblock:
                atoms = row.split()
                natoms = len(atoms)
                assert natoms in (0, 7)
                if natoms:
                    gene, taxa = atoms[0], atoms[5]
                    component.append((gene, taxa))
                else:
                    inblock = False
                    components.append(tuple(component))

            if row.strip().startswith("---"):
                inblock = True
                component = []

        if inblock:
            components.append(tuple(component))
        self.components = components

    def best(self):
        bb = set()
        for component in self.components:
            size = len(component)
            if size > 1:
                bb.add(component)
        return bb


def main():

    actions = (
        ("tandem", "identify tandem gene groups within certain distance"),
        ("ortholog", "run a combined synteny and RBH pipeline to call orthologs"),
        ("group", "cluster the anchors into ortho-groups"),
        ("omgprepare", "prepare weights file to run Sankoff OMG algorithm"),
        ("omg", "generate a series of Sankoff OMG algorithm inputs"),
        ("omgparse", "parse the OMG outputs to get gene lists"),
        ("enrich", "enrich OMG output by pulling genes missed by OMG"),
        ("layout", "layout the gene lists"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_weights(weightsfiles=None):
    if weightsfiles is None:
        weightsfiles = glob("*.weights")

    weights = defaultdict(list)
    for row in must_open(weightsfiles):
        a, b, c = row.split()
        weights[a].append((a, b, c))
    return weights


def get_edges(weightsfiles=None):
    if weightsfiles is None:
        weightsfiles = glob("*.weights")

    edges = {}
    for row in must_open(weightsfiles):
        a, b, c = row.split()
        c = int(c)
        edges[(a, b)] = c
        edges[(b, a)] = c
    return edges


def get_info():
    infofiles = glob("*.info")
    info = {}
    for row in must_open(infofiles):
        a = row.split()[0]
        info[a] = row.rstrip()
    return info


def enrich(args):
    """
    %prog enrich omgfile groups ntaxa > enriched.omg

    Enrich OMG output by pulling genes misses by OMG.
    """
    p = OptionParser(enrich.__doc__)
    p.add_option(
        "--ghost",
        default=False,
        action="store_true",
        help="Add ghost homologs already used",
    )
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    omgfile, groupsfile, ntaxa = args
    ntaxa = int(ntaxa)
    ghost = opts.ghost

    # Get gene pair => weight mapping
    weights = get_edges()
    info = get_info()
    # Get gene => taxon mapping
    info = dict((k, v.split()[5]) for k, v in info.items())

    groups = Grouper()

    fp = open(groupsfile)
    for row in fp:
        members = row.strip().split(",")
        groups.join(*members)

    logging.debug(
        "Imported {0} families with {1} members.".format(
            len(groups), groups.num_members
        )
    )

    seen = set()
    omggroups = Grouper()
    fp = open(omgfile)
    for row in fp:
        genes, idxs = row.split()
        genes = genes.split(",")
        seen.update(genes)
        omggroups.join(*genes)

    nmembers = omggroups.num_members
    logging.debug(
        "Imported {0} OMG families with {1} members.".format(len(omggroups), nmembers)
    )
    assert nmembers == len(seen)

    alltaxa = set(str(x) for x in range(ntaxa))
    recruited = []
    fp = open(omgfile)
    for row in fp:
        genes, idxs = row.split()
        genes = genes.split(",")
        a = genes[0]

        idxs = set(idxs.split(","))
        missing_taxa = alltaxa - idxs
        if not missing_taxa:
            print(row.rstrip())
            continue

        leftover = groups[a]
        if not ghost:
            leftover = set(leftover) - seen

        if not leftover:
            print(row.rstrip())
            continue

        leftover_sorted_by_taxa = dict(
            (k, [x for x in leftover if info[x] == k]) for k in missing_taxa
        )

        # print genes, leftover
        # print leftover_sorted_by_taxa
        solutions = []
        for solution in product(*leftover_sorted_by_taxa.values()):
            score = sum(weights.get((a, b), 0) for a in solution for b in genes)
            if score == 0:
                continue
            score += sum(weights.get((a, b), 0) for a, b in combinations(solution, 2))
            solutions.append((score, solution))
            # print solution, score

        best_solution = max(solutions) if solutions else None
        if best_solution is None:
            print(row.rstrip())
            continue

        # print "best ==>", best_solution
        best_score, best_addition = best_solution
        genes.extend(best_addition)
        recruited.extend(best_addition)

        genes = sorted([(info[x], x) for x in genes])
        idxs, genes = zip(*genes)

        if ghost:  # decorate additions so it's clear that they were added
            pgenes = []
            for g in genes:
                if g in recruited and g in seen:
                    pgenes.append("|{0}|".format(g))
                else:
                    pgenes.append(g)
            genes = pgenes

        print("\t".join((",".join(genes), ",".join(idxs))))
        if not ghost:
            seen.update(best_addition)

    logging.debug("Recruited {0} new genes.".format(len(recruited)))


def pairwise_distance(a, b, threadorder):
    d = 0
    for x, y in zip(a, b)[:-1]:  # Last column not used
        x, y = x.strip("|"), y.strip("|")
        if "." in (x, y):
            dd = 50
        else:
            xi, x = threadorder[x]
            yi, y = threadorder[y]
            dd = min(abs(xi - yi), 50)
        d += dd
    return d


def insert_into_threaded(atoms, threaded, threadorder):
    min_idx, min_d = 0, 1000
    for i, t in enumerate(threaded):
        # calculate distance
        d = pairwise_distance(atoms, t, threadorder)
        if d < min_d:
            min_idx = i
            min_d = d

    i = min_idx
    t = threaded[i]
    threaded.insert(i, atoms)
    logging.debug("Insert {0} before {1} (d={2})".format(atoms, t, min_d))


def sort_layout(thread, listfile, column=0):
    """
    Sort the syntelog table according to chromomomal positions. First orient the
    contents against threadbed, then for contents not in threadbed, insert to
    the nearest neighbor.
    """
    from jcvi.formats.base import DictFile

    outfile = listfile.rsplit(".", 1)[0] + ".sorted.list"
    threadorder = thread.order
    fw = open(outfile, "w")
    lt = DictFile(listfile, keypos=column, valuepos=None)
    threaded = []
    imported = set()
    for t in thread:
        accn = t.accn
        if accn not in lt:
            continue

        imported.add(accn)
        atoms = lt[accn]
        threaded.append(atoms)

    assert len(threaded) == len(imported)

    total = sum(1 for x in open(listfile))
    logging.debug("Total: {0}, currently threaded: {1}".format(total, len(threaded)))
    fp = open(listfile)
    for row in fp:
        atoms = row.split()
        accn = atoms[0]
        if accn in imported:
            continue
        insert_into_threaded(atoms, threaded, threadorder)

    for atoms in threaded:
        print("\t".join(atoms), file=fw)

    fw.close()
    logging.debug("File `{0}` sorted to `{1}`.".format(outfile, thread.filename))


def layout(args):
    """
    %prog layout omgfile taxa

    Build column formatted gene lists after omgparse(). Use species list
    separated by comma in place of taxa, e.g. "BR,BO,AN,CN"
    """
    p = OptionParser(layout.__doc__)
    p.add_option("--sort", help="Sort layout file based on bedfile")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    omgfile, taxa = args
    listfile = omgfile.rsplit(".", 1)[0] + ".list"
    taxa = taxa.split(",")
    ntaxa = len(taxa)
    fw = open(listfile, "w")

    data = []
    fp = open(omgfile)
    for row in fp:
        genes, idxs = row.split()
        row = ["."] * ntaxa
        genes = genes.split(",")
        ixs = [int(x) for x in idxs.split(",")]
        for gene, idx in zip(genes, ixs):
            row[idx] = gene
        txs = ",".join(taxa[x] for x in ixs)
        print("\t".join(("\t".join(row), txs)), file=fw)
        data.append(row)

    coldata = zip(*data)
    ngenes = []
    for i, tx in enumerate(taxa):
        genes = [x for x in coldata[i] if x != "."]
        genes = set(x.strip("|") for x in genes)
        ngenes.append((len(genes), tx))

    details = ", ".join("{0} {1}".format(a, b) for a, b in ngenes)
    total = sum(a for a, b in ngenes)
    s = "A list of {0} orthologous families that collectively".format(len(data))
    s += " contain a total of {0} genes ({1})".format(total, details)
    print(s, file=sys.stderr)

    fw.close()
    lastcolumn = ntaxa + 1
    cmd = "sort -k{0},{0} {1} -o {1}".format(lastcolumn, listfile)
    sh(cmd)

    logging.debug("List file written to `{0}`.".format(listfile))
    sort = opts.sort
    if sort:
        thread = Bed(sort)
        sort_layout(thread, listfile)


def omgparse(args):
    """
    %prog omgparse work

    Parse the OMG outputs to get gene lists.
    """
    p = OptionParser(omgparse.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (work,) = args
    omgfiles = glob(op.join(work, "gf*.out"))
    for omgfile in omgfiles:
        omg = OMGFile(omgfile)
        best = omg.best()
        for bb in best:
            genes, taxa = zip(*bb)
            print("\t".join((",".join(genes), ",".join(taxa))))


def group(args):
    """
    %prog group anchorfiles

    Group the anchors into ortho-groups. Can input multiple anchor files.
    """
    p = OptionParser(group.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    anchorfiles = args
    groups = Grouper()
    for anchorfile in anchorfiles:
        ac = AnchorFile(anchorfile)
        for a, b, idx in ac.iter_pairs():
            groups.join(a, b)

    logging.debug(
        "Created {0} groups with {1} members.".format(len(groups), groups.num_members)
    )

    outfile = opts.outfile
    fw = must_open(outfile, "w")
    for g in groups:
        print(",".join(sorted(g)), file=fw)
    fw.close()

    return outfile


def omg(args):
    """
    %prog omg weightsfile

    Run Sankoff's OMG algorithm to get orthologs. Download OMG code at:
    <http://137.122.149.195/IsbraSoftware/OMGMec.html>

    This script only writes the partitions, but not launch OMGMec. You may need to:

    $ parallel "java -cp ~/code/OMGMec TestOMGMec {} 4 > {}.out" ::: work/gf?????

    Then followed by omgparse() to get the gene lists.
    """
    p = OptionParser(omg.__doc__)

    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    weightsfiles = args
    groupfile = group(weightsfiles + ["--outfile=groups"])

    weights = get_weights(weightsfiles)
    info = get_info()

    fp = open(groupfile)

    work = "work"
    mkdir(work)
    for i, row in enumerate(fp):
        gf = op.join(work, "gf{0:05d}".format(i))
        genes = row.rstrip().split(",")

        fw = open(gf, "w")
        contents = ""
        npairs = 0
        for gene in genes:
            gene_pairs = weights[gene]
            for a, b, c in gene_pairs:
                if b not in genes:
                    continue

                contents += "weight {0}".format(c) + "\n"
                contents += info[a] + "\n"
                contents += info[b] + "\n\n"
                npairs += 1

        header = "a group of genes  :length ={0}".format(npairs)
        print(header, file=fw)
        print(contents, file=fw)

        fw.close()


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
        try:
            chr = int(chr)
        except ValueError:
            chr = "0"

        print(
            "\t".join(str(x) for x in (s.accn, chr, s.start, s.end, s.strand, idx, pd)),
            file=fwinfo,
        )
    fwinfo.close()

    logging.debug("Update info file `{0}`.".format(infofile))

    return infofile


def omgprepare(args):
    """
    %prog omgprepare ploidy anchorsfile blastfile

    Prepare to run Sankoff's OMG algorithm to get orthologs.
    """
    from jcvi.formats.blast import cscore
    from jcvi.formats.base import DictFile

    p = OptionParser(omgprepare.__doc__)
    p.add_option("--norbh", action="store_true", help="Disable RBH hits")
    p.add_option(
        "--pctid", default=0, type="int", help="Percent id cutoff for RBH hits"
    )
    p.add_option("--cscore", default=90, type="int", help="C-score cutoff for RBH hits")
    p.set_stripnames()
    p.set_beds()

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    ploidy, anchorfile, blastfile = args
    norbh = opts.norbh
    pctid = opts.pctid
    cs = opts.cscore
    qbed, sbed, qorder, sorder, is_self = check_beds(anchorfile, p, opts)

    fp = open(ploidy)
    genomeidx = dict((x.split()[0], i) for i, x in enumerate(fp))
    fp.close()

    ploidy = DictFile(ploidy)

    geneinfo(qbed, qorder, genomeidx, ploidy)
    geneinfo(sbed, sorder, genomeidx, ploidy)

    pf = blastfile.rsplit(".", 1)[0]
    cscorefile = pf + ".cscore"
    cscore([blastfile, "-o", cscorefile, "--cutoff=0", "--pct"])
    ac = AnchorFile(anchorfile)
    pairs = set((a, b) for a, b, i in ac.iter_pairs())
    logging.debug("Imported {0} pairs from `{1}`.".format(len(pairs), anchorfile))

    weightsfile = pf + ".weights"
    fp = open(cscorefile)
    fw = open(weightsfile, "w")
    npairs = 0
    for row in fp:
        a, b, c, pct = row.split()
        c, pct = float(c), float(pct)
        c = int(c * 100)
        if (a, b) not in pairs:
            if norbh:
                continue
            if c < cs:
                continue
            if pct < pctid:
                continue
            c /= 10  # This severely penalizes RBH against synteny

        print("\t".join((a, b, str(c))), file=fw)
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
        if b == ".":
            if a in adict:
                b = adict[a]
                nrecruited += 1
                b += "'"
        print("\t".join((a, b)), file=fw)

    logging.debug("Recruited {0} pairs from RBH.".format(nrecruited))
    fp.close()
    fw.close()


def ortholog(args):
    """
    %prog ortholog species_a species_b

    Run a sensitive pipeline to find orthologs between two species a and b.
    The pipeline runs LAST and generate .lifted.anchors.

    `--full` mode would assume 1-to-1 quota synteny blocks as the backbone of
    such predictions. Extra orthologs will be recruited from reciprocal best
    match (RBH).
    """
    from jcvi.apps.align import last as last_main
    from jcvi.compara.blastfilter import main as blastfilter_main
    from jcvi.compara.quota import main as quota_main
    from jcvi.compara.synteny import scan, mcscan, liftover
    from jcvi.formats.blast import cscore, filter

    p = OptionParser(ortholog.__doc__)
    p.add_option(
        "--dbtype",
        default="nucl",
        choices=("nucl", "prot"),
        help="Molecule type of subject database",
    )
    p.add_option(
        "--full",
        default=False,
        action="store_true",
        help="Run in full 1x1 mode, including blocks and RBH",
    )
    p.add_option("--cscore", default=0.7, type="float", help="C-score cutoff")
    p.add_option(
        "--dist", default=20, type="int", help="Extent of flanking regions to search"
    )
    p.add_option(
        "-n",
        "--min_size",
        dest="n",
        type="int",
        default=4,
        help="minimum number of anchors in a cluster",
    )
    p.add_option("--quota", help="Quota align parameter")
    p.add_option("--exclude", help="Remove anchors from a previous run")
    p.add_option(
        "--no_strip_names",
        default=False,
        action="store_true",
        help="Do not strip alternative splicing (e.g. At5g06540.1 -> At5g06540)",
    )
    p.add_option(
        "--liftover_dist",
        type="int",
        help="Distance to extend from liftover. Defaults to half of --dist",
    )
    p.set_cpus()
    p.set_dotplot_opts()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    a, b = args
    dbtype = opts.dbtype
    suffix = ".cds" if dbtype == "nucl" else ".pep"
    abed, afasta = a + ".bed", a + suffix
    bbed, bfasta = b + ".bed", b + suffix
    ccscore = opts.cscore
    quota = opts.quota
    exclude = opts.exclude
    dist = "--dist={0}".format(opts.dist)
    minsize_flag = "--min_size={}".format(opts.n)
    cpus_flag = "--cpus={}".format(opts.cpus)

    aprefix = afasta.split(".")[0]
    bprefix = bfasta.split(".")[0]
    pprefix = ".".join((aprefix, bprefix))
    qprefix = ".".join((bprefix, aprefix))
    last = pprefix + ".last"
    if need_update((afasta, bfasta), last):
        last_main([bfasta, afasta, cpus_flag], dbtype)

    if a == b:
        lastself = last + ".P98L0.inverse"
        if need_update(last, lastself):
            filter([last, "--hitlen=0", "--pctid=98", "--inverse", "--noself"])
        last = lastself

    filtered_last = last + ".filtered"
    if need_update(last, filtered_last):
        # If we are doing filtering based on another file then we don't run cscore anymore
        dargs = [last, "--cscore={}".format(ccscore)]
        if exclude:
            dargs += ["--exclude={}".format(exclude)]
        if opts.no_strip_names:
            dargs += ["--no_strip_names"]
        blastfilter_main(dargs)

    anchors = pprefix + ".anchors"
    lifted_anchors = pprefix + ".lifted.anchors"
    pdf = pprefix + ".pdf"
    if not opts.full:
        if need_update(filtered_last, lifted_anchors):
            dargs = [
                filtered_last,
                anchors,
                minsize_flag,
                dist,
                "--liftover={0}".format(last),
            ]
            if opts.no_strip_names:
                dargs += ["--no_strip_names"]
            if opts.liftover_dist:
                dargs += ["--liftover_dist={}".format(opts.liftover_dist)]
            scan(dargs)
        if quota:
            quota_main([lifted_anchors, "--quota={0}".format(quota), "--screen"])
        if need_update(anchors, pdf):
            from jcvi.graphics.dotplot import dotplot_main

            dargs = [anchors]
            if opts.nostdpf:
                dargs += ["--nostdpf"]
            if opts.nochpf:
                dargs += ["--nochpf"]
            if opts.skipempty:
                dargs += ["--skipempty"]
            if opts.genomenames:
                dargs += ["--genomenames", opts.genomenames]
            if opts.theme:
                dargs += ["--theme", opts.theme]
            dotplot_main(dargs)
        return

    if need_update(filtered_last, anchors):
        if opts.no_strip_names:
            scan([filtered_last, anchors, dist, "--no_strip_names"])
        else:
            scan([filtered_last, anchors, dist])

    ooanchors = pprefix + ".1x1.anchors"
    if need_update(anchors, ooanchors):
        quota_main([anchors, "--quota=1:1", "--screen"])

    lifted_anchors = pprefix + ".1x1.lifted.anchors"
    if need_update((last, ooanchors), lifted_anchors):
        if opts.no_strip_names:
            liftover([last, ooanchors, dist, "--no_strip_names"])
        else:
            liftover([last, ooanchors, dist])

    pblocks = pprefix + ".1x1.blocks"
    qblocks = qprefix + ".1x1.blocks"
    if need_update(lifted_anchors, [pblocks, qblocks]):
        mcscan([abed, lifted_anchors, "--iter=1", "-o", pblocks])
        mcscan([bbed, lifted_anchors, "--iter=1", "-o", qblocks])

    rbh = pprefix + ".rbh"
    if need_update(last, rbh):
        cscore([last, "-o", rbh])

    portho = pprefix + ".ortholog"
    qortho = qprefix + ".ortholog"
    if need_update([pblocks, qblocks, rbh], [portho, qortho]):
        make_ortholog(pblocks, rbh, portho)
        make_ortholog(qblocks, rbh, qortho)


def tandem_main(
    blast_file,
    cds_file,
    bed_file,
    N=3,
    P=50,
    is_self=True,
    evalue=0.01,
    strip_name=".",
    ofile=sys.stderr,
    genefam=False,
):

    if genefam:
        N = 1e5

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
            if b.hitlen < min(query_len, subject_len) * P / 100.0:
                continue

            query = gene_name(b.query, sep=strip_name)
            subject = gene_name(b.subject, sep=strip_name)
            qi, q = order[query]
            si, s = order[subject]

            if abs(qi - si) <= N and b.evalue <= evalue:
                if genefam:
                    g.join(query, subject)
                elif q.seqid == s.seqid:
                    g.join(query, subject)

    else:
        homologs = Grouper()
        fp = open(blast_file)
        for row in fp:
            b = BlastLine(row)
            query_len = sizes[b.query]
            subject_len = sizes[b.subject]
            if b.hitlen < min(query_len, subject_len) * P / 100.0:
                continue
            if b.evalue > evalue:
                continue

            query = gene_name(b.query, sep=strip_name)
            subject = gene_name(b.subject, sep=strip_name)
            homologs.join(query, subject)

        if genefam:
            g = homologs
        else:
            g = Grouper()
            for i, atom in enumerate(bed):
                for x in range(1, N + 1):
                    if all(
                        [
                            i - x >= 0,
                            bed[i - x].seqid == atom.seqid,
                            homologs.joined(bed[i - x].accn, atom.accn),
                        ]
                    ):
                        leni = sizes[bed[i].accn]
                        lenx = sizes[bed[i - x].accn]
                        if abs(leni - lenx) > max(leni, lenx) * (1 - P / 100.0):
                            continue
                        g.join(bed[i - x].accn, atom.accn)

    # dump the grouper
    fw = must_open(ofile, "w")
    ngenes, nfamilies = 0, 0
    families = []
    for group in sorted(g):
        if len(group) >= 2:
            print(",".join(sorted(group)), file=fw)
            ngenes += len(group)
            nfamilies += 1
            families.append(sorted(group))

    longest_family = max(families, key=lambda x: len(x))

    # generate reports
    print("Proximal paralogues (dist=%d):" % N, file=sys.stderr)
    print("Total %d genes in %d families" % (ngenes, nfamilies), file=sys.stderr)
    print(
        "Longest families (%d): %s" % (len(longest_family), ",".join(longest_family)),
        file=sys.stderr,
    )

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
    p.add_option(
        "--tandem_Nmax",
        dest="tandem_Nmax",
        type="int",
        default=3,
        help="merge tandem genes within distance",
    )
    p.add_option(
        "--percent_overlap",
        type="int",
        default=50,
        help="tandem genes have >=x% aligned sequence, x=0-100",
    )
    p.set_align(evalue=0.01)
    p.add_option(
        "--not_self",
        default=False,
        action="store_true",
        help="provided is not self blast file",
    )
    p.add_option(
        "--strip_gene_name",
        dest="sep",
        type="string",
        default=".",
        help="strip alternative splicing. Use None for no stripping.",
    )
    p.add_option(
        "--genefamily",
        dest="genefam",
        action="store_true",
        help="compile gene families based on similarity",
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    blast_file, cds_file, bed_file = args
    N = opts.tandem_Nmax
    P = opts.percent_overlap
    is_self = not opts.not_self
    sep = opts.sep
    ofile = opts.outfile

    tandem_main(
        blast_file,
        cds_file,
        bed_file,
        N=N,
        P=P,
        is_self=is_self,
        evalue=opts.evalue,
        strip_name=sep,
        ofile=ofile,
        genefam=opts.genefam,
    )


if __name__ == "__main__":
    main()
