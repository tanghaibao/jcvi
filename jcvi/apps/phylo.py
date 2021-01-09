#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Construct and visualize phylogenetic trees from:
1.  MCSCAN output
2.  CDS sequences in FASTA format

Options are provided for each step:
1.  sequence alignment:
    ClustalW2 or MUSCLE (wrapped on Biopython)
2.  alignment editting:
    GBlocks (optional)
3.  build trees:
    NJ: PHYLIP
    ML: RAxML or PHYML

Optional steps:
- reroot tree
- alternative topology test (SH test)
- TreeFix

The external software needs be installed first.
"""
from __future__ import print_function

import sys
import os
import os.path as op
import logging
import re
import warnings

from math import ceil
from itertools import chain
from functools import partial

import numpy as np
from ete3 import Tree
from Bio import SeqIO, AlignIO
from Bio.Data import CodonTable
from Bio.Emboss.Applications import (
    FSeqBootCommandline,
    FDNADistCommandline,
    FNeighborCommandline,
    FConsenseCommandline,
)
from Bio.Phylo.Applications import PhymlCommandline, RaxmlCommandline

from jcvi.apps.ks import (
    AbstractCommandline,
    find_first_isoform,
    run_mrtrans,
    clustal_align_protein,
    muscle_align_protein,
)
from jcvi.formats.base import must_open, DictFile, LineFile
from jcvi.formats.fasta import Fasta
from jcvi.utils.orderedcollections import OrderedDict
from jcvi.graphics.base import plt, savefig
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, sh, getpath


GBLOCKS_BIN = partial(getpath, name="GBLOCKS", warn="warn")
PHYML_BIN = partial(getpath, name="PHYML", warn="warn")
RAXML_BIN = partial(getpath, name="RAXML", warn="warn")
FPHYLIP_BIN = partial(getpath, name="FPHYLIP", warn="warn")
TREEFIX_BIN = partial(getpath, name="TREEFIX", warn="warn")


class GblocksCommandline(AbstractCommandline):
    """Little commandline for Gblocks
    (http://molevol.cmima.csic.es/castresana/Gblocks.html).

    Accepts alignment in FASTA or NBRF/PIR format.
    """

    def __init__(
        self, aln_file, aln_type="c", command=GBLOCKS_BIN("Gblocks"), **kwargs
    ):

        self.aln_file = aln_file
        self.aln_type = aln_type
        self.command = command

        params = {"b4": 5, "b5": "h", "p": "n"}
        params.update(kwargs)
        self.parameters = ["-{0}={1}".format(k, v) for k, v in params.items()]

    def __str__(self):
        return (
            self.command
            + " %s -t=%s " % (self.aln_file, self.aln_type)
            + " ".join(self.parameters)
        )


class FfitchCommandline(AbstractCommandline):
    """Little commandline for ffitch in EMBOSS
    (http://www.molgen.mpg.de/~beck/embassy/phylipnew/ffitch.html).

    Infer branch lengths of tree.
    """

    def __init__(
        self,
        datafile,
        outtreefile,
        command=FPHYLIP_BIN("ffitch"),
        intreefile=None,
        **kwargs
    ):

        self.datafile = datafile
        self.outtreefile = outtreefile
        self.outfile = datafile.rsplit(".", 1)[0] + ".ffitch"
        self.command = command
        self.intreefile = intreefile if intreefile else '""'

        self.parameters = ["-{0} {1}".format(k, v) for k, v in kwargs.items()]

    def __str__(self):
        return (
            self.command
            + " -datafile %s -intreefile %s -outfile %s -outtreefile %s "
            % (
                self.datafile,
                self.intreefile,
                self.outfile,
                self.outtreefile,
            )
            + " ".join(self.parameters)
        )


class TreeFixCommandline(AbstractCommandline):
    """Little commandline for TreeFix
    (http://compbio.mit.edu/treefix/).
    """

    def __init__(
        self,
        input,
        stree_file,
        smap_file,
        a_ext,
        command=TREEFIX_BIN("treefix"),
        r=False,
        **kwargs
    ):

        self.input = input
        self.s = stree_file
        self.S = smap_file
        self.A = a_ext
        self.command = command

        params = {"V": 1, "l": input.rsplit(".", 1)[0] + ".treefix.log"}
        params.update(kwargs)
        self.parameters = ["-{0} {1}".format(k, v) for k, v in params.items()]
        if r:
            self.parameters.append("-r")

    def __str__(self):
        return (
            self.command
            + " -s %s -S %s -A %s " % (self.s, self.S, self.A)
            + " ".join(self.parameters)
            + " %s" % self.input
        )


def run_treefix(
    input,
    stree_file,
    smap_file,
    a_ext=".fasta",
    o_ext=".dnd",
    n_ext=".treefix.dnd",
    **kwargs
):
    """
    get the ML tree closest to the species tree
    """
    cl = TreeFixCommandline(
        input=input,
        stree_file=stree_file,
        smap_file=smap_file,
        a_ext=a_ext,
        o=o_ext,
        n=n_ext,
        **kwargs
    )
    outtreefile = input.rsplit(o_ext, 1)[0] + n_ext
    print("TreeFix:", cl, file=sys.stderr)
    r, e = cl.run()

    if e:
        print("***TreeFix could not run", file=sys.stderr)
        return None
    else:
        logging.debug("new tree written to {0}".format(outtreefile))
        return outtreefile


def run_gblocks(align_fasta_file, **kwargs):
    """
    remove poorly aligned positions and divergent regions with Gblocks
    """
    cl = GblocksCommandline(aln_file=align_fasta_file, **kwargs)
    r, e = cl.run()

    print("Gblocks:", cl, file=sys.stderr)

    if e:
        print("***Gblocks could not run", file=sys.stderr)
        return None
    else:
        print(r, file=sys.stderr)
        alignp = re.sub(
            r".*Gblocks alignment:.*\(([0-9]{1,3}) %\).*", r"\1", r, flags=re.DOTALL
        )
        alignp = int(alignp)
        if alignp <= 10:
            print(
                "** WARNING ** Only %s %% positions retained by Gblocks. "
                "Results aborted. Using original alignment instead.\n" % alignp,
                file=sys.stderr,
            )
            return None
        else:
            return align_fasta_file + "-gb"


def run_ffitch(distfile, outtreefile, intreefile=None, **kwargs):
    """
    Infer tree branch lengths using ffitch in EMBOSS PHYLIP
    """
    cl = FfitchCommandline(
        datafile=distfile, outtreefile=outtreefile, intreefile=intreefile, **kwargs
    )
    r, e = cl.run()

    if e:
        print("***ffitch could not run", file=sys.stderr)
        return None
    else:
        print("ffitch:", cl, file=sys.stderr)
        return outtreefile


def smart_reroot(treefile, outgroupfile, outfile, format=0):
    """
    simple function to reroot Newick format tree using ete2

    Tree reading format options see here:
    http://packages.python.org/ete2/tutorial/tutorial_trees.html#reading-newick-trees
    """
    tree = Tree(treefile, format=format)
    leaves = [t.name for t in tree.get_leaves()][::-1]
    outgroup = []
    for o in must_open(outgroupfile):
        o = o.strip()
        for leaf in leaves:
            if leaf[: len(o)] == o:
                outgroup.append(leaf)
        if outgroup:
            break

    if not outgroup:
        print(
            "Outgroup not found. Tree {0} cannot be rerooted.".format(treefile),
            file=sys.stderr,
        )
        return treefile

    try:
        tree.set_outgroup(tree.get_common_ancestor(*outgroup))
    except ValueError:
        assert type(outgroup) == list
        outgroup = outgroup[0]
        tree.set_outgroup(outgroup)
    tree.write(outfile=outfile, format=format)

    logging.debug("Rerooted tree printed to {0}".format(outfile))
    return outfile


def build_nj_phylip(alignment, outfile, outgroup, work_dir="."):
    """
    build neighbor joining tree of DNA seqs with PHYLIP in EMBOSS

    PHYLIP manual
    http://evolution.genetics.washington.edu/phylip/doc/
    """

    phy_file = op.join(work_dir, "work", "aln.phy")
    try:
        AlignIO.write(alignment, file(phy_file, "w"), "phylip")
    except ValueError:
        print(
            "Repeated seq name, possibly due to truncation. NJ tree not built.",
            file=sys.stderr,
        )
        return None

    seqboot_out = phy_file.rsplit(".", 1)[0] + ".fseqboot"
    seqboot_cl = FSeqBootCommandline(
        FPHYLIP_BIN("fseqboot"),
        sequence=phy_file,
        outfile=seqboot_out,
        seqtype="d",
        reps=100,
        seed=12345,
    )
    stdout, stderr = seqboot_cl()
    logging.debug("Resampling alignment: %s" % seqboot_cl)

    dnadist_out = phy_file.rsplit(".", 1)[0] + ".fdnadist"
    dnadist_cl = FDNADistCommandline(
        FPHYLIP_BIN("fdnadist"), sequence=seqboot_out, outfile=dnadist_out, method="f"
    )
    stdout, stderr = dnadist_cl()
    logging.debug("Calculating distance for bootstrapped alignments: %s" % dnadist_cl)

    neighbor_out = phy_file.rsplit(".", 1)[0] + ".njtree"
    e = phy_file.rsplit(".", 1)[0] + ".fneighbor"
    neighbor_cl = FNeighborCommandline(
        FPHYLIP_BIN("fneighbor"),
        datafile=dnadist_out,
        outfile=e,
        outtreefile=neighbor_out,
    )
    stdout, stderr = neighbor_cl()
    logging.debug("Building Neighbor Joining tree: %s" % neighbor_cl)

    consense_out = phy_file.rsplit(".", 1)[0] + ".consensustree.nodesupport"
    e = phy_file.rsplit(".", 1)[0] + ".fconsense"
    consense_cl = FConsenseCommandline(
        FPHYLIP_BIN("fconsense"),
        intreefile=neighbor_out,
        outfile=e,
        outtreefile=consense_out,
    )
    stdout, stderr = consense_cl()
    logging.debug("Building consensus tree: %s" % consense_cl)

    # distance without bootstrapping
    dnadist_out0 = phy_file.rsplit(".", 1)[0] + ".fdnadist0"
    dnadist_cl0 = FDNADistCommandline(
        FPHYLIP_BIN("fdnadist"), sequence=phy_file, outfile=dnadist_out0, method="f"
    )
    stdout, stderr = dnadist_cl0()
    logging.debug("Calculating distance for original alignment: %s" % dnadist_cl0)

    # infer branch length on consensus tree
    consensustree1 = phy_file.rsplit(".", 1)[0] + ".consensustree.branchlength"
    run_ffitch(
        distfile=dnadist_out0, outtreefile=consensustree1, intreefile=consense_out
    )

    # write final tree
    ct_s = Tree(consense_out)

    if outgroup:
        t1 = consensustree1 + ".rooted"
        t2 = smart_reroot(consensustree1, outgroup, t1)
        if t2 == t1:
            outfile = outfile.replace(".unrooted", "")
        ct_b = Tree(t2)
    else:
        ct_b = Tree(consensustree1)

    nodesupport = {}
    for node in ct_s.traverse("postorder"):
        node_children = tuple(sorted([f.name for f in node]))
        if len(node_children) > 1:
            nodesupport[node_children] = node.dist / 100.0

    for k, v in nodesupport.items():
        ct_b.get_common_ancestor(*k).support = v
    print(ct_b)
    ct_b.write(format=0, outfile=outfile)

    try:
        s = op.getsize(outfile)
    except OSError:
        s = 0
    if s:
        logging.debug("NJ tree printed to %s" % outfile)
        return outfile, phy_file
    else:
        logging.debug("Something was wrong. NJ tree was not built.")
        return None


def build_ml_phyml(alignment, outfile, work_dir=".", **kwargs):
    """
    build maximum likelihood tree of DNA seqs with PhyML
    """
    phy_file = op.join(work_dir, "work", "aln.phy")
    AlignIO.write(alignment, file(phy_file, "w"), "phylip-relaxed")

    phyml_cl = PhymlCommandline(cmd=PHYML_BIN("phyml"), input=phy_file, **kwargs)
    logging.debug("Building ML tree using PhyML: %s" % phyml_cl)
    stdout, stderr = phyml_cl()

    tree_file = phy_file + "_phyml_tree.txt"
    if not op.exists(tree_file):
        print("***PhyML failed.", file=sys.stderr)
        return None
    sh("cp {0} {1}".format(tree_file, outfile), log=False)

    logging.debug("ML tree printed to %s" % outfile)

    return outfile, phy_file


def build_ml_raxml(alignment, outfile, work_dir=".", **kwargs):
    """
    build maximum likelihood tree of DNA seqs with RAxML
    """
    work_dir = op.join(work_dir, "work")
    mkdir(work_dir)
    phy_file = op.join(work_dir, "aln.phy")
    AlignIO.write(alignment, file(phy_file, "w"), "phylip-relaxed")

    raxml_work = op.abspath(op.join(op.dirname(phy_file), "raxml_work"))
    mkdir(raxml_work)
    raxml_cl = RaxmlCommandline(
        cmd=RAXML_BIN("raxmlHPC"),
        sequences=phy_file,
        algorithm="a",
        model="GTRGAMMA",
        parsimony_seed=12345,
        rapid_bootstrap_seed=12345,
        num_replicates=100,
        name="aln",
        working_dir=raxml_work,
        **kwargs
    )

    logging.debug("Building ML tree using RAxML: %s" % raxml_cl)
    stdout, stderr = raxml_cl()

    tree_file = "{0}/RAxML_bipartitions.aln".format(raxml_work)
    if not op.exists(tree_file):
        print("***RAxML failed.", file=sys.stderr)
        sh("rm -rf %s" % raxml_work, log=False)
        return None
    sh("cp {0} {1}".format(tree_file, outfile), log=False)

    logging.debug("ML tree printed to %s" % outfile)
    sh("rm -rf %s" % raxml_work)

    return outfile, phy_file


def SH_raxml(reftree, querytree, phy_file, shout="SH_out.txt"):
    """
    SH test using RAxML

    querytree can be a single tree or a bunch of trees (eg. from bootstrapping)
    """
    assert op.isfile(reftree)
    shout = must_open(shout, "a")

    raxml_work = op.abspath(op.join(op.dirname(phy_file), "raxml_work"))
    mkdir(raxml_work)
    raxml_cl = RaxmlCommandline(
        cmd=RAXML_BIN("raxmlHPC"),
        sequences=phy_file,
        algorithm="h",
        model="GTRGAMMA",
        name="SH",
        starting_tree=reftree,
        bipartition_filename=querytree,
        working_dir=raxml_work,
    )

    logging.debug("Running SH test in RAxML: %s" % raxml_cl)
    o, stderr = raxml_cl()
    # hard coded
    try:
        pval = re.search("(Significantly.*:.*)", o).group(0)
    except:
        print("SH test failed.", file=sys.stderr)
    else:
        pval = pval.strip().replace("\t", " ").replace("%", "\%")
        print("{0}\t{1}".format(op.basename(querytree), pval), file=shout)
        logging.debug("SH p-value appended to %s" % shout.name)

    shout.close()
    return shout.name


CODON_TRANSLATION = CodonTable.standard_dna_table.forward_table
FOURFOLD = {
    "CTT": "L",
    "ACA": "T",
    "ACG": "T",
    "CCT": "P",
    "CTG": "L",
    "CTA": "L",
    "ACT": "T",
    "CCG": "P",
    "CCA": "P",
    "CCC": "P",
    "GGT": "G",
    "CGA": "R",
    "CGC": "R",
    "CGG": "R",
    "GGG": "G",
    "GGA": "G",
    "GGC": "G",
    "CGT": "R",
    "GTA": "V",
    "GTC": "V",
    "GTG": "V",
    "GTT": "V",
    "CTC": "L",
    "TCT": "S",
    "TCG": "S",
    "TCC": "S",
    "ACC": "T",
    "TCA": "S",
    "GCA": "A",
    "GCC": "A",
    "GCG": "A",
    "GCT": "A",
}


def subalignment(alnfle, subtype, alntype="fasta"):
    """
    Subset synonymous or fourfold degenerate sites from an alignment

    input should be a codon alignment
    """
    aln = AlignIO.read(alnfle, alntype)
    alnlen = aln.get_alignment_length()
    nseq = len(aln)
    subaln = None
    subalnfile = alnfle.rsplit(".", 1)[0] + "_{0}.{1}".format(subtype, alntype)

    if subtype == "synonymous":
        for j in range(0, alnlen, 3):
            aa = None
            for i in range(nseq):
                codon = str(aln[i, j : j + 3].seq)
                if codon not in CODON_TRANSLATION:
                    break
                if aa and CODON_TRANSLATION[codon] != aa:
                    break
                else:
                    aa = CODON_TRANSLATION[codon]
            else:
                if subaln is None:
                    subaln = aln[:, j : j + 3]
                else:
                    subaln += aln[:, j : j + 3]

    if subtype == "fourfold":
        for j in range(0, alnlen, 3):
            for i in range(nseq):
                codon = str(aln[i, j : j + 3].seq)
                if codon not in FOURFOLD:
                    break
            else:
                if subaln is None:
                    subaln = aln[:, j : j + 3]
                else:
                    subaln += aln[:, j : j + 3]

    if subaln:
        AlignIO.write(subaln, subalnfile, alntype)
        return subalnfile
    else:
        print("No sites {0} selected.".format(subtype), file=sys.stderr)
        return None


def merge_rows_local(
    filename, ignore=".", colsep="\t", local=10, fieldcheck=True, fsep=","
):
    """
    merge overlapping rows within given row count distance
    """
    fw = must_open(filename + ".merged", "w")
    rows = file(filename).readlines()
    rows = [row.strip().split(colsep) for row in rows]
    l = len(rows[0])

    for rowi, row in enumerate(rows):
        n = len(rows)
        i = rowi + 1
        while i <= min(rowi + local, n - 1):
            merge = 1
            row2 = rows[i]
            for j in range(l):
                a = row[j]
                b = row2[j]
                if fieldcheck:
                    a = set(a.split(fsep))
                    a = fsep.join(sorted(list(a)))
                    b = set(b.split(fsep))
                    b = fsep.join(sorted(list(b)))

                if all([a != ignore, b != ignore, a not in b, b not in a]):
                    merge = 0
                    i += 1
                    break

            if merge:
                for x in range(l):
                    if row[x] == ignore:
                        rows[rowi][x] = row2[x]
                    elif row[x] in row2[x]:
                        rows[rowi][x] = row2[x]
                    else:
                        rows[rowi][x] = row[x]
                row = rows[rowi]
                rows.remove(row2)

        print(colsep.join(row), file=fw)
    fw.close()

    return fw.name


def add_tandems(mcscanfile, tandemfile):
    """
    add tandem genes to anchor genes in mcscan file
    """
    tandems = [f.strip().split(",") for f in file(tandemfile)]
    fw = must_open(mcscanfile + ".withtandems", "w")
    fp = must_open(mcscanfile)
    seen = set()
    for i, row in enumerate(fp):
        if row[0] == "#":
            continue
        anchorslist = row.strip().split("\t")
        anchors = set([a.split(",")[0] for a in anchorslist])
        anchors.remove(".")
        if anchors & seen == anchors:
            continue

        newanchors = []
        for a in anchorslist:
            if a == ".":
                newanchors.append(a)
                continue
            for t in tandems:
                if a in t:
                    newanchors.append(",".join(t))
                    seen.update(t)
                    break
            else:
                newanchors.append(a)
                seen.add(a)
        print("\t".join(newanchors), file=fw)

    fw.close()
    newmcscanfile = merge_rows_local(fw.name)

    logging.debug(
        "Tandems added to `{0}`. Results in `{1}`".format(mcscanfile, newmcscanfile)
    )
    fp.seek(0)
    logging.debug(
        "{0} rows merged to {1} rows".format(
            len(fp.readlines()), len(file(newmcscanfile).readlines())
        )
    )
    sh("rm %s" % fw.name)

    return newmcscanfile


def main():

    actions = (
        ("prepare", "prepare cds sequences from .mcscan"),
        ("build", "build NJ and ML trees from cds"),
        ("draw", "draw Newick formatted trees"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare mcscanfile cdsfile [options]

    Pick sequences from cdsfile to form fasta files, according to multiple
    alignment in the mcscanfile.
    The fasta sequences can then be used to construct phylogenetic tree.

    Use --addtandem=tandemfile to collapse tandems of anchors into single row.
    The tandemfile must be provided with *ALL* genomes involved, otherwise
    result will be incomplete and redundant.
    """
    from jcvi.graphics.base import discrete_rainbow

    p = OptionParser(prepare.__doc__)
    p.add_option("--addtandem", help="path to tandemfile")
    p.add_option(
        "--writecolors",
        default=False,
        action="store_true",
        help="generate a gene_name to color mapping file which will be taken "
        "by jcvi.apps.phylo.draw",
    )
    p.set_outdir(outdir="sequences")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mcscanfile, cdsfile = args

    if opts.addtandem:
        tandemfile = opts.addtandem
        mcscanfile_with_tandems = add_tandems(mcscanfile, tandemfile)
        mcscanfile = mcscanfile_with_tandems

    seqdir = opts.outdir
    mkdir(seqdir)
    f = Fasta(cdsfile)
    fp = must_open(mcscanfile)
    if opts.writecolors:
        fc = must_open("leafcolors.txt", "w")

    n = 0
    for i, row in enumerate(fp):
        row = row.strip().split("\t")
        if i == 0:
            l = len(row)
            if l <= 20:
                colors = discrete_rainbow(l, shuffle=False)[1]
            else:
                colors = discrete_rainbow(l, usepreset=False, shuffle=False)[1]
                warnings.warn(
                    "*** WARNING ***\n"
                    "Too many columns. Colors may not be all distinctive."
                )

        assert len(row) == l, "All rows should have same number of fields."

        anchors = set()
        for j, atom in enumerate(row):
            color = "%s,%s,%s" % colors[j]
            if atom == ".":
                continue
            elif "," in atom:
                atom = atom.split(",")
                for a in atom:
                    fc.write("{0}\t{1}\n".format(a, color))
                    anchors.add(a)
            else:
                fc.write("{0}\t{1}\n".format(atom, color))
                anchors.add(atom)

        if len(anchors) <= 3:
            print(
                "Not enough seqs to build trees for {0}".format(anchors),
                file=sys.stderr,
            )
            continue

        pivot = row[0]
        fw = must_open("%s/%s.cds" % (seqdir, pivot), "w")
        for a in anchors:
            if a not in f:
                print(a)
                a = find_first_isoform(a, f)
                assert a, a
            arec = f[a]
            SeqIO.write((arec), fw, "fasta")
        fw.close()
        n += 1

    if opts.writecolors:
        fc.close()
        logging.debug("leaf colors written to `{0}`".format(fc.name))

    logging.debug("cds of {0} syntelog groups written to {1}/".format(n, seqdir))

    return seqdir


def build(args):
    """
    %prog build [prot.fasta] cds.fasta [options] --outdir=outdir

    This function wraps on the following steps:
    1. msa using ClustalW2 or MUSCLE(default)
    2. (optional) alignment editing using Gblocks
    3. build NJ tree using PHYLIP in EMBOSS package
       seq names should be unique by first 10 chars (restriction of PHYLIP)
    4. build ML tree using RAxML(default) or PHYML, use keywords raxml or phyml,
       *WARNING* maybe slow with large dataset

    If an outgroup file is provided, the result tree will be rooted on the
    outgroup according to order in the file, i.e. the name in row1 will be
    tried first. If not found, row2 will be used, etc.
    Tail truncated names can be provided so long as it is unique among the seqs.
    If not uniq, the first occurrence will be used. For example, if you have
    two moss sequences in your input, then the tree will be rooted on the
    first moss sequence encountered by the program, unless they are monophylic,
     in which case the root will be their common ancestor.

    --stree and --smap are required if --treefix is set.

    Trees can be edited again using an editor such as Dendroscope. This
    is the recommended way to get highly customized trees.

    Newick format trees will be deposited into outdir (. by default).
    """
    from jcvi.formats.fasta import translate

    p = OptionParser(build.__doc__)
    p.add_option(
        "--longest",
        action="store_true",
        help="Get longest ORF, only works if no pep file, e.g. ESTs",
    )
    p.add_option(
        "--nogblocks",
        action="store_true",
        help="don't use Gblocks to edit alignment",
    )
    p.add_option(
        "--synonymous",
        action="store_true",
        help="extract synonymous sites of the alignment",
    )
    p.add_option(
        "--fourfold",
        action="store_true",
        help="extract fourfold degenerate sites of the alignment",
    )
    p.add_option(
        "--msa",
        default="muscle",
        choices=("clustalw", "muscle"),
        help="software used to align the proteins",
    )
    p.add_option(
        "--noneighbor",
        action="store_true",
        help="don't build NJ tree",
    )
    p.add_option(
        "--ml",
        default=None,
        choices=("raxml", "phyml"),
        help="software used to build ML tree",
    )
    p.add_option("--outgroup", help="path to file containing outgroup orders")
    p.add_option("--SH", help="path to reference Newick tree")
    p.add_option("--shout", default="SH_out.txt", help="SH output file name")
    p.add_option(
        "--treefix",
        action="store_true",
        help="use TreeFix to rearrange ML tree",
    )
    p.add_option("--stree", help="path to species Newick tree")
    p.add_option(
        "--smap",
        help="path to smap file: gene_name_pattern<tab>species_name",
    )
    p.set_outdir()

    opts, args = p.parse_args(args)
    gblocks = not opts.nogblocks
    synonymous = opts.synonymous
    fourfold = opts.fourfold
    neighbor = not opts.noneighbor
    outgroup = opts.outgroup
    outdir = opts.outdir

    if len(args) == 1:
        protein_file, dna_file = None, args[0]
    elif len(args) == 2:
        protein_file, dna_file = args
    else:
        print("Incorrect arguments", file=sys.stderr)
        sys.exit(not p.print_help())

    if opts.treefix:
        stree = opts.stree
        smap = opts.smap
        assert stree and smap, "TreeFix requires stree and smap files."
        opts.ml = "raxml"

    treedir = op.join(outdir, "tree")
    mkdir(treedir)

    if not protein_file:
        protein_file = dna_file + ".pep"
        translate_args = [dna_file, "--outfile=" + protein_file]
        if opts.longest:
            translate_args += ["--longest"]
        dna_file, protein_file = translate(translate_args)

    work_dir = op.join(outdir, "alignment")
    mkdir(work_dir)
    p_recs = list(SeqIO.parse(open(protein_file), "fasta"))
    if opts.msa == "clustalw":
        align_fasta = clustal_align_protein(p_recs, work_dir)
    elif opts.msa == "muscle":
        align_fasta = muscle_align_protein(p_recs, work_dir)

    n_recs = list(SeqIO.parse(open(dna_file), "fasta"))
    mrtrans_fasta = run_mrtrans(align_fasta, n_recs, work_dir, outfmt="fasta")

    if not mrtrans_fasta:
        logging.debug("pal2nal aborted. Cannot reliably build tree for %s", dna_file)
        return

    codon_aln_fasta = mrtrans_fasta
    if gblocks:
        gb_fasta = run_gblocks(mrtrans_fasta)
        codon_aln_fasta = gb_fasta if gb_fasta else codon_aln_fasta

    else:
        if synonymous:
            codon_aln_fasta = subalignment(mrtrans_fasta, "synonymous")

        if fourfold:
            codon_aln_fasta = subalignment(mrtrans_fasta, "fourfold")

    if not neighbor and not opts.ml:
        return codon_aln_fasta

    alignment = AlignIO.read(codon_aln_fasta, "fasta")
    if len(alignment) <= 3:
        raise ValueError("Too few seqs to build tree.")

    mkdir(op.join(treedir, "work"))
    if neighbor:
        out_file = op.join(
            treedir, op.basename(dna_file).rsplit(".", 1)[0] + ".NJ.unrooted.dnd"
        )
        try:
            outfile, phy_file = build_nj_phylip(
                alignment, outfile=out_file, outgroup=outgroup, work_dir=treedir
            )
        except:
            print("NJ tree cannot be built for {0}".format(dna_file))

        if opts.SH:
            reftree = opts.SH
            querytree = outfile
            SH_raxml(reftree, querytree, phy_file, shout=opts.shout)

    if opts.ml:
        out_file = op.join(
            treedir, op.basename(dna_file).rsplit(".", 1)[0] + ".ML.unrooted.dnd"
        )

        if opts.ml == "phyml":
            try:
                outfile, phy_file = build_ml_phyml(
                    alignment, outfile=out_file, work_dir=treedir
                )
            except:
                print("ML tree cannot be built for {0}".format(dna_file))

        elif opts.ml == "raxml":
            try:
                outfile, phy_file = build_ml_raxml(
                    alignment, outfile=out_file, work_dir=treedir
                )
            except:
                print("ML tree cannot be built for {0}".format(dna_file))

        if outgroup:
            new_out_file = out_file.replace(".unrooted", "")
            t = smart_reroot(
                treefile=out_file, outgroupfile=outgroup, outfile=new_out_file
            )
            if t == new_out_file:
                sh("rm %s" % out_file)
                outfile = new_out_file

        if opts.SH:
            reftree = opts.SH
            querytree = outfile
            SH_raxml(reftree, querytree, phy_file, shout=opts.shout)

        if opts.treefix:
            treefix_dir = op.join(treedir, "treefix")
            assert mkdir(treefix_dir, overwrite=True)

            sh("cp {0} {1}/".format(outfile, treefix_dir))
            input = op.join(treefix_dir, op.basename(outfile))
            aln_file = input.rsplit(".", 1)[0] + ".fasta"
            SeqIO.write(alignment, aln_file, "fasta")

            outfile = run_treefix(
                input=input,
                stree_file=stree,
                smap_file=smap,
                a_ext=".fasta",
                o_ext=".dnd",
                n_ext=".treefix.dnd",
            )

    return outfile


def _draw_trees(
    trees, nrow=1, ncol=1, rmargin=0.3, iopts=None, outdir=".", shfile=None, **kwargs
):
    """
    Draw one or multiple trees on one plot.
    """
    from jcvi.graphics.tree import draw_tree

    if shfile:
        SHs = DictFile(shfile, delimiter="\t")

    ntrees = len(trees)
    n = nrow * ncol
    for x in range(int(ceil(float(ntrees) / n))):
        fig = plt.figure(1, (iopts.w, iopts.h)) if iopts else plt.figure(1, (5, 5))
        root = fig.add_axes([0, 0, 1, 1])

        xiv = 1.0 / ncol
        yiv = 1.0 / nrow
        xstart = list(np.arange(0, 1, xiv)) * nrow
        ystart = list(chain(*zip(*[list(np.arange(0, 1, yiv))[::-1]] * ncol)))
        for i in range(n * x, n * (x + 1)):
            if i == ntrees:
                break
            ax = fig.add_axes([xstart[i % n], ystart[i % n], xiv, yiv])
            f = trees.keys()[i]
            tree = trees[f]
            try:
                SH = SHs[f]
            except:
                SH = None
            draw_tree(
                ax,
                tree,
                rmargin=rmargin,
                reroot=False,
                supportcolor="r",
                SH=SH,
                **kwargs
            )

        root.set_xlim(0, 1)
        root.set_ylim(0, 1)
        root.set_axis_off()

        format = iopts.format if iopts else "pdf"
        dpi = iopts.dpi if iopts else 300
        if n == 1:
            image_name = f.rsplit(".", 1)[0] + "." + format
        else:
            image_name = "trees{0}.{1}".format(x, format)
        image_name = op.join(outdir, image_name)
        savefig(image_name, dpi=dpi, iopts=iopts)
        plt.clf()


def draw(args):
    """
    %prog draw --input newicktrees [options]

    Draw phylogenetic trees into single or combined plots.
    Input trees should be one of the following:
    1.  single Newick format tree file
    2.  a dir containing *ONLY* the tree files to be drawn

    Newick format:
    http://evolution.genetics.washington.edu/phylip/newicktree.html

    This function wraps on jcvi.graphics.tree
    This function is better used for trees generated by jcvi.apps.phylo (rooted
    if possible). For drawing general Newick trees from external sources invoke
    jcvi.graphics.tree directly, which also gives more drawing options.
    """
    trunc_name_options = ["headn", "oheadn", "tailn", "otailn"]
    p = OptionParser(draw.__doc__)
    p.add_option(
        "--input",
        help="path to single input tree file or a dir "
        "containing ONLY the input tree files",
    )
    p.add_option(
        "--combine",
        type="string",
        default="1x1",
        help="combine multiple trees into one plot in nrowxncol",
    )
    p.add_option(
        "--trunc_name",
        default=None,
        help="Options are: {0}. "
        "truncate first n chars, retains only first n chars, "
        "truncate last n chars, retain only last chars. "
        "n=1~99.".format(trunc_name_options),
    )
    p.add_option(
        "--SH",
        default=None,
        help="path to a file containing SH test p-values in format:"
        "tree_file_name<tab>p-values "
        "This file can be generated with jcvi.apps.phylo build",
    )
    p.add_option(
        "--scutoff",
        default=50,
        type="int",
        help="cutoff for displaying node support, 0-100",
    )
    p.add_option(
        "--barcode",
        default=None,
        help="path to seq/taxon name barcode mapping file: "
        "barcode<tab>new_name "
        "This option is downstream of `--trunc_name`",
    )
    p.add_option(
        "--leafcolorfile",
        default=None,
        help="path to a mapping file containing font colors "
        "for the OTUs: leafname<tab>color",
    )
    p.set_outdir()
    opts, args, iopts = p.set_image_options(figsize="8x6")
    input = opts.input
    outdir = opts.outdir
    combine = opts.combine.split("x")
    trunc_name = opts.trunc_name
    SH = opts.SH

    mkdir(outdir)
    if not input:
        sys.exit(not p.print_help())
    elif op.isfile(input):
        trees_file = input
        treenames = [op.basename(input)]
    elif op.isdir(input):
        trees_file = op.join(outdir, "alltrees.dnd")
        treenames = []
        for f in sorted(os.listdir(input)):
            sh("cat {0}/{1} >> {2}".format(input, f, trees_file), log=False)
            treenames.append(f)
    else:
        sys.exit(not p.print_help())

    trees = OrderedDict()
    tree = ""
    i = 0
    for row in LineFile(trees_file, comment="#", load=True).lines:
        if i == len(treenames):
            break
        if not len(row):
            continue

        if ";" in row:
            # sanity check
            if row.index(";") != len(row) - 1:
                ts = row.split(";")
                for ii in range(len(ts) - 1):
                    ts[ii] += ";"
            else:
                ts = [row]
            for t in ts:
                if ";" in t:
                    tree += t
                    if tree:
                        trees[treenames[i]] = tree
                        tree = ""
                        i += 1
                else:
                    tree += t
        else:
            tree += row

    logging.debug("A total of {0} trees imported.".format(len(trees)))
    sh("rm {0}".format(op.join(outdir, "alltrees.dnd")))

    _draw_trees(
        trees,
        nrow=int(combine[0]),
        ncol=int(combine[1]),
        rmargin=0.3,
        iopts=iopts,
        outdir=outdir,
        shfile=SH,
        trunc_name=trunc_name,
        scutoff=opts.scutoff,
        barcodefile=opts.barcode,
        leafcolorfile=opts.leafcolorfile,
    )


if __name__ == "__main__":
    main()
