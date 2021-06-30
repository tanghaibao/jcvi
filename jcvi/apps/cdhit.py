#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using CD-HIT to remove duplicate reads.
"""
import os.path as op
import sys
import logging

from collections import defaultdict

from jcvi.formats.base import LineFile, read_block, must_open
from jcvi.formats.fastq import fasta
from jcvi.utils.cbook import percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


class ClstrLine(object):
    """
    Lines like these:
    0       12067nt, >LAP012517... at -/99.85%
    1       15532nt, >MOL158919... *
    2       15515nt, >SES069071... at +/99.85%
    """

    def __init__(self, row):
        a, b = row.split(">", 1)
        a = a.split("nt")[0]
        sid, size = a.split()
        self.size = int(size)
        self.name = b.split("...")[0]
        self.rep = row.rstrip()[-1] == "*"


class ClstrFile(LineFile):
    def __init__(self, filename):
        super(ClstrFile, self).__init__(filename)
        assert filename.endswith(".clstr")

        fp = open(filename)
        for clstr, members in read_block(fp, ">"):
            self.append([ClstrLine(x) for x in members])

    def iter_sizes(self):
        for members in self:
            yield len(members)

    def iter_reps(self):
        for i, members in enumerate(self):
            for b in members:
                if b.rep:
                    yield i, b.name

    def iter_reps_prefix(self, prefix=3):
        for i, members in enumerate(self):
            d = defaultdict(list)
            for b in members:
                pp = b.name[:prefix]
                d[pp].append(b)

            for pp, members_with_same_pp in sorted(d.items()):
                yield i, max(members_with_same_pp, key=lambda x: x.size).name


def main():

    actions = (
        ("ids", "get the representative ids from clstr file"),
        ("deduplicate", "use `cd-hit-est` to remove duplicate reads"),
        ("filter", "filter consensus sequence with min cluster size"),
        ("summary", "parse cdhit.clstr file to get distribution of cluster sizes"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def filter(args):
    """
    %prog filter *.consensus.fasta

    Filter consensus sequence with min cluster size.
    """
    from jcvi.formats.fasta import Fasta, SeqIO

    p = OptionParser(filter.__doc__)
    p.add_option("--minsize", default=2, type="int", help="Minimum cluster size")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastafiles = args
    minsize = opts.minsize
    totalreads = totalassembled = 0
    fw = must_open(opts.outfile, "w")
    for i, fastafile in enumerate(fastafiles):
        f = Fasta(fastafile, lazy=True)
        pf = "s{0:03d}".format(i)
        nreads = nsingletons = nclusters = 0
        for desc, rec in f.iterdescriptions_ordered():
            nclusters += 1
            if desc.startswith("singleton"):
                nsingletons += 1
                nreads += 1
                continue
            # consensus_for_cluster_0 with 63 sequences
            name, w, size, seqs = desc.split()
            assert w == "with"
            size = int(size)
            nreads += size
            if size < minsize:
                continue
            rec.description = rec.description.split(None, 1)[-1]
            rec.id = pf + "_" + rec.id
            SeqIO.write(rec, fw, "fasta")
        logging.debug(
            "Scanned {0} clusters with {1} reads ..".format(nclusters, nreads)
        )
        cclusters, creads = nclusters - nsingletons, nreads - nsingletons
        logging.debug(
            "Saved {0} clusters (min={1}) with {2} reads (avg:{3}) [{4}]".format(
                cclusters, minsize, creads, creads / cclusters, pf
            )
        )
        totalreads += nreads
        totalassembled += nreads - nsingletons
    logging.debug("Total assembled: {0}".format(percentage(totalassembled, totalreads)))


def ids(args):
    """
    %prog ids cdhit.clstr

    Get the representative ids from clstr file.
    """
    p = OptionParser(ids.__doc__)
    p.add_option("--prefix", type="int", help="Find rep id for prefix of len")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (clstrfile,) = args
    cf = ClstrFile(clstrfile)
    prefix = opts.prefix
    if prefix:
        reads = list(cf.iter_reps_prefix(prefix=prefix))
    else:
        reads = list(cf.iter_reps())

    nreads = len(reads)
    idsfile = clstrfile.replace(".clstr", ".ids")
    fw = open(idsfile, "w")
    for i, name in reads:
        print("\t".join(str(x) for x in (i, name)), file=fw)

    logging.debug(
        "A total of {0} unique reads written to `{1}`.".format(nreads, idsfile)
    )
    fw.close()

    return idsfile


def summary(args):
    """
    %prog summary cdhit.clstr

    Parse cdhit.clstr file to get distribution of cluster sizes.
    """
    from jcvi.graphics.histogram import loghistogram

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (clstrfile,) = args
    cf = ClstrFile(clstrfile)
    data = list(cf.iter_sizes())
    loghistogram(data, summary=True)


def deduplicate(args):
    """
    %prog deduplicate fastafile

    Wraps `cd-hit-est` to remove duplicate sequences.
    """
    p = OptionParser(deduplicate.__doc__)
    p.set_align(pctid=96, pctcov=0)
    p.add_option(
        "--fast",
        default=False,
        action="store_true",
        help="Place sequence in the first cluster",
    )
    p.add_option(
        "--consensus",
        default=False,
        action="store_true",
        help="Compute consensus sequences",
    )
    p.add_option(
        "--reads",
        default=False,
        action="store_true",
        help="Use `cd-hit-454` to deduplicate",
    )
    p.add_option(
        "--samestrand",
        default=False,
        action="store_true",
        help="Enforce same strand alignment",
    )
    p.set_home("cdhit")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    identity = opts.pctid / 100.0
    fastafile, qualfile = fasta([fastafile, "--seqtk"])

    ocmd = "cd-hit-454" if opts.reads else "cd-hit-est"
    cmd = op.join(opts.cdhit_home, ocmd)
    cmd += " -c {0}".format(identity)
    if ocmd == "cd-hit-est":
        cmd += " -d 0"  # include complete defline
        if opts.samestrand:
            cmd += " -r 0"
    if not opts.fast:
        cmd += " -g 1"
    if opts.pctcov != 0:
        cmd += " -aL {0} -aS {0}".format(opts.pctcov / 100.0)

    dd = fastafile + ".P{0}.cdhit".format(opts.pctid)
    clstr = dd + ".clstr"

    cmd += " -M 0 -T {0} -i {1} -o {2}".format(opts.cpus, fastafile, dd)
    if need_update(fastafile, (dd, clstr)):
        sh(cmd)

    if opts.consensus:
        cons = dd + ".consensus"
        cmd = op.join(opts.cdhit_home, "cdhit-cluster-consensus")
        cmd += " clustfile={0} fastafile={1} output={2} maxlen=1".format(
            clstr, fastafile, cons
        )
        if need_update((clstr, fastafile), cons):
            sh(cmd)

    return dd


if __name__ == "__main__":
    main()
