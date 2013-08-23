#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utility script for annotations based on MAKER.

Many of the routines in this script is to select among a set of conflicting
models, either through accuracy (batcheval) or simply the length (longest).
"""

import os
import os.path as op
import sys
import logging

from collections import defaultdict
from optparse import OptionParser

from jcvi.formats.base import LineFile, write_file
from jcvi.apps.softlink import get_abs_path
from jcvi.apps.grid import GridProcess
from jcvi.apps.base import ActionDispatcher, need_update, popen, debug, sh, \
        mkdir, set_grid_opts
debug()


class CTLine (object):

    def __init__(self, row):
        row = row.strip()
        tag = value = real = comment = ""
        if "#" in row:
            real, comment = row.split("#", 1)
        if "=" in real:
            tag, value = real.split("=", 1)

        self.tag = tag.strip()
        self.value = value.strip()
        self.comment = comment.strip()

    def __str__(self):
        tag = self.tag
        value = self.value
        comment = self.comment

        s = "=".join((tag, value)) if tag else ""
        if s:
            if comment:
                s += "  # " + comment
        else:
            if comment:
                s += "# " + comment
        return s


class CTLFile (LineFile):

    def __init__(self, filename):
        super(CTLFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            self.append(CTLine(row))
        fp.close()

    def update_abs_path(self):
        for r in self:
            path = r.value
            if path and op.exists(path):
                npath = get_abs_path(path)
                logging.debug("{0}={1} => {2}".format(r.tag, path, npath))
                r.value = npath

    def update_tag(self, key, value):
        for r in self:
            if r.tag == key:
                logging.debug("{0}={1} => {2}".format(r.tag, r.value, value))
                r.value = value
                break

    def write_file(self, filename):
        fw = open(filename, "w")
        for r in self:
            print >> fw, r
        fw.close()
        logging.debug("File written to `{0}`".format(filename))


def main():

    actions = (
        ('parallel', 'partition the genome into parts and run separately'),
        ('merge', 'generate the gff files after parallel'),
        ('datastore', 'generate a list of gff filenames to merge'),
        ('split', 'split MAKER models by checking against evidences'),
        ('batcheval', 'calls bed.evaluate() in batch'),
        ('longest', 'pick the longest model per group'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


arraysh = """
DIR=`awk "NR==$SGE_TASK_ID" {0}`
cd $DIR
{1} --ignore_nfs_tmp"""


def parallel(args):
    """
    %prog parallel genome.fasta N

    Partition the genome into parts and run separately. This is useful if MAKER
    is to be run on the grid.
    """
    from jcvi.formats.base import split

    p = OptionParser(parallel.__doc__)
    p.add_option("--maker_home", default="~/htang/export/maker",
                 help="Home directory for MAKER [default: %default]")
    set_grid_opts(p)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    genome, NN = args
    threaded = opts.threaded
    N = int(NN)
    assert 1 < N < 1000, "Required: 1 < N < 1000!"

    outdir = "outdir"
    fs = split([genome, outdir, "--cycle", "-n", NN])

    c = CTLFile("maker_opts.ctl")
    c.update_abs_path()
    if threaded > 1:
        c.update_tag("cpus", threaded)

    cwd = os.getcwd()
    dirs = []
    for name in fs.names:
        fn = get_abs_path(name)
        bn = op.basename(name)
        dirs.append(bn)
        c.update_tag("genome", fn)
        mkdir(bn)
        sh("cp *.ctl {0}".format(bn))

        os.chdir(bn)
        c.write_file("maker_opts.ctl")
        os.chdir(cwd)

    jobs = "jobs"
    fw = open(jobs, "w")
    print >> fw, "\n".join(dirs)
    fw.close()

    # Submit to grid
    ncmds = len(dirs)
    runfile = "array.sh"
    cmd = op.join(opts.maker_home, "bin/maker")
    contents = arraysh.format(jobs, cmd)
    write_file(runfile, contents, meta="run script")

    # qsub script
    outfile = "maker.\$TASK_ID.out"
    p = GridProcess(runfile, outfile=outfile, errfile=outfile,
                    queue=opts.queue, threaded=threaded,
                    arr=ncmds)
    qsubfile = "qsub.sh"
    qsub = p.build()
    write_file(qsubfile, qsub, meta="run script")


mergesh = """
BASE=$1
cd $1.fa/$1.maker.output
{0} -n -d $1_master_datastore_index.log
mv $1.all.gff ../../
"""


def merge(args):
    """
    %prog merge outdir output.gff

    Follow-up command after grid jobs are completed after parallel().
    """
    from glob import glob

    p = OptionParser(merge.__doc__)
    p.add_option("--maker_home", default="~/htang/export/maker",
                 help="Home directory for MAKER [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    outdir, outputgff = args
    fsnames = [op.basename(x).replace(".fa", "") for x in glob(op.join(outdir, "*.fa"))]
    cmd = op.join(opts.maker_home, "bin/gff3_merge")

    outfile = "merge.sh"
    write_file(outfile, mergesh.format(cmd), meta="run script")

    # Generate per split directory
    # Note that gff3_merge write to /tmp, so I limit processes here to avoid
    # filling up disk space
    sh("parallel -j 8 merge.sh {} ::: " + " ".join(fsnames))

    # One final output
    gffnames = glob("*.all.gff")
    assert len(gffnames) == len(fsnames)
    # Again, DO NOT USE gff3_merge to merge with a smallist /tmp/ area
    print >> sys.stderr, "Use the following command to merge:\n" \
                         "$ python -m jcvi.formats.gff merge gfflist"


def longest(args):
    """
    %prog longest pile.txt cds.fasta

    Pick the longest model per group. `pile.txt` can be generated by
    formats.bed.pile().
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(longest.__doc__)
    p.add_option("--samesize", default=False, action="store_true",
                 help="Only report where the group has same size "\
                      "[default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pilefile, cdsfasta = args
    sizes = Sizes(cdsfasta).mapping

    fp = open(pilefile)
    fw = open("Problems.ids", "w")
    for row in fp:
        models, = row.split()
        all_models = models.split("|")
        all_lengths = [(x, sizes[x]) for x in all_models]
        max_model = max(all_lengths, key=lambda x: x[-1])[0]

        if opts.samesize:
            mms, lengths = zip(*all_lengths)
            if len(set(lengths)) != 1:
                continue

        modelmsg = "|".join("{0}({1})".format(a, b) for a, b in all_lengths)
        print "\t".join((max_model, modelmsg))

        problems = [x for x in all_models if x != max_model]
        print >> fw, "\n".join(problems)

    fw.close()


def batcheval(args):
    """
    %prog batcheval model.ids gff_file evidences.bed fastafile

    Get the accuracy for a list of models against evidences in the range of the
    genes. For example:

    $ %prog batcheval all.gff3 isoforms.ids proteins.bed scaffolds.fasta

    Outfile contains the scores for the models can be found in models.scores
    """
    from jcvi.formats.bed import evaluate
    from jcvi.formats.gff import make_index

    p = OptionParser(evaluate.__doc__)
    p.add_option("--type", default="CDS",
            help="list of features to extract, use comma to separate (e.g."
            "'five_prime_UTR,CDS,three_prime_UTR') [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    model_ids, gff_file, evidences_bed, fastafile = args
    type = set(opts.type.split(","))

    g = make_index(gff_file)
    fp = open(model_ids)
    prefix = model_ids.rsplit(".", 1)[0]
    fwscores = open(prefix + ".scores", "w")

    for row in fp:
        cid = row.strip()
        b = g.parents(cid, 1).next()
        query = "{0}:{1}-{2}".format(b.chrom, b.start, b.stop)
        children = [c for c in g.children(cid, 1)]

        cidbed = prefix + ".bed"
        fw = open(cidbed, "w")
        for c in children:
            if c.featuretype not in type:
                continue

            fw.write(c.to_bed())

        fw.close()

        b = evaluate([cidbed, evidences_bed, fastafile, "--query={0}".format(query)])
        print >> fwscores, "\t".join((cid, b.score))
        fwscores.flush()


def get_bed_file(gff_file, stype, key):

    from jcvi.formats.gff import bed

    opr = stype.replace(",", "") + ".bed"
    bed_opts = ["--type=" + stype, "--key=" + key]
    bed_file = ".".join((gff_file.split(".")[0], opr))

    if need_update(gff_file, bed_file):
        bed([gff_file, "--outfile={0}".format(bed_file)] + bed_opts)

    return bed_file


def get_splits(split_bed, gff_file, stype, key):
    """
    Use intersectBed to find the fused gene => split genes mappings.
    """
    bed_file = get_bed_file(gff_file, stype, key)
    cmd = "intersectBed -a {0} -b {1} -wao".format(split_bed, bed_file)
    cmd += " | cut -f4,10"
    p = popen(cmd)
    splits = defaultdict(set)
    for row in p:
        a, b = row.split()
        splits[a].add(b)

    return splits


def get_accuracy(query, gff_file, evidences_bed, sizesfile, type, key):
    """
    Get sensitivity, specificity and accuracy given gff_file, and a query range
    that look like "chr1:1-10000".
    """
    from jcvi.formats.bed import evaluate

    bed_file = get_bed_file(gff_file, type, key)
    b = evaluate([bed_file, evidences_bed, sizesfile, "--query={0}".format(query)])

    return b


def split(args):
    """
    %prog split split.bed evidences.bed predictor1.gff predictor2.gff fastafile

    Split MAKER models by checking against predictors (such as AUGUSTUS and
    FGENESH). For each region covered by a working model. Find out the
    combination of predictors that gives the best accuracy against evidences
    (such as PASA).

    `split.bed` can be generated by pulling out subset from a list of ids
    $ python -m jcvi.formats.base join split.ids working.bed
        --column=0,3 --noheader | cut -f2-7 > split.bed
    """
    from jcvi.formats.bed import Bed

    p = OptionParser(split.__doc__)
    p.add_option("--key", default="Name",
            help="Key in the attributes to extract predictor.gff [default: %default]")
    p.add_option("--parents", default="match",
            help="list of features to extract, use comma to separate (e.g."
            "'gene,mRNA') [default: %default]")
    p.add_option("--children", default="match_part",
            help="list of features to extract, use comma to separate (e.g."
            "'five_prime_UTR,CDS,three_prime_UTR') [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 5:
        sys.exit(not p.print_help())

    split_bed, evidences_bed, p1_gff, p2_gff, fastafile = args
    parents = opts.parents
    children = opts.children
    key = opts.key

    bed = Bed(split_bed)

    s1 = get_splits(split_bed, p1_gff, parents, key)
    s2 = get_splits(split_bed, p2_gff, parents, key)

    for b in bed:
        query = "{0}:{1}-{2}".format(b.seqid, b.start, b.end)
        b1 = get_accuracy(query, p1_gff, evidences_bed, fastafile, children, key)
        b2 = get_accuracy(query, p2_gff, evidences_bed, fastafile, children, key)
        accn = b.accn
        c1 = "|".join(s1[accn])
        c2 = "|".join(s2[accn])
        ac1 = b1.accuracy
        ac2 = b2.accuracy
        tag = p1_gff if ac1 >= ac2 else p2_gff
        tag = tag.split(".")[0]

        ac1 = "{0:.3f}".format(ac1)
        ac2 = "{0:.3f}".format(ac2)

        print "\t".join((accn, tag, ac1, ac2, c1, c2))


def datastore(args):
    """
    %prog datastore datastore.log > gfflist.log

    Generate a list of gff filenames to merge. The `datastore.log` file can be
    generated by something like:

    $ find
    /usr/local/scratch/htang/EVM_test/gannotation/maker/1132350111853_default/i1/
    -maxdepth 4 -name "*datastore*.log" > datastore.log
    """
    p = OptionParser(datastore.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    ds, = args
    fp = open(ds)
    for row in fp:
        fn = row.strip()
        assert op.exists(fn)
        pp, logfile = op.split(fn)
        flog = open(fn)
        for row in flog:
            ctg, folder, status = row.split()
            if status != "FINISHED":
                continue

            gff_file = op.join(pp, folder, ctg + ".gff")
            assert op.exists(gff_file)
            print gff_file


if __name__ == '__main__':
    main()
