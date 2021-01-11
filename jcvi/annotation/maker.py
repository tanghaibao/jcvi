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

from collections import Counter, defaultdict

from jcvi.formats.base import BaseFile, LineFile, write_file
from jcvi.apps.grid import GridProcess, get_grid_engine, PBS_STANZA
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    need_update,
    popen,
    sh,
    mkdir,
    glob,
    get_abs_path,
)


class CTLine(object):
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

        s = "=".join(str(x) for x in (tag, value)) if tag else ""
        if s:
            if comment:
                s += "  # " + comment
        else:
            if comment:
                s += "# " + comment
        return s


class CTLFile(LineFile):
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
            print(r, file=fw)
        fw.close()
        logging.debug("File written to `%s`.", filename)


class DatastoreIndexFile(BaseFile):
    def __init__(self, filename):
        super(DatastoreIndexFile, self).__init__(filename)
        scaffold_status = {}
        failed = []

        fp = open(filename)
        for row in fp:
            scaffold, dir, status = row.strip().split("\t")
            scaffold_status[scaffold] = status
        for scaffold, status in scaffold_status.items():
            if status != "FINISHED":
                failed.append(scaffold)

        self.scaffold_status = scaffold_status
        self.failed = failed


def main():

    actions = (
        ("parallel", "partition the genome into parts and run separately"),
        ("merge", "generate the gff files after parallel"),
        ("validate", "validate after MAKER run to check for failures"),
        ("datastore", "generate a list of gff filenames to merge"),
        ("split", "split MAKER models by checking against evidences"),
        ("batcheval", "calls bed.evaluate() in batch"),
        ("longest", "pick the longest model per group"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


arraysh = """
DIR=`awk "NR==$SGE_TASK_ID" {0}`
cd $DIR
{1} --ignore_nfs_tmp"""

arraysh_ua = (
    PBS_STANZA
    + """
cd $PBS_O_WORKDIR
DIR=`awk "NR==$PBS_ARRAY_INDEX" {2}`
cd $DIR
{3} --ignore_nfs_tmp > ../maker.$PBS_ARRAY_INDEX.out 2>&1
"""
)


def parallel(args):
    """
    %prog parallel genome.fasta N

    Partition the genome into parts and run separately. This is useful if MAKER
    is to be run on the grid.
    """
    from jcvi.formats.base import split

    p = OptionParser(parallel.__doc__)
    p.set_home("maker")
    p.set_tmpdir(tmpdir="tmp")
    p.set_grid_opts(array=True)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    genome, NN = args
    threaded = opts.threaded or 1
    tmpdir = opts.tmpdir

    mkdir(tmpdir)
    tmpdir = get_abs_path(tmpdir)

    N = int(NN)
    assert 1 <= N < 1000, "Required: 1 < N < 1000!"

    outdir = "outdir"
    fs = split([genome, outdir, NN])

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
    print("\n".join(dirs), file=fw)
    fw.close()

    # Submit to grid
    ncmds = len(dirs)
    runfile = "array.sh"
    cmd = op.join(opts.maker_home, "bin/maker")
    if tmpdir:
        cmd += " -TMP {0}".format(tmpdir)

    engine = get_grid_engine()
    contents = (
        arraysh.format(jobs, cmd)
        if engine == "SGE"
        else arraysh_ua.format(N, threaded, jobs, cmd)
    )
    write_file(runfile, contents)

    if engine == "PBS":
        return

    # qsub script
    outfile = "maker.\$TASK_ID.out"
    p = GridProcess(
        runfile, outfile=outfile, errfile=outfile, arr=ncmds, grid_opts=opts
    )
    qsubfile = "qsub.sh"
    qsub = p.build()
    write_file(qsubfile, qsub)


mergesh = """
BASE=$1
cd $1{0}/$1.maker.output
{1} -n -d $1_master_datastore_index.log
mv $1.all.gff ../../
"""


def get_fsnames(outdir):
    fnames = glob(op.join(outdir, "*.fa*"))
    suffix = "." + fnames[0].split(".")[-1]
    fsnames = [op.basename(x).rsplit(".", 1)[0] for x in fnames]

    return fsnames, suffix


def merge(args):
    """
    %prog merge outdir output.gff

    Follow-up command after grid jobs are completed after parallel().
    """
    from jcvi.formats.gff import merge as gmerge

    p = OptionParser(merge.__doc__)
    p.set_home("maker")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    outdir, outputgff = args
    fsnames, suffix = get_fsnames(outdir)
    nfs = len(fsnames)
    cmd = op.join(opts.maker_home, "bin/gff3_merge")

    outfile = "merge.sh"
    write_file(outfile, mergesh.format(suffix, cmd))

    # Generate per split directory
    # Note that gff3_merge write to /tmp, so I limit processes here to avoid
    # filling up disk space
    sh("parallel -j 8 merge.sh {} ::: " + " ".join(fsnames))

    # One final output
    gffnames = glob("*.all.gff")
    assert len(gffnames) == nfs

    # Again, DO NOT USE gff3_merge to merge with a smallish /tmp/ area
    gfflist = "gfflist"
    fw = open(gfflist, "w")
    print("\n".join(gffnames), file=fw)
    fw.close()

    nlines = sum(1 for x in open(gfflist))
    assert nlines == nfs  # Be extra, extra careful to include all results
    gmerge([gfflist, "-o", outputgff])
    logging.debug("Merged GFF file written to `{0}`".format(outputgff))


def validate(args):
    """
    %prog validate outdir genome.fasta

    Validate current folder after MAKER run and check for failures. Failed batch
    will be written to a directory for additional work.
    """
    p = OptionParser(validate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    outdir, genome = args
    counter = Counter()

    fsnames, suffix = get_fsnames(outdir)
    dsfile = "{0}{1}/{0}.maker.output/{0}_master_datastore_index.log"
    dslogs = [dsfile.format(x, suffix) for x in fsnames]
    all_failed = []
    for f, d in zip(fsnames, dslogs):
        dslog = DatastoreIndexFile(d)
        counter.update(dslog.scaffold_status.values())
        all_failed.extend([(f, x) for x in dslog.failed])

    cmd = 'tail maker.*.out | grep -c "now finished"'
    n = int(popen(cmd).read())
    assert len(fsnames) == n
    print("ALL jobs have been finished", file=sys.stderr)

    nfailed = len(all_failed)
    if nfailed == 0:
        print("ALL scaffolds are completed with no errors", file=sys.stderr)
        return

    print("Scaffold status:", file=sys.stderr)
    print(counter, file=sys.stderr)
    failed = "FAILED"
    fw = open(failed, "w")
    print("\n".join(["\t".join((f, x)) for f, x in all_failed]), file=fw)
    fw.close()

    nlines = sum(1 for x in open("FAILED"))
    assert nlines == nfailed
    print("FAILED !! {0} instances.".format(nfailed), file=sys.stderr)

    # Rebuild the failed batch
    failed_ids = failed + ".ids"
    failed_fasta = failed + ".fasta"
    cmd = "cut -f2 {0}".format(failed)
    sh(cmd, outfile=failed_ids)
    if need_update((genome, failed_ids), failed_fasta):
        cmd = "faSomeRecords {} {} {}".format(genome, failed_ids, failed_fasta)
        sh(cmd)


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
    p.add_option(
        "--type",
        default="CDS",
        help="list of features to extract, use comma to separate (e.g."
        "'five_prime_UTR,CDS,three_prime_UTR')",
    )
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
        b = next(g.parents(cid, 1))
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
        print("\t".join((cid, b.score)), file=fwscores)
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
    p.add_option(
        "--key",
        default="Name",
        help="Key in the attributes to extract predictor.gff",
    )
    p.add_option(
        "--parents",
        default="match",
        help="list of features to extract, use comma to separate (e.g.'gene,mRNA')",
    )
    p.add_option(
        "--children",
        default="match_part",
        help="list of features to extract, use comma to separate (e.g."
        "'five_prime_UTR,CDS,three_prime_UTR')",
    )
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

        print("\t".join((accn, tag, ac1, ac2, c1, c2)))


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

    (ds,) = args
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
            print(gff_file)


if __name__ == "__main__":
    main()
