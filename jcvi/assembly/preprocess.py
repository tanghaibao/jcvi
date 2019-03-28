#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to trim and correct sequence data.
"""
from __future__ import print_function

import os
import os.path as op
import sys
import logging

from jcvi.formats.base import BaseFile, write_file, must_open
from jcvi.formats.fastq import guessoffset
from jcvi.utils.cbook import depends, human_size
from jcvi.apps.base import OptionParser, ActionDispatcher, download, \
            sh, mkdir, need_update, datadir


class FastQCdata (BaseFile, dict):

    def __init__(self, filename, human=False):
        super(FastQCdata, self).__init__(filename)
        if not op.exists(filename):
            logging.debug("File `{0}` not found.".format(filename))
            # Sample_RF37-1/RF37-1_GATCAG_L008_R2_fastqc =>
            # RF37-1_GATCAG_L008_R2
            self["Filename"] = op.basename(\
                    op.split(filename)[0]).rsplit("_", 1)[0]
            self["Total Sequences"] = self["Sequence length"] = \
                self["Total Bases"] = "na"
            return

        fp = open(filename)
        for row in fp:
            atoms = row.rstrip().split("\t")
            if atoms[0] in ("#", ">"):
                continue
            if len(atoms) != 2:
                continue

            a, b = atoms
            self[a] = b

        ts = self["Total Sequences"]
        sl = self["Sequence length"]
        if "-" in sl:
            a, b = sl.split("-")
            sl = (int(a) + int(b)) / 2
            if a == "30":
                sl = int(b)

        ts, sl = int(ts), int(sl)
        tb = ts * sl

        self["Total Sequences"] = human_size(ts).rstrip("b") if human else ts
        self["Total Bases"] = human_size(tb).rstrip("b") if human else tb


def main():

    actions = (
        ('count', 'count reads based on FASTQC results'),
        ('trim', 'trim reads using TRIMMOMATIC'),
        ('correct', 'correct reads using ALLPATHS-LG'),
        ('hetsmooth', 'reduce K-mer diversity using het-smooth'),
        ('alignextend', 'increase read length by extending based on alignments'),
        ('contamination', 'check reads contamination against Ecoli'),
        ('diginorm', 'run K-mer based normalization'),
        ('expand', 'expand sequences using short reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def diginorm(args):
    """
    %prog diginorm fastqfile

    Run K-mer based normalization. Based on tutorial:
    <http://ged.msu.edu/angus/diginorm-2012/tutorial.html>

    Assume input is either an interleaved pairs file, or two separate files.

    To set up khmer:
    $ git clone git://github.com/ged-lab/screed.git
    $ git clone git://github.com/ged-lab/khmer.git
    $ cd screed
    $ python setup.py install
    $ cd ../khmer
    $ make test
    $ export PYTHONPATH=~/export/khmer
    """
    from jcvi.formats.fastq import shuffle, pairinplace, split
    from jcvi.apps.base import getfilesize

    p = OptionParser(diginorm.__doc__)
    p.add_option("--single", default=False, action="store_true",
                 help="Single end reads")
    p.add_option("--tablesize", help="Memory size")
    p.add_option("--npass", default="1", choices=("1", "2"),
                 help="How many passes of normalization")
    p.set_depth(depth=50)
    p.set_home("khmer", default="/usr/local/bin/")
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    if len(args) == 2:
        fastq = shuffle(args + ["--tag"])
    else:
        fastq, = args

    kh = opts.khmer_home
    depth = opts.depth
    PE = not opts.single
    sys.path.insert(0, op.join(kh, "python"))

    pf = fastq.rsplit(".", 1)[0]
    keepfile = fastq + ".keep"
    hashfile = pf + ".kh"
    mints = 10000000
    ts = opts.tablesize or ((getfilesize(fastq) / 16 / mints + 1) * mints)

    norm_cmd = op.join(kh, "normalize-by-median.py")
    filt_cmd = op.join(kh, "filter-abund.py")
    if need_update(fastq, (hashfile, keepfile)):
        cmd = norm_cmd
        cmd += " -C {0} -k 20 -N 4 -x {1}".format(depth, ts)
        if PE:
            cmd += " -p"
        cmd += " -s {0} {1}".format(hashfile, fastq)
        sh(cmd)

    abundfiltfile = keepfile + ".abundfilt"
    if need_update((hashfile, keepfile), abundfiltfile):
        cmd = filt_cmd
        cmd += " {0} {1}".format(hashfile, keepfile)
        sh(cmd)

    if opts.npass == "1":
        seckeepfile = abundfiltfile
    else:
        seckeepfile = abundfiltfile + ".keep"
        if need_update(abundfiltfile, seckeepfile):
            cmd = norm_cmd
            cmd += " -C {0} -k 20 -N 4 -x {1}".format(depth - 10, ts / 2)
            cmd += " {0}".format(abundfiltfile)
            sh(cmd)

    if PE:
        pairsfile = pairinplace([seckeepfile,
                                "--base={0}".format(pf + "_norm"), "--rclip=2"])
        split([pairsfile])


def expand(args):
    """
    %prog expand bes.fasta reads.fastq

    Expand sequences using short reads. Useful, for example for getting BAC-end
    sequences. The template to use, in `bes.fasta` may just contain the junction
    sequences, then align the reads to get the 'flanks' for such sequences.
    """
    import math

    from jcvi.formats.fasta import Fasta, SeqIO
    from jcvi.formats.fastq import readlen, first, fasta
    from jcvi.formats.blast import Blast
    from jcvi.formats.base import FileShredder
    from jcvi.apps.bowtie import align, get_samfile
    from jcvi.apps.align import blast

    p = OptionParser(expand.__doc__)
    p.set_depth(depth=200)
    p.set_firstN()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bes, reads = args
    size = Fasta(bes).totalsize
    rl = readlen([reads])
    expected_size = size + 2 * rl
    nreads = expected_size * opts.depth / rl
    nreads = int(math.ceil(nreads / 1000.)) * 1000

    # Attract reads
    samfile, logfile = align([bes, reads, "--reorder", "--mapped",
           "--firstN={0}".format(opts.firstN)])

    samfile, mapped, _ = get_samfile(reads, bes, bowtie=True, mapped=True)
    logging.debug("Extract first {0} reads from `{1}`.".format(nreads, mapped))

    pf = mapped.split(".")[0]
    pf = pf.split("-")[0]
    bespf = bes.split(".")[0]
    reads = pf + ".expand.fastq"
    first([str(nreads), mapped, "-o", reads])

    # Perform mini-assembly
    fastafile = reads.rsplit(".", 1)[0] + ".fasta"
    qualfile = ""
    if need_update(reads, fastafile):
        fastafile, qualfile = fasta([reads])

    contigs = op.join(pf, "454LargeContigs.fna")
    if need_update(fastafile, contigs):
        cmd = "runAssembly -o {0} -cpu 8 {1}".format(pf, fastafile)
        sh(cmd)
    assert op.exists(contigs)

    # Annotate contigs
    blastfile = blast([bes, contigs])
    mapping = {}
    for query, b in Blast(blastfile).iter_best_hit():
        mapping[query] = b

    f = Fasta(contigs, lazy=True)
    annotatedfasta = ".".join((pf, bespf, "fasta"))
    fw = open(annotatedfasta, "w")
    keys = list(Fasta(bes).iterkeys_ordered())  # keep an ordered list
    recs = []
    for key, v in f.iteritems_ordered():
        vid = v.id
        if vid not in mapping:
            continue
        b = mapping[vid]
        subject = b.subject
        rec = v.reverse_complement() if b.orientation == '-' else v
        rec.id = rid = "_".join((pf, vid, subject))
        rec.description = ""
        recs.append((keys.index(subject), rid, rec))

    recs = [x[-1] for x in sorted(recs)]
    SeqIO.write(recs, fw, "fasta")
    fw.close()

    FileShredder([samfile, logfile, mapped, reads, fastafile, qualfile, blastfile, pf])
    logging.debug("Annotated seqs (n={0}) written to `{1}`.".\
                    format(len(recs), annotatedfasta))

    return annotatedfasta


def contamination(args):
    """
    %prog contamination Ecoli.fasta genome.fasta read.fastq

    Check read contamination on a folder of paired reads. Use bowtie2 to compare
    the reads against:
    1. Ecoli.fsata - this will tell us the lower bound of contamination
    2. genome.fasta - this will tell us the upper bound of contamination
    """
    from jcvi.apps.bowtie import BowtieLogFile, align

    p = OptionParser(contamination.__doc__)
    p.set_firstN()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    ecoli, genome, fq = args
    firstN_opt = "--firstN={0}".format(opts.firstN)
    samfile, logfile = align([ecoli, fq, firstN_opt])
    bl = BowtieLogFile(logfile)
    lowerbound = bl.rate
    samfile, logfile = align([genome, fq, firstN_opt])
    bl = BowtieLogFile(logfile)
    upperbound = 100 - bl.rate

    median = (lowerbound + upperbound) / 2

    clogfile = fq + ".Ecoli"
    fw = open(clogfile, "w")
    lowerbound = "{0:.1f}".format(lowerbound)
    upperbound = "{0:.1f}".format(upperbound)
    median = "{0:.1f}".format(median)

    print("\t".join((fq, lowerbound, median, upperbound)), file=fw)
    print("{0}: Ecoli contamination rate {1}-{2}".\
                        format(fq, lowerbound, upperbound), file=sys.stderr)
    fw.close()


def alignextend(args):
    """
    %prog alignextend ref.fasta read.1.fastq read.2.fastq

    Wrapper around AMOS alignextend.
    """
    choices = "prepare,align,filter,rmdup,genreads".split(",")
    p = OptionParser(alignextend.__doc__)
    p.add_option("--nosuffix", default=False, action="store_true",
                 help="Do not add /1/2 suffix to the read [default: %default]")
    p.add_option("--rc", default=False, action="store_true",
                 help="Reverse complement the reads before alignment")
    p.add_option("--len", default=100, type="int",
                 help="Extend to this length")
    p.add_option("--stage", default="prepare", choices=choices,
                 help="Start from certain stage")
    p.add_option("--dup", default=10, type="int",
                 help="Filter duplicates with coordinates within this distance")
    p.add_option("--maxdiff", default=1, type="int",
                 help="Maximum number of differences")
    p.set_home("amos")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    ref, r1, r2 = args
    pf = op.basename(r1).split(".")[0]
    cmd = op.join(opts.amos_home, "src/Experimental/alignextend.pl")
    if not opts.nosuffix:
        cmd += " -suffix"
    bwa_idx = "{0}.ref.fa.sa".format(pf)
    if not need_update(ref, bwa_idx):
        cmd += " -noindex"
    cmd += " -threads {0}".format(opts.cpus)
    offset = guessoffset([r1])
    if offset == 64:
        cmd += " -I"
    if opts.rc:
        cmd += " -rc"
    cmd += " -allow -len {0} -dup {1}".format(opts.len, opts.dup)
    cmd += " -min {0} -max {1}".format(2 * opts.len, 20 * opts.len)
    cmd += " -maxdiff {0}".format(opts.maxdiff)
    cmd += " -stage {0}".format(opts.stage)
    cmd += " ".join(("", pf, ref, r1, r2))
    sh(cmd)


def count(args):
    """
    %prog count *.gz

    Count reads based on FASTQC results. FASTQC needs to be run on all the input
    data given before running this command.
    """
    from jcvi.utils.table import loadtable, write_csv

    p = OptionParser(count.__doc__)
    p.add_option("--dir",
                help="Sub-directory where FASTQC was run [default: %default]")
    p.add_option("--human", default=False, action="store_true",
                help="Human friendly numbers [default: %default]")
    p.set_table()
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    filenames = args
    subdir = opts.dir
    header = "Filename|Total Sequences|Sequence length|Total Bases".split("|")
    rows = []
    human = opts.human
    for f in filenames:
        folder = f.replace(".gz", "").rsplit(".", 1)[0] + "_fastqc"
        if subdir:
            folder = op.join(subdir, folder)
        summaryfile = op.join(folder, "fastqc_data.txt")

        fqcdata = FastQCdata(summaryfile, human=human)
        row = [fqcdata[x] for x in header]
        rows.append(row)

    print(loadtable(header, rows), file=sys.stderr)
    write_csv(header, rows, sep=opts.sep,
              filename=opts.outfile, align=opts.align)


def hetsmooth(args):
    """
    %prog hetsmooth reads_1.fq reads_2.fq jf-23_0

    Wrapper against het-smooth. Below is the command used in het-smooth manual.

    $ het-smooth --kmer-len=23 --bottom-threshold=38 --top-threshold=220
           --no-multibase-replacements --jellyfish-hash-file=23-mers.jf
               reads_1.fq reads_2.fq
    """
    p = OptionParser(hetsmooth.__doc__)
    p.add_option("-K", default=23, type="int",
                 help="K-mer size [default: %default]")
    p.add_option("-L", type="int",
                 help="Bottom threshold, first min [default: %default]")
    p.add_option("-U", type="int",
                 help="Top threshold, second min [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    reads1fq, reads2fq, jfdb = args
    K = opts.K
    L = opts.L
    U = opts.U

    assert L is not None and U is not None, "Please specify -L and -U"

    cmd = "het-smooth --kmer-len={0}".format(K)
    cmd += " --bottom-threshold={0} --top-threshold={1}".format(L, U)
    cmd += " --no-multibase-replacements --jellyfish-hash-file={0}".format(jfdb)
    cmd += " --no-reads-log"
    cmd += " " + " ".join((reads1fq, reads2fq))

    sh(cmd)


def trim(args):
    """
    %prog trim fastqfiles

    Trim reads using TRIMMOMATIC. If two fastqfiles are given, then it invokes
    the paired reads mode. See manual:

    <http://www.usadellab.org/cms/index.php?page=trimmomatic>
    """
    tv = "0.32"
    TrimJar = "trimmomatic-{0}.jar".format(tv)
    p = OptionParser(trim.__doc__)
    p.add_option("--path", default=op.join("~/bin", TrimJar),
            help="Path to trimmomatic jar file [default: %default]")
    p.set_phred()
    p.add_option("--nofrags", default=False, action="store_true",
            help="Discard frags file in PE mode [default: %default]")
    p.add_option("--minqv", default=15, type="int",
            help="Average qv after trimming [default: %default]")
    p.add_option("--minlen", default=36, type="int",
            help="Minimum length after trimming [default: %default]")
    p.add_option("--adapteronly", default=False, action="store_true",
            help="Only trim adapters with no qv trimming [default: %default]")
    p.add_option("--nogz", default=False, action="store_true",
            help="Do not write to gzipped files [default: %default]")
    p.add_option("--log", default=None, dest="trimlog",
            help="Specify a `trimlog` file [default: %default]")
    p.set_cpus(cpus=4)
    opts, args = p.parse_args(args)

    if len(args) not in (1, 2):
        sys.exit(not p.print_help())

    path = op.expanduser(opts.path)
    url = \
    "http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-{0}.zip"\
    .format(tv)

    if not op.exists(path):
        path = download(url)
        TrimUnzipped = "Trimmomatic-" + tv
        if not op.exists(TrimUnzipped):
            sh("unzip " + path)
        os.remove(path)
        path = op.join(TrimUnzipped, TrimJar)

    assert op.exists(path), \
        "Couldn't find Trimmomatic jar file at `{0}`".\
        format(path)

    adaptersfile = "adapters.fasta"
    Adapters = must_open(op.join(datadir, adaptersfile)).read()
    write_file(adaptersfile, Adapters, skipcheck=True)

    assert op.exists(adaptersfile), \
        "Please place the illumina adapter sequence in `{0}`".\
        format(adaptersfile)

    if opts.phred is None:
        offset = guessoffset([args[0]])
    else:
        offset = int(opts.phred)

    phredflag = " -phred{0}".format(offset)
    threadsflag = " -threads {0}".format(opts.cpus)
    if opts.trimlog:
        trimlog = " -trimlog {0}".format(opts.trimlog)

    cmd = "java -Xmx4g -jar {0}".format(path)
    frags = ".frags.fastq"
    pairs = ".pairs.fastq"
    if not opts.nogz:
        frags += ".gz"
        pairs += ".gz"

    get_prefix = lambda x: op.basename(x).replace(".gz", "").rsplit(".", 1)[0]
    get_dirname = lambda x: "{0}/".format(op.dirname(x)) if op.dirname(x) else ''
    if len(args) == 1:
        cmd += " SE"
        cmd += phredflag
        cmd += threadsflag
        if opts.trimlog:
            cmd += trimlog
        fastqfile, = args
        prefix = get_prefix(fastqfile)
        dirname = get_dirname(fastqfile)
        frags1 = dirname + prefix + frags
        cmd += " {0}".format(" ".join((fastqfile, frags1)))
    else:
        cmd += " PE"
        cmd += phredflag
        cmd += threadsflag
        if opts.trimlog:
            cmd += trimlog
        fastqfile1, fastqfile2 = args
        prefix1 = get_prefix(fastqfile1)
        dirname1 = get_dirname(fastqfile1)
        prefix2 = get_prefix(fastqfile2)
        dirname2 = get_dirname(fastqfile2)
        pairs1 = dirname1 + prefix1 + pairs
        pairs2 = dirname2 + prefix2 + pairs
        frags1 = dirname1 + prefix1 + frags
        frags2 = dirname2 + prefix2 + frags
        if opts.nofrags:
            frags1 = "/dev/null"
            frags2 = "/dev/null"
        cmd += " {0}".format(" ".join((fastqfile1, fastqfile2, \
                pairs1, frags1, pairs2, frags2)))

    cmd += " ILLUMINACLIP:{0}:2:30:10".format(adaptersfile)

    if not opts.adapteronly:
        cmd += " LEADING:3 TRAILING:3"
        cmd += " SLIDINGWINDOW:4:{0}".format(opts.minqv)

    cmd += " MINLEN:{0}".format(opts.minlen)

    if offset != 33:
        cmd += " TOPHRED33"
    sh(cmd)


@depends
def run_RemoveDodgyReads(infile=None, outfile=None, workdir=None,
        removeDuplicates=True, rc=False, nthreads=32):
    # orig.fastb => filt.fastb
    assert op.exists(infile)
    orig = infile.rsplit(".", 1)[0]
    filt = outfile.rsplit(".", 1)[0]

    cmd = "RemoveDodgyReads IN_HEAD={0} OUT_HEAD={1}".format(orig, filt)
    if not removeDuplicates:
        cmd += " REMOVE_DUPLICATES=False"
    if rc:
        cmd += " RC=True"
    cmd += nthreads
    sh(cmd)


@depends
def run_FastbAndQualb2Fastq(infile=None, outfile=None, rc=False):
    corr = op.basename(infile).rsplit(".", 1)[0]
    cmd = "FastbQualbToFastq HEAD_IN={0} HEAD_OUT={0}".format(corr)
    cmd += " PAIRED=False PHRED_OFFSET=33"
    if rc:
        cmd += " FLIP=True"
    sh(cmd)


@depends
def run_pairs(infile=None, outfile=None, suffix=False):
    from jcvi.assembly.allpaths import pairs
    args = infile
    if suffix:
        args.append("--suffix")
    pairs(args)


def correct(args):
    """
    %prog correct *.fastq

    Correct the fastqfile and generated corrected fastqfiles. This calls
    assembly.allpaths.prepare() to generate input files for ALLPATHS-LG. The
    naming convention for your fastqfiles are important, and are listed below.

    By default, this will correct all PE reads, and remove duplicates of all MP
    reads, and results will be placed in `frag_reads.corr.{pairs,frags}.fastq`
    and `jump_reads.corr.{pairs,frags}.fastq`.
    """
    from jcvi.assembly.allpaths import prepare
    from jcvi.assembly.base import FastqNamings

    p = OptionParser(correct.__doc__ + FastqNamings)
    p.add_option("--dir", default="data",
                help="Working directory [default: %default]")
    p.add_option("--fragsdedup", default=False, action="store_true",
                 help="Don't deduplicate the fragment reads [default: %default]")
    p.add_option("--ploidy", default="2", choices=("1", "2"),
                 help="Ploidy [default: %default]")
    p.add_option("--haploidify", default=False, action="store_true",
                 help="Set HAPLOIDIFY=True [default: %default]")
    p.add_option("--suffix", default=False, action="store_true",
                 help="Add suffix /1, /2 to read names")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fastq = args
    tag, tagj, taglj = "frag_reads", "jump_reads", "long_jump_reads"

    ploidy = opts.ploidy
    haploidify = opts.haploidify
    suffix = opts.suffix
    assert (not haploidify) or (haploidify and ploidy == '2')

    prepare(["Unknown"] + fastq + ["--norun"])

    datadir = opts.dir
    mkdir(datadir)
    fullpath = op.join(os.getcwd(), datadir)
    nthreads = " NUM_THREADS={0}".format(opts.cpus)
    phred64 = (guessoffset([args[0]]) == 64)

    orig = datadir + "/{0}_orig".format(tag)
    origfastb = orig + ".fastb"
    if need_update(fastq, origfastb):
        cmd = "PrepareAllPathsInputs.pl DATA_DIR={0} HOSTS='{1}' PLOIDY={2}".\
                format(fullpath, opts.cpus, ploidy)
        if phred64:
            cmd += " PHRED_64=True"
        sh(cmd)

    if op.exists(origfastb):
        correct_frag(datadir, tag, origfastb, nthreads, dedup=opts.fragsdedup,
                     haploidify=haploidify, suffix=suffix)

    origj = datadir + "/{0}_orig".format(tagj)
    origjfastb = origj + ".fastb"
    if op.exists(origjfastb):
        correct_jump(datadir, tagj, origjfastb, nthreads, suffix=suffix)

    origlj = datadir + "/{0}_orig".format(taglj)
    origljfastb = origlj + ".fastb"
    if op.exists(origljfastb):
        correct_jump(datadir, taglj, origljfastb, nthreads, suffix=suffix)


def export_fastq(datadir, corrfastb, rc=False, suffix=False):
    pf = op.basename(corrfastb.rsplit(".", 1)[0])

    cwd = os.getcwd()
    os.chdir(datadir)
    corrfastq = pf + ".fastq"
    run_FastbAndQualb2Fastq(infile=op.basename(corrfastb), \
                            outfile=corrfastq, rc=rc)
    os.chdir(cwd)

    pairsfile = pf + ".pairs"
    fragsfastq = pf + ".corr.fastq"
    run_pairs(infile=[op.join(datadir, pairsfile), op.join(datadir, corrfastq)],
                      outfile=fragsfastq, suffix=suffix)


def correct_frag(datadir, tag, origfastb, nthreads,
                 dedup=False, haploidify=False, suffix=False):
    filt = datadir + "/{0}_filt".format(tag)
    filtfastb = filt + ".fastb"
    run_RemoveDodgyReads(infile=origfastb, outfile=filtfastb,
                         removeDuplicates=dedup, rc=False, nthreads=nthreads)

    filtpairs = filt + ".pairs"
    edit = datadir + "/{0}_edit".format(tag)
    editpairs = edit + ".pairs"
    if need_update(filtpairs, editpairs):
        cmd = "ln -sf {0} {1}.pairs".format(op.basename(filtpairs), edit)
        sh(cmd)

    editfastb = edit + ".fastb"
    if need_update(filtfastb, editfastb):
        cmd = "FindErrors HEAD_IN={0} HEAD_OUT={1}".format(filt, edit)
        cmd += " PLOIDY_FILE=data/ploidy"
        cmd += nthreads
        sh(cmd)

    corr = datadir + "/{0}_corr".format(tag)
    corrfastb = corr + ".fastb"
    if need_update(editfastb, corrfastb):
        cmd = "CleanCorrectedReads DELETE=True"
        cmd += " HEAD_IN={0} HEAD_OUT={1}".format(edit, corr)
        cmd += " PLOIDY_FILE={0}/ploidy".format(datadir)
        if haploidify:
            cmd += " HAPLOIDIFY=True"
        cmd += nthreads
        sh(cmd)

    export_fastq(datadir, corrfastb, suffix=suffix)


def correct_jump(datadir, tagj, origjfastb, nthreads, suffix=False):
    # Pipeline for jump reads does not involve correction
    filt = datadir + "/{0}_filt".format(tagj)
    filtfastb = filt + ".fastb"
    run_RemoveDodgyReads(infile=origjfastb, outfile=filtfastb, \
                         removeDuplicates=True, rc=True, nthreads=nthreads)

    export_fastq(datadir, filtfastb, rc=True, suffix=suffix)


if __name__ == '__main__':
    main()
