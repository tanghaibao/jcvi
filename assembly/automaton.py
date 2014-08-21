#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Automate genome assembly by iterating assembly on a set of files, individually.
This is useful for BAC pool assembly where you want to assemble per pool.
"""

import os
import os.path as op
import sys
import logging

from jcvi.utils.iter import grouper
from jcvi.formats.base import LineFile, write_file
from jcvi.formats.fastq import first, pairspf
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, \
            mkdir, sh, glob, get_abs_path


class Meta (object):

    def __init__(self, fastq, guess=True):
        # Note the guesswork is largely based on JIRA LIMS naming convention
        self.fastq = fastq.strip()
        self.suffix = op.splitext(fastq)[-1]
        if ".1." in fastq or ".2." in fastq:
            paired = ".1" if ".1." in fastq else ".2"
        elif "_R1_" in fastq or "_R2_" in fastq:
            paired = ".1" if "_R1_" in fastq else ".2"
        else:
            paired = ""
        self.paired = paired
        if guess:
            self.guess()

    def __str__(self):
        return "\t".join((self.genome, self.tag, self.fastq))

    @property
    def link(self):
        linkname = "{0}{1}{2}".format(self.tag, self.paired, self.suffix)
        return op.join(self.genome, linkname)

    def make_link(self, firstN=0):
        mkdir(self.genome)
        if firstN > 0:
            first([str(firstN), self.fastq, "--outfile={0}".format(self.link)])
            return

        if op.islink(self.link):
            os.unlink(self.link)
        os.symlink(get_abs_path(self.fastq), self.link)

    def guess(self):
        # Try to guess library info based on file name
        # SUBAC47-MP-IL73-1_CGGAAT_L001_R1_filtered.fastq
        basename = op.basename(self.fastq)
        baseparts = basename.split("-")
        self.genome = baseparts[0]
        self.tag = baseparts[1]

        if self.genome.endswith("BP"):
            self.genome, bp = self.genome[:-5], self.genome[-5:-2]
            self.tag = "-".join((self.tag, bp))  # 500BP


class MetaFile (LineFile):

    def __init__(self, filename):
        super(MetaFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            genome, tag, fastq = row.split()
            m = Meta(fastq, guess=False)
            m.genome, m.tag = genome, tag
            self.append(m)

        self.sort(key=lambda x: (x.genome, x.tag, x.fastq))


def main():

    actions = (
        ('tophat', 'run tophat on a list of inputs'),
        ('cufflinks', 'run cufflinks following tophat'),
        ('prepare', 'parse list of FASTQ files and prepare input'),
        ('pairs', 'estimate insert sizes for input files'),
        ('contamination', 'remove contaminated reads'),
        ('allpaths', 'run automated ALLPATHS'),
        ('spades', 'run automated SPADES assembly'),
        ('allpathsX', 'run automated ALLPATHS on list of files'),
        ('soapX', 'run automated SOAP on list of files'),
        ('correctX', 'run automated ALLPATHS correction on list of files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def spades(args):
    """
    %prog spades folder

    Run automated SPADES.
    """
    from jcvi.formats.fastq import readlen

    p = OptionParser(spades.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    folder, = args
    for p, pf in iter_project(folder, 2):
        rl = readlen([p[0], "--silent"])

        # <http://spades.bioinf.spbau.ru/release3.1.0/manual.html#sec3.4>
        kmers = None
        if rl >= 150:
            kmers = "21,33,55,77"
        elif rl >= 250:
            kmers = "21,33,55,77,99,127"

        cmd = "spades.py"
        if kmers:
            cmd += " -k {0}".format(kmers)
        cmd += " --careful"
        cmd += " --pe1-1 {0} --pe1-2 {1}".format(*p)
        cmd += " -o {0}_spades".format(pf)
        print cmd


def contamination(args):
    """
    %prog contamination folder Ecoli.fasta

    Remove contaminated reads. The FASTQ files in the folder will automatically
    pair and filtered against Ecoli.fasta to remove contaminants using BOWTIE2.
    """
    from jcvi.apps.bowtie import align

    p = OptionParser(contamination.__doc__)
    p.add_option("--mapped", default=False, action="store_true",
                 help="Retain contaminated reads instead [default: %default]")
    p.set_cutoff(cutoff=800)
    p.set_mateorientation(mateorientation="+-")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, ecoli = args
    ecoli = get_abs_path(ecoli)
    tag = "--mapped" if opts.mapped else "--unmapped"
    for p, pf in iter_project(folder, 2):
        align_opts = [ecoli] + p + [tag]
        align_opts += ["--cutoff={0}".format(opts.cutoff), "--null"]
        if opts.mateorientation:
            align_opts += ["--mateorientation={0}".format(opts.mateorientation)]
        samfile, logfile = align(align_opts)


def cufflinks(args):
    """
    %prog cufflinks folder reference

    Run cufflinks on a folder containing tophat results.
    """
    p = OptionParser(cufflinks.__doc__)
    p.add_option("--gtf", help="Reference annotation [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, reference = args
    os.chdir(folder)
    bams = glob("*tophat/accepted_hits.bam")
    for bam in bams:
        pf, ab = op.split(bam)
        outdir = op.join(pf, "cufflinks")
        if op.exists(outdir):
            logging.debug("Directory {0} found. Skipping.".format(outdir))
            continue
        cmd = "cufflinks"
        cmd += " -o {0}".format(outdir)
        cmd += " -p {0}".format(opts.cpus)
        if opts.gtf:
            cmd += " -g {0}".format(opts.gtf)
        cmd += " --frag-bias-correct {0}".format(reference)
        cmd += " --multi-read-correct"
        cmd += " {0}".format(bam)
        sh(cmd)


def tophat(args):
    """
    %prog tophat folder reference

    Run tophat on a folder of reads.
    """
    p = OptionParser(tophat.__doc__)
    p.add_option("--gtf", help="Reference annotation [default: %default]")
    p.add_option("--single", default=False, action="store_true",
                 help="Single end mapping")
    p.add_option("--intron", default=15000, type="int",
                 help="Max intron size [default: %default]")
    p.add_option("--dist", default=-50, type="int",
                 help="Mate inner distance [default: %default]")
    p.add_option("--stdev", default=50, type="int",
                 help="Mate standard deviation [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    num = 1 if opts.single else 2
    folder, reference = args
    for p, prefix in iter_project(folder, num):
        outdir = "{0}_tophat".format(prefix)
        outfile = op.join(outdir, "accepted_hits.bam")
        if op.exists(outfile):
            logging.debug("File `{0}` found. Skipping.".format(outfile))
            continue

        cmd = "tophat -p {0}".format(opts.cpus)
        if opts.gtf:
            cmd += " -G {0}".format(opts.gtf)
        cmd += " -o {0}".format(outdir)

        if num == 1:  # Single-end
            a, = p
            cmd += " {0} {1}".format(reference, a)
        else:  # Paired-end
            a, b = p
            cmd += " --max-intron-length {0}".format(opts.intron)
            cmd += " --mate-inner-dist {0}".format(opts.dist)
            cmd += " --mate-std-dev {0}".format(opts.stdev)
            cmd += " {0} {1} {2}".format(reference, a, b)

        sh(cmd)


def pairs(args):
    """
    %prog pairs folder reference.fasta

    Estimate insert size distribution. Compatible with a variety of aligners,
    including CLC, BOWTIE and BWA.
    """
    p = OptionParser(pairs.__doc__)
    p.set_firstN()
    p.set_mates()
    p.set_aligner()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    cwd = os.getcwd()
    aligner = opts.aligner
    work = "-".join(("pairs", aligner))
    mkdir(work)

    if aligner == "clc":
        from jcvi.apps.clc import align
        from jcvi.formats.cas import pairs as ps
    else:
        from jcvi.formats.sam import pairs as ps

    if aligner == "bowtie":
        from jcvi.apps.bowtie import align
    elif aligner == "bwa":
        from jcvi.apps.bwa import align

    folder, ref = args
    ref = get_abs_path(ref)
    messages = []
    for p, prefix in iter_project(folder, 2):
        samplefq = op.join(work, prefix + ".first.fastq")
        first([str(opts.firstN)] + p + ["-o", samplefq])

        os.chdir(work)
        align_args = [ref, op.basename(samplefq)]
        outfile, logfile = align(align_args)
        bedfile, stats = ps([outfile, "--rclip={0}".format(opts.rclip)])
        os.chdir(cwd)

        median = stats.median
        tag = "MP" if median > 1000 else "PE"
        median = str(median)
        pf, sf = median[:2], median[2:]
        if sf and int(sf) != 0:
            pf = str(int(pf) + 1)  # Get the first two effective digits
        lib = "{0}-{1}".format(tag, pf.ljust(len(median), '0'))
        for i, xp in enumerate(p):
            suffix = "fastq.gz" if xp.endswith(".gz") else "fastq"
            link = "{0}-{1}.{2}.{3}".format(lib, prefix.replace("-", ""),
                                            i + 1, suffix)
            m = "\t".join(str(x) for x in (xp, link))
            messages.append(m)

    messages = "\n".join(messages)
    write_file("f.meta", messages, tee=True)


def allpaths(args):
    """
    %prog allpaths folder1 folder2 ...

    Run automated ALLPATHS on list of dirs.
    """
    p = OptionParser(allpaths.__doc__)
    p.add_option("--ploidy", default="1", choices=("1", "2"),
                 help="Ploidy [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(not p.print_help())

    folders = args
    for pf in folders:
        if not op.isdir(pf):
            continue
        assemble_dir(pf, target=["final.contigs.fasta", "final.assembly.fasta"],
                     ploidy=opts.ploidy)


def prepare(args):
    """
    %prog prepare jira.txt

    Parse JIRA report and prepare input. Look for all FASTQ files in the report
    and get the prefix. Assign fastq to a folder and a new file name indicating
    the library type (e.g. PE-500, MP-5000, etc.).

    Note that JIRA report can also be a list of FASTQ files.
    """
    p = OptionParser(prepare.__doc__)
    p.add_option("--first", default=0, type="int",
                 help="Use only first N reads [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    jfile, = args
    metafile = jfile + ".meta"

    if need_update(jfile, metafile):
        fp = open(jfile)
        fastqfiles = [x.strip() for x in fp if ".fastq" in x]
        metas = [Meta(x) for x in fastqfiles]

        fw = open(metafile, "w")
        print >> fw, "\n".join(str(x) for x in metas)
        print >> sys.stderr, "Now modify `{0}`, and restart this script.".\
                            format(metafile)
        print >> sys.stderr, "Each line is : genome library fastqfile"
        fw.close()
        return

    mf = MetaFile(metafile)
    for m in mf:
        m.make_link(firstN=opts.first)


def slink(p, pf, tag, extra=None):

    mkdir(pf, overwrite=True)
    cwd = os.getcwd()
    os.chdir(pf)

    # Create sym-links for the input files
    i = 1
    for f in sorted(p):
        gz = ".gz" if f.endswith(".gz") else ""
        if "PE-0" in f:
            sh("ln -sf ../{0} PE-0.fastq{1}".format(f, gz))
            continue
        for t in tag:
            sh("ln -sf ../{0} {1}.{2}.fastq{3}".format(f, t, i, gz))
        i += 1

    if extra:
        for e in extra:
            sh("ln -sf {0}".format(e))

    os.chdir(cwd)


def assemble_pairs(p, pf, tag, target=["final.contigs.fasta"]):
    """
    Take one pair of reads and assemble to contigs.fasta.
    """
    slink(p, pf, tag)
    assemble_dir(pf, target)


def assemble_dir(pf, target, ploidy="1"):
    from jcvi.assembly.allpaths import prepare

    logging.debug("Work on {0}".format(pf))
    asm = [x.replace("final", pf) for x in target]
    if not need_update(pf, asm):
        logging.debug("Assembly found: {0}. Skipped.".format(asm))
        return

    cwd = os.getcwd()
    os.chdir(pf)
    prepare([pf] + sorted(glob("*.fastq") + glob("*.fastq.gz")) + \
            ["--ploidy={0}".format(ploidy)])
    sh("./run.sh")

    for a, t in zip(asm, target):
        sh("cp allpaths/ASSEMBLIES/run/{0} ../{1}".format(t, a))

    logging.debug("Assembly finished: {0}".format(asm))
    os.chdir(cwd)


def correct_pairs(p, pf, tag):
    """
    Take one pair of reads and correct to generate *.corr.fastq.
    """
    from jcvi.assembly.preprocess import correct as cr

    logging.debug("Work on {0} ({1})".format(pf, ','.join(p)))
    itag = tag[0]
    cm = ".".join((pf, itag))
    targets = (cm + ".1.corr.fastq", cm + ".2.corr.fastq", \
                pf + ".PE-0.corr.fastq")
    if not need_update(p, targets):
        logging.debug("Corrected reads found: {0}. Skipped.".format(targets))
        return

    slink(p, pf, tag)

    cwd = os.getcwd()
    os.chdir(pf)
    cr(sorted(glob("*.fastq") + glob("*.fastq.gz")) + ["--nofragsdedup"])
    sh("mv {0}.1.corr.fastq ../{1}".format(itag, targets[0]))
    sh("mv {0}.2.corr.fastq ../{1}".format(itag, targets[1]))
    sh("mv frag_reads_corr.corr.fastq ../{0}".format(targets[2]))

    logging.debug("Correction finished: {0}".format(targets))
    os.chdir(cwd)


def soap_trios(p, pf, tag, extra):
    """
    Take one pair of reads and 'widow' reads after correction and run SOAP.
    """
    from jcvi.assembly.soap import prepare

    logging.debug("Work on {0} ({1})".format(pf, ','.join(p)))
    asm = "{0}.closed.scafSeq".format(pf)
    if not need_update(p, asm):
        logging.debug("Assembly found: {0}. Skipped.".format(asm))
        return

    slink(p, pf, tag, extra)

    cwd = os.getcwd()
    os.chdir(pf)
    prepare(sorted(glob("*.fastq") + glob("*.fastq.gz")) + \
            ["--assemble_1st_rank_only", "-K 31"])
    sh("./run.sh")
    sh("cp asm31.closed.scafSeq ../{0}".format(asm))

    logging.debug("Assembly finished: {0}".format(asm))
    os.chdir(cwd)


def iter_project(folder, n=2):
    # Check for paired reads and extract project id
    filelist = [x for x in glob(folder + "/*.*") \
                    if x.rsplit(".", 1)[-1] in ("fq", "fastq", "txt", "gz")]
    for p in grouper(filelist, n):
        if len(p) != n:
            continue

        pp = [op.basename(x) for x in p]
        pf = pairspf(pp)
        yield list(p), pf


def soapX(args):
    """
    %prog soapX folder tag [*.fastq]

    Run SOAP on a folder of paired reads and apply tag before assembly.
    Optional *.fastq in the argument list will be symlinked in each folder and
    co-assembled.
    """
    p = OptionParser(soapX.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    folder, tag = args[:2]
    extra = args[2:]
    extra = [get_abs_path(x) for x in extra]
    tag = tag.split(",")
    for p, pf in iter_project(folder, n=3):
        soap_trios(p, pf, tag, extra)


def correctX(args):
    """
    %prog correctX folder tag

    Run ALLPATHS correction on a folder of paired reads and apply tag.
    """
    p = OptionParser(correctX.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, tag = args
    tag = tag.split(",")
    for p, pf in iter_project(folder):
        correct_pairs(p, pf, tag)


def allpathsX(args):
    """
    %prog allpathsX folder tag

    Run assembly on a folder of paired reads and apply tag (PE-200, PE-500).
    Allow multiple tags separated by comma, e.g. PE-350,TT-1050
    """
    p = OptionParser(allpathsX.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    folder, tag = args
    tag = tag.split(",")
    for p, pf in iter_project(folder):
        assemble_pairs(p, pf, tag)


if __name__ == '__main__':
    main()
