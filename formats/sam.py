"""
SAM alignment format. There are other tools that handles better SAM and BAM.
This script simply parses the lines in SAM into human readable fields.

http://samtools.sourceforge.net/SAM1.pdf
"""

import os.path as op
import sys
import logging

from itertools import groupby

from jcvi.formats.base import LineFile
from jcvi.formats.fasta import Fasta
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import fill
from jcvi.assembly.base import Astat
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, \
            mkdir, glob, get_abs_path


class SamLine (object):

    def __init__(self, row):

        args = row.strip().split("\t")
        self.qname = args[0]
        self.flag = args[1]
        self.rname = args[2]
        self.pos = args[3]
        self.mapq = args[4]
        self.cigar = args[5]
        self.mrnm = args[6]
        self.mpos = args[7]
        self.isize = args[8]
        self.seq = args[9]
        self.qual = args[10]

    @property
    def pairline(self):
        qpos = self.cigar.split('H', 1)[0]
        return "%s:%s\t%s:%s" % (self.qname, qpos, self.rname, self.pos)


class Sam (LineFile):

    def __init__(self, filename, callback=None):

        fp = open(filename)
        for row in fp:
            if row[0] == '@':
                continue
            s = SamLine(row)
            if callback:
                callback(s)


def output_bam(cmd, outfile):
    bam = outfile.endswith(".bam")
    if not bam:
        return cmd + " > {0}".format(outfile)

    outcmd, mflag = ("samtools view -bS", "-F 4")
    cmd += " | {0} {1} - > {2}".format(outcmd, mflag, outfile)

    return cmd


class GenomeCoverageLine (object):

    def __init__(self, row):
        args = row.split()
        self.seqid = args[0]
        self.depth = int(args[1])
        self.positions = int(args[2])
        self.length = int(args[3])
        self.freq = float(args[4])


class GenomeCoverageFile (LineFile):

    def __init__(self, filename):
        super(GenomeCoverageFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            self.append(GenomeCoverageLine(row))

    def iter_coverage_seqid(self):
        for seqid, lines in groupby(self, key=lambda x: x.seqid):
            lines = list(lines)
            length = lines[0].length
            counts = 0
            for r in lines:
                counts += r.depth * r.positions
            yield seqid, counts * 1. / length


def get_prefix(readfile, dbfile):
    rdpf = op.basename(readfile).replace(".gz", "").rsplit(".", 1)[0]
    dbpf = op.basename(dbfile).split(".")[0]
    return ".".join((rdpf, dbpf))


def get_samfile(readfile, dbfile, bam=False, mapped=False,
                unmapped=False, bowtie=False):
    prefix = get_prefix(readfile, dbfile)
    ext = ".bam" if bam else ".sam"
    samfile = prefix + ext
    ext = ".fastq" if bowtie else ext
    mapped = (prefix + ".mapped" + ext) if mapped else None
    unmapped = (prefix + ".unmapped" + ext) if unmapped else None
    return samfile, mapped, unmapped


def main():

    actions = (
        ('bed', 'convert bam files to bed'),
        ('pair', 'parse sam file and get pairs'),
        ('pairs', 'print paired-end reads from BAM file'),
        ('chimera', 'parse sam file from `bwasw` and list multi-hit reads'),
        ('ace', 'convert sam file to ace'),
        ('index', 'convert to bam, sort and then index'),
        ('consensus', 'convert bam alignments to consensus FASTA'),
        ('fpkm', 'calculate FPKM values from BAM file'),
        ('coverage', 'calculate depth for BAM file'),
        ('vcf', 'call SNPs on a set of bam files'),
        ('mapped', 'extract mapped/unmapped reads from samfile'),
        ('count', 'count the number of reads mapped using htseq'),
        ('merge', 'merge bam files'),
            )

    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed bedfile bamfiles

    Convert bam files to bed.
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    bedfile = args[0]
    bamfiles = args[1:]
    for bamfile in bamfiles:
        cmd = "bamToBed -i {0}".format(bamfile)
        sh(cmd, outfile=bedfile, append=True)


def merge(args):
    """
    %prog merge merged_bams bams1_dir bams2_dir ...

    Merge BAM files. Treat the bams with the same prefix as a set.
    Output the commands first.
    """
    from jcvi.apps.grid import MakeManager

    p = OptionParser(merge.__doc__)
    p.set_sep(sep="_", help="Separator to group per prefix")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    merged_bams = args[0]
    bamdirs = args[1:]

    mkdir(merged_bams)
    bams = []
    for x in bamdirs:
        bams += glob(op.join(x, "*.bam"))
    bams = [x for x in bams if "nsorted" not in x]

    logging.debug("Found a total of {0} BAM files.".format(len(bams)))

    sep = opts.sep
    key = lambda x: op.basename(x).split(sep)[0]
    bams.sort(key=key)
    mm = MakeManager()
    for prefix, files in groupby(bams, key=key):
        files = sorted(list(files))
        nfiles = len(files)
        source = " ".join(files)
        target = op.join(merged_bams, op.basename(files[0]))
        if nfiles == 1:
            source = get_abs_path(source)
            cmd = "ln -s {0} {1}".format(source, target)
            mm.add("", target, cmd)
        else:
            cmd = "samtools merge {0} {1}".format(target, source)
            mm.add(files, target, cmd, remove=True)
    mm.write()


def count(args):
    """
    %prog count bamfile gtf

    Count the number of reads mapped using `htseq-count`.
    """
    p = OptionParser(count.__doc__)
    p.add_option("--type", default="exon",
                 help="Only count feature type")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, gtf = args
    pf = bamfile.split(".")[0]
    countfile = pf + ".count"
    if not need_update(bamfile, countfile):
        return

    nsorted = pf + "_nsorted"
    nsortedbam, nsortedsam = nsorted + ".bam", nsorted + ".sam"
    if need_update(bamfile, nsortedsam):
        cmd = "samtools sort -n {0} {1}".format(bamfile, nsorted)
        sh(cmd)
        cmd = "samtools view -h {0}".format(nsortedbam)
        sh(cmd, outfile=nsortedsam)

    if need_update(nsortedsam, countfile):
        cmd = "htseq-count --stranded=no --minaqual=10"
        cmd += " -t {0}".format(opts.type)
        cmd += " {0} {1}".format(nsortedsam, gtf)
        sh(cmd, outfile=countfile)


def coverage(args):
    """
    %prog coverage fastafile bamfile

    Calculate coverage for BAM file. BAM file must be sorted.
    """
    p = OptionParser(coverage.__doc__)
    p.add_option("--format", default="bigwig",
                 choices=("bedgraph", "bigwig", "coverage"),
                 help="Output format")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, bamfile = args
    format = opts.format
    pf = bamfile.rsplit(".", 2)[0]
    sizesfile = Sizes(fastafile).filename
    cmd = "genomeCoverageBed -ibam {0} -g {1}".format(bamfile, sizesfile)
    if format in ("bedgraph", "bigwig"):
        cmd += " -bg"
        bedgraphfile = pf + ".bedgraph"
        sh(cmd, outfile=bedgraphfile)

        if format == "bedgraph":
            return bedgraphfile

        bigwigfile = pf + ".bigwig"
        cmd = "bedGraphToBigWig {0} {1} {2}".\
                    format(bedgraphfile, sizesfile, bigwigfile)
        sh(cmd)
        return bigwigfile

    coveragefile = pf + ".coverage"
    if need_update(fastafile, coveragefile):
        sh(cmd, outfile=coveragefile)

    gcf = GenomeCoverageFile(coveragefile)
    for seqid, cov in gcf.iter_coverage_seqid():
        print "\t".join((seqid, "{0:.1f}".format(cov)))


def fpkm(args):
    """
    %prog fpkm fastafile *.bam

    Calculate FPKM values from BAM file.
    """
    p = OptionParser(fpkm.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    fastafile = args[0]
    bamfiles = args[1:]
    # Create a DUMMY gff file for cuffdiff
    gffile = fastafile.rsplit(".", 1)[0] + ".gff"
    if need_update(fastafile, gffile):
        fw = open(gffile, "w")
        f = Fasta(fastafile, lazy=True)
        for key, size in f.itersizes_ordered():
            print >> fw, "\t".join(str(x) for x in (key, "dummy", "transcript",\
                1, size, ".", ".", ".", "ID=" + key))
        fw.close()
        logging.debug("Dummy GFF created: {0}".format(gffile))

    cmd = "cuffdiff {0} {1}".format(gffile, " ".join(bamfiles))
    sh(cmd)


def pairs(args):
    """
    See __doc__ for OptionParser.set_pairs().
    """
    import jcvi.formats.bed

    p = OptionParser(pairs.__doc__)
    p.set_pairs()
    opts, targs = p.parse_args(args)

    if len(targs) != 1:
        sys.exit(not p.print_help())

    samfile, = targs
    bedfile = samfile.rsplit(".", 1)[0] + ".bed"
    if need_update(samfile, bedfile):
        cmd = "bamToBed -i {0}".format(samfile)
        sh(cmd, outfile=bedfile)

    args[args.index(samfile)] = bedfile

    return jcvi.formats.bed.pairs(args)


def consensus(args):
    """
    %prog consensus fastafile bamfile

    Convert bam alignments to consensus FASTQ/FASTA.
    """
    p = OptionParser(consensus.__doc__)
    p.add_option("--fasta", default=False, action="store_true",
            help="Generate consensus FASTA sequences [default: %default]")
    p.add_option("--mask", default=0, type="int",
            help="Mask bases with quality lower than")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    fastafile, bamfile = args
    fasta = opts.fasta
    suffix = "fasta" if fasta else "fastq"
    pf = bamfile.rsplit(".", 1)[0]
    cnsfile = pf + ".cns.{0}".format(suffix)
    vcfgzfile = pf + ".vcf.gz"
    vcf([fastafile, bamfile, "-o", vcfgzfile])
    cmd += "zcat {0} | vcfutils.pl vcf2fq".format(vcfgzfile)
    if fasta:
        cmd += " | seqtk seq -q {0} -A -".format(opts.mask)

    sh(cmd, outfile=cnsfile)


def vcf(args):
    """
    %prog vcf fastafile bamfiles > out.vcf.gz

    Call SNPs on bam files.
    """
    from jcvi.apps.grid import Jobs

    valid_callers = ("mpileup", "freebayes")
    p = OptionParser(vcf.__doc__)
    p.set_outfile(outfile="out.vcf.gz")
    p.add_option("--nosort", default=False, action="store_true",
                 help="Do not sort the BAM files")
    p.add_option("--caller", default="mpileup", choices=valid_callers,
                 help="Use variant caller [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    fastafile = args[0]
    bamfiles = args[1:]
    caller = opts.caller

    unsorted = [x for x in bamfiles if ".sorted." not in x]
    if opts.nosort:
        bamfiles = unsorted
    else:
        jargs = [[[x, "--unique"]] for x in unsorted]
        jobs = Jobs(index, args=jargs)
        jobs.run()
        bamfiles = [x.replace(".sorted.bam", ".bam") for x in bamfiles]
        bamfiles = [x.replace(".bam", ".sorted.bam") for x in bamfiles]

    if caller == "mpileup":
        cmd = "samtools mpileup -E -uf"
        cmd += " {0} {1}".format(fastafile, " ".join(bamfiles))
        cmd += " | bcftools call -vmO v"
    elif caller == "freebayes":
        cmd = "freebayes -f"
        cmd += " {0} {1}".format(fastafile, " ".join(bamfiles))
    sh(cmd, outfile=opts.outfile)


def chimera(args):
    """
    %prog index samfile

    Parse SAM file from `bwasw` and list multi-hit reads.
    """
    p = OptionParser(chimera.__doc__)
    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    samfile, = args
    fp = open(samfile)

    def key_fun(x):
        if x[0] == '@':
            return x[0]
        s = SamLine(x)
        return s.qname

    for read, samlines in groupby(fp, key=key_fun):
        if read == '@':
            continue

        samlines = [SamLine(x) for x in samlines]
        if len(samlines) == 1:
            continue
        print read, [x.rname for x in samlines]


def index(args):
    """
    %prog index samfile/bamfile

    If SAM file, convert to BAM, sort and then index, using SAMTOOLS
    """
    p = OptionParser(index.__doc__)
    p.add_option("--fasta", dest="fasta", default=None,
            help="add @SQ header to the BAM file [default: %default]")
    p.add_option("--unique", default=False, action="store_true",
            help="only retain uniquely mapped reads [default: %default]")

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    samfile, = args
    fastafile = opts.fasta
    if fastafile:
        assert op.exists(fastafile)

    bamfile = samfile.replace(".sam", ".bam")
    if fastafile:
        faifile = fastafile + ".fai"
        if need_update(fastafile, faifile):
            sh("samtools faidx {0}".format(fastafile))
        cmd = "samtools view -bt {0} {1} -F 4 -o {2}".\
                format(faifile, samfile, bamfile)
    else:
        cmd = "samtools view -bS {0} -F 4 -o {1}".\
                format(samfile, bamfile)

    if opts.unique:
        cmd += " -q 1"

    if samfile.endswith(".sam"):
        sh(cmd)

    # Already sorted?
    if bamfile.endswith(".sorted.bam"):
        sortedbamfile = bamfile
    else:
        prefix = bamfile.replace(".bam", "")
        sortedbamfile = prefix + ".sorted.bam"

    if need_update(bamfile, sortedbamfile):
        sh("samtools sort {0} {1}.sorted".format(bamfile, prefix))

    baifile = sortedbamfile + ".bai"
    if need_update(sortedbamfile, baifile):
        sh("samtools index {0}".format(sortedbamfile))

    return sortedbamfile


def mapped(args):
    """
    %prog mapped sam/bamfile

    Given an input sam/bam file, output a sam/bam file containing only the mapped reads.
    Optionally, extract the unmapped reads into a separate file
    """
    import pysam
    from jcvi.apps.grid import Jobs

    p = OptionParser(mapped.__doc__)
    p.set_sam_options(extra=False)

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    samfile, = args

    view_opts = []
    oext, mopts = (".sam", ["-S"]) \
            if samfile.endswith(".sam") else (".bam", [])

    flag, ext = ("-b", ".bam") if opts.bam else ("-h", ".sam")
    mopts.append(flag)

    if opts.uniq:
        mopts.append("-q1")
        ext = ".uniq{0}".format(ext)

    if opts.unmapped:
        uopts = [x for x in mopts]
        uoutfile = samfile.replace(oext, ".unmapped{0}".format(ext))
        uopts.extend(["-f4", samfile, "-o{0}".format(uoutfile)])
        view_opts.append(uopts)

    outfile = samfile.replace(oext, ".mapped{0}".format(ext))
    mopts.extend(["-F4", samfile, "-o{0}".format(outfile)])
    view_opts.append(mopts)

    for vo in view_opts:
        logging.debug('samtools view {0}'.format(" ".join(vo)))

    jobs = Jobs(pysam.view, [(z for z in x) for x in view_opts])
    jobs.run()


def pair(args):
    """
    %prog pair samfile

    Parses the sam file and retrieve in pairs format,
    query:pos ref:pos
    """
    p = OptionParser(pair.__doc__)

    opts, args = p.parse_args(args)
    if len(args) != 1:
        sys.exit(p.print_help())

    def callback(s):
        print s.pairline
    Sam(args[0], callback=callback)


def cigar_to_seq(a, gap='*'):
    """
    Accepts a pysam row.

    cigar alignment is presented as a list of tuples (operation,length). For
    example, the tuple [ (0,3), (1,5), (0,2) ] refers to an alignment with 3
    matches, 5 insertions and another 2 matches.

    Op BAM Description
    M 0 alignment match (can be a sequence match or mismatch)
    I 1 insertion to the reference
    D 2 deletion from the reference
    N 3 skipped region from the reference
    S 4 soft clipping (clipped sequences present in SEQ)
    H 5 hard clipping (clipped sequences NOT present in SEQ)
    P 6 padding (silent deletion from padded reference)
    = 7 sequence match
    X 8 sequence mismatch

    convert the sequence based on the cigar string. For example:
    """
    seq, cigar = a.seq, a.cigar
    start = 0
    subseqs = []
    npadded = 0
    if cigar is None:
        return None, npadded

    for operation, length in cigar:
        end = start if operation == 2 else start + length

        if operation == 0:  # match
            subseq = seq[start:end]
        elif operation == 1:  # insertion
            subseq = ""
        elif operation == 2:  # deletion
            subseq = gap * length
            npadded += length
        elif operation == 3:  # skipped
            subseq = 'N' * length
        elif operation in (4, 5):  # clip
            subseq = ""
        else:
            raise NotImplementedError

        subseqs.append(subseq)
        start = end

    return "".join(subseqs), npadded


def ace(args):
    """
    %prog ace bamfile fastafile

    convert bam format to ace format. This often allows the remapping to be
    assessed as a denovo assembly format. bam file needs to be indexed. also
    creates a .mates file to be used in amos/bambus, and .astat file to mark
    whether the contig is unique or repetitive based on A-statistics in Celera
    assembler.
    """
    p = OptionParser(ace.__doc__)
    p.add_option("--splitdir", dest="splitdir", default="outRoot",
            help="split the ace per contig to dir [default: %default]")
    p.add_option("--unpaired", dest="unpaired", default=False,
            help="remove read pairs on the same contig [default: %default]")
    p.add_option("--minreadno", dest="minreadno", default=3, type="int",
            help="minimum read numbers per contig [default: %default]")
    p.add_option("--minctgsize", dest="minctgsize", default=100, type="int",
            help="minimum contig size per contig [default: %default]")
    p.add_option("--astat", default=False, action="store_true",
            help="create .astat to list repetitiveness [default: %default]")
    p.add_option("--readids", default=False, action="store_true",
            help="create file of mapped and unmapped ids [default: %default]")

    from pysam import Samfile

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, fastafile = args
    astat = opts.astat
    readids = opts.readids

    f = Fasta(fastafile)
    prefix = bamfile.split(".")[0]
    acefile = prefix + ".ace"
    readsfile = prefix + ".reads"
    astatfile = prefix + ".astat"

    logging.debug("Load {0}".format(bamfile))
    s = Samfile(bamfile, "rb")

    ncontigs = s.nreferences
    genomesize = sum(x for a, x in f.itersizes())
    logging.debug("Total {0} contigs with size {1} base".format(ncontigs,
        genomesize))
    qual = "20"  # default qual

    totalreads = sum(s.count(x) for x in s.references)
    logging.debug("Total {0} reads mapped".format(totalreads))

    fw = open(acefile, "w")
    if astat:
        astatfw = open(astatfile, "w")
    if readids:
        readsfw = open(readsfile, "w")

    print >> fw, "AS {0} {1}".format(ncontigs, totalreads)
    print >> fw

    for i, contig in enumerate(s.references):
        cseq = f[contig]
        nbases = len(cseq)

        mapped_reads = [x for x in s.fetch(contig) if not x.is_unmapped]
        nreads = len(mapped_reads)

        nsegments = 0
        print >> fw, "CO {0} {1} {2} {3} U".format(contig, nbases, nreads,
                nsegments)
        print >> fw, fill(str(cseq.seq))
        print >> fw

        if astat:
            astat = Astat(nbases, nreads, genomesize, totalreads)
            print >> astatfw, "{0}\t{1:.1f}".format(contig, astat)

        text = fill([qual] * nbases, delimiter=" ", width=30)
        print >> fw, "BQ\n{0}".format(text)
        print >> fw

        rnames = []
        for a in mapped_reads:
            readname = a.qname
            rname = readname

            if readids:
                print >> readsfw, readname
            rnames.append(rname)

            strand = "C" if a.is_reverse else "U"
            paddedstart = a.pos + 1  # 0-based to 1-based
            af = "AF {0} {1} {2}".format(rname, strand, paddedstart)
            print >> fw, af

        print >> fw

        for a, rname in zip(mapped_reads, rnames):
            aseq, npadded = cigar_to_seq(a)
            if aseq is None:
                continue

            ninfos = 0
            ntags = 0
            alen = len(aseq)
            rd = "RD {0} {1} {2} {3}\n{4}".format(rname, alen, ninfos, ntags,
                    fill(aseq))
            qs = "QA 1 {0} 1 {0}".format(alen)

            print >> fw, rd
            print >> fw
            print >> fw, qs
            print >> fw


if __name__ == '__main__':
    main()
