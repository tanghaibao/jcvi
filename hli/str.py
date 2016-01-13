#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Identify repeat numbers in STR repeats.
"""

import re
import os.path as op
import sys
import logging

import pyfasta

from math import log, ceil
from collections import defaultdict
from jcvi.graphics.histogram import stem_leaf_plot
from jcvi.utils.cbook import percentage
from jcvi.formats.bed import natsorted
from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, sh


READLEN = 150
MINSCORE = 30
YSEARCH_HAPLOTYPE = """
DYS393  DYS390 DYS19/DYS394  DYS19b        DYS391        DYS385a       DYS385b DYS426  DYS388  DYS439
DYS389I DYS392 DYS389B       DYS458        DYS459a/b     DYS459a/b     DYS455  DYS454  DYS447  DYS437
DYS448  DYS449 DYS464a/b/c/d DYS464a/b/c/d DYS464a/b/c/d DYS464a/b/c/d DYS464e DYS464f DYS464g DYS460
GATA-H4 YCAIIa YCAIIb        DYS456        DYS607        DYS576        DYS570  CDYa    CDYb    DYS442
DYS438  DYS531 DYS578        DYS395S1a/b   DYS395S1a/b   DYS590        DYS537  DYS641  DYS472  DYS406S1
DYS511  DYS425 DYS413a       DYS413b       DYS557        DYS594        DYS436  DYS490  DYS534  DYS450
DYS444  DYS481 DYS520        DYS446        DYS617        DYS568        DYS487  DYS572  DYS640  DYS492
DYS565  DYS461 DYS462        GATA-A10      DYS635        GAAT1B07      DYS441  DYS445  DYS452  DYS463
DYS434  DYS435 DYS485        DYS494        DYS495        DYS505        DYS522  DYS533  DYS549  DYS556
DYS575  DYS589 DYS636        DYS638        DYS643        DYS714        DYS716  DYS717  DYS726  DXYS156-Y
""".split()
YSEARCH_LL = """
L1  L2  L3  L4  L5  L6  L7  L8  L9  L10
L11 L12 L13 L14 L15 L16 L17 L18 L19 L20
L21 L22 L23 L24 L25 L26 L27 L28 L29 L30
L31 L32 L33 L34 L35 L36 L37 L38 L39 L40
L41 L54 L55 L56 L57 L58 L59 L60 L61 L62
L63 L42 L64 L65 L66 L67 L68 L69 L70 L71
L49 L72 L73 L51 L74 L75 L76 L77 L78 L79
L80 L43 L44 L45 L46 L47 L48 L50 L52 L53
L81 L82 L83 L84 L85 L86 L87 L88 L89 L90
L91 L92 L93 L94 L95 L96 L97 L98 L99 L100
""".split()
YHRD_YFILER = """
DYS456 DYS389I DYS390 DYS389B DYS458 DYS19/DYS394 DYS385
DYS393 DYS391 DYS439 DYS635 DYS392 GATA-H4 DYS437 DYS438 DYS448
""".split()
YHRD_YFILERPLUS = """
DYS576 DYS389I DYS635 DYS389B DYS627 DYS460 DYS458 DYS19/DYS394 GATA-H4 DYS448 DYS391
DYS456 DYS390 DYS438 DYS392 DYS518 DYS570 DYS437 DYS385a DYS449
DYS393 DYS439 DYS481 DYF387S1 DYS533
""".split()
USYSTR_ALL = """
DYF387S1 DYS19/DYS394 DYS385 DYS389I
DYS389B DYS390 DYS391 DYS392
DYS393 DYS437 DYS438 DYS439
DYS448 DYS449 DYS456 DYS458
DYS460 DYS481 DYS518 DYS533
DYS549 DYS570 DYS576 DYS627
DYS635 DYS643 GATA-H4
""".split()


class STRLine(object):

    def __init__(self, line):
        args = line.split()
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.period = int(args[3])
        self.copynum = args[4]
        self.consensusSize = int(args[5])
        self.pctmatch = int(args[6])
        self.pctindel = int(args[7])
        self.score = args[8]
        self.A = args[9]
        self.C = args[10]
        self.G = args[11]
        self.T = args[12]
        self.entropy = float(args[13])
        self.motif = args[14]
        assert self.period == len(self.motif)
        self.name = args[15] if len(args) == 16 else None

    def __str__(self):
        fields = [self.seqid, self.start, self.end,
                  self.period, self.copynum, self.consensusSize,
                  self.pctmatch, self.pctindel, self.score,
                  self.A, self.C, self.G, self.T,
                  "{0:.2f}".format(self.entropy), self.motif]
        if self.name is not None:
            fields += [self.name]
        return "\t".join(str(x) for x in fields)

    def is_valid(self, maxperiod=6, maxlength=READLEN, minscore=MINSCORE):
        return 2 <= self.period <= maxperiod and \
               (self.end - self.start + 1) <= maxlength and \
               self.score >= minscore

    def calc_entropy(self):
        total = self.A + self.C + self.G + self.T
        if total == 0:  # Perhaps they are all Ns - might crash in lobstrindex()
            return 0
        fractions = [x * 1.0 / total for x in [self.A, self.C, self.G, self.T]]
        entropy = sum([-1.0 * x * log(x, 2) for x in fractions if x != 0])
        return entropy

    def iter_exact_str(self, genome):
        pat = re.compile("(({0}){{2,}})".format(self.motif))
        start = self.start
        s = genome[self.seqid][self.start - 1: self.end].upper()
        for m in re.finditer(pat, s):
            self.start = start + m.start()
            length = m.end() - m.start()
            subseq = m.group(0)
            assert length % self.period == 0
            assert subseq.startswith(self.motif)

            self.end = self.start - 1 + length
            self.copynum = length / self.period
            self.pctmatch = 100
            self.pctindel = 0
            self.score = 2 * length
            self.fix_counts(subseq)
            yield self

    def fix_counts(self, subseq):
        length = int(ceil(self.period * float(self.copynum)))
        # Sanity check for length, otherwise lobSTR misses units
        self.end = max(self.end, self.start + length - 1)
        self.A = subseq.count('A')
        self.C = subseq.count('C')
        self.G = subseq.count('G')
        self.T = subseq.count('T')
        self.entropy = self.calc_entropy()


def main():

    actions = (
        ('pe', 'infer paired-end reads spanning a certain region'),
        ('htt', 'extract HTT region and run lobSTR'),
        ('liftover', 'liftOver CODIS/Y-STR markers'),
        ('lobstr', 'run lobSTR on a big BAM'),
        ('lobstrindex', 'make lobSTR index'),
        ('trf', 'run TRF on FASTA files'),
        ('ystr', 'print out Y-STR info given VCF'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def build_ysearch_link(r, ban=["DYS520", "DYS413a", "DYS413b"]):
    template = \
    "http://www.ysearch.org/search_search.asp?fail=2&uid=&freeentry=true&"
    markers = []
    for i, marker in zip(YSEARCH_LL, YSEARCH_HAPLOTYPE):
        z = r.get(marker, "null")
        if "a/b" in marker or marker in ban:
            z = "null"
        m = "{0}={1}".format(i, z)
        markers.append(m)
    print template + "&".join(markers)


def build_yhrd_link(r, panel, ban=["DYS385"]):
    L = []
    for marker in panel:
        z = r.get(marker, "--")
        if marker in ban:
            z = "--"
        L.append(z)
    print " ".join(str(x) for x in L)


def ystr(args):
    """
    %prog ystr chrY.vcf

    Print out Y-STR info given VCF. Marker name extracted from tabfile.
    """
    import vcf
    from jcvi.utils.table import write_csv

    p = OptionParser(ystr.__doc__)
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    vcffile, = args
    db = vcffile.split(".")[-2]
    tabfile = op.join(opts.lobstr_home, "{0}/index.info".format(db))
    reader = vcf.Reader(filename=vcffile)
    register = {}
    fp = open(tabfile)
    for row in fp:
        s = STRLine(row)
        register[(s.seqid, s.start)] = s.name

    header = "Marker|Reads|Ref|Genotype|Motif".split("|")
    contents = []

    simple_register = {}
    for record in reader:
        name = register[(record.CHROM, record.POS)]
        info = record.INFO
        ref = int(float(info["REF"]))
        rpa = info.get("RPA", ref)
        if isinstance(rpa, list):
            rpa = "|".join(str(int(float(x))) for x in rpa)
        ru = info["RU"]
        simple_register[name] = int(rpa)
        for sample in record.samples:
            contents.append((name, sample["ALLREADS"], ref, rpa, ru))

    # Multi-part markers
    a, b, c = "DYS389I", "DYS389B.1", "DYS389B"
    if a in simple_register and b in simple_register:
        simple_register[c] = simple_register[a] + simple_register[b]

    # Multi-copy markers
    mm = ["DYS385", "DYS413", "YCAII"]
    for m in mm:
        ma, mb = m + 'a', m + 'b'
        if ma not in simple_register or mb not in simple_register:
            simple_register[ma] = simple_register[mb] = None
            del simple_register[ma]
            del simple_register[mb]
            continue
        if simple_register[ma] > simple_register[mb]:
            simple_register[ma], simple_register[mb] = \
                    simple_register[mb], simple_register[ma]

    write_csv(header, contents, sep=" ")
    print "[YSEARCH]"
    build_ysearch_link(simple_register)
    print "[YFILER]"
    build_yhrd_link(simple_register, panel=YHRD_YFILER)
    print "[YFILERPLUS]"
    build_yhrd_link(simple_register, panel=YHRD_YFILERPLUS)
    print "[YSTR-ALL]"
    build_yhrd_link(simple_register, panel=USYSTR_ALL)


def get_motif(s, motif_length):
    sl = len(s)
    kmers = set()
    # Get all kmers
    for i in xrange(sl - motif_length):
        ss = s[i: i + motif_length]
        kmers.add(ss)

    kmer_counts = []
    for kmer in kmers:
        kmer_counts.append((s.count(kmer), -s.index(kmer), kmer))

    return sorted(kmer_counts, reverse=True)[0][-1]


def liftover(args):
    """
    %prog liftover lobstr_v3.0.2_hg38_ref.bed hg38.upper.fa

    LiftOver CODIS/Y-STR markers.
    """
    p = OptionParser(liftover.__doc__)
    p.add_option("--checkvalid", default=False, action="store_true",
                help="Do not check minscore, period and length")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    refbed, fastafile = args
    genome = pyfasta.Fasta(fastafile)
    edits = []
    fp = open(refbed)
    for i, row in enumerate(fp):
        s = STRLine(row)
        seq = genome[s.seqid][s.start - 1: s.end].upper()
        s.motif = get_motif(seq, len(s.motif))
        s.fix_counts(seq)
        if opts.checkvalid and not s.is_valid():
            continue
        edits.append(s)
        if i % 10000 == 0:
            print >> sys.stderr, i, "lines read"

    edits = natsorted(edits, key=lambda x: (x.seqid, x.start))
    for e in edits:
        print str(e)


def trf(args):
    """
    %prog trf outdir

    Run TRF on FASTA files.
    """
    from jcvi.apps.base import iglob

    p = OptionParser(trf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    outdir, = args
    mm = MakeManager()

    params = "2 31 31 80 10 {0} 6".format(MINSCORE).split()
    bedfiles = []
    for fastafile in natsorted(iglob(outdir, "*.fa,*.fasta")):
        pf = op.basename(fastafile).split(".")[0]
        cmd1 = "trf {0} {1} -d -h".format(fastafile, " ".join(params))
        datfile = op.basename(fastafile) + "." + ".".join(params) + ".dat"
        bedfile = "{0}.trf.bed".format(pf)
        cmd2 = "cat {0} | awk '($9 >0)' | sed 's/ /\\t/g'".format(datfile)
        cmd2 += " | awk '{{print \"{0}\\t\" $0}}' > {1}".format(pf, bedfile)
        mm.add(fastafile, datfile, cmd1)
        mm.add(datfile, bedfile, cmd2)
        bedfiles.append(bedfile)

    bedfile = "trf.bed"
    cmd = "cat {0} > {1}".format(" ".join(natsorted(bedfiles)), bedfile)
    mm.add(bedfiles, bedfile, cmd)

    mm.write()


def lobstr(args):
    """
    %prog lobstr bamfile lobstr_index

    Run lobSTR on a big BAM file.
    """
    p = OptionParser(lobstr.__doc__)
    p.add_option("--chr", help="Run only this chromosome")
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, lbidx = args
    lhome = opts.lobstr_home

    mm = MakeManager()
    vcffiles = []
    chrs = [opts.chr] or (range(1, 23) + ["X", "Y"])
    for chr in chrs:
        cmd, vcffile = allelotype_on_chr(bamfile, chr, lhome, lbidx)
        mm.add(bamfile, vcffile, cmd)
        vcffiles.append(vcffile)

    gzfile = bamfile.split(".")[0] + ".vcf.gz"
    cmd = "vcf-concat {0} | vcf-sort".format(" ".join(vcffiles))
    cmd += " | bgzip -c > {0}".format(gzfile)
    mm.add(vcffiles, gzfile, cmd)

    tbifile = gzfile + ".tbi"
    cmd = "tabix -p vcf {0}".format(gzfile)
    mm.add(gzfile, tbifile, cmd)

    mm.write()


def allelotype_on_chr(bamfile, chr, lhome, lbidx):
    outfile = "{0}.chr{1}".format(bamfile.split(".")[0], chr)
    cmd = "{0}/bin/allelotype --command classify --bam {1}".format(lhome, bamfile)
    cmd += " --noise_model {0}/models/illumina_v3.pcrfree".format(lhome)
    cmd += " --strinfo {0}/{1}/index.tab".format(lhome, lbidx)
    cmd += " --index-prefix {0}/{1}/lobSTR_".format(lhome, lbidx)
    cmd += " --chrom chr{0} --out {1}.{2}".format(chr, outfile, lbidx)
    cmd += " --max-diff-ref {0}".format(READLEN)
    cmd += " --haploid chrX,chrY"
    return cmd, ".".join((outfile, lbidx, "vcf"))


def htt(args):
    """
    %prog htt bamfile

    Extract HTT region and run lobSTR.
    """
    p = OptionParser(htt.__doc__)
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bamfile, = args
    lhome = opts.lobstr_home

    minibamfile = bamfile.split("/")[-1]
    cmd = "samtools view {0} chr4:3070000-3080000 -b".format(bamfile)
    cmd += " -o {0}".format(minibamfile)
    sh(cmd)

    sh("rm {0}.bai".format(minibamfile))
    sh("samtools index {0}".format(minibamfile))

    cmd = allelotype_on_chr(minibamfile, 4, lhome)
    sh(cmd)


def lobstrindex(args):
    """
    %prog lobstrindex hg38.trf.bed hg38.upper.fa hg38

    Make lobSTR index. Make sure the FASTA contain only upper case (so use
    fasta.format --upper to convert from UCSC fasta). The bed file is generated
    by str().
    """
    p = OptionParser(lobstrindex.__doc__)
    p.add_option("--fixseq", action="store_true", default=False,
                 help="Scan sequences to extract perfect STRs")
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    trfbed, fastafile, pf = args
    lhome = opts.lobstr_home
    mkdir(pf)

    if opts.fixseq:
        genome = pyfasta.Fasta(fastafile)
        newbedfile = trfbed + ".new"
        newbed = open(newbedfile, "w")
        fp = open(trfbed)
        retained = total = 0
        for row in fp:
            s = STRLine(row)
            total += 1
            for ns in s.iter_exact_str(genome):
                if not ns.is_valid():
                    continue
                print >> newbed, ns
                retained += 1
        newbed.close()
        logging.debug("Retained: {0}".format(percentage(retained, total)))
    else:
        newbedfile = trfbed

    mm = MakeManager()
    cmd = "python {0}/scripts/lobstr_index.py".format(lhome)
    cmd += " --str {0} --ref {1} --out {2}".format(newbedfile, fastafile, pf)
    mm.add((newbedfile, fastafile), op.join(pf, "lobSTR_ref.fasta.rsa"), cmd)

    tabfile = "{0}/index.tab".format(pf)
    cmd = "python {0}/scripts/GetSTRInfo.py".format(lhome)
    cmd += " {0} {1} > {2}".format(newbedfile, fastafile, tabfile)
    mm.add((newbedfile, fastafile), tabfile, cmd)

    infofile = "{0}/index.info".format(pf)
    cmd = "cp {0} {1}".format(trfbed, infofile)
    mm.add(trfbed, infofile, cmd)
    mm.write()


def pe(args):
    """
    %prog pe bam chr start end

    Infer distance paired-end reads spanning a certain region.
    """
    import pysam

    p = OptionParser(pe.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    bam, chr, start, end = args
    start, end = int(start), int(end)
    target_len = end - start + 1
    pad = 1000  # How far do we look beyond the target
    pstart = start - pad
    pend = end + pad
    logging.debug("Target length={0}".format(target_len))

    samfile = pysam.AlignmentFile(bam, "rb")
    iter = samfile.fetch(chr, pstart, pend)
    cache = defaultdict(list)
    for x in iter:
        if not x.is_paired:
            continue
        cache[x.query_name].append(x)

    tlens = []
    for name, reads in cache.iteritems():
        if len(reads) < 2:
            continue
        a, b = reads[:2]
        # Get all pairs where read1 is on left flank and read2 is on right flank
        if a.reference_start >= start or b.reference_end <= end:
            continue
        for x in (a, b):
            print x.query_name, x.is_reverse, \
                    x.query_alignment_start, x.query_alignment_end, \
                    x.reference_start, x.reference_end, \
                    x.tlen
        print '=' * 60
        tlen = abs(x.tlen)
        tlens.append(tlen)

    stem_leaf_plot(tlens, 300, 600, 20)


if __name__ == '__main__':
    main()
