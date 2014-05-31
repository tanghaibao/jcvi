#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Analyze SNPs in resequencing panels.
"""

import sys
import logging

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update


def main():

    actions = (
        ('frommaf', 'convert to four-column tabular format from MAF'),
        ('freq', 'call snp frequencies and keep AO and RO'),
        ('rmdup', 'remove PCR duplicates from BAM files'),
        ('freebayes', 'call snps using freebayes'),
        ('mpileup', 'call snps using samtools-mpileup'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def rmdup(args):
    """
    %prog rmdup *.bam > rmdup.cmds

    Remove PCR duplicates from BAM files, generate a list of commands.
    """
    p = OptionParser(rmdup.__doc__)
    p.add_option("-S", default=False, action="store_true",
                 help="Treat PE reads as SE in rmdup")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    bams = args
    cmd = "samtools rmdup"
    if opts.S:
        cmd += " -S"
    for b in bams:
        if "rmdup" in b:
            continue
        rb = b.rsplit(".", 1)[0] + ".rmdup.bam"
        if not need_update(b, rb):
            continue
        print " ".join((cmd, b, rb))


def mpileup(args):
    """
    %prog mpileup prefix ref.fa *.bam

    Call SNPs using samtools mpileup.
    """
    p = OptionParser(mpileup.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    prefix, ref = args[0:2]
    bams = args[2:]
    cmd = "samtools mpileup -P ILLUMINA -E -ugD -r {0}"
    cmd += " -f {0} {1}".format(ref, " ".join(bams))
    fmd = "bcftools view -cvg -"
    seqids = list(Fasta(ref).iterkeys_ordered())
    for s in seqids:
        outfile = prefix + ".{0}.vcf".format(s)
        print cmd.format(s), "|", fmd, ">", outfile


def freebayes(args):
    """
    %prog freebayes prefix ref.fa *.bam

    Call SNPs using freebayes.
    """
    p = OptionParser(freebayes.__doc__)
    p.add_option("--mindepth", default=3, type="int",
                 help="Minimum depth [default: %default]")
    p.add_option("--minqual", default=20, type="int",
                 help="Minimum quality [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    prefix, ref = args[0:2]
    bams = args[2:]
    cmd = "bamaddrg -R {0}"
    cmd += " " + " ".join("-b {0}".format(x) for x in bams)
    fmd = "freebayes --stdin -C {0} -f {1}".format(opts.mindepth, ref)
    seqids = list(Fasta(ref).iterkeys_ordered())
    for s in seqids:
        outfile = prefix + ".{0}.vcf".format(s)
        print cmd.format(s), "|", fmd + " -r {0} -v {1}".format(s, outfile)


def freq(args):
    """
    %prog freq fastafile bamfile

    Call SNP frequencies and generate GFF file.
    """
    p = OptionParser(freq.__doc__)
    p.add_option("--mindepth", default=3, type="int",
                 help="Minimum depth [default: %default]")
    p.add_option("--minqual", default=20, type="int",
                 help="Minimum quality [default: %default]")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, bamfile = args
    cmd = "freebayes -f {0} --pooled-continuous {1}".format(fastafile, bamfile)
    cmd += " -F 0 -C {0}".format(opts.mindepth)
    cmd += ' | vcffilter -f "QUAL > {0}"'.format(opts.minqual)
    cmd += " | vcfkeepinfo - AO RO TYPE"
    sh(cmd, outfile=opts.outfile)


def frommaf(args):
    """
    %prog frommaf maffile

    Convert to four-column tabular format from MAF.
    """
    p = OptionParser(frommaf.__doc__)
    p.add_option("--validate",
                 help="Validate coordinates against FASTA [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    maf, = args
    snpfile = maf.rsplit(".", 1)[0] + ".vcf"
    fp = open(maf)
    fw = open(snpfile, "w")
    total = 0
    id = "."
    qual = 20
    filter = "PASS"
    info = "DP=20"
    print >> fw, "##fileformat=VCFv4.0"
    print >> fw, "#CHROM POS ID REF ALT QUAL FILTER INFO".replace(" ", "\t")
    for row in fp:
        atoms = row.split()
        c, pos, ref, alt = atoms[:4]
        try:
            c = int(c)
        except:
            continue
        c = "chr{0:02d}".format(c)
        pos = int(pos)
        print >> fw, "\t".join(str(x) for x in \
                (c, pos, id, ref, alt, qual, filter, info))
        total += 1
    fw.close()

    validate = opts.validate
    if not validate:
        return

    from jcvi.utils.cbook import percentage

    f = Fasta(validate)
    fp = open(snpfile)
    nsnps = 0
    for row in fp:
        if row[0] == '#':
            continue

        c, pos, id, ref, alt, qual, filter, info = row.split("\t")
        pos = int(pos)
        feat = dict(chr=c, start=pos, stop=pos)
        s = f.sequence(feat)
        s = str(s)
        assert s == ref, "Validation error: {0} is {1} (expect: {2})".\
                        format(feat, s, ref)
        nsnps += 1
        if nsnps % 50000 == 0:
            logging.debug("SNPs parsed: {0}".format(percentage(nsnps, total)))
    logging.debug("A total of {0} SNPs validated and written to `{1}`.".\
                        format(nsnps, snpfile))


if __name__ == '__main__':
    main()
