#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert delly output to BED format.
"""

import os.path as op
import sys
import logging

from jcvi.formats.base import BaseFile, read_until, must_open
from jcvi.utils.aws import ls_s3, push_to_s3
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update


class DelLine (object):

    def __init__(self, line):
        args = line.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        self.size = int(args[3])
        assert self.size == self.end - self.start + 1
        self.supporting_pairs = int(args[4])
        self.avg_mapping_quality = float(args[5])
        self.accn = args[6]

    @property
    def bedline(self):
        return "\t".join(str(x) for x in (self.seqid,
                          self.start - 1, self.end, self.accn,
                          self.supporting_pairs, "+"))


class Delly (BaseFile):

    def __init__(self, filename):
        super(Delly, self).__init__(filename)

    def __iter__(self):
        fp = must_open(self.filename)
        while True:
            read_until(fp, "-----")
            nextline = fp.readline()
            nextline = fp.readline()
            if not nextline.strip():
                break
            d = DelLine(nextline)
            yield d

    def write_bed(self, bedfile="stdout"):
        fw = must_open(bedfile, "w")
        for d in self:
            print >> fw, d.bedline
        logging.debug("File written to `{0}`.".format(bedfile))


def main():

    actions = (
        ('bed', 'Convert del.txt to del.bed'),
        ('mito', 'Find mito deletions in BAM'),
        ('mitocompile', 'Compile mito deletions from multiple VCF files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed del.txt

    Convert `del.txt` to BED format. DELLY manual here:
    <http://www.embl.de/~rausch/delly.html>

    Deletion:
    chr, start, end, size, #supporting_pairs, avg._mapping_quality, deletion_id
    chr1, 10180, 10509, 329, 75, 15.8667, Deletion_Sample_00000000
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    delt, = args
    dt = Delly(delt)
    dt.write_bed("del.bed")


def mitocompile(args):
    """
    %prog mitcompile *.vcf.gz

    Extract information about deletions in vcf file.
    """
    from jcvi.formats.vcf import VcfLine
    from urlparse import parse_qsl

    p = OptionParser(mitocompile.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    vcfs = args
    print "\t".join("vcf seqid pos alt svlen pe sr".split())
    for vcf in vcfs:
        fp = must_open(vcf)
        for row in fp:
            if row[0] == '#':
                continue
            v = VcfLine(row)
            info = dict(parse_qsl(v.info))
            print "\t".join(str(x) for x in (vcf, v.seqid, v.pos, v.alt,
                        info.get("SVLEN"), info["PE"], info["SR"]))


def mito(args):
    """
    %prog mito chrM.fa input.bam

    Identify mitochondrial deletions.
    """
    p = OptionParser(mito.__doc__)
    p.set_aws_opts(store="hli-mv-data-science/htang/mito-deletions")
    p.add_option("--realignonly", default=False, action="store_true",
                 help="Realign only")
    p.add_option("--svonly", default=False, action="store_true",
                 help="Run Realign => SV calls only")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    chrMfa, bamfile = args
    store = opts.output_path
    cleanup = not opts.nocleanup

    if not op.exists(chrMfa):
        logging.debug("File `{}` missing. Exiting.".format(chrMfa))
        return

    chrMfai = chrMfa + ".fai"
    if not op.exists(chrMfai):
        cmd = "samtools index {}".format(chrMfa)
        sh(cmd)

    if not bamfile.endswith(".bam"):
        bamfiles = [x.strip() for x in open(bamfile)]
    else:
        bamfiles = [bamfile]

    if store:
        computed = ls_s3(store)
        computed = [op.basename(x).split('.')[0] for x in computed if \
                        x.endswith(".sv.vcf.gz")]
        remaining_samples = [x for x in bamfiles \
                    if op.basename(x).split(".")[0] not in computed]

        logging.debug("Already computed on `{}`: {}".\
                        format(store, len(bamfiles) - len(remaining_samples)))
        bamfiles = remaining_samples

    logging.debug("Total samples: {}".format(len(bamfiles)))

    for bamfile in bamfiles:
        run_mito(chrMfa, bamfile, opts,
                 realignonly=opts.realignonly,
                 svonly=opts.svonly,
                 store=store, cleanup=cleanup)


def run_mito(chrMfa, bamfile, opts, realignonly=False, svonly=False,
             store=None, cleanup=False):
    from jcvi.formats.sam import get_minibam
    region = "chrM"
    minibam = op.basename(bamfile).replace(".bam", ".{}.bam".format(region))
    if not op.exists(minibam):
        get_minibam(bamfile, region)
    else:
        logging.debug("{} found. Skipped.".format(minibam))

    realign = minibam.rsplit(".", 1)[0] + ".realign"
    realignbam = realign + ".bam"
    margs = " -v -t {} -o {}".format(opts.cpus, realign)
    if need_update(minibam, realign + ".bam"):
        cmd = "speedseq realign"
        cmd += margs
        cmd += " {} {}".format(chrMfa, minibam)
        sh(cmd)
    else:
        logging.debug("{} found. Skipped.".format(realignbam))

    if realignonly:
        return

    vcffile = realign + ".sv.vcf.gz"
    if need_update(realignbam, vcffile):
        cmd = "speedseq sv"
        cmd += margs
        cmd += " -R {}".format(chrMfa)
        cmd += " -B {} -D {} -S {}".format(realignbam,
                        realign + ".discordants.bam", realign + ".splitters.bam")
        sh(cmd)
    else:
        logging.debug("{} found. Skipped.".format(vcffile))

    if store:
        push_to_s3(store, vcffile)

    if svonly:
        if cleanup:
            do_cleanup(minibam, realignbam)
        return

    piledriver = realign + ".piledriver"
    if need_update(realignbam, piledriver):
        cmd = "bamtools piledriver -fasta {}".format(chrMfa)
        cmd += " -in {}".format(realignbam)
        sh(cmd, outfile=piledriver)

    if store:
        push_to_s3(store, piledriver)

    if cleanup:
        do_cleanup(minibam, realignbam)


def do_cleanup(minibam, realignbam):
    sh("rm -f {}* {}*".format(minibam, realignbam))


if __name__ == '__main__':
    main()
