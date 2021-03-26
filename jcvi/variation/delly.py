#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert delly output to BED format.
"""

import os.path as op
import sys
import logging

from jcvi.formats.base import BaseFile, read_until, must_open
from jcvi.formats.sam import coverage
from jcvi.utils.cbook import percentage
from jcvi.utils.aws import ls_s3, push_to_s3
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update


class DelLine(object):
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
        return "\t".join(
            str(x)
            for x in (
                self.seqid,
                self.start - 1,
                self.end,
                self.accn,
                self.supporting_pairs,
                "+",
            )
        )


class Delly(BaseFile):
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
            print(d.bedline, file=fw)
        logging.debug("File written to `%s`.", bedfile)


def main():

    actions = (
        ("bed", "Convert del.txt to del.bed"),
        ("mito", "Find mito deletions in BAM"),
        ("mitosomatic", "Find mito mosaic somatic mutations in piledriver results"),
        ("mitocompile", "Compile mito deletions from multiple VCF files"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mitosomatic(args):
    """
    %prog mitosomatic t.piledriver

    Find mito mosaic somatic mutations in piledriver results.
    """
    import pandas as pd

    p = OptionParser(mitosomatic.__doc__)
    p.add_option("--minaf", default=0.005, type="float", help="Minimum allele fraction")
    p.add_option("--maxaf", default=0.1, type="float", help="Maximum allele fraction")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (df,) = args
    af_file = df.rsplit(".", 1)[0] + ".af"
    fw = open(af_file, "w")
    df = pd.read_csv(df, sep="\t")
    for i, row in df.iterrows():
        na = row["num_A"]
        nt = row["num_T"]
        nc = row["num_C"]
        ng = row["num_G"]
        nd = row["num_D"]
        ni = row["num_I"]
        depth = row["depth"]
        # major, minor = sorted([na, nt, nc, ng], reverse=True)[:2]
        # af = minor * 1. / (major + minor)
        af = (nd + ni) * 1.0 / depth
        if not (opts.minaf <= af <= opts.maxaf):
            continue
        print("{}\t{}\t{:.6f}".format(row["chrom"], row["start"], af), file=fw)
    fw.close()

    logging.debug("Allele freq written to `{}`".format(af_file))


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

    (delt,) = args
    dt = Delly(delt)
    dt.write_bed("del.bed")


def mitocompile(args):
    """
    %prog mitcompile *.vcf.gz

    Extract information about deletions in vcf file.
    """
    from urllib.parse import parse_qsl
    from jcvi.formats.vcf import VcfLine

    p = OptionParser(mitocompile.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    vcfs = args
    print("\t".join("vcf samplekey depth seqid pos alt svlen pe sr".split()))
    for i, vcf in enumerate(vcfs):
        if (i + 1) % 100 == 0:
            logging.debug("Process `{}` [{}]".format(vcf, percentage(i + 1, len(vcfs))))
        depthfile = vcf.replace(".sv.vcf.gz", ".depth")
        fp = must_open(depthfile)
        _, depth = next(fp).split()
        depth = int(float(depth))
        samplekey = op.basename(vcf).split("_")[0]

        fp = must_open(vcf)
        for row in fp:
            if row[0] == "#":
                continue
            v = VcfLine(row)
            info = dict(parse_qsl(v.info))
            print(
                "\t".join(
                    str(x)
                    for x in (
                        vcf,
                        samplekey,
                        depth,
                        v.seqid,
                        v.pos,
                        v.alt,
                        info.get("SVLEN"),
                        info["PE"],
                        info["SR"],
                    )
                )
            )


def mito(args):
    """
    %prog mito chrM.fa input.bam

    Identify mitochondrial deletions.
    """
    p = OptionParser(mito.__doc__)
    p.set_aws_opts(store="hli-mv-data-science/htang/mito-deletions")
    p.add_option(
        "--realignonly", default=False, action="store_true", help="Realign only"
    )
    p.add_option(
        "--svonly",
        default=False,
        action="store_true",
        help="Run Realign => SV calls only",
    )
    p.add_option(
        "--support", default=1, type="int", help="Minimum number of supporting reads"
    )
    p.set_home("speedseq", default="/mnt/software/speedseq/bin")
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
        computed = [
            op.basename(x).split(".")[0] for x in computed if x.endswith(".depth")
        ]
        remaining_samples = [
            x for x in bamfiles if op.basename(x).split(".")[0] not in computed
        ]

        logging.debug(
            "Already computed on `{}`: {}".format(
                store, len(bamfiles) - len(remaining_samples)
            )
        )
        bamfiles = remaining_samples

    logging.debug("Total samples: {}".format(len(bamfiles)))

    for bamfile in bamfiles:
        run_mito(
            chrMfa,
            bamfile,
            opts,
            realignonly=opts.realignonly,
            svonly=opts.svonly,
            store=store,
            cleanup=cleanup,
        )


def run_mito(
    chrMfa, bamfile, opts, realignonly=False, svonly=False, store=None, cleanup=False
):
    from jcvi.formats.sam import get_minibam

    region = "chrM"
    minibam = op.basename(bamfile).replace(".bam", ".{}.bam".format(region))
    if not op.exists(minibam):
        get_minibam(bamfile, region)
    else:
        logging.debug("{} found. Skipped.".format(minibam))

    speedseq_bin = op.join(opts.speedseq_home, "speedseq")

    realign = minibam.rsplit(".", 1)[0] + ".realign"
    realignbam = realign + ".bam"
    margs = " -v -t {} -o {}".format(opts.cpus, realign)
    if need_update(minibam, realign + ".bam"):
        cmd = speedseq_bin + " realign"
        cmd += margs
        cmd += " {} {}".format(chrMfa, minibam)
        sh(cmd)
    else:
        logging.debug("{} found. Skipped.".format(realignbam))

    if realignonly:
        return

    depthfile = realign + ".depth"
    if need_update(realignbam, depthfile):
        coverage(
            [
                chrMfa,
                realignbam,
                "--nosort",
                "--format=coverage",
                "--outfile={}".format(depthfile),
            ]
        )

    if store:
        push_to_s3(store, depthfile)

    vcffile = realign + ".sv.vcf.gz"
    if need_update(realignbam, vcffile):
        cmd = speedseq_bin + " sv"
        cmd += margs
        cmd += " -R {}".format(chrMfa)
        cmd += " -m {}".format(opts.support)
        cmd += " -B {} -D {} -S {}".format(
            realignbam, realign + ".discordants.bam", realign + ".splitters.bam"
        )
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


if __name__ == "__main__":
    main()
