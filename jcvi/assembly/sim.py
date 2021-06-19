#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Simulate Illumina sequencing reads.
"""
import os
import os.path as op
import random
import shutil
import sys
import logging
import math

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


def main():

    actions = (
        ("wgsim", "sample paired end reads using dwgsim"),
        ("eagle", "simulate Illumina reads using EAGLE"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add_sim_options(p):
    """
    Add options shared by eagle or wgsim.
    """
    p.add_option(
        "--distance",
        default=500,
        type="int",
        help="Outer distance between the two ends",
    )
    p.add_option("--readlen", default=150, type="int", help="Length of the read")
    p.set_depth(depth=10)
    p.set_outfile(outfile=None)


def eagle(args):
    """
    %prog eagle fastafile

    """
    p = OptionParser(eagle.__doc__)
    p.add_option(
        "--share", default="/usr/local/share/EAGLE/", help="Default EAGLE share path"
    )
    add_sim_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    share = opts.share
    depth = opts.depth
    readlen = opts.readlen
    distance = opts.distance
    pf = op.basename(fastafile).split(".")[0]

    # Since EAGLE does not natively support read length other than 100bp and
    # 250bp - for an arbitrary read length we need to generate a bunch of
    # support files

    # First file is the Runinfo
    runinfo_readlen = "RunInfo_PairedReads2x{}Cycles1x1Tiles.xml".format(readlen)
    if not op.exists(runinfo_readlen):
        runinfo = op.join(share, "RunInfo/RunInfo_PairedReads2x251Cycles1x1Tiles.xml")
        runinfo_xml = open(runinfo).read()
        runinfo_xml = (
            runinfo_xml.replace("251", str(readlen))
            .replace("252", str(readlen + 1))
            .replace("502", str(2 * readlen))
        )
        fw = open(runinfo_readlen, "w")
        print(runinfo_xml.strip(), file=fw)
        fw.close()

    # Generate quality profiles
    quality_file1 = "QualityTable.read1.length{}.qval".format(readlen)
    quality_file2 = "QualityTable.read2.length{}.qval".format(readlen)
    if not (op.exists(quality_file1) and op.exists(quality_file2)):
        for i, qq in enumerate([quality_file1, quality_file2]):
            cmd = "/usr/local/libexec/EAGLE/scaleQualityTable.pl"
            cmd += " --input {}".format(
                op.join(
                    share,
                    "QualityTables/DefaultQualityTable.read{}.length101.qval".format(
                        i + 1
                    ),
                )
            )
            cmd += " --cycles {}".format(readlen)
            cmd += " --output {}".format(qq)
            sh(cmd, silent=True)

    # Since distance is different from the default distribution which is
    # centered around 319, we shift our peak to the new peak
    template_lengths = op.join(
        share, "TemplateLengthTables/DefaultTemplateLengthTable.tsv"
    )
    template_distance = "TemplateLengthTable{}.tsv".format(distance)
    shift = distance - 319
    if not op.exists(template_distance):
        fp = open(template_lengths)
        fw = open(template_distance, "w")
        for row in fp:
            size, counts = row.split()
            size = int(size)
            counts = int(counts)
            size += shift
            if size < readlen:
                continue
            print("\t".join(str(x) for x in (size, counts)), file=fw)
        fw.close()

    # All done, let's simulate!
    cmd = "configureEAGLE.pl"
    cmd += " --reference-genome {}".format(fastafile)
    cmd += " --coverage-depth {}".format(depth)
    cmd += " --gc-coverage-fit-table {}".format(
        op.join(share, "GcCoverageFitTables/Homo_sapiens.example1.tsv")
    )
    cmd += " --run-info {}".format(runinfo_readlen)
    cmd += " --quality-table {}".format(quality_file1)
    cmd += " --quality-table {}".format(quality_file2)
    cmd += " --template-length-table {}".format(template_distance)
    cmd += " --random-seed {}".format(random.randint(1, 65535))
    sh(cmd, silent=True)

    # Retrieve results
    outpf = opts.outfile or "{0}.{1}bp.{2}x".format(pf, distance, depth)
    outpf += ".bwa"
    cwd = os.getcwd()
    eagle_dir = "EAGLE"
    os.chdir(eagle_dir)
    sh("make bam", silent=True)

    # Convert BAM to FASTQ
    from jcvi.formats.sam import fastq

    a, b = fastq(["eagle.bam", outpf])
    sh("mv {} {} ../".format(a, b))
    os.chdir(cwd)

    # Clean-up
    shutil.rmtree(eagle_dir)


def wgsim(args):
    """
    %prog wgsim fastafile

    Run dwgsim on fastafile.
    """
    p = OptionParser(wgsim.__doc__)
    p.add_option(
        "--erate",
        default=0.01,
        type="float",
        help="Base error rate of the read",
    )
    p.add_option(
        "--noerrors",
        default=False,
        action="store_true",
        help="Simulate reads with no errors",
    )
    p.add_option(
        "--genomesize",
        type="int",
        help="Genome size in Mb [default: estimate from data]",
    )
    add_sim_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    pf = op.basename(fastafile).split(".")[0]

    genomesize = opts.genomesize
    size = genomesize * 1000000 if genomesize else Fasta(fastafile).totalsize
    depth = opts.depth
    readlen = opts.readlen
    readnum = int(math.ceil(size * depth / (2 * readlen)))

    distance = opts.distance
    stdev = distance / 10

    outpf = opts.outfile or "{0}.{1}bp.{2}x".format(pf, distance, depth)

    logging.debug("Total genome size: {0} bp".format(size))
    logging.debug("Target depth: {0}x".format(depth))
    logging.debug("Number of read pairs (2x{0}): {1}".format(readlen, readnum))

    if opts.noerrors:
        opts.erate = 0

    cmd = "dwgsim -e {0} -E {0}".format(opts.erate)
    if opts.noerrors:
        cmd += " -r 0 -R 0 -X 0 -y 0"

    cmd += " -d {0} -s {1}".format(distance, stdev)
    cmd += " -N {0} -1 {1} -2 {1}".format(readnum, readlen)
    cmd += " {0} {1}".format(fastafile, outpf)
    sh(cmd)


if __name__ == "__main__":
    main()
