#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Impute unknown variations given an input vcf file.
"""

import os.path as op
import logging
import sys

from jcvi.utils.cbook import percentage
from jcvi.apps.grid import MakeManager
from jcvi.formats.base import must_open
from jcvi.formats.vcf import VcfLine
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('beagle', 'use BEAGLE4.1 to impute vcf'),
        ('impute', 'use IMPUTE2 to impute vcf'),
        ('minimac', 'use MINIMAC3 to impute vcf'),
        ('validate', 'validate imputation against withheld variants'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def validate(args):
    """
    %prog validate imputed.vcf withheld.vcf

    Validate imputation against withheld variants.
    """
    p = OptionParser(validate.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    imputed, withheld = args
    register = {}
    fp = open(withheld)
    for row in fp:
        if row[0] == "#":
            continue
        v = VcfLine(row)
        register[(v.seqid, v.pos)] = v.genotype

    logging.debug("Imported {0} records from `{1}`".\
                    format(len(register), withheld))

    fp = must_open(imputed)
    hit = concordant = 0
    seen = set()
    for row in fp:
        if row[0] == "#":
            continue
        chr, pos, genotype = v.seqid, v.pos, v.gentoype
        v = VcfLine(row)
        if (v.seqid, v.pos) in seen:
            continue
        seen.add((chr, pos))
        if (chr, pos) not in register:
            continue
        truth = register[(chr, pos)]
        imputed = genotype.split(":")[0]
        if "|" in imputed:
            probs = [float(x) for x in genotype.split(":")[-1].split(",")]
            imputed = max(zip(probs, ["0/0", "0/1", "1/1"]))[-1]
        hit += 1
        if truth == imputed:
            concordant += 1

    logging.debug("Total concordant: {0}".\
            format(percentage(concordant, hit)))


def minimac(args):
    """
    %prog batchminimac input.txt

    Use MINIMAC3 to impute vcf on all chromosomes.
    """
    p = OptionParser(minimac.__doc__)
    p.set_home("shapeit")
    p.set_home("minimac")
    p.set_outfile()
    p.set_chr()
    p.set_ref()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    txtfile, = args
    ref = opts.ref
    mm = MakeManager()
    pf = txtfile.split(".")[0]
    allrawvcf = []
    alloutvcf = []
    chrs = opts.chr.split(",")
    for x in chrs:
        chrvcf = pf + ".chr{0}.vcf".format(x)
        cmd = "python -m jcvi.formats.vcf from23andme {0} {1}".format(txtfile, x)
        cmd += " --ref {0}".format(ref)
        mm.add(txtfile, chrvcf, cmd)

        chrvcf_hg38 = pf + ".chr{0}.23andme.hg38.vcf".format(x)
        minimac_liftover(mm, chrvcf, chrvcf_hg38, opts)
        allrawvcf.append(chrvcf_hg38)

        if x == "X":
            continue
        elif x in ["Y", "MT"]:
            continue
        else:
            minimac_autosome(mm, x, chrvcf, opts)

        minimacvcf = "{0}.chr{1}.minimac.dose.vcf".format(pf, x)
        minimacvcf_hg38 = "{0}.chr{1}.minimac.hg38.vcf".format(pf, x)
        minimac_liftover(mm, minimacvcf, minimacvcf_hg38, opts)
        alloutvcf.append(minimacvcf)

    rawhg38vcfgz = pf + ".all.23andme.hg38.vcf.gz"
    cmd = "vcf-concat {0} | bgzip > {1}".format(" ".join(allrawvcf), rawhg38vcfgz)
    mm.add(allrawvcf, rawhg38vcfgz, cmd)

    outhg38vcfgz = pf + ".all.minimac.hg38.vcf.gz"
    cmd = "vcf-concat {0} | bgzip > {1}".format(" ".join(alloutvcf), outhg38vcfgz)
    mm.add(alloutvcf, outhg38vcfgz, cmd)

    mm.write()


def minimac_liftover(mm, chrvcf, chrvcf_hg38, opts):
    cmd = "python -m jcvi.formats.vcf liftover {0} {1}/hg19ToHg38.over.chain.gz {2}".\
            format(chrvcf, opts.ref, chrvcf_hg38)
    mm.add(chrvcf, chrvcf_hg38, cmd)


def minimac_X(mm, chr, vcffile, opts):
    minimac_autosome(mm, chr, vcffile, opts)


def minimac_autosome(mm, chr, vcffile, opts):
    pf = vcffile.rsplit(".", 1)[0]
    hapsfile = pf + ".haps"
    kg = op.join(opts.ref, "1000GP_Phase3")
    shapeit_cmd = op.join(opts.shapeit_home, "shapeit")
    cmd = get_phasing_cmd(shapeit_cmd, vcffile, chr, opts.cpus, kg)
    mm.add(vcffile, hapsfile, cmd)

    phasedfile = pf + ".phased.vcf"
    cmd = shapeit_cmd + " -convert --input-haps {0}".format(pf)
    cmd += " --output-vcf {0}".format(phasedfile)
    mm.add(hapsfile, phasedfile, cmd)

    opf = pf + ".minimac"
    minimac_cmd = op.join(opts.minimac_home, "Minimac3-omp")
    cmd = minimac_cmd + " --chr {0} --cpus {1}".format(chr, opts.cpus)
    cmd += " --refHaps {0}/{1}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz".format(kg, chr)
    cmd += " --haps {0} --prefix {1}".format(phasedfile, opf)
    cmd += " --format GT,GP --nobgzip"
    outvcf = opf + ".dose.vcf"
    mm.add(phasedfile, outvcf, cmd)


def beagle(args):
    """
    %prog beagle input.vcf 1

    Use BEAGLE4.1 to impute vcf on chromosome 1.
    """
    p = OptionParser(beagle.__doc__)
    p.set_home("beagle")
    p.set_ref()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    vcffile, chr = args
    pf = vcffile.rsplit(".", 1)[0]
    outpf = pf + ".beagle"
    outfile = outpf + ".vcf.gz"

    mm = MakeManager()
    beagle_cmd = opts.beagle_home
    kg = op.join(opts.ref, "1000GP_Phase3")
    cmd = beagle_cmd + " gt={0}".format(vcffile)
    cmd += " ref={0}/chr{1}.1kg.phase3.v5a.bref".format(kg, chr)
    cmd += " map={0}/plink.chr{1}.GRCh37.map".format(kg, chr)
    cmd += " out={0}".format(outpf)
    cmd += " nthreads=16 gprobs=true"
    mm.add(vcffile, outfile, cmd)

    mm.write()


def get_phasing_cmd(shapeit_cmd, vcffile, chr, cpus, kg):
    rpf = "{0}/1000GP_Phase3_chr{1}".format(kg, chr)
    pf = vcffile.rsplit(".", 1)[0]

    mapfile = "{0}/genetic_map_chr{1}_combined_b37.txt".format(kg, chr)
    cmd = shapeit_cmd + " --input-vcf {0}".format(vcffile)
    cmd += " --input-map {0}".format(mapfile)
    cmd += " --thread {0} --effective-size 20000".format(cpus)
    cmd += " --output-max {0}.haps {0}.sample".format(pf)
    cmd += " --input-ref {0}.hap.gz {0}.legend.gz".format(rpf)
    cmd += " {0}.sample".format(rpf.rsplit("_", 1)[0])
    cmd += " --output-log {0}.log".format(pf)
    return cmd


def impute(args):
    """
    %prog impute input.vcf hs37d5.fa 1

    Use IMPUTE2 to impute vcf on chromosome 1.
    """
    from pyfaidx import Fasta

    p = OptionParser(impute.__doc__)
    p.set_home("shapeit")
    p.set_home("impute")
    p.set_ref()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    vcffile, fastafile, chr = args
    mm = MakeManager()
    pf = vcffile.rsplit(".", 1)[0]
    hapsfile = pf + ".haps"
    shapeit_cmd = op.join(opts.shapeit_home, "shapeit")
    kg = op.join(opts.ref, "1000GP_Phase3")
    cmd = get_phasing_cmd(shapeit_cmd, vcffile, chr, opts.cpus, kg)
    mm.add(vcffile, hapsfile, cmd)

    fasta = Fasta(fastafile)
    size = len(fasta[chr])
    binsize = 5000000
    bins = size / binsize  # 5Mb bins
    if size % binsize:
        bins += 1
    impute_cmd = op.join(opts.impute_home, "impute2")
    chunks = []
    for x in xrange(bins + 1):
        chunk_start = x * binsize + 1
        chunk_end = min(chunk_start + binsize - 1, size)
        outfile = pf + ".chunk{0:02d}.impute2".format(x)
        mapfile = "{0}/genetic_map_chr{1}_combined_b37.txt".format(kg, chr)
        rpf = "{0}/1000GP_Phase3_chr{1}".format(kg, chr)
        cmd = impute_cmd + " -m {0}".format(mapfile)
        cmd += " -known_haps_g {0}".format(hapsfile)
        cmd += " -h {0}.hap.gz -l {0}.legend.gz".format(rpf)
        cmd += " -Ne 20000 -int {0} {1}".format(chunk_start, chunk_end)
        cmd += " -o {0} -allow_large_regions -seed 367946".format(outfile)
        cmd += " && touch {0}".format(outfile)
        mm.add(hapsfile, outfile, cmd)
        chunks.append(outfile)

    # Combine all the files
    imputefile = pf + ".impute2"
    cmd = "cat {0} > {1}".format(" ".join(chunks), imputefile)
    mm.add(chunks, imputefile, cmd)

    # Convert to vcf
    vcffile = pf + ".impute2.vcf"
    cmd = "python -m jcvi.formats.vcf fromimpute2 {0} {1} {2} > {3}".\
                format(imputefile, fastafile, chr, vcffile)
    mm.add(imputefile, vcffile, cmd)
    mm.write()


if __name__ == '__main__':
    main()
