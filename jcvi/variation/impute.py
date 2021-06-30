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
from jcvi.formats.vcf import VcfLine, CM
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ("beagle", "use BEAGLE4.1 to impute vcf"),
        ("impute", "use IMPUTE2 to impute vcf"),
        ("minimac", "use MINIMAC3 to impute vcf"),
        ("passthrough", "pass through Y and MT vcf"),
        ("validate", "validate imputation against withheld variants"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def passthrough(args):
    """
    %prog passthrough chrY.vcf chrY.new.vcf

    Pass through Y and MT vcf.
    """
    p = OptionParser(passthrough.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    vcffile, newvcffile = args
    fp = open(vcffile)
    fw = open(newvcffile, "w")
    gg = ["0/0", "0/1", "1/1"]
    for row in fp:
        if row[0] == "#":
            print(row.strip(), file=fw)
            continue

        v = VcfLine(row)
        v.filter = "PASS"
        v.format = "GT:GP"
        probs = [0] * 3
        probs[gg.index(v.genotype)] = 1
        v.genotype = v.genotype.replace("/", "|") + ":{0}".format(
            ",".join("{0:.3f}".format(x) for x in probs)
        )
        print(v, file=fw)
    fw.close()


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

    logging.debug("Imported {0} records from `{1}`".format(len(register), withheld))

    fp = must_open(imputed)
    hit = concordant = 0
    seen = set()
    for row in fp:
        if row[0] == "#":
            continue
        v = VcfLine(row)
        chr, pos, genotype = v.seqid, v.pos, v.genotype
        if (chr, pos) in seen:
            continue
        seen.add((chr, pos))
        if (chr, pos) not in register:
            continue
        truth = register[(chr, pos)]
        imputed = genotype.split(":")[0]
        if "|" in imputed:
            imputed = "/".join(sorted(genotype.split(":")[0].split("|")))
            # probs = [float(x) for x in genotype.split(":")[-1].split(",")]
            # imputed = max(zip(probs, ["0/0", "0/1", "1/1"]))[-1]
        hit += 1
        if truth == imputed:
            concordant += 1
        else:
            print(row.strip(), "truth={0}".format(truth), file=sys.stderr)

    logging.debug("Total concordant: {0}".format(percentage(concordant, hit)))


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

    (txtfile,) = args
    ref = opts.ref
    mm = MakeManager()
    pf = txtfile.split(".")[0]
    allrawvcf = []
    alloutvcf = []
    chrs = opts.chr.split(",")
    for x in chrs:
        px = CM[x]
        chrvcf = pf + ".{0}.vcf".format(px)
        if txtfile.endswith(".vcf"):
            cmd = "vcftools --vcf {0} --chr {1}".format(txtfile, x)
            cmd += " --out {0}.{1} --recode".format(pf, px)
            cmd += " && mv {0}.{1}.recode.vcf {2}".format(pf, px, chrvcf)
        else:  # 23andme
            cmd = "python -m jcvi.formats.vcf from23andme {0} {1}".format(txtfile, x)
            cmd += " --ref {0}".format(ref)
        mm.add(txtfile, chrvcf, cmd)

        chrvcf_hg38 = pf + ".{0}.23andme.hg38.vcf".format(px)
        minimac_liftover(mm, chrvcf, chrvcf_hg38, opts)
        allrawvcf.append(chrvcf_hg38)

        minimacvcf = "{0}.{1}.minimac.dose.vcf".format(pf, px)
        if x == "X":
            minimac_X(mm, x, chrvcf, opts)
        elif x in ["Y", "MT"]:
            cmd = "python -m jcvi.variation.impute passthrough"
            cmd += " {0} {1}".format(chrvcf, minimacvcf)
            mm.add(chrvcf, minimacvcf, cmd)
        else:
            minimac_autosome(mm, x, chrvcf, opts)

        # keep the best line for multi-allelic markers
        uniqvcf = "{0}.{1}.minimac.uniq.vcf".format(pf, px)
        cmd = "python -m jcvi.formats.vcf uniq {0} > {1}".format(minimacvcf, uniqvcf)
        mm.add(minimacvcf, uniqvcf, cmd)

        minimacvcf_hg38 = "{0}.{1}.minimac.hg38.vcf".format(pf, px)
        minimac_liftover(mm, uniqvcf, minimacvcf_hg38, opts)
        alloutvcf.append(minimacvcf_hg38)

    if len(allrawvcf) > 1:
        rawhg38vcfgz = pf + ".all.23andme.hg38.vcf.gz"
        cmd = "vcf-concat {0} | bgzip > {1}".format(" ".join(allrawvcf), rawhg38vcfgz)
        mm.add(allrawvcf, rawhg38vcfgz, cmd)

    if len(alloutvcf) > 1:
        outhg38vcfgz = pf + ".all.minimac.hg38.vcf.gz"
        cmd = "vcf-concat {0} | bgzip > {1}".format(" ".join(alloutvcf), outhg38vcfgz)
        mm.add(alloutvcf, outhg38vcfgz, cmd)

    mm.write()


def minimac_liftover(mm, chrvcf, chrvcf_hg38, opts):
    cmd = "python -m jcvi.formats.vcf liftover {0} {1}/hg19ToHg38.over.chain.gz {2}".format(
        chrvcf, opts.ref, chrvcf_hg38
    )
    mm.add(chrvcf, chrvcf_hg38, cmd)


def minimac_X(mm, chr, vcffile, opts):
    """See details here:
    http://genome.sph.umich.edu/wiki/Minimac3_Cookbook_:_Chromosome_X_Imputation
    """
    pf = vcffile.rsplit(".", 1)[0]
    ranges = [(1, 2699519), (2699520, 154931043), (154931044, 155270560)]
    tags = ["PAR1", "NONPAR", "PAR2"]
    Xvcf = []
    phasedfiles = []
    for tag, (start, end) in zip(tags, ranges):
        recodefile = pf + "_{0}.recode.vcf".format(tag)
        cmd = "vcftools --vcf {0} --out {1}_{2}".format(vcffile, pf, tag)
        cmd += " --chr X --from-bp {0} --to-bp {1} --recode".format(start, end)
        mm.add(vcffile, recodefile, cmd)

        phasedfile = shapeit_phasing(mm, chr + "_{0}".format(tag), recodefile, opts)
        phasedfiles.append(phasedfile)

    pars = [x for x in phasedfiles if "_PAR" in x]
    parfile = pf + "_PAR.recode.phased.vcf"
    nonparfile = pf + "_NONPAR.recode.phased.vcf"
    cmd = "vcf-concat {0} > {1}".format(" ".join(pars), parfile)
    mm.add(pars, parfile, cmd)

    for phasedfile in (parfile, nonparfile):
        outvcf = minimac_autosome(mm, chr, phasedfile, opts, phasing=False)
        Xvcf.append(outvcf)

    minimacvcf = pf + ".minimac.dose.vcf"
    cmd = "vcf-concat {0} | vcf-sort -c > {1}".format(" ".join(Xvcf), minimacvcf)
    mm.add(Xvcf, minimacvcf, cmd)


def minimac_autosome(mm, chr, vcffile, opts, phasing=True):
    pf = vcffile.rsplit(".", 1)[0]
    kg = op.join(opts.ref, "1000GP_Phase3")
    if phasing:
        shapeit_phasing(mm, chr, vcffile, opts)
        phasedfile = pf + ".phased.vcf"
    else:
        phasedfile = vcffile

    chrtag = chr
    if chr == "X":
        chrtag = "X.Non.Pseudo.Auto" if "NONPAR" in vcffile else "X.Pseudo.Auto"

    opf = pf + ".minimac"
    minimac_cmd = op.join(opts.minimac_home, "Minimac3")

    cmd = minimac_cmd + " --chr {0}".format(chr)
    cmd += (
        " --refHaps {0}/{1}.1000g.Phase3.v5.With.Parameter.Estimates.m3vcf.gz".format(
            kg, chrtag
        )
    )
    cmd += " --haps {0} --prefix {1}".format(phasedfile, opf)
    cmd += " --format GT,GP --nobgzip"
    outvcf = opf + ".dose.vcf"
    mm.add(phasedfile, outvcf, cmd)

    return outvcf


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


def shapeit_phasing(mm, chr, vcffile, opts, vcf=True):
    kg = op.join(opts.ref, "1000GP_Phase3")
    shapeit_cmd = op.join(opts.shapeit_home, "shapeit")

    rpf = "{0}/1000GP_Phase3_chr{1}".format(kg, chr)
    pf = vcffile.rsplit(".", 1)[0]
    mapfile = "{0}/genetic_map_chr{1}_combined_b37.txt".format(kg, chr)
    mapfile = mapfile.replace("NONPAR", "nonPAR")

    hapsfile = pf + ".haps"
    cmd = shapeit_cmd + " --input-vcf {0}".format(vcffile)
    cmd += " --input-map {0}".format(mapfile)
    cmd += " --effective-size 11418"
    cmd += " --output-max {0}.haps {0}.sample".format(pf)
    cmd += " --input-ref {0}.hap.gz {0}.legend.gz".format(rpf)
    cmd += " {0}/1000GP_Phase3.sample --output-log {1}.log".format(kg, pf)
    if chr == "X":
        cmd += " --chrX"
    mm.add(vcffile, hapsfile, cmd)

    if not vcf:
        return

    phasedfile = pf + ".phased.vcf"
    cmd = shapeit_cmd + " -convert --input-haps {0}".format(pf)
    cmd += " --output-vcf {0}".format(phasedfile)
    mm.add(hapsfile, phasedfile, cmd)

    return phasedfile


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
    kg = op.join(opts.ref, "1000GP_Phase3")
    shapeit_phasing(mm, chr, vcffile, opts)

    fasta = Fasta(fastafile)
    size = len(fasta[chr])
    binsize = 5000000
    bins = size / binsize  # 5Mb bins
    if size % binsize:
        bins += 1
    impute_cmd = op.join(opts.impute_home, "impute2")
    chunks = []
    for x in range(bins + 1):
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
    cmd = "python -m jcvi.formats.vcf fromimpute2 {0} {1} {2} > {3}".format(
        imputefile, fastafile, chr, vcffile
    )
    mm.add(imputefile, vcffile, cmd)
    mm.write()


if __name__ == "__main__":
    main()
