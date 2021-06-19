#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Analyze SNPs in re-sequencing panels.
"""
import sys
import logging

from jcvi.formats.fasta import Fasta
from jcvi.formats.base import is_number, write_file
from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update


def main():

    actions = (
        ("frommaf", "convert to four-column tabular format from MAF"),
        ("freq", "call snp frequencies and keep AO and RO"),
        ("rmdup", "remove PCR duplicates from BAM files"),
        ("freebayes", "call snps using freebayes"),
        ("mpileup", "call snps using samtools-mpileup"),
        ("gatk", "call snps using GATK"),
        ("somatic", "generate series of SPEEDSESQ-somatic commands"),
        ("mappability", "generate 50mer mappability for reference genome"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mappability(args):
    """
    %prog mappability reference.fasta

    Generate 50mer mappability for reference genome. Commands are based on gem
    mapper. See instructions:
    <https://github.com/xuefzhao/Reference.Mappability>
    """
    p = OptionParser(mappability.__doc__)
    p.add_option("--mer", default=50, type="int", help="User mer size")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (ref,) = args
    K = opts.mer
    pf = ref.rsplit(".", 1)[0]
    mm = MakeManager()

    gem = pf + ".gem"
    cmd = "gem-indexer -i {} -o {}".format(ref, pf)
    mm.add(ref, gem, cmd)

    mer = pf + ".{}mer".format(K)
    mapb = mer + ".mappability"
    cmd = "gem-mappability -I {} -l {} -o {} -T {}".format(gem, K, mer, opts.cpus)
    mm.add(gem, mapb, cmd)

    wig = mer + ".wig"
    cmd = "gem-2-wig -I {} -i {} -o {}".format(gem, mapb, mer)
    mm.add(mapb, wig, cmd)

    bw = mer + ".bw"
    cmd = "wigToBigWig {} {}.sizes {}".format(wig, mer, bw)
    mm.add(wig, bw, cmd)

    bg = mer + ".bedGraph"
    cmd = "bigWigToBedGraph {} {}".format(bw, bg)
    mm.add(bw, bg, cmd)

    merged = mer + ".filtered-1.merge.bed"
    cmd = "python -m jcvi.formats.bed filterbedgraph {} 1".format(bg)
    mm.add(bg, merged, cmd)

    mm.write()


def gatk(args):
    """
    %prog gatk bamfile reference.fasta

    Call SNPs based on GATK best practices.
    """
    p = OptionParser(gatk.__doc__)
    p.add_option(
        "--indelrealign",
        default=False,
        action="store_true",
        help="Perform indel realignment",
    )
    p.set_home("gatk")
    p.set_home("picard")
    p.set_phred()
    p.set_cpus(cpus=24)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, ref = args
    pf = bamfile.rsplit(".", 1)[0]
    mm = MakeManager()
    picard = "java -Xmx32g -jar {0}/picard.jar".format(opts.picard_home)
    tk = "java -Xmx32g -jar {0}/GenomeAnalysisTK.jar".format(opts.gatk_home)
    tk += " -R {0}".format(ref)

    # Step 0 - build reference
    dictfile = ref.rsplit(".", 1)[0] + ".dict"
    cmd1 = picard + " CreateSequenceDictionary"
    cmd1 += " R={0} O={1}".format(ref, dictfile)
    cmd2 = "samtools faidx {0}".format(ref)
    mm.add(ref, dictfile, (cmd1, cmd2))

    # Step 1 - sort bam
    sortedbamfile = pf + ".sorted.bam"
    cmd = picard + " SortSam"
    cmd += " INPUT={0} OUTPUT={1}".format(bamfile, sortedbamfile)
    cmd += " SORT_ORDER=coordinate CREATE_INDEX=true"
    mm.add(bamfile, sortedbamfile, cmd)

    # Step 2 - mark duplicates
    dedupbamfile = pf + ".dedup.bam"
    cmd = picard + " MarkDuplicates"
    cmd += " INPUT={0} OUTPUT={1}".format(sortedbamfile, dedupbamfile)
    cmd += " METRICS_FILE=dedup.log CREATE_INDEX=true"
    mm.add(sortedbamfile, dedupbamfile, cmd)

    if opts.indelrealign:
        # Step 3 - create indel realignment targets
        intervals = pf + ".intervals"
        cmd = tk + " -T RealignerTargetCreator"
        cmd += " -I {0} -o {1}".format(dedupbamfile, intervals)
        mm.add(dedupbamfile, intervals, cmd)

        # Step 4 - indel realignment
        realignedbamfile = pf + ".realigned.bam"
        cmd = tk + " -T IndelRealigner"
        cmd += " -targetIntervals {0}".format(intervals)
        cmd += " -I {0} -o {1}".format(dedupbamfile, realignedbamfile)
        mm.add((dictfile, intervals), realignedbamfile, cmd)
    else:
        realignedbamfile = dedupbamfile

    # Step 5 - SNP calling
    vcf = pf + ".vcf"
    cmd = tk + " -T HaplotypeCaller"
    cmd += " -I {0}".format(realignedbamfile)
    cmd += " --genotyping_mode DISCOVERY"
    cmd += " -stand_emit_conf 10 -stand_call_conf 30"
    cmd += " -nct {0}".format(opts.cpus)
    cmd += " -o {0}".format(vcf)
    if opts.phred == "64":
        cmd += " --fix_misencoded_quality_scores"
    mm.add(realignedbamfile, vcf, cmd)

    # Step 6 - SNP filtering
    filtered_vcf = pf + ".filtered.vcf"
    cmd = tk + " -T VariantFiltration"
    cmd += " -V {0}".format(vcf)
    cmd += ' --filterExpression "DP < 10 || DP > 300 || QD < 2.0 || FS > 60.0 || MQ < 40.0"'
    cmd += ' --filterName "LOWQUAL"'
    cmd += ' --genotypeFilterExpression "isHomVar == 1"'
    cmd += ' --genotypeFilterName "HOMOVAR"'
    cmd += ' --genotypeFilterExpression "isHet == 1"'
    cmd += ' --genotypeFilterName "HET"'
    cmd += " -o {0}".format(filtered_vcf)
    mm.add(vcf, filtered_vcf, cmd)

    mm.write()


def somatic(args):
    """
    %prog somatic ref.fasta *.bam > somatic.sh

    Useful to identify somatic mutations in each sample compared to all other
    samples. Script using SPEEDSEQ-somatic will be written to stdout.
    """
    p = OptionParser(somatic.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 3:
        sys.exit(not p.print_help())

    ref, bams = args[0], args[1:]
    tcmd = "~/export/speedseq/bin/speedseq somatic"
    tcmd += " -t 32 -F .2 -C 3 -q 30"
    cmds = []
    for b in bams:
        pf = b.split(".")[0]
        cmd = tcmd
        cmd += " -o {0}".format(pf)
        others = ",".join(sorted(set(bams) - {b}))
        cmd += " {0} {1} {2}".format(ref, others, b)
        cmds.append(cmd)

    write_file("somatic.sh", "\n".join(cmds))


def rmdup(args):
    """
    %prog rmdup *.bam > rmdup.cmds

    Remove PCR duplicates from BAM files, generate a list of commands.
    """
    p = OptionParser(rmdup.__doc__)
    p.add_option(
        "-S", default=False, action="store_true", help="Treat PE reads as SE in rmdup"
    )
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
        print(" ".join((cmd, b, rb)))


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
        print(cmd.format(s), "|", fmd, ">", outfile)


def freebayes(args):
    """
    %prog freebayes prefix ref.fa *.bam

    Call SNPs using freebayes.
    """
    p = OptionParser(freebayes.__doc__)
    p.add_option("--mindepth", default=3, type="int", help="Minimum depth")
    p.add_option("--minqual", default=20, type="int", help="Minimum quality")
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
        print(cmd.format(s), "|", fmd + " -r {0} -v {1}".format(s, outfile))


def freq(args):
    """
    %prog freq fastafile bamfile

    Call SNP frequencies and generate GFF file.
    """
    p = OptionParser(freq.__doc__)
    p.add_option("--mindepth", default=3, type="int", help="Minimum depth")
    p.add_option("--minqual", default=20, type="int", help="Minimum quality")
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
    p.add_option("--validate", help="Validate coordinates against FASTA")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (maf,) = args
    snpfile = maf.rsplit(".", 1)[0] + ".vcf"
    fp = open(maf)
    fw = open(snpfile, "w")
    total = 0
    id = "."
    qual = 20
    filter = "PASS"
    info = "DP=20"
    print("##fileformat=VCFv4.0", file=fw)
    print("#CHROM POS ID REF ALT QUAL FILTER INFO".replace(" ", "\t"), file=fw)
    for row in fp:
        atoms = row.split()
        c, pos, ref, alt = atoms[:4]
        if is_number(c, int):
            c = int(c)
        else:
            continue
        c = "chr{0:02d}".format(c)
        pos = int(pos)
        print(
            "\t".join(str(x) for x in (c, pos, id, ref, alt, qual, filter, info)),
            file=fw,
        )
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
        if row[0] == "#":
            continue

        c, pos, id, ref, alt, qual, filter, info = row.split("\t")
        pos = int(pos)
        feat = dict(chr=c, start=pos, stop=pos)
        s = f.sequence(feat)
        s = str(s)
        assert s == ref, "Validation error: {0} is {1} (expect: {2})".format(
            feat, s, ref
        )
        nsnps += 1
        if nsnps % 50000 == 0:
            logging.debug("SNPs parsed: {0}".format(percentage(nsnps, total)))
    logging.debug(
        "A total of {0} SNPs validated and written to `{1}`.".format(nsnps, snpfile)
    )


if __name__ == "__main__":
    main()
