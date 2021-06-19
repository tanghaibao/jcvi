#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Variant call format.
"""
import os.path as op
import sys
import logging

from collections import defaultdict
from itertools import groupby
from pyfaidx import Fasta
from pyliftover import LiftOver

from jcvi.formats.base import must_open
from jcvi.formats.sizes import Sizes
from jcvi.utils.cbook import percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


class VcfLine:
    def __init__(self, row):
        args = row.strip().split("\t")
        self.seqid = args[0]
        self.pos = int(args[1])
        self.rsid = args[2]
        self.ref = args[3]
        self.alt = args[4]
        self.qual = args[5]
        self.filter = args[6]
        self.info = args[7]
        self.format = args[8]
        self.genotype = args[9]

    def __str__(self):
        return "\t".join(
            str(x)
            for x in (
                self.seqid,
                self.pos,
                self.rsid,
                self.ref,
                self.alt,
                self.qual,
                self.filter,
                self.info,
                self.format,
                self.genotype,
            )
        )


class UniqueLiftover(object):
    def __init__(self, chainfile):
        """
        This object will perform unique single positional liftovers - it will only lift over chromosome positions that
        map unique to the new genome and if the strand hasn't changed.
        Note: You should run a VCF Normalization sweep on all lifted ofer CPRAs to check for variants that need to be
        re-normalized, and to remove variants where the REF now doesn't match after a liftover.
        The combination of these steps will ensure high quality liftovers. However, it should be noted that this won't
        prevent the situation where multiple positions in the old genome pile up uniquely in the new genome, so one
        needs to check for this.
        It's organised as an object rather than a collection of functions  so that the LiftOver chainfile
        only gets opened/passed once and not for every position to be lifted over.
        :param chainfile: A string containing the path to the local UCSC .gzipped chainfile
        :return:
        """

        self.liftover = LiftOver(chainfile)

    def liftover_cpra(self, chromosome, position, verbose=False):
        """
        Given chromosome, position in 1-based co-ordinates,
        This will use pyliftover to liftover a CPRA, will return a (c,p) tuple or raise NonUniqueLiftover if no unique
        and strand maintaining liftover is possible
        :param chromosome: string with the chromosome as it's represented in the from_genome
        :param position: position on chromosome (will be cast to int)
        :param verbose: print verbose information for debugging
        :return: ((str) chromosome, (int) position) or None if no liftover
        """

        chromosome = str(chromosome)
        position = int(position)

        # Perform the liftover lookup, shift the position by 1 as pyliftover deals in 0-based co-ords
        new = self.liftover.convert_coordinate(chromosome, position - 1)
        # This has to be here as new will be NoneType when the chromosome doesn't exist in the chainfile
        if new:
            # If the liftover is unique
            if len(new) == 1:
                # If the liftover hasn't changed strand
                if new[0][2] == "+":
                    # Set the co-ordinates to the lifted-over ones and write out
                    new_chromosome = str(new[0][0])
                    # Shift the position forward by one to convert back to a 1-based co-ords
                    new_position = int(new[0][1]) + 1
                    return new_chromosome, new_position
                else:
                    exception_string = (
                        "{},{} has a flipped strand in liftover: {}".format(
                            chromosome, position, new
                        )
                    )
            else:
                exception_string = "{},{} lifts over to multiple positions: {}".format(
                    chromosome, position, new
                )
        elif new is None:
            exception_string = "Chromosome '{}' provided not in chain file".format(
                chromosome
            )

        if verbose:
            logging.error(exception_string)
        return None, None


CM = dict(
    list(
        zip([str(x) for x in range(1, 23)], ["chr{0}".format(x) for x in range(1, 23)])
    )
    + [("X", "chrX"), ("Y", "chrY"), ("MT", "chrM")]
)


def main():

    actions = (
        ("from23andme", "convert 23andme file to vcf file"),
        ("fromimpute2", "convert impute2 output to vcf file"),
        ("liftover", "lift over coordinates in vcf file"),
        ("location", "given SNP locations characterize the locations"),
        ("mstmap", "convert vcf format to mstmap input"),
        ("refallele", "make refAllele file"),
        ("sample", "sample subset of vcf file"),
        ("summary", "summarize the genotype calls in table"),
        ("uniq", "retain only the first entry in vcf file"),
        ("validate", "fast validation of vcf file"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def validate(args):
    """
    %prog validate input.vcf genome.fasta

    Fasta validation of vcf file.
    """
    import pyfasta

    p = OptionParser(validate.__doc__)
    p.add_option("--prefix", help="Add prefix to seqid")
    opts, args = p.parse_args(args)

    vcffile, fastafile = args
    pf = opts.prefix
    genome = pyfasta.Fasta(fastafile, record_class=pyfasta.MemoryRecord)
    fp = must_open(vcffile)
    match_ref = match_alt = total = 0
    for row in fp:
        if row[0] == "#":
            continue
        seqid, pos, id, ref, alt = row.split()[:5]
        total += 1
        if pf:
            seqid = pf + seqid
        pos = int(pos)
        if seqid not in genome:
            continue
        true_ref = genome[seqid][pos - 1]
        if total % 100000 == 0:
            print(total, "sites parsed", file=sys.stderr)
        if ref == true_ref:
            match_ref += 1
        elif alt == true_ref:
            match_alt += 1

    logging.debug("Match REF: {}".format(percentage(match_ref, total)))
    logging.debug("Match ALT: {}".format(percentage(match_alt, total)))


def uniq(args):
    """
    %prog uniq vcffile

    Retain only the first entry in vcf file.
    """
    from urllib.parse import parse_qs

    p = OptionParser(uniq.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (vcffile,) = args
    fp = must_open(vcffile)
    data = []
    for row in fp:
        if row[0] == "#":
            print(row.strip())
            continue
        v = VcfLine(row)
        data.append(v)

    for pos, vv in groupby(data, lambda x: x.pos):
        vv = list(vv)
        if len(vv) == 1:
            print(vv[0])
            continue
        bestv = max(vv, key=lambda x: float(parse_qs(x.info)["R2"][0]))
        print(bestv)


def sample(args):
    """
    %prog sample vcffile 0.9

    Sample subset of vcf file.
    """
    from random import random

    p = OptionParser(sample.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    vcffile, ratio = args
    ratio = float(ratio)
    fp = open(vcffile)
    pf = vcffile.rsplit(".", 1)[0]
    kept = pf + ".kept.vcf"
    withheld = pf + ".withheld.vcf"
    fwk = open(kept, "w")
    fww = open(withheld, "w")
    nkept = nwithheld = 0
    for row in fp:
        if row[0] == "#":
            print(row.strip(), file=fwk)
            continue
        if random() < ratio:
            nkept += 1
            print(row.strip(), file=fwk)
        else:
            nwithheld += 1
            print(row.strip(), file=fww)
    logging.debug("{0} records kept to `{1}`".format(nkept, kept))
    logging.debug("{0} records withheld to `{1}`".format(nwithheld, withheld))


def get_vcfstanza(fastafile, sampleid="SAMP_001"):
    from jcvi.formats.base import timestamp

    # VCF spec
    m = "##fileformat=VCFv4.1\n"
    m += "##fileDate={0}\n".format(timestamp())
    m += "##source={0}\n".format(__file__)
    m += "##reference=file://{0}\n".format(op.abspath(fastafile).strip("/"))
    m += '##INFO=<ID=PR,Number=0,Type=Flag,Description="Provisional genotype">\n'
    m += '##INFO=<ID=IM,Number=0,Type=Flag,Description="Imputed genotype">\n'
    m += '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    m += '##FORMAT=<ID=GP,Number=3,Type=Float,Description="Estimated Genotype Probability">\n'
    header = "CHROM POS ID REF ALT QUAL FILTER INFO FORMAT\n".split() + [sampleid]
    m += "#" + "\t".join(header)
    return m


def fromimpute2(args):
    """
    %prog fromimpute2 impute2file fastafile 1

    Convert impute2 output to vcf file. Imputed file looks like:

    --- 1:10177:A:AC 10177 A AC 0.451 0.547 0.002
    """
    p = OptionParser(fromimpute2.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    impute2file, fastafile, chr = args
    fasta = Fasta(fastafile)
    print(get_vcfstanza(fastafile))
    fp = open(impute2file)
    seen = set()
    for row in fp:
        snp_id, rsid, pos, ref, alt, aa, ab, bb = row.split()
        pos = int(pos)
        if pos in seen:
            continue
        seen.add(pos)
        code = max((float(aa), "0/0"), (float(ab), "0/1"), (float(bb), "1/1"))[-1]
        tag = "PR" if snp_id == chr else "IM"
        print(
            "\t".join(
                str(x)
                for x in (
                    chr,
                    pos,
                    rsid,
                    ref,
                    alt,
                    ".",
                    ".",
                    tag,
                    "GT:GP",
                    code + ":" + ",".join((aa, ab, bb)),
                )
            )
        )


def read_rsid(seqid, legend):
    if seqid in ["Y", "MT"]:
        return {}
    # Read rsid
    fp = open(legend)
    # rs145072688:10352:T:TA
    register = {}
    for row in fp:
        atoms = row.strip().split(":")
        if len(atoms) == 4:
            rsid, pos, ref, alt = atoms
        else:
            continue
        pos = int(pos)
        # Use position for non-rsid
        rsids = [pos] if rsid == seqid else [rsid, pos]
        for rsid in rsids:
            if rsid in register:
                pos1, ref1, alt1 = register[rsid]
                if alt not in alt1:
                    register[rsid][-1].append(alt)
            else:
                register[rsid] = (pos, ref, [alt])
    logging.debug(
        "A total of {0} sites imported from `{1}`".format(len(register), legend)
    )
    return register


def from23andme(args):
    """
    %prog from23andme txtfile 1

    Convert from23andme file to vcf file.

    --ref points to the folder that contains chr1.rsids

    $ zcat 1000GP_Phase3/1000GP_Phase3_chr1.legend.gz \\
            | cut -d" " -f1 | grep ":" > chr1.rsids
    """
    p = OptionParser(from23andme.__doc__)
    p.set_ref()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    txtfile, seqid = args
    ref_dir = opts.ref
    fastafile = op.join(ref_dir, "hs37d5.fa")
    fasta = Fasta(fastafile)

    pf = txtfile.rsplit(".", 1)[0]
    px = CM[seqid]
    chrvcf = pf + ".{0}.vcf".format(px)
    legend = op.join(ref_dir, "1000GP_Phase3/{0}.rsids".format(px))
    register = read_rsid(seqid, legend)

    fw = open(chrvcf, "w")
    print(get_vcfstanza(fastafile, txtfile), file=fw)

    fp = open(txtfile)
    seen = set()
    duplicates = skipped = missing = 0
    for row in fp:
        if row[0] == "#":
            continue
        rsid, chr, pos, genotype = row.split()
        if chr != seqid:
            continue
        pos = int(pos)
        if (chr, pos) in seen:
            duplicates += 1
            continue
        seen.add((chr, pos))
        genotype = list(genotype)
        if "-" in genotype:  # missing daa
            missing += 1
            continue

        # Y or MT
        if not register:
            assert len(genotype) == 1
            ref = fasta[chr][pos - 1].seq.upper()
            if "D" in genotype or "I" in genotype:
                skipped += 1
                continue
            genotype = genotype[0]
            code = "0/0" if ref == genotype else "1/1"
            alt = "." if ref == genotype else genotype
            print(
                "\t".join(
                    str(x)
                    for x in (chr, pos, rsid, ref, alt, ".", ".", "PR", "GT", code)
                ),
                file=fw,
            )
            continue

        # If rsid is seen in the db, use that
        if rsid in register:
            pos, ref, alt = register[rsid]
        elif pos in register:
            pos, ref, alt = register[pos]
        else:
            skipped += 1  # Not in reference panel
            continue

        assert fasta[chr][pos - 1 : pos + len(ref) - 1].seq.upper() == ref
        # Keep it bi-allelic
        not_seen = [x for x in alt if x not in genotype]
        while len(alt) > 1 and not_seen:
            alt.remove(not_seen.pop())
        if len(alt) > 1:
            alt = [alt[0]]
        alleles = [ref] + alt

        if len(genotype) == 1:
            genotype = [genotype[0]] * 2

        alt = ",".join(alt) or "."
        if "D" in genotype or "I" in genotype:
            max_allele = max((len(x), x) for x in alleles)[1]
            alleles = [("I" if x == max_allele else "D") for x in alleles]
            assert "I" in alleles and "D" in alleles
        a, b = genotype
        try:
            ia, ib = alleles.index(a), alleles.index(b)
        except ValueError:  # alleles not seen
            logging.error(
                "{0}: alleles={1}, genotype={2}".format(rsid, alleles, genotype)
            )
            skipped += 1
            continue
        code = "/".join(str(x) for x in sorted((ia, ib)))

        print(
            "\t".join(
                str(x) for x in (chr, pos, rsid, ref, alt, ".", ".", "PR", "GT", code)
            ),
            file=fw,
        )

    logging.debug(
        "duplicates={0} skipped={1} missing={2}".format(duplicates, skipped, missing)
    )


def refallele(args):
    """
    %prog refallele vcffile > out.refAllele

    Make refAllele file which can be used to convert PLINK file to VCF file.
    """
    p = OptionParser(refallele.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (vcffile,) = args
    fp = open(vcffile)
    for row in fp:
        if row[0] == "#":
            continue
        atoms = row.split()
        marker = "{0}:{1}".format(*atoms[:2])
        ref = atoms[3]
        print("\t".join((marker, ref)))


def location(args):
    """
    %prog location bedfile fastafile

    Given SNP locations, summarize the locations in the sequences. For example,
    find out if there are more 3`-SNPs than 5`-SNPs.
    """
    from jcvi.formats.bed import BedLine
    from jcvi.graphics.histogram import stem_leaf_plot

    p = OptionParser(location.__doc__)
    p.add_option(
        "--dist",
        default=100,
        type="int",
        help="Distance cutoff to call 5` and 3`",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, fastafile = args
    dist = opts.dist
    sizes = Sizes(fastafile).mapping
    fp = open(bedfile)
    fiveprime = threeprime = total = 0
    percentages = []
    for row in fp:
        b = BedLine(row)
        pos = b.start
        size = sizes[b.seqid]
        if pos < dist:
            fiveprime += 1
        if size - pos < dist:
            threeprime += 1
        total += 1
        percentages.append(100 * pos / size)

    m = "Five prime (within {0}bp of start codon): {1}\n".format(dist, fiveprime)
    m += "Three prime (within {0}bp of stop codon): {1}\n".format(dist, threeprime)
    m += "Total: {0}".format(total)
    print(m, file=sys.stderr)

    bins = 10
    title = "Locations within the gene [0=Five-prime, 100=Three-prime]"
    stem_leaf_plot(percentages, 0, 100, bins, title=title)


def summary(args):
    """
    %prog summary txtfile fastafile

    The txtfile can be generated by: %prog mstmap --noheader --freq=0

    Tabulate on all possible combinations of genotypes and provide results
    in a nicely-formatted table. Give a fastafile for SNP rate (average
    # of SNPs per Kb).

    Only three-column file is supported:
    locus_id    intra- genotype    inter- genotype
    """
    from jcvi.utils.cbook import thousands
    from jcvi.utils.table import tabulate

    p = OptionParser(summary.__doc__)
    p.add_option("--counts", help="Print SNP counts in a txt file")
    p.add_option("--bed", help="Print SNPs locations in a bed file")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    txtfile, fastafile = args
    bedfw = open(opts.bed, "w") if opts.bed else None

    fp = open(txtfile)
    header = next(fp).split()  # Header
    snps = defaultdict(list)  # contig => list of loci
    combinations = defaultdict(int)
    intraSNPs = interSNPs = 0
    distinctSet = set()  # set of genes that show A-B pattern
    ref, alt = header[1:3]
    snpcounts, goodsnpcounts = defaultdict(int), defaultdict(int)
    for row in fp:
        atoms = row.split()
        assert len(atoms) == 3, "Only three-column file is supported"
        locus, intra, inter = atoms
        ctg, pos = locus.rsplit(".", 1)
        pos = int(pos)
        snps[ctg].append(pos)
        snpcounts[ctg] += 1

        if intra == "X":
            intraSNPs += 1
        if inter in ("B", "X"):
            interSNPs += 1
        if intra == "A" and inter == "B":
            distinctSet.add(ctg)
            goodsnpcounts[ctg] += 1
        # Tabulate all possible combinations
        intra = ref + "-" + intra
        inter = alt + "-" + inter
        combinations[(intra, inter)] += 1

        if bedfw:
            print("\t".join(str(x) for x in (ctg, pos - 1, pos, locus)), file=bedfw)

    if bedfw:
        logging.debug("SNP locations written to `{0}`.".format(opts.bed))
        bedfw.close()

    nsites = sum(len(x) for x in snps.values())
    sizes = Sizes(fastafile)
    bpsize = sizes.totalsize
    snprate = lambda a: a * 1000.0 / bpsize
    m = "Dataset `{0}` contains {1} contigs ({2} bp).\n".format(
        fastafile, len(sizes), thousands(bpsize)
    )
    m += "A total of {0} SNPs within {1} contigs ({2} bp).\n".format(
        nsites, len(snps), thousands(sum(sizes.mapping[x] for x in snps.keys()))
    )
    m += "SNP rate: {0:.1f}/Kb, ".format(snprate(nsites))
    m += "IntraSNPs: {0} ({1:.1f}/Kb), InterSNPs: {2} ({3:.1f}/Kb)".format(
        intraSNPs, snprate(intraSNPs), interSNPs, snprate(interSNPs)
    )
    print(m, file=sys.stderr)
    print(tabulate(combinations), file=sys.stderr)

    leg = "Legend: A - homozygous same, B - homozygous different, X - heterozygous"
    print(leg, file=sys.stderr)

    tag = (ref + "-A", alt + "-B")
    distinctSNPs = combinations[tag]
    tag = str(tag).replace("'", "")
    print(
        "A total of {0} disparate {1} SNPs in {2} contigs.".format(
            distinctSNPs, tag, len(distinctSet)
        ),
        file=sys.stderr,
    )

    if not opts.counts:
        return

    snpcountsfile = opts.counts
    fw = open(snpcountsfile, "w")
    header = "\t".join(("Contig", "#_SNPs", "#_AB_SNP"))
    print(header, file=fw)

    assert sum(snpcounts.values()) == nsites
    assert sum(goodsnpcounts.values()) == distinctSNPs

    for ctg in sorted(snps.keys()):
        snpcount = snpcounts[ctg]
        goodsnpcount = goodsnpcounts[ctg]
        print("\t".join(str(x) for x in (ctg, snpcount, goodsnpcount)), file=fw)

    fw.close()
    logging.debug("SNP counts per contig is written to `{0}`.".format(snpcountsfile))


g2x = {"0/0": "A", "0/1": "X", "1/1": "B", "./.": "-", ".": "-"}


def encode_genotype(s, mindepth=3, depth_index=2, nohet=False):
    """
    >>> encode_genotype("1/1:128,18,0:6:18")  # homozygote B
    'B'
    >>> encode_genotype("0/1:0,0,0:0:3")      # missing data
    '-'
    >>> encode_genotype("0/1:128,0,26:7:22")  # heterozygous A/B
    'X'
    """
    atoms = s.split(":")
    if len(atoms) < 3:
        return g2x[atoms[0]]

    inferred = atoms[0]
    depth = int(atoms[depth_index])
    if depth < mindepth:
        return "-"
    if inferred == "0/0":
        return "A"
    if inferred == "0/1":
        return "-" if nohet else "X"
    if inferred == "1/1":
        return "B"
    return "-"


def mstmap(args):
    """
    %prog mstmap bcffile/vcffile > matrixfile

    Convert bcf/vcf format to mstmap input.
    """
    from jcvi.assembly.geneticmap import MSTMatrix

    p = OptionParser(mstmap.__doc__)
    p.add_option(
        "--dh",
        default=False,
        action="store_true",
        help="Double haploid population, no het",
    )
    p.add_option(
        "--freq",
        default=0.2,
        type="float",
        help="Allele must be above frequency",
    )
    p.add_option(
        "--mindepth",
        default=3,
        type="int",
        help="Only trust genotype calls with depth",
    )
    p.add_option(
        "--missing_threshold",
        default=0.25,
        type="float",
        help="Fraction missing must be below",
    )
    p.add_option(
        "--noheader",
        default=False,
        action="store_true",
        help="Do not print MSTmap run parameters",
    )
    p.add_option(
        "--pv4",
        default=False,
        action="store_true",
        help="Enable filtering strand-bias, tail distance bias, etc.",
    )
    p.add_option(
        "--freebayes",
        default=False,
        action="store_true",
        help="VCF output from freebayes",
    )
    p.set_sep(sep=".", help="Use separator to simplify individual names")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (vcffile,) = args
    if vcffile.endswith(".bcf"):
        bcffile = vcffile
        vcffile = bcffile.rsplit(".", 1)[0] + ".vcf"
        cmd = "bcftools view {0}".format(bcffile)
        cmd += " | vcfutils.pl varFilter"
        if not opts.pv4:
            cmd += " -1 0 -2 0 -3 0 -4 0 -e 0"
        if need_update(bcffile, vcffile):
            sh(cmd, outfile=vcffile)

    freq = opts.freq
    sep = opts.sep
    depth_index = 1 if opts.freebayes else 2

    ptype = "DH" if opts.dh else "RIL6"
    nohet = ptype == "DH"
    fp = open(vcffile)
    genotypes = []
    for row in fp:
        if row[:2] == "##":
            continue
        atoms = row.split()
        if row[0] == "#":
            ind = [x.split(sep)[0] for x in atoms[9:]]
            nind = len(ind)
            mh = ["locus_name"] + ind
            continue

        marker = "{0}.{1}".format(*atoms[:2])

        geno = atoms[9:]
        geno = [
            encode_genotype(
                x, mindepth=opts.mindepth, depth_index=depth_index, nohet=nohet
            )
            for x in geno
        ]
        assert len(geno) == nind
        f = 1.0 / nind

        if geno.count("A") * f < freq:
            continue
        if geno.count("B") * f < freq:
            continue
        if geno.count("-") * f > opts.missing_threshold:
            continue

        genotype = [marker] + geno
        genotypes.append(genotype)

    mm = MSTMatrix(genotypes, mh, ptype, opts.missing_threshold)
    mm.write(opts.outfile, header=(not opts.noheader))


def liftover(args):
    """
    %prog liftover old.vcf hg19ToHg38.over.chain.gz new.vcf

    Lift over coordinates in vcf file.
    """
    p = OptionParser(liftover.__doc__)
    p.add_option(
        "--newid", default=False, action="store_true", help="Make new identifiers"
    )
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    oldvcf, chainfile, newvcf = args
    ul = UniqueLiftover(chainfile)
    num_excluded = 0
    fp = open(oldvcf)
    fw = open(newvcf, "w")
    for row in fp:
        row = row.strip()
        if row[0] == "#":
            if row.startswith("##source="):
                row = "##source={0}".format(__file__)
            elif row.startswith("##reference="):
                row = "##reference=hg38"
            elif row.startswith("##contig="):
                continue
            print(row.strip(), file=fw)
            continue

        v = VcfLine(row)
        # GRCh37.p2 has the same MT sequence as hg38 (but hg19 is different)
        if v.seqid == "MT":
            v.seqid = "chrM"
            print(v, file=fw)
            continue

        try:
            new_chrom, new_pos = ul.liftover_cpra(CM[v.seqid], v.pos)
        except:
            num_excluded += 1
            continue

        if new_chrom is not None and new_pos is not None:
            v.seqid, v.pos = new_chrom, new_pos
            if opts.newid:
                v.rsid = "{0}:{1}".format(new_chrom.replace("chr", ""), new_pos)
            print(v, file=fw)
        else:
            num_excluded += 1

    logging.debug("Excluded {0}".format(num_excluded))


if __name__ == "__main__":
    main()
