#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Read-based phasing.
"""
import sys
import logging
import vcf
import pysam

from jcvi.apps.base import OptionParser, ActionDispatcher


class CPRA:
    def __init__(self, vcf_record):
        r = vcf_record
        self.chr = r.CHROM
        self.pos = r.POS
        self.ref = r.REF
        self.alt = r.ALT
        self.alleles = [self.ref] + self.alt

    @property
    def is_valid(self):
        """Only retain SNPs or single indels, and are bi-allelic"""
        return len(self.ref) == 1 and len(self.alt) == 1 and len(self.alt[0]) == 1

    def __str__(self):
        return "_".join(str(x) for x in (self.chr, self.pos, self.ref, self.alt[0]))

    __repr__ = __str__


def main():

    actions = (
        ("prepare", "convert vcf and bam to variant list"),
        ("counts", "collect allele counts from RO/AO fields"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def counts(args):
    """
    %prog counts vcffile

    Collect allele counts from RO and AO fields.
    """
    p = OptionParser(counts.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (vcffile,) = args
    vcf_reader = vcf.Reader(open(vcffile))
    for r in vcf_reader:
        v = CPRA(r)
        if not v.is_valid:
            continue
        for sample in r.samples:
            ro = sample["RO"]
            ao = sample["AO"]
            print("\t".join(str(x) for x in (v, ro, ao)))


def prepare(args):
    """
    %prog prepare vcffile bamfile

    Convert vcf and bam to variant list. Inputs are:
    - vcffile: contains the positions of variants
    - bamfile: contains the reads that hold the variants

    Outputs:
    - reads_to_phase: phasing for each read
    - variants_to_phase: in format of phased vcf
    """
    p = OptionParser(prepare.__doc__)
    p.add_option("--accuracy", default=0.85, help="Sequencing per-base accuracy")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    vcffile, bamfile = args
    right = "{:.2f}".format(opts.accuracy)
    wrong = "{:.2f}".format(1 - opts.accuracy)
    vcf_reader = vcf.Reader(open(vcffile))
    variants = []
    for r in vcf_reader:
        v = CPRA(r)
        if not v.is_valid:
            continue
        variants.append(v)

    logging.debug(
        "A total of {} bi-allelic SNVs imported from `{}`".format(
            len(variants), vcffile
        )
    )

    bamfile = pysam.AlignmentFile(bamfile, "rb")
    for v in variants:
        pos = v.pos - 1
        for column in bamfile.pileup(v.chr, pos, pos + 1, truncate=True):
            for read in column.pileups:
                query_position = read.query_position
                if query_position is None:
                    continue
                read_name = read.alignment.query_name
                query_base = read.alignment.query_sequence[query_position]
                a, b = v.alleles
                if query_base == a:
                    other_base = b
                elif query_base == b:
                    other_base = a
                else:
                    continue
                print(
                    " ".join(
                        str(x)
                        for x in (v, read_name, query_base, right, other_base, wrong)
                    )
                )


if __name__ == "__main__":
    main()
