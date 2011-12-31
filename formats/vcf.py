#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Variant call format.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('mstmap', 'convert vcf format to mstmap input'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def encode_genotype(s):
    """
    >>> encode_genotype("1/1:128,18,0:6:18")  # homozygote B
    'B'
    >>> encode_genotype("0/1:0,0,0:0:3")      # missing data
    '-'
    >>> encode_genotype("0/1:128,0,26:7:22")  # heterozygous A/B
    'X'
    """
    atoms = s.split(":")
    inferred, likelihood, depth = atoms[:3]
    depth = int(depth)
    if depth == 0:
        return '-'
    if inferred == '0/0':
        return 'A'
    if inferred == '0/1':
        return 'X'
    if inferred == '1/1':
        return 'B'


def mstmap(args):
    """
    %prog mstmap vcffile

    Convert vcf format to mstmap input.
    """
    from urlparse import parse_qs

    p = OptionParser(mstmap.__doc__)
    p.add_option("--freq", default=.25, type="float",
                 help="alle must be above frequencey [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    vcffile, = args

    header = """population_type RIL6
population_name LG
distance_function kosambi
cut_off_p_value 2.0
no_map_dist 15.0
no_map_size 0
missing_threshold 0.25
estimation_before_clustering no
detect_bad_data yes
objective_function COUNT
number_of_loci {0}
number_of_individual {1}
    """

    fp = open(vcffile)
    genotypes = []
    for row in fp:
        if row[:2] == "##":
            continue
        atoms = row.split()
        if row[0] == '#':
            ind = [x.split(".")[0] for x in atoms[9:]]
            nind = len(ind)
            mh = "\t".join(["locus_name"] + ind)
            continue

        marker = "{0}.{1}".format(*atoms[:2])
        score = atoms[7]
        score = parse_qs(score)
        dp4= score["DP4"][0]

        a, b, c, d = dp4.split(",")
        ref = int(a) + int(b)
        alt = int(c) + int(d)
        total = ref + alt
        if ref * 1. / total < opts.freq:
            continue
        if alt * 1. / total < opts.freq:
            continue

        geno = atoms[9:]
        geno = [encode_genotype(x) for x in geno]
        assert len(geno) == nind

        genotype = "\t".join([marker] + geno)
        genotypes.append(genotype)

    ngenotypes = len(genotypes)
    logging.debug("Imported {0} markers and {1} individuals.".\
                  format(ngenotypes, nind))

    print header.format(ngenotypes, nind)
    print mh
    print "\n".join(genotypes)


if __name__ == '__main__':
    main()
