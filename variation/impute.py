#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Impute unknown variations given an input vcf file.
"""

import os.path as op
import sys

from jcvi.apps.grid import MakeManager
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('impute', 'use IMPUTE2 to impute vcf'),
        ('beagle', 'use BEAGLE4.1 to impute vcf'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def impute(args):
    """
    %prog impute input.vcf hs37d5.fa 1

    Use IMPUTE2 to impute vcf on chromosome 1.
    """
    from pyfaidx import Fasta

    p = OptionParser(impute.__doc__)
    p.set_home("shapeit")
    p.set_home("impute")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    vcffile, fastafile, chr = args
    mm = MakeManager()
    pf = vcffile.rsplit(".", 1)[0]
    rpf = "1000GP_Phase3/1000GP_Phase3_chr{0}".format(chr)
    shapeit_cmd = op.join(opts.shapeit_home, "shapeit")
    mapfile = "1000GP_Phase3/genetic_map_chr{0}_combined_b37.txt".format(chr)
    cmd = shapeit_cmd + " --input-vcf {0}".format(vcffile)
    cmd += " --input-map {0}".format(mapfile)
    cmd += " --threads {0} --effective-size 20000".format(opts.cpus)
    cmd += " --output-max {0}.haps {0}.sample".format(pf)
    cmd += " --input-ref {0}.hap.gz {0}.legend.gz {0}.sample".format(rpf)
    cmd += " --output-log {0}.log".format(pf)
    hapsfile = pf + ".haps"
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
        cmd = impute_cmd + " -m {0}".format(mapfile)
        cmd += " -known_haps_g {0}".format(hapsfile)
        cmd += " -h {0}.hap.gz -l {0}.legend.gz".format(rpf)
        cmd += " -Ne 20000 -int {0} {1}".format(chunk_start, chunk_end)
        cmd += " -o {0} -allow_large_regions -seed 367946".format(outfile)
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
