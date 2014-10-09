#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Reference-free tGBS related functions.
"""

import os.path as op
import sys

from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


def main():

    actions = (
        ('snp', 'run SNP calling on GSNAP output'),
        ('bam', 'convert GSNAP output to BAM'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bam(args):
    """
    %prog snp input.gsnap ref.fasta

    Convert GSNAP output to BAM.
    """
    from jcvi.formats.sizes import Sizes
    from jcvi.formats.sam import index

    p = OptionParser(bam.__doc__)
    p.set_home("eddyyeh")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gsnapfile, fastafile = args
    EYHOME = opts.eddyyeh_home
    pf = gsnapfile.rsplit(".", 1)[0]
    uniqsam = pf + ".unique.sam"
    if need_update((gsnapfile, fastafile), uniqsam):
        cmd = op.join(EYHOME, "gsnap2gff3.pl")
        sizesfile = Sizes(fastafile).filename
        cmd += " --format sam -i {0} -o {1}".format(gsnapfile, uniqsam)
        cmd += " -u -l {0} -p {1}".format(sizesfile, opts.cpus)
        sh(cmd)

    index([uniqsam])


def snp(args):
    """
    %prog snp input.gsnap

    Run SNP calling on GSNAP output after apps.gsnap.align().
    """
    p = OptionParser(snp.__doc__)
    p.set_home("eddyyeh")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    gsnapfile, = args
    EYHOME = opts.eddyyeh_home
    pf = gsnapfile.rsplit(".", 1)[0]
    nativefile = pf + ".native"
    if need_update(gsnapfile, nativefile):
        cmd = op.join(EYHOME, "convert2native.pl")
        cmd += " --gsnap {0} -o {1}".format(gsnapfile, nativefile)
        cmd += " -proc {0}".format(opts.cpus)
        sh(cmd)

    snpfile = pf + ".snp"
    if need_update(nativefile, snpfile):
        cmd = op.join(EYHOME, "SNPs/SNP_Discovery-short.pl")
        cmd += " --native {0} -o {1}".format(nativefile, snpfile)
        cmd += " -a 2 -ac 0.3 -c 0.8"
        sh(cmd)


if __name__ == '__main__':
    main()
