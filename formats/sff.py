#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Procedures to handle Roche 454 sff data, mostly use sff-tools supplied by Roche.
"""

import os.path as op
import sys
import logging

from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, sh, glob


def main():

    actions = (
        ('mid', 'produce a MID configuration file from 2-column mapping'),
        ('assemble', 'assemble each BAC separately using newbler'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def assemble(args):
    """
    %prog assemble sffdir

    Assemble each BAC separately using newbler.
    """
    from jcvi.formats.fasta import join

    p = OptionParser(assemble.__doc__)
    p.add_option("--overwrite", default=False, action="store_true",
            help="Overwrite the separate BAC assembly [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    sffdir, = args
    asmdir = "newbler"
    fastadir = "fasta"
    mkdir(asmdir, overwrite=opts.overwrite)
    mkdir(fastadir, overwrite=opts.overwrite)
    cmd = "runAssembly -cpu 8 -o {0} {1}"
    for sffile in glob("{0}/*.sff".format(sffdir)):
        pf = op.basename(sffile).split(".")[1]
        pf = pf.lower()
        outdir = op.join(asmdir, pf)
        if op.exists(outdir):
            logging.debug("`{0}` exists. Ignored.".format(outdir))
            continue

        acmd = cmd.format(outdir, sffile)
        sh(acmd)

        ctgfile = op.join(outdir, "454LargeContigs.fna")
        if not op.exists(ctgfile):  # newbler failure
            logging.error("File `{0}` not found (newbler failure).".\
                    format(ctgfile))
            continue
        outfile = op.join(fastadir, "{0}.fasta".format(pf))
        newidopt = "--newid={0}".format(pf)
        minctgsizeopt = "--minctgsize=200"
        join([ctgfile, outfile, newidopt, minctgsizeopt])


def mid(args):
    """
    %prog mid mappingfile

    Produce a MID configuration file from primer mapping. The primer mapping can
    be two or three columns, the third column is the optional 3`-primer.
    """
    from jcvi.formats.fasta import rc

    p = OptionParser(mid.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    mappingfile, = args

    fp = open(mappingfile)
    data = [row.split() for row in fp]
    templatefile = "/usr/local/seq454-64_v2.6/config/MIDConfig.parse"
    midfile = op.basename(templatefile)
    fw = open(midfile, "w")

    lines = open(templatefile).readlines()

    insertline = 45
    first, second = lines[:insertline], lines[insertline:]
    for row in first:
        fw.write(row)

    # The inserted block
    print >> fw, "MYMIDs\n{"
    for atoms in data:
        natoms = len(atoms)
        assert natoms in (2, 3)
        name, p5seq = atoms[:2]
        line = '        mid = "{0}", "{1}", 1'.format(name, p5seq)

        # Since most often I find p3 primers by parsing newbler progress file
        # The string I grepped needs to be reverse-complemented
        p3seq = atoms[2]
        p3seq = rc(p3seq)

        if natoms == 3:
            line += ', "{0}"'.format(p3seq)
        line += ';'
        print >> fw, line

    print >> fw, "}"

    for row in second:
        fw.write(row)

    logging.debug("Barcodes written to `{0}`.".format(midfile))


if __name__ == '__main__':
    main()
