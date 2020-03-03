#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run cap3 command
The whole pipeline is following CAP3 documentation at
<http://deepc2.psi.iastate.edu/aat/cap/capdoc.html>

The action `prepare` is specific to the JCVI Legume FL-cDNA project.
"""
from __future__ import print_function

import re
import sys
import logging
import os
import os.path as op

from jcvi.apps.base import OptionParser, OptionGroup, ActionDispatcher, sh


def main():

    actions = (
        ("prepare", "prepares data for per-clone assembly"),
        ("assemble", "wraps the cap3 assembly command"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare --rearray_lib=<rearraylibrary> --orig_lib_file=<origlibfile>

    Inferred file names
    ---------------------------------------------
    `lookuptblfile` : rearraylibrary.lookup
    `rearraylibfile`: rearraylibrary.fasta

    Pick sequences from the original library file and the rearrayed library file
    based on the mapping information provided in the `lookuptblfile`.

    # lookuptblfile format: column number (index)
    # 1 (0)          2 (1)          3 (2)         4 (3)        5 (4)        6 (5)
    # source_clone   source_plate   source_well   dest_clone   dest_plate   dest_well

    The 1st and 4th column in the `lookuptblfile` form the pair of clones which
    constitute the elements used for the per-clone assembly.
    """
    from operator import itemgetter
    from jcvi.formats.fasta import Fasta, SeqIO

    p = OptionParser(prepare.__doc__)
    p.add_option(
        "--rearray_lib", default=None, help="name of the rearrayed library",
    )
    p.add_option(
        "--orig_lib_file",
        help="fasta file containing reads from the original libraries",
    )

    g = OptionGroup(p, "Optional parameters")
    g.add_option(
        "--output_folder",
        default="to_assemble",
        help="output folder to write the FASTA files to",
    )
    p.add_option_group(g)

    opts, args = p.parse_args(args)

    if not opts.rearray_lib or not opts.orig_lib_file:
        logging.error("Please specify the required parameters")
        sys.exit(not p.print_help())

    rearraylib, origlibfile = opts.rearray_lib, opts.orig_lib_file

    if not op.isfile(origlibfile):
        logging.error(
            "Original library reads file `{0}` does not exist!".format(origlibfile)
        )
        sys.exit()

    lookuptblfile = rearraylib + ".lookup"
    logging.debug(lookuptblfile)
    if not op.isfile(lookuptblfile):
        logging.error("Lookup table file `{0}` does not exist!".format(lookuptblfile))
        sys.exit()

    rearraylibfile = rearraylib + ".fasta"
    logging.debug(rearraylibfile)
    if not op.isfile(rearraylibfile):
        logging.error(
            "Rearrayed library reads file `{0}` does not exist!".format(rearraylibfile)
        )
        sys.exit()

    origlibFasta = Fasta(origlibfile)
    rearraylibFasta = Fasta(rearraylibfile)

    origlibids = [o for o in origlibFasta.iterkeys_ordered()]
    rearraylibids = [r for r in rearraylibFasta.iterkeys_ordered()]

    if not op.isdir(opts.output_folder):
        logging.warning(
            "Output directory `{0}` missing. Creating it now...".format(
                opts.output_folder
            )
        )
        os.makedirs(opts.output_folder)

    logfile = rearraylib + ".log"
    log = open(logfile, "w")

    fp = open(lookuptblfile, "r")
    for row in fp:
        origprefix, rearrayprefix = itemgetter(0, 3)(row.split("\t"))
        libpair = origprefix + "_" + rearrayprefix
        outfile = opts.output_folder + "/" + libpair + ".fasta"
        ofp = open(outfile, "w")

        for o in origlibids:
            if re.match(origprefix, o):
                SeqIO.write(origlibFasta[o], ofp, "fasta")

        for r in rearraylibids:
            if re.match(rearrayprefix, r):
                SeqIO.write(rearraylibFasta[r], ofp, "fasta")

        ofp.close()
        print(outfile, file=log)

    log.close()
    logging.debug("Wrote log file `{0}`".format(logfile))


def assemble(args):
    """
    Run `cap3` on a single multi FASTA file containing reads or a folder containing several
    multi FASTA files. Allows for tweaking of `cap3` parameters max_gap_len, ovl_pct_id, etc.
    """
    p = OptionParser(assemble.__doc__)
    g1 = OptionGroup(
        p,
        "Input file options (required)",
        "Note: Please choose from and provide values for one of the following parameters",
    )
    g1.add_option("--input_file", default=None, help="input file of reads")
    g1.add_option(
        "--input_folder",
        default=None,
        help="input folder containing multi FASTA files of reads",
    )
    g1.add_option(
        "--input_file_list",
        default=None,
        help="list file containing paths to multi FASTA files of reads",
    )
    p.add_option_group(g1)

    g2 = OptionGroup(
        p, "Optional parameters", "Note: If not specified, `cap3` defaults will be used"
    )
    g2.add_option(
        "-f",
        "--max_gap_len",
        default=20,
        type="int",
        help="maximum gap length in any overlap\n" + "Same as cap3 `-f` parameter.",
    )
    g2.add_option(
        "-p",
        "--ovl_pct_id",
        default=90,
        type="int",
        help="overlap percent identity cutoff\n" + "Same as cap3 `-p` parameter.",
    )
    g2.add_option(
        "-s",
        "--ovl_sim_score",
        default=900,
        type="int",
        help="overlap similarity score cutoff\n" + "Same as cap3 `-s` parameter.",
    )
    g2.add_option(
        "-x",
        "--prefix",
        dest="prefix",
        default="cap3",
        help="prefix string for output file name",
    )
    p.add_option_group(g2)

    p.set_params()

    opts, args = p.parse_args(args)

    if opts.max_gap_len and opts.max_gap_len <= 1:
        logging.error("--max_gap_len should be > 1")
        sys.exit()
    elif opts.ovl_pct_id and opts.ovl_pct_id <= 65:
        logging.error("--ovl_pct_id should be > 65")
        sys.exit()
    elif opts.ovl_sim_score and opts.ovl_sim_score <= 250:
        logging.error("--ovl_sim_score should be > 250")
        sys.exit()

    file_list = []
    if opts.input_file_list:
        if not op.isfile(opts.input_file_list):
            logging.error(
                "Input file list {0} does not exist".format(opts.input_file_list)
            )
            sys.exit()
        with open(opts.input_file_list, "r") as f:
            file_list = f.read().splitlines()
    elif opts.input_folder:
        if not op.isdir(opts.input_folder):
            logging.error("Input folder {0} does not exist".format(opts.input_folder))
            sys.exit()

        file_list = [
            file
            for file in os.listdir(opts.input_folder)
            if file.lower().endswith((".fa", ".fasta"))
        ]
        folder = opts.input_folder
        folder = folder.rstrip("/")
        for i in range(len(file_list)):
            file_list[i] = folder + "/" + file_list[i]
    elif opts.input_file:
        file_list.append(opts.input_file)
    else:
        logging.error("Please specify one of the options for input files")
        sys.exit(not p.print_help())

    if len(file_list) == 0:
        logging.warning("List of files to process is empty. Please check your input!")
        sys.exit()

    for file in file_list:
        if not op.isfile(file):
            logging.warning("Input file {0} does not exist".format(file))
        else:
            cmd = "cap3 {0} -f {1} -p {2} -s {3} -x {4}".format(
                file, opts.max_gap_len, opts.ovl_pct_id, opts.ovl_sim_score, opts.prefix
            )
            if opts.extra:
                cmd += " {0}".format(opts.extra)
            logfile = "{0}.{1}.log".format(file, opts.prefix)

            sh(cmd, outfile=logfile)


if __name__ == "__main__":
    main()
