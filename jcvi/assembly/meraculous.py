#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Subroutines to aid MERACULOUS assembly.
"""

import os.path as op
import sys

from jcvi.formats.fastq import readlen
from jcvi.formats.base import write_file
from jcvi.assembly.base import FastqNamings, get_libs
from jcvi.utils.table import comment_banner, load_csv
from jcvi.apps.base import OptionParser, ActionDispatcher


header = """
lib_seq [ wildcard ][ prefix ][ insAvg ][ insSdev ][ avgReadLen ][
hasInnieArtifact ][ isRevComped ][ useForContigging ][ scaffRound ][
useForGapClosing ][ 5pWiggleRoom ][3pWiggleRoom]
""".translate(None, "[]\n")


def main():

    actions = (
        ('prepare', 'prepare meraculous input'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare genomesize *.fastq

    Prepare MERACULOUS configuation file. Genome size should be entered in Mb.
    """
    p = OptionParser(prepare.__doc__ + FastqNamings)
    p.add_option("-K", default=51, type="int", help="K-mer size")
    p.set_cpus(cpus=32)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    genomesize = float(args[0]) / 1000
    fnames = args[1:]
    for x in fnames:
        assert op.exists(x), "File `{0}` not found.".format(x)

    s = comment_banner("Meraculous params file") + "\n"
    s += comment_banner("Basic parameters") + "\n"
    s += "# Describe the libraries ( one line per library )\n"
    s += "# " + " ".join(header.split()) + "\n"

    libs = get_libs(fnames)
    lib_seqs = []
    rank = 0
    for lib, fs in libs:
        size = lib.size
        if size == 0:
            continue
        rank += 1
        library_name = lib.library_name
        name = library_name.replace("-", "")
        wildcard = "{0}*.1.*,{0}*.2.*".format(library_name)
        rl = max(readlen([x]) for x in fs)
        lib_seq = lib.get_lib_seq(wildcard, name, rl, rank)
        lib_seqs.append(lib_seq)

    s += "\n" + "\n".join(load_csv(None, lib_seqs, sep=" ")) + "\n"
    params = [("genome_size", genomesize),
              ("is_diploid", 0),
              ("mer_size", opts.K),
              ("num_prefix_blocks", 1),
              ("no_read_validation", 0),
              ("local_num_procs", opts.cpus)]
    s += "\n" + "\n".join(load_csv(None, params, sep=" ")) + "\n"

    cfgfile = "meraculous.config"
    write_file(cfgfile, s, tee=True)

    s = "~/export/meraculous/bin/run_meraculous.sh -c {0}"\
                .format(cfgfile)
    runsh = "run.sh"
    write_file(runsh, s)


if __name__ == '__main__':
    main()
