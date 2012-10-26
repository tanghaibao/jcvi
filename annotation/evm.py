#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for running series of EVM commands. There are two flavors of running
EVM - TIGR only mode which communicates with the Sybase db; evm mode which
communicates with GFF file.
"""

import sys
import logging

from optparse import OptionParser

from jcvi.formats.fasta import ids
from jcvi.formats.base import check_exists
from jcvi.apps.base import ActionDispatcher, need_update, debug, sh
debug()


EVMRUN = r"""
#!/bin/bash

W=`pwd`/weights.txt

$EVM/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights $W \
    --gene_predictions gene_predictions.gff3 \
    --protein_alignments protein_alignments.gff3 \
    --transcript_alignments transcript_alignments.fixed.gff3 \
    --terminalExons pasa.terminal_exons.gff3 \
    --output_file_name evm.out --partitions partitions_list.out > commands.list

$EGC_SCRIPTS/run_cmds_on_grid.pl commands.list 04048"""


def main():

    actions = (
        ('pasa', 'extract terminal exons'),
        ('tigrprepare', 'run EVM in TIGR-only mode'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pasa(args):
    """
    %prog pasa pasa_db

    Run EVM in TIGR-only mode.
    """
    p = OptionParser(pasa.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    pasa_db, = args

    termexons = "pasa.terminal_exons.gff3"
    if need_update(weightsfile, termexons):
        cmd = "$ANNOT_DEVEL/PASA2/scripts/pasa_asmbls_to_training_set.dbi"
        cmd += ' -M "{0}:mysql.tigr.org" -p "access:access"'.format(pasa_db)
        cmd += ' -g {0}'.format(fastafile)
        sh(cmd)

        cmd = "$EVM/PasaUtils/retrieve_terminal_CDS_exons.pl"
        cmd += " trainingSetCandidates.fasta trainingSetCandidates.gff"
        sh(cmd, outfile=termexons)

    return termexons


def fix_transcript():
    # Fix `transcript_alignments.gff3`
    transcript = "transcript_alignments.gff3"
    fixedtranscript = "transcript_alignments.fixed.gff3"
    if need_update(transcript, fixedtranscript):
        fp = open(transcript)
        fw = open(fixedtranscript, "w")
        stack = ""
        for row in fp:
            row = row.rstrip()
            goodline = len(row.split()) == 9
            if goodline:
                if stack:
                    print >> fw, stack
                stack = row
            else:
                print >> fw, stack + row
                stack = ""

        fw.close()

    return fixedtranscript


def tigrprepare(args):
    """
    %prog tigrprepare asmbl.fasta [all|asmbl_id] db pasa.terminal_exons.gff3

    Run EVM in TIGR-only mode.
    """
    p = OptionParser(tigrprepare.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    fastafile, asmbl_id, db, pasa_db = args
    if asmbl_id == 'all':
        idsfile = fastafile + ".ids"
        if need_update(fastafile, idsfile):
            ids([fastafile, "-o", idsfile])
    else:
        idsfile = asmbl_id + ".ids"
        if need_update(fastafile, idsfile):
            fw = open(idsfile, "w")
            print >> fw, asmbl_id
            fw.close()

    oneid = open(idsfile).next().strip()

    weightsfile = "weights.txt"
    if need_update(idsfile, weightsfile):
        cmd = "$EVM/TIGR-only/create_sample_weights_file.dbi"
        cmd += " mta4 {0} | tee weights.txt".format(oneid)
        sh(cmd)

    evs = ["gene_predictions.gff3", "transcript_alignments.gff3",
           "protein_alignments.gff3"]
    if need_update(weightsfile, evs):
        cmd = "$EVM/TIGR-only/write_GFF3_files.dbi"
        cmd += " --db {0} --asmbl_id {1} --weights {2}".\
                format(db, idsfile, weightsfile)
        sh(cmd)

    transcript = fix_transcript()
    partition_list = "partitions_list.out"
    evs[1] = transcript
    evs += ["pasa.terminal_exons.gff3"]
    if need_update(evs, partition_list):
        cmd = "$EVM/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta"
        cmd += " --gene_predictions gene_predictions.gff3"
        cmd += " --protein_alignments protein_alignments.gff3"
        cmd += " --transcript_alignments transcript_alignments.fixed.gff3"
        cmd += " --pasaTerminalExons pasa.terminal_exons.gff3"
        cmd += " --segmentSize 500000 --overlapSize 10000 "
        cmd += " --partition_listing partitions_list.out"
        sh(cmd)

    runfile = "run.sh"
    if check_exists(runfile):
        fw = open(runfile, "w")
        print >> fw, EVMRUN
        logging.debug("Run script written to `{0}`.".format(runfile))
        fw.close()


if __name__ == '__main__':
    main()
