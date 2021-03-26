#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper for running series of EVM commands. There are two flavors of running
EVM - TIGR only mode which communicates with the Sybase db; evm mode which
communicates with GFF file.
"""
from __future__ import print_function

import os.path as op
import sys

from collections import defaultdict

from jcvi.formats.fasta import ids
from jcvi.formats.base import write_file
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


EVMRUN = r"""
W=`pwd`/weights.txt

$EVM/EvmUtils/write_EVM_commands.pl --genome genome.fasta --weights $W \
    --gene_predictions {0} \
    --transcript_alignments {1} \
    --protein_alignments {2} \
    --terminalExons pasa.terminal_exons.gff3 \
    --output_file_name evm.out --partitions partitions_list.out > commands.list

$EGC_SCRIPTS/run_cmds_on_grid.pl commands.list 0372

#$EVM/EvmUtils/execute_EVM_commands.pl commands.list
"""

EVMLOAD = r"""
$EVM/EvmUtils/recombine_EVM_partial_outputs.pl  \
    --partitions partitions_list.out \
    --output_file_name evm.out

$EVM/TIGR-only/TIGR_EVM_loader.pl --db {0} \
    --partitions partitions_list.out \
    --output_file_name evm.out \
    --ev_type {1}

#$EVM/EvmUtils/convert_EVM_outputs_to_GFF3.pl \
#    --partitions partitions_list.out \
#    --output evm.out
"""


def main():

    actions = (
        ("pasa", "extract terminal exons"),
        ("tigrprepare", "run EVM in TIGR-only mode"),
        ("tigrload", "load EVM results into TIGR db"),
        ("maker", "run EVM based on MAKER output"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def partition(evs):
    partition_list = "partitions_list.out"
    A, T, P = evs
    if not need_update(evs, partition_list):
        return

    cmd = "$EVM/EvmUtils/partition_EVM_inputs.pl --genome genome.fasta"
    cmd += " --gene_predictions {0}".format(A)
    cmd += " --transcript_alignments {0}".format(T)
    cmd += " --protein_alignments {0}".format(P)
    cmd += " --segmentSize 500000 --overlapSize 10000 "
    cmd += " --partition_listing partitions_list.out"

    termexons = "pasa.terminal_exons.gff3"
    if op.exists(termexons):
        cmd += " --pasaTerminalExons {0}".format(termexons)

    sh(cmd)


def maker(args):
    """
    %prog maker maker.gff3 genome.fasta

    Prepare EVM inputs by separating tracks from MAKER.
    """
    from jcvi.formats.base import SetFile, FileShredder

    A, T, P = "ABINITIO_PREDICTION", "TRANSCRIPT", "PROTEIN"
    # Stores default weights and types
    Registry = {
        "maker": (A, 5),
        "augustus_masked": (A, 1),
        "snap_masked": (A, 1),
        "genemark": (A, 1),
        "est2genome": (T, 5),
        "est_gff": (T, 5),
        "protein2genome": (P, 5),
        "blastx": (P, 1),
    }

    p = OptionParser(maker.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gffile, fastafile = args

    types = "type.ids"
    if need_update(gffile, types):
        cmd = "cut -f2 -s {0} | sort -u".format(gffile)
        sh(cmd, outfile=types)

    types = SetFile(types)
    reg = defaultdict(list)
    weightsfile = "weights.txt"
    contents = []
    for s in types:
        rs = s.split(":")[0]
        if rs not in Registry:
            continue

        type, weight = Registry[rs]
        reg[type].append(s)
        contents.append("\t".join(str(x) for x in (type, s, weight)))

    contents = "\n".join(sorted(contents))
    write_file(weightsfile, contents)

    evs = [x + ".gff" for x in (A, T, P)]
    FileShredder(evs)

    for type, tracks in reg.items():
        for t in tracks:
            cmd = "grep '\t{0}' {1} | grep -v '_match\t' >> {2}.gff".format(
                t, gffile, type
            )
            sh(cmd)

    partition(evs)
    runfile = "run.sh"
    contents = EVMRUN.format(*evs)
    write_file(runfile, contents)


def tigrload(args):
    """
    %prog tigrload db ev_type

    Load EVM results into TIGR db. Actually, just write a load.sh script. The
    ev_type should be set, e.g. "EVM1", "EVM2", etc.
    """
    p = OptionParser(tigrload.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    db, ev_type = args

    runfile = "load.sh"
    contents = EVMLOAD.format(db, ev_type)
    write_file(runfile, contents)


def pasa(args):
    """
    %prog pasa pasa_db fastafile

    Run EVM in TIGR-only mode.
    """
    p = OptionParser(pasa.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pasa_db, fastafile = args

    termexons = "pasa.terminal_exons.gff3"
    if need_update(fastafile, termexons):
        cmd = "$ANNOT_DEVEL/PASA2/scripts/pasa_asmbls_to_training_set.dbi"
        cmd += ' -M "{0}:mysql.tigr.org" -p "access:access"'.format(pasa_db)
        cmd += " -g {0}".format(fastafile)
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
                    print(stack, file=fw)
                stack = row
            else:
                print(stack + row, file=fw)
                stack = ""

        fw.close()

    return fixedtranscript


def tigrprepare(args):
    """
    %prog tigrprepare asmbl.fasta asmbl.ids db pasa.terminal_exons.gff3

    Run EVM in TIGR-only mode.
    """
    p = OptionParser(tigrprepare.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    fastafile, asmbl_id, db, pasa_db = args
    if asmbl_id == "all":
        idsfile = fastafile + ".ids"
        if need_update(fastafile, idsfile):
            ids([fastafile, "-o", idsfile])
    else:
        idsfile = asmbl_id

    oneid = next(open(idsfile)).strip()

    weightsfile = "weights.txt"
    if need_update(idsfile, weightsfile):
        cmd = "$EVM/TIGR-only/create_sample_weights_file.dbi"
        cmd += " {0} {1} | tee weights.txt".format(db, oneid)
        sh(cmd)

    evs = [
        "gene_predictions.gff3",
        "transcript_alignments.gff3",
        "protein_alignments.gff3",
    ]
    if need_update(weightsfile, evs):
        cmd = "$EVM/TIGR-only/write_GFF3_files.dbi"
        cmd += " --db {0} --asmbl_id {1} --weights {2}".format(db, idsfile, weightsfile)
        sh(cmd)

    evs[1] = fix_transcript()

    partition(evs)
    runfile = "run.sh"
    contents = EVMRUN.format(*evs)
    write_file(runfile, contents)


if __name__ == "__main__":
    main()
