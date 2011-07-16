#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to call tigStore and utgcns, for debugging failed utgcns runs.

See full commands:
http://sf.net/apps/mediawiki/wgs-assembler/index.php?title=Unitig_Consensus_Failures_in_CA_6

It is expected to be executed within 5-consensus/ folder.
"""

import os
import os.path as op
import sys
import shutil
import logging

from glob import glob
from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, sh, mkdir, debug
from jcvi.assembly.base import CAPATH
debug()


def main():

    actions = (
        ('error', 'find all errors in ../5-consensus/*.err'),
        ('pull', 'pull unitig from tigStore'),
        ('trace', 'find the error messages with the unitig'),
        ('test', 'test the modified unitig layout'),
        ('push', 'push the modified unitig into tigStore'),
        ('delete', 'delete specified unitig'),
        ('cut', 'cut the unitig at a given fragment ID'),
        ('shred', 'shred the unitig as a desperate way of forcing a fix'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_prefix(dir="../"):
    """
    Look for prefix.gkpStore in the upper directory.
    """
    prefix = glob(dir + "*.gkpStore")[0]
    prefix = op.basename(prefix).rsplit(".", 1)[0]

    return prefix


def get_ID(s):
    """
    unitig3.16057 => ('3', '16057')
    """
    return s.replace("unitig", "").split(".")


working = "Working"
failed = "failed"


def error(args):
    """
    %prog error backup_folder

    Find all errors in ../5-consensus/*.err and pull the error unitigs into
    backup/ folder.
    """
    p = OptionParser(error.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    backup_folder, = args
    mkdir(backup_folder)

    fw = open("errors.log", "w")

    for g in sorted(glob("../5-consensus/*.err")):
        if "partitioned" in g:
            continue

        fp = open(g)
        partID = op.basename(g).rsplit(".err", 1)[0]
        partID = int(partID.split("_")[-1])

        for row in fp:
            if row.startswith(working):
                unitigID = row.split("(")[0].split()[-1]
                continue

            if not failed in row:
                continue

            print >> fw, "\t".join(str(x) for x in (partID, unitigID))

            cmd = "{0} {1}".format(partID, unitigID)
            unitigfile = pull(cmd.split())
            cmd = "mv {0} {1}".format(unitigfile, backup_folder)
            sh(cmd)

        fp.close()


def trace(args):
    """
    %prog trace unitig{partID}.{unitigID}

    Call `grep` to get the erroneous fragment placement.
    """
    p = OptionParser(trace.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    s, = args
    partID, unitigID = get_ID(s)

    flist = glob("../5-consensus/*_{0:03d}.err".format(int(partID)))
    assert len(flist) == 1
    fp = open(flist[0])

    instate = False
    for row in fp:
        if working in row and unitigID in row:
            rows = []
            instate = True
        if instate:
            rows.append(row)
        if failed in row:
            instate = False
            if len(rows) > 20:
                ignore_line = "... ({0} lines skipped)\n".format(len(rows) - 20)
                rows = rows[:10] + [ignore_line] + rows[-10:]

    print >> sys.stderr, "".join(rows)


def cut(args):
    """
    %prog cut unitigfile fragID

    Cut the unitig at a given fragment. Run `%prog trace unitigfile` first to
    see which fragment breaks the unitig.
    """
    p = OptionParser(cut.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    s, fragID = args
    fixfile = s + ".fix"
    fp = open(s)
    fw = open(fixfile, "w")
    lines = fp.readlines()
    unitigline = lines[0]
    newunitig = "unitig -1\n"
    next10lines = lines[1:11]
    fraglines = lines[11:]
    unitigheader = [unitigline] + next10lines[:]
    altheader = [newunitig] + next10lines[:]

    # The line looks like:
    # FRG type R ident  41655818 container ...
    found = False
    for i, fline in enumerate(fraglines):
        if fline.split()[4] == fragID:
            found = True
            break

    assert i != 0, "FragID {0} is the first fragment. Something is wrong.".\
            format(fragID)
    assert found, "FragID {0} not found.".format(fragID)

    afrags = fraglines[:i]
    bfrags = fraglines[i:]
    for i, line in enumerate(unitigheader):
        if line.startswith("data.num_frags"):
            unitigheader[i] = "data.num_frags            {0}\n".\
                    format(len(afrags))
            altheader[i] = "data.num_frags            {0}\n".\
                    format(len(bfrags))

    unitigheader = "".join(unitigheader)
    altheader = "".join(altheader)

    fw.write(unitigheader)
    fw.write("".join(afrags))
    fw.write(altheader)
    fw.write("".join(bfrags))
    fw.close()

    nfrags = len(fraglines)
    nlines = sum(1 for x in open(fixfile))
    assert nfrags + 11 * 2 == nlines

    shutil.move(fixfile, s)


def shred(args):
    """
    %prog shred unitigfile

    Shred the unitig into one fragment per unitig to fix. This is the last
    resort as a desperate fix.
    """
    p = OptionParser(shred.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    s, =  args
    partID, unitigID = get_ID(s)

    # Read the header
    fixfile = s + ".fix"
    fp = open(s)
    fw = open(fixfile, "w")
    lines = fp.readlines()
    for i, line in enumerate(lines[:11]):
        if line.startswith("data.num_frags"):
            lines[i] = "data.num_frags            1\n"
    unitigline = lines[0]
    newunitig = "unitig -1\n"
    next10lines = lines[1:11]
    fraglines = lines[11:]
    unitigheader = "".join([unitigline] + next10lines)
    altheader = "".join([newunitig] + next10lines)
    frag = fraglines[0]
    fw.write(unitigheader)
    fw.write(frag)
    for frag in fraglines[1:]:
        fw.write(altheader)
        fw.write(frag)
    fw.close()
    fp.close()

    nfrags = len(fraglines)
    nlines = sum(1 for x in open(fixfile))
    assert nfrags * 12 == nlines

    shutil.move(fixfile, s)


def pull(args):
    """
    %prog pull partID unitigID

    For example, `%prog pull 5 530` will pull the utg530 from partition 5
    The layout is written to `unitig530`
    """
    p = OptionParser(pull.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    prefix = get_prefix()
    partID, unitigID = args

    cmd = CAPATH + "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "-up {0} -d layout -u {1} > unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)
    unitigfile = "unitig{0}.{1}".format(partID, unitigID)
    return unitigfile


def test(args):
    """
    %prog test unitig{partID}.{unitigID}

    For example, `%prog test unitig5.530` will test the modified `unitig530`
    """
    p = OptionParser(test.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    prefix = get_prefix()
    s, = args
    partID, unitigID = get_ID(s)

    cmd = CAPATH + "utgcns -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "{0} -T unitig{0}.{1} -V -V -V -v 2> unitig{0}.{1}.log".\
            format(partID, unitigID)

    sh(cmd)


def push(args):
    """
    %prog push unitig{partID}.{unitigID}

    For example, `%prog push unitig5.530` will push the modified `unitig530`
    and replace the one in the tigStore
    """
    p = OptionParser(push.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    prefix = get_prefix()
    s, = args
    partID, unitigID = get_ID(s)

    cmd = CAPATH + "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "-up {0} -R unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)


def delete(args):
    """
    %prog delete unitig{partID}.{unitigID}

    For example, `%prog delete unitig5.530` will delete unitig 530
    NOTICE: DELETE DOES NOT WORK!!
    """
    p = OptionParser(delete.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    prefix = get_prefix()
    s, = args
    partID, unitigID = get_ID(s)

    cmd = CAPATH + "tigStore -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += "{0} -D -u {1}".format(partID, unitigID)

    sh(cmd)


if __name__ == '__main__':
    main()
