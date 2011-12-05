#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Wrapper to call tigStore and utgcns, for debugging failed utgcns runs.

See full commands:
http://sf.net/apps/mediawiki/wgs-assembler/index.php?title=Unitig_Consensus_Failures_in_CA_6

It is expected to be executed within 5-consensus/ folder, or (like me) executed
in a fix_unitig folder in the assembly folder. In any case, from the current
dir, it needs to get access to ../5-consensus.
"""

import os
import os.path as op
import sys
import shutil
import logging

from glob import glob
from optparse import OptionParser

from jcvi.formats.base import BaseFile
from jcvi.apps.command import CAPATH
from jcvi.apps.base import ActionDispatcher, sh, mkdir, debug
debug()


class UnitigLayout (BaseFile):
    """
    Reads in a unitig layout file.
    """
    def __init__(self, filename):

        super(UnitigLayout, self).__init__(filename)

        fp = open(filename)
        lines = fp.readlines()
        unitigline = lines[0]
        assert unitigline.startswith("unitig")

        tag, self.unitig = unitigline.rsplit(None, 1)

        next10lines = lines[1:11]
        fraglines = lines[11:]
        # Sanity check to make sure every line is frag
        for fl in fraglines:
            assert fl.startswith("FRG")

        self.num_frags = len(fraglines)
        self.fraglines = fraglines
        self.header = [unitigline] + next10lines
        self.parts = [self.fraglines]
        fp.close()

    def cut(self, fragID):

        fraglines = self.fraglines

        # The line looks like:
        # FRG type R ident  41655818 container ...
        found = False
        for i, fline in enumerate(fraglines):
            if fline.split()[4] == fragID:
                found = True
                break

        assert i, "FragID {0} is the first fragment. Something is wrong.".\
                format(fragID)
        assert found, "FragID {0} not found.".format(fragID)

        afrags = fraglines[:i]
        bfrags = fraglines[i:]
        self.parts = [afrags, bfrags]

    def shred(self):
        self.parts = [[x] for x in self.fraglines]

    def shredafter(self, fragID):
        self.cut(fragID)
        self.parts = [self.parts[0]] + [[x] for x in self.parts[1]]

    def get_header(self, unitig=-1, num_frags=1):
        header = self.header[:]
        for i, line in enumerate(header):
            if line.startswith("unitig"):
                header[i] = "unitig {0}\n".format(unitig)
            if line.startswith("data.num_frags"):
                header[i] = "data.num_frags            {0}\n".format(num_frags)
        return header

    def print_to_file(self, inplace=False):
        fixfile = self.filename + ".fix"
        fw = open(fixfile, "w")

        part = self.parts[0]
        header = self.get_header(unitig=self.unitig, num_frags=len(part))
        fw.write("".join(header))
        fw.write("".join(part))
        for part in self.parts[1:]:
            header = self.get_header(unitig=-1, num_frags=len(part))
            fw.write("".join(header))
            fw.write("".join(part))
        fw.close()

        logging.debug("Writing {0} parts.".format(len(self.parts)))

        nlines = sum(1 for x in open(fixfile))
        expectlines = self.num_frags + 11 * len(self.parts)
        assert expectlines == nlines, \
                "Expecting {0} line, you have {1}.".format(expectlines, nlines)

        if inplace:
            shutil.move(fixfile, self.filename)


def main():

    actions = (
        ('error', 'find all errors in ../5-consensus/*.err'),
        ('pull', 'pull unitig from tigStore'),
        ('trace', 'find the error messages with the unitig'),
        ('test', 'test the modified unitig layout'),
        ('push', 'push the modified unitig into tigStore'),
        ('pushall', 'push all the modified unitigs into tigStore'),
        ('delete', 'delete specified unitig'),
        ('cut', 'cut the unitig at a given fragment ID'),
        ('shred', 'shred the unitig as a desperate way of forcing a fix'),
        ('shredafter', 'shred the unitig after a given fragment ID'),
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

            if not failed.upper() in row.upper():
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
    p.add_option("-s", dest="shredafter", default=False, action="store_true",
                 help="Shred fragments after the given fragID [default: %default]")
    p.add_option("--notest", default=False, action="store_true",
                 help="Do not test the unitigfile after edits [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    s, fragID = args
    u = UnitigLayout(s)
    if opts.shredafter:
        u.shredafter(fragID)
    else:
        u.cut(fragID)
    u.print_to_file(inplace=True)

    if not opts.notest:
        test([s])


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
    u = UnitigLayout(s)
    u.shred()
    u.print_to_file(inplace=True)


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

    cmd = CAPATH("tigStore")
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore 1".format(prefix)
    cmd += " -up {0} -d layout -u {1} > unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)
    unitigfile = "unitig{0}.{1}".format(partID, unitigID)
    return unitigfile


def test(args):
    """
    %prog test unitig{partID}.{unitigID}

    For example, `%prog test unitig5.530` will test the modified `unitig530`
    """
    p = OptionParser(test.__doc__)
    p.add_option("--verbose", default=False, action="store_true",
                 help="Turn on verbose debugging [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    prefix = get_prefix()
    s, = args
    partID, unitigID = get_ID(s)

    cmd = CAPATH("utgcns")
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore 1".format(prefix)
    cmd += " {0} -T {1}".format(partID, s)
    if opts.verbose:
        cmd += " -V -V"
    cmd += " -V -v 2> {0}.log".format(s)

    sh(cmd)

    # Show log
    cmd = "tail {0}.log".format(s)
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

    cmd = CAPATH("tigStore")
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore 1".format(prefix)
    cmd += " -up {0} -R unitig{0}.{1}".format(partID, unitigID)

    sh(cmd)


def pushall(args):
    """
    %prog pushall .

    Push a bunch of unitig layout changes.
    """
    p = OptionParser(pushall.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    cwd, = args
    assert cwd == ".", "Must execute the command in current folder"

    flist = glob("unitig*")
    for f in flist:
        if f.endswith(".log"):
            continue

        push([f])


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

    cmd = CAPATH("tigStore")
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore 1 ".format(prefix)
    cmd += " {0} -D -u {1}".format(partID, unitigID)

    sh(cmd)


if __name__ == '__main__':
    main()
