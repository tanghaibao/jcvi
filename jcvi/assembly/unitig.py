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
from __future__ import print_function

import os.path as op
import sys
import shutil
import logging

from jcvi.formats.base import BaseFile
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, mkdir, glob


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
            fid = fline.split()[4]
            if fid == fragID:
                found = True
                break

        assert i, "FragID {0} is the first fragment. Something is wrong.".\
                format(fragID)
        assert found, "FragID {0} not found.".format(fragID)

        afrags = fraglines[:i]
        bfrags = fraglines[i:]
        self.parts = [afrags, bfrags]

    def pop(self, black):
        fraglines = self.fraglines
        good, save = [], []
        for fline in fraglines:
            fid = fline.split()[4]
            if fid in black:
                save.append(fline)
            else:
                good.append(fline)
        assert good, "good empty"
        assert save, "save empty"
        self.parts = [good, save]

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
        ('cut', 'cut the unitig at a given fragment ID'),
        ('shred', 'shred the unitig as a desperate way of forcing a fix'),
        ('cnsfix', 'parse consensus-fix.out to fix unitigs'),
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
    unitig1.3.16057 => ('1', '3', '16057')
    """
    return s.replace("unitig", "").split(".")


def cnsfix(args):
    """
    %prog cnsfix consensus-fix.out.FAILED > blacklist.ids

    Parse consensus-fix.out to extract layouts for fixed unitigs. This will
    mark all the failed fragments detected by utgcnsfix and pop them out of the
    existing unitigs.
    """
    from jcvi.formats.base import read_block

    p = OptionParser(cnsfix.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    cnsfixout, = args
    fp = open(cnsfixout)
    utgs = []
    saves = []
    for header, contents in read_block(fp, "Evaluating"):
        contents = list(contents)
        utg = header.split()[2]
        utgs.append(utg)
        # Look for this line:
        #   save fragment idx=388 ident=206054426 for next pass
        for c in contents:
            if not c.startswith("save"):
                continue
            ident = c.split()[3].split("=")[-1]
            saves.append(ident)
    print("\n".join(saves))


working = "Working"
failed = "failed"


def error(args):
    """
    %prog error version backup_folder

    Find all errors in ../5-consensus/*.err and pull the error unitigs into
    backup/ folder.
    """
    p = OptionParser(error.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    version, backup_folder = args
    mkdir(backup_folder)

    fw = open("errors.log", "w")

    seen = set()
    for g in glob("../5-consensus/*.err"):
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

            uu = (version, partID, unitigID)
            if uu in seen:
                continue
            seen.add(uu)

            print("\t".join(str(x) for x in (partID, unitigID)), file=fw)

            s = [str(x) for x in uu]
            unitigfile = pull(s)
            cmd = "mv {0} {1}".format(unitigfile, backup_folder)
            sh(cmd)

        fp.close()

    logging.debug("A total of {0} unitigs saved to {1}.".\
                 format(len(seen), backup_folder))


def trace(args):
    """
    %prog trace unitig{version}.{partID}.{unitigID}

    Call `grep` to get the erroneous fragment placement.
    """
    p = OptionParser(trace.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    s, = args
    version, partID, unitigID = get_ID(s)

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

    print("".join(rows), file=sys.stderr)


def cut(args):
    """
    %prog cut unitigfile fragID

    Cut the unitig at a given fragment. Run `%prog trace unitigfile` first to
    see which fragment breaks the unitig.
    """
    from jcvi.formats.base import SetFile

    p = OptionParser(cut.__doc__)
    p.add_option("-s", dest="shredafter", default=False, action="store_true",
                 help="Shred fragments after the given fragID [default: %default]")
    p.add_option("--notest", default=False, action="store_true",
                 help="Do not test the unitigfile after edits [default: %default]")
    p.add_option("--blacklist",
                 help="File that contains blacklisted fragments to be popped "
                      "[default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    s, fragID = args
    u = UnitigLayout(s)
    blacklist = opts.blacklist
    black = SetFile(blacklist) if blacklist else None

    if opts.shredafter:
        u.shredafter(fragID)
    elif black:
        assert fragID == "0", "Must set fragID to 0 when --blacklist is on"
        u.pop(black)
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
        sys.exit(not p.print_help())

    s, = args
    u = UnitigLayout(s)
    u.shred()
    u.print_to_file(inplace=True)


def pull(args):
    """
    %prog pull version partID unitigID

    For example, `%prog pull 5 530` will pull the utg530 from partition 5
    The layout is written to `unitig530`
    """
    p = OptionParser(pull.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    prefix = get_prefix()
    version, partID, unitigID = args
    s = ".".join(args)

    cmd = "tigStore"
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore".format(prefix)
    cmd += " {0} -up {1} -d layout -u {2}".format(version, partID, unitigID)

    unitigfile = "unitig" + s
    sh(cmd, outfile=unitigfile)

    return unitigfile


def test(args):
    """
    %prog test unitig{version}.{partID}.{unitigID}

    For example, `%prog test unitig1.5.530` will test the modified `unitig530`
    """
    p = OptionParser(test.__doc__)
    p.set_verbose(help="Turn on verbose debugging")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    prefix = get_prefix()
    s, = args
    version, partID, unitigID = get_ID(s)

    cmd = "utgcns"
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore".format(prefix)
    cmd += " {0} {1} -T {2}".format(version, partID, s)
    if opts.verbose:
        cmd += " -V -V"
    cmd += " -V -v 2> {0}.log".format(s)

    sh(cmd)

    # Show log
    cmd = "tail {0}.log".format(s)
    sh(cmd)


def push(args):
    """
    %prog push unitig{version}.{partID}.{unitigID}

    For example, `%prog push unitig5.530` will push the modified `unitig530`
    and replace the one in the tigStore
    """
    p = OptionParser(push.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    s, = args
    prefix = get_prefix()
    version, partID, unitigID = get_ID(s)

    cmd = "tigStore"
    cmd += " -g ../{0}.gkpStore -t ../{0}.tigStore".format(prefix)
    cmd += " {0} -up {1} -R {2}".format(version, partID, s)

    sh(cmd)


def pushall(args):
    """
    %prog pushall unitig1.*

    Push a bunch of unitig layout changes.
    """
    p = OptionParser(pushall.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    flist = args
    for f in flist:
        if f.endswith(".log"):
            continue

        push([f])


if __name__ == '__main__':
    main()
