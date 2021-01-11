#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
POSMAP (POSitional MAPping) files are part of the Celera Assembler output.

Specs:
http://sourceforge.net/apps/mediawiki/wgs-assembler/index.php?title=POSMAP
"""

import os.path as op
import sys
import csv
import logging

from collections import namedtuple, defaultdict
from itertools import groupby

from jcvi.formats.base import BaseFile, LineFile
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


class MateLine(object):
    def __init__(self, line):
        args = line.split()
        self.read1 = args[0]
        self.read2 = args[1]
        self.library = args[2]

    def __str__(self):
        return "\t".join((self.read1, self.read2, self.library))


LibraryTag = "library"


class LibraryLine(object):
    def __init__(self, line):
        args = line.split()
        assert args[0] == LibraryTag

        self.library = args[1]
        self.min = int(args[2])
        self.max = int(args[3])
        assert self.min <= self.max

        self.sd = (self.max - self.min) / 4
        self.mean = (self.max + self.min) / 2

    def zscore(self, distance):
        return (distance - self.mean) / self.sd


class MatesFile(BaseFile, dict):
    def __init__(self, filename):
        super(MatesFile, self).__init__(filename)

        libraries = {}

        fp = open(filename)
        for row in fp:
            atoms = row.split()
            if atoms[0] == LibraryTag:
                lib = LibraryLine(row)
                libraries[lib.library] = lib
            else:
                mate = MateLine(row)
                lib = libraries[mate.library]
                r1, r2 = mate.read1, mate.read2
                self[r1] = (r2, lib)
                self[r2] = (r1, lib)

        logging.debug("Libraries: {0}, Mates: {1}".format(len(libraries), len(self)))

        self.libraries = libraries


MatesLine = namedtuple("MatesLine", "firstReadID secondReadID mateStatus")


class Mates(object):
    def __init__(self, filename):
        fp = csv.reader(open(filename), delimiter="\t")
        for row in fp:
            MatesLine._make(row)


class FrgScfLine(object):
    def __init__(self, row):
        atoms = row.split()
        self.fragmentID = atoms[0]
        self.scaffoldID = atoms[1]
        self.begin = int(atoms[2]) + 1  # convert to 1-based
        self.end = int(atoms[3])
        self.orientation = "+" if atoms[4] == "f" else "-"

    @property
    def bedline(self):
        s = "\t".join(
            str(x)
            for x in (
                self.scaffoldID,
                self.begin - 1,
                self.end,
                self.fragmentID,
                "na",
                self.orientation,
            )
        )
        return s


class FrgScf(object):
    def __init__(self, filename):
        fp = csv.reader(open(filename), delimiter="\t")
        for row in fp:
            b = FrgScfLine(row)


class Posmap(LineFile):

    # dispatch based on filename
    mapping = {
        "mates": Mates,
        "frgscf": FrgScf,
    }

    def __init__(self, filename):
        super(Posmap, self).__init__(filename)

    def parse(self):
        filename = self.filename
        suffix = filename.rsplit(".", 1)[-1]
        assert suffix in self.mapping, "`{0}` unknown format".format(filename)

        # dispatch to the proper handler
        klass = self.mapping[suffix]
        return klass(filename)


def main():

    actions = (
        ("bed", "convert to bed format"),
        ("index", "index the frgscf.posmap.sorted file"),
        ("query", "query region from frgscf index"),
        ("dup", "estimate level of redundancy based on position collision"),
        ("reads", "report read counts per scaffold (based on frgscf)"),
        ("pairs", "report insert statistics for read pairs"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def index(args):
    """
    %prog index frgscf.sorted

    Compress frgscffile.sorted and index it using `tabix`.
    """
    p = OptionParser(index.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (frgscffile,) = args
    gzfile = frgscffile + ".gz"
    cmd = "bgzip -c {0}".format(frgscffile)

    if not op.exists(gzfile):
        sh(cmd, outfile=gzfile)

    tbifile = gzfile + ".tbi"
    # Sequence, begin, end in 2, 3, 4-th column, respectively
    cmd = "tabix -s 2 -b 3 -e 4 {0}".format(gzfile)

    if not op.exists(tbifile):
        sh(cmd)


def query(args):
    """
    %prog query frgscf.sorted scfID:start-end

    Query certain region to get frg placement, using random access. Build index
    if not present.
    """
    p = OptionParser(query.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    frgscffile, region = args
    gzfile = frgscffile + ".gz"
    tbifile = gzfile + ".tbi"
    outfile = region + ".posmap"

    if not (op.exists(gzfile) and op.exists(tbifile)):
        index(frgscffile)

    assert op.exists(gzfile) and op.exists(tbifile)

    cmd = "tabix {0} {1}".format(gzfile, region)
    sh(cmd, outfile=outfile)


def reads(args):
    """
    %prog reads frgscffile

    Report read counts per scaffold (based on frgscf).
    """
    p = OptionParser(reads.__doc__)
    p.add_option(
        "-p",
        dest="prefix_length",
        default=4,
        type="int",
        help="group the reads based on the first N chars [default: %default]",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (frgscffile,) = args
    prefix_length = opts.prefix_length

    fp = open(frgscffile)
    keyfn = lambda: defaultdict(int)
    counts = defaultdict(keyfn)
    for row in fp:
        f = FrgScfLine(row)
        fi = f.fragmentID[:prefix_length]
        counts[f.scaffoldID][fi] += 1

    for scf, count in sorted(counts.items()):
        print(
            "{0}\t{1}".format(
                scf, ", ".join("{0}:{1}".format(*x) for x in sorted(count.items()))
            )
        )


def bed(args):
    """
    %prog bed frgscffile

    Convert the frgscf posmap file to bed format.
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (frgscffile,) = args
    bedfile = frgscffile.rsplit(".", 1)[0] + ".bed"
    fw = open(bedfile, "w")

    fp = open(frgscffile)
    for row in fp:
        f = FrgScfLine(row)
        print(f.bedline, file=fw)

    logging.debug("File written to `%s`.", bedfile)

    return bedfile


def dup(args):
    """
    %prog dup frgscffile

    Use the frgscf posmap file as an indication of the coverage of the library.
    Large insert libraries are frequently victims of high levels of redundancy.
    """
    p = OptionParser(dup.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    (frgscffile,) = args

    fp = open(frgscffile)
    data = [FrgScfLine(row) for row in fp]
    # we need to separate forward and reverse reads, because the position
    # collisions are handled differently
    forward_data = [x for x in data if x.orientation == "+"]
    reverse_data = [x for x in data if x.orientation == "-"]

    counts = defaultdict(int)
    key = lambda x: (x.scaffoldID, x.begin)
    forward_data.sort(key=key)
    for k, data in groupby(forward_data, key=key):
        data = list(data)
        count = len(data)
        counts[count] += 1

    key = lambda x: (x.scaffoldID, x.end)
    reverse_data.sort(key=key)
    for k, data in groupby(forward_data, key=key):
        data = list(data)
        count = len(data)
        counts[count] += 1

    prefix = frgscffile.split(".")[0]
    print("Duplication level in `{0}`".format(prefix), file=sys.stderr)
    print("=" * 40, file=sys.stderr)
    for c, v in sorted(counts.items()):
        if c > 10:
            break
        label = "unique" if c == 1 else "{0} copies".format(c)
        print("{0}: {1}".format(label, v), file=sys.stderr)


if __name__ == "__main__":
    main()
