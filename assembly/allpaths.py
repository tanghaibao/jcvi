#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Subroutines to aid ALLPATHS-LG assembly.
"""

import os.path as op
import sys
import logging

from glob import glob
from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('prepare', 'prepare ALLPATHS csv files and run script'),
        ('log', 'prepare a log of created files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare "B. oleracea" *.fastq

    Scans the current folder looking for input fastq files that look like:
    PE-376.fastq (paired end)
    MP-3000.fastq (mate pairs)
    TT-3000.fastq (mate pairs, but from 454 data)

    Then create "in_groups.csv" and "in_libs.csv". Finally write a generic run
    script.
    """
    from jcvi.utils.table import write_csv

    p = OptionParser(prepare.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    organism_name = args[0]
    project_name = "".join(x[0] for x in organism_name.split()).upper()
    fnames = sorted(glob("*.fastq") if len(args) == 1 else args[1:])

    groupheader = "group_name library_name file_name".split()
    libheader = "library_name project_name organism_name type paired "\
        "frag_size frag_stddev insert_size insert_stddev read_orientation "\
        "genomic_start genomic_end".split()
    groupcontents = []
    libs = []
    for file_name in fnames:
        group_name, ext = op.splitext(file_name)
        library_name = "-".join(group_name.split("-")[:2])
        groupcontents.append((group_name, library_name, file_name))
        if library_name not in libs:
            libs.append(library_name)

    libcontents = []
    for library_name in libs:
        pf, size = library_name.split("-")
        size = int(size)
        stddev = size / 5
        type = "fragment" if pf == "PE" else "jumping"
        paired = 1
        frag_size = size if type == "fragment" else ""
        frag_stddev = stddev if type == "fragment" else ""
        insert_size = size if type != "fragment" else ""
        insert_stddev = stddev if type != "fragment" else ""
        read_orientation = "outward" if pf == "MP" else "inward"
        genomic_start, genomic_end = 0, 0
        libcontents.append((library_name, project_name, organism_name, type, \
            paired, frag_size, frag_stddev, insert_size, insert_stddev, \
            read_orientation, genomic_start, genomic_end))

    write_csv(groupheader, groupcontents, "in_groups.csv", tee=True)
    logging.debug("`in_group.csv` created (# of groups = {0}).".\
        format(len(groupcontents)))

    write_csv(libheader, libcontents, "in_libs.csv", tee=True)
    logging.debug("`in_libs.csv` created (# of libs = {0}).".\
        format(len(libcontents)))


def log(args):
    """
    %prog log logfile

    Prepare a log of created files, ordered by their creation data. The purpose
    for this script is to touch these files sequentially to reflect their build
    order. On the JCVI scratch area, the files are touched regularly to avoid
    getting deleted, losing their respective timestamps. However, this created a
    problem for the make system adopted by ALLPATHS.

    An example block to be extracted ==>
    [PC] Calling PreCorrect to create 2 file(s):
    [PC]
    [PC] $(RUN)/frag_reads_prec.fastb
    [PC] $(RUN)/frag_reads_prec.qualb
    [PC]
    [PC] from 2 file(s):
    [PC]
    [PC] $(RUN)/frag_reads_filt.fastb
    [PC] $(RUN)/frag_reads_filt.qualb
    """
    from jcvi.algorithms.graph import nx, topological_sort

    p = OptionParser(log.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    g = nx.DiGraph()

    logfile, = args
    fp = open(logfile)
    row = fp.readline()
    incalling = False
    basedb = {}

    while row:
        atoms = row.split()
        if len(atoms) < 3:
            row = fp.readline()
            continue

        tag, token, trailing = atoms[0], atoms[1], atoms[-1]
        if trailing == 'file(s):':
            numfiles = int(atoms[-2])
            row = fp.readline()
            assert row.strip() == tag

            if token == "Calling" and not incalling:
                createfiles = []
                for i in xrange(numfiles):
                    row = fp.readline()
                    createfiles.append(row.split()[-1])
                incalling = True

            if token == "from" and incalling:
                fromfiles = []
                for i in xrange(numfiles):
                    row = fp.readline()
                    fromfiles.append(row.split()[-1])

                for a in fromfiles:
                    for b in createfiles:
                        ba, bb = op.basename(a), op.basename(b)
                        basedb[ba] = a
                        basedb[bb] = b
                        g.add_edge(ba, bb)

                incalling = False

        if token == "ln":
            fromfile, createfile = atoms[-2:]
            ba, bb = op.basename(fromfile), op.basename(createfile)
            #print ba, "-->", bb
            if ba != bb:
                g.add_edge(ba, bb)

        row = fp.readline()

    ts = [basedb[x] for x in topological_sort(g) if x in basedb]
    print "\n".join(ts)


if __name__ == '__main__':
    main()
