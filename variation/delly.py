#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert delly output to BED format.
"""

import sys
import logging

from jcvi.formats.base import BaseFile, read_until, must_open
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


class DelLine (object):

    def __init__(self, line):
        args = line.strip().split("\t")
        self.seqid = args[0]
        self.start = int(args[1]) + 1
        self.end = int(args[2])
        self.size = int(args[3])
        assert self.size == self.end - self.start + 1
        self.supporting_pairs = int(args[4])
        self.avg_mapping_quality = float(args[5])
        self.accn = args[6]

    @property
    def bedline(self):
        return "\t".join(str(x) for x in (self.seqid,
                          self.start - 1, self.end, self.accn,
                          self.supporting_pairs, "+"))


class Delly (BaseFile):

    def __init__(self, filename):
        super(Delly, self).__init__(filename)

    def __iter__(self):
        fp = must_open(self.filename)
        while True:
            read_until(fp, "-----")
            nextline = fp.readline()
            nextline = fp.readline()
            if not nextline.strip():
                break
            d = DelLine(nextline)
            yield d

    def write_bed(self, bedfile="stdout"):
        fw = must_open(bedfile, "w")
        for d in self:
            print >> fw, d.bedline
        logging.debug("File written to `{0}`.".format(bedfile))


def main():

    actions = (
        ('bed', 'Convert del.txt to del.bed'),
        ('mito', 'Find mito deletions in BAM'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed del.txt

    Convert `del.txt` to BED format. DELLY manual here:
    <http://www.embl.de/~rausch/delly.html>

    Deletion:
    chr, start, end, size, #supporting_pairs, avg._mapping_quality, deletion_id
    chr1, 10180, 10509, 329, 75, 15.8667, Deletion_Sample_00000000
    """
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    delt, = args
    dt = Delly(delt)
    dt.write_bed("del.bed")


def mito(args):
    """
    %prog mito input.bam chrM chrM.fa

    Identify mitochondrial deletions.
    """
    from jcvi.formats.sam import get_minibam

    p = OptionParser(mito.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    inputbam, region, chrMfa = args
    minibam = get_minibam(inputbam, region)

    realign = minibam.rsplit(".", 1)[0] + ".realign"
    cmd = "speedseq realign -v -t {}".format(opts.cpus)
    cmd += " -o {} {} {}".format(realign, chrMfa, minibam)
    sh(cmd)


if __name__ == '__main__':
    main()
