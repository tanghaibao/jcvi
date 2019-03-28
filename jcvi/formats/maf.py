#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
MAF format specification:
<http://genome.ucsc.edu/FAQ/FAQformat#format5>
"""
from __future__ import print_function

import sys

from bx import interval_index_file
from bx.align import maf

from jcvi.formats.base import BaseFile
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update
from jcvi.apps.lastz import blastz_score_to_ncbi_expectation, \
            blastz_score_to_ncbi_bits


class Maf (BaseFile, dict):

    def __init__(self, filename, index=False):
        super(Maf, self).__init__(filename)

        indexfile = filename + ".idx"
        if index:
            if need_update(filename, indexfile):
                self.build_index(filename, indexfile)

            self.index = maf.Index(filename, indexfile)

        fp = open(filename)
        self.reader = maf.Reader(fp)

    def build_index(self, filename, indexfile):
        """
        Recipe from Brad Chapman's blog
        <http://bcbio.wordpress.com/2009/07/26/sorting-genomic-alignments-using-python/>
        """
        indexes = interval_index_file.Indexes()
        in_handle = open(filename)

        reader = maf.Reader(in_handle)
        while True:
            pos = reader.file.tell()
            rec = next(reader)
            if rec is None:
                break
            for c in rec.components:
                indexes.add(c.src, c.forward_strand_start,
                        c.forward_strand_end, pos, max=c.src_size )

        index_handle = open(indexfile, "w")
        indexes.write(index_handle)
        index_handle.close()


def main():

    actions = (
        ('bed', 'convert MAF to BED format'),
        ('blast', 'convert MAF to BLAST tabular format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed maffiles > out.bed

    Convert a folder of maf alignments to the bed features
    then useful to check coverage, etc.
    """
    p = OptionParser(bed.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    flist = args
    prefix = flist[0].split(".")[0]

    j = 0
    for f in flist:
        reader = Maf(f).reader
        for rec in reader:
            a, b = rec.components

            for a, tag in zip((a, b), "ab"):
                name = "{0}_{1:07d}{2}".format(prefix, j, tag)
                print("\t".join(str(x) for x in (a.src, a.forward_strand_start, \
                        a.forward_strand_end, name)))

            j += 1


def alignment_details(a, b):
    nmatch = 0
    nmismatch = 0
    ngaps = 0

    assert len(a) == len(b)
    l = len(a)

    for i in range(l):
        if a[i] == b[i]:
            nmatch += 1
        elif a[i] == "-" or b[i] == "-":
            ngaps += 1
        else:
            nmismatch += 1

    pctid = 100. * nmatch / l
    return pctid, nmismatch, ngaps


def maf_to_blast8(f):
    reader = Maf(f).reader
    for rec in reader:
        a, b = rec.components
        query = a.src
        subject = b.src
        qstart = a.forward_strand_start
        qstop = a.forward_strand_end
        sstart = b.forward_strand_start
        sstop = b.forward_strand_end
        score = rec.score

        evalue = blastz_score_to_ncbi_expectation(score)
        score = blastz_score_to_ncbi_bits(score)
        evalue, score = "{0:.2g}".format(evalue), "{0:.1f}".format(score)
        hitlen = len(a.text)

        pctid, nmismatch, ngaps = alignment_details(a.text, b.text)
        print("\t".join(str(x) for x in (query, subject, pctid, hitlen,
            nmismatch, ngaps, qstart, qstop, sstart, sstop,
            evalue, score)))


def blast(args):
    '''
    %prog blast maffiles > out.blast

    From a folder of .maf files, generate .blast file with tabular format.
    '''
    p = OptionParser(blast.__doc__)
    opts, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    flist = args

    for f in flist:
        maf_to_blast8(f)


if __name__ == '__main__':
    main()
