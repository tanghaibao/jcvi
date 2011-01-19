#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import sys

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.base import LineFile
from jcvi.apps.base import ActionDispatcher, debug
from jcvi.utils.range import range_distance
debug()


class CasTabLine (LineFile):
    """
    The table generate by command `assembly_table -n -s -p` 
    from clcbio assembly cell
    """
    def __init__(self, line):
        args = line.split()
        self.readnum = int(args[0])
        self.readname = args[1]
        self.readlen = int(args[2])
        # 0-based indexing
        self.readstart = int(args[3])
        if self.readstart >=0: self.readstart += 1
        
        self.readstop = int(args[4])
        self.refnum = int(args[5])
        
        self.refstart = int(args[6])
        if self.refstart >=0: self.refstart += 1

        self.refstop = int(args[7])
        self.is_reversed = (int(args[8])==1)
        
        self.orientation = '-' if self.is_reversed else '+'

        self.nummatches = int(args[9])
        self.is_paired = (int(args[10])==1)
        self.score = int(args[11])

    def __str__(self):
        return "\t".join(str(x) for x in (self.readname, self.refnum, 
            self.refstart-1, self.refstop, self.score, self.orientation))

def main():
    
    actions = (
        ('bed', 'convert cas tabular output to bed format'),
        ('summary', 'print summary of cas tabular output'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def bed(args):
    """
    %prog bed cas_tabbed

    convert the format into bed format
    """
    
    p = OptionParser(bed.__doc__)
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    castabfile = args[0]
    fp = open(castabfile)
    for row in fp:
        b = CasTabLine(row)
        if b.readstart!=-1:
            print b


def summary(args):
    """
    %prog bed cas_tabbed
    
    report summary of the cas tabular results, how many paired ends mapped, avg
    distance between paired ends, etc
    """
    p = OptionParser(summary.__doc__)
    p.add_option("--cutoff", dest="cutoff", default=300000, type="int",
            help="distance to call valid links between PE [default: %default]")
    opts, args = p.parse_args(args)

    if len(args)!=1:
        sys.exit(p.print_help())

    cutoff = opts.cutoff

    castabfile = args[0]
    fp = open(castabfile)
    data = [CasTabLine(row) for row in fp]
    num_fragments, num_pairs = 0, 0
    num_linked = 0
    for pe, lines in groupby(data, key=lambda x: x.readname.split("/")[0]):   
        lines = list(lines)
        if len(lines)==2: 
            a, b = lines
            dist = range_distance((a.refnum, a.refstart, a.refstop), (b.refnum,
                b.refstart, b.refstop))
            num_pairs += 1
            if 0 <= dist <= cutoff:
                num_linked += 1
        else:
            num_fragments += 1
    
    print "%d fragments, %d pairs" % (num_fragments, num_pairs)
    print "%d pairs are linked (cutoff=%d)" % (num_linked, cutoff)


if __name__ == '__main__':
    main()
