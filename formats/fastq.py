#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Processing fastq files
"""

import os.path as op
import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


class FastqRecord (object):
    def __init__(self, fh, offset=0, key=None):
        self.name = fh.readline().split()[0]
        if not self.name: return
        self.seq = fh.readline().rstrip()
        self.l3 = fh.readline().rstrip()
        self.qual = fh.readline().rstrip()
        if offset != 0:
            self.qual = "".join(chr(ord(x) + offset) for x in self.qual)
        self.id = key(self.name) if key else self.name

    def __str__(self):
        return "\n".join((self.name, self.seq, "+", self.qual))
    

def main():

    actions = (
        ('pair', 'pair up two fastq files and combine pairs'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pair(args):
    """
    %prog pair 1.fastq 2.fastq frags.fastq pairs.fastq 

    pair up the records in 1.fastq and 2.fastq, pairs are indicated by trailing
    "/1" and "/2". If using raw sequences, this is trivial, since we can just
    iterate one at a time for both files; however if two files are not matching,
    (e.g. due to trimming), we need a third fastq that provides the order.
    """
    p = OptionParser(pair.__doc__)
    p.add_option("-Q", dest="infastq", default="illumina", 
            help="input fastq [default: %default]")
    p.add_option("-q", dest="outfastq", default="sanger",
            help="output fastq format [default: %default]")
    p.add_option("-r", dest="ref", default=None,
            help="a reference fastq that provides order")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(p.print_help())

    afastq, bfastq, frags, pairs = args

    assert op.exists(afastq) and op.exists(bfastq)
    assert not op.exists(frags) and not op.exists(pairs)

    qual_offset = lambda x: 33 if x=="sanger" else 64 
    in_offset = qual_offset(opts.infastq)
    out_offset = qual_offset(opts.outfastq)
    offset = out_offset - in_offset
    ref = opts.ref

    ah = open(afastq)
    bh = open(bfastq)
    if ref: 
        rh = open(ref) 
        totalsize = op.getsize(ref) 
    else:
        totalsize = op.getsize(afastq) 

    from progressbar import ProgressBar
    bar = ProgressBar(maxval=totalsize).start()
    strip_name = lambda x: x.rsplit("/", 1)[0]

    j = 0
    if ref:
        while True:
            rec = FastqRecord(rh, offset=0, key=strip_name)
            pos = rh.tell()
            bar.update(pos)
    else:
        # TODO: unimplemented
        pass
    
    bar.finish()
    sys.stdout.write("\n")


if __name__ == '__main__':
    main()
