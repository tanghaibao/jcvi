#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using CD-HIT (CD-HIT-454) in particular to remove duplicate reads.
"""

import os.path as op
import sys
import logging

from collections import defaultdict

from jcvi.formats.base import LineFile, read_block
from jcvi.apps.base import OptionParser, ActionDispatcher, sh


class ClstrLine (object):
    """
    Lines like these:
    0       12067nt, >LAP012517... at -/99.85%
    1       15532nt, >MOL158919... *
    2       15515nt, >SES069071... at +/99.85%
    """
    def __init__(self, row):
        a, b = row.split('>', 1)
        a = a.split("nt")[0]
        sid, size = a.split()
        self.size = int(size)
        self.name = b.split("...")[0]
        self.rep = (row.rstrip()[-1] == '*')


class ClstrFile (LineFile):

    def __init__(self, filename):
        super(ClstrFile, self).__init__(filename)
        assert filename.endswith(".clstr")

        fp = open(filename)
        for clstr, members in read_block(fp, ">"):
            self.append([ClstrLine(x) for x in members])

    def iter_sizes(self):
        for members in self:
            yield len(members)

    def iter_reps(self):
        for i, members in enumerate(self):
            for b in members:
                if b.rep:
                    yield i, b.name

    def iter_reps_prefix(self, prefix=3):
        for i, members in enumerate(self):
            d = defaultdict(list)
            for b in members:
                pp = b.name[:prefix]
                d[pp].append(b)

            for pp, members_with_same_pp in sorted(d.items()):
                yield i, max(members_with_same_pp, \
                             key=lambda x: x.size).name


def main():

    actions = (
        ('ids', 'get the representative ids from clstr file'),
        ('deduplicate', 'use `cd-hit-est` to remove duplicate reads'),
        ('summary', 'parse cdhit.clstr file to get distribution of cluster sizes'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def ids(args):
    """
    %prog ids cdhit.clstr

    Get the representative ids from clstr file.
    """
    p = OptionParser(ids.__doc__)
    p.add_option("--prefix", type="int",
                 help="Find rep id for prefix of len [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clstrfile, = args
    cf = ClstrFile(clstrfile)
    prefix = opts.prefix
    if prefix:
        reads = list(cf.iter_reps_prefix(prefix=prefix))
    else:
        reads = list(cf.iter_reps())

    nreads = len(reads)
    idsfile = clstrfile.replace(".clstr", ".ids")
    fw = open(idsfile, "w")
    for i, name in reads:
        print >> fw, "\t".join(str(x) for x in (i, name))

    logging.debug("A total of {0} unique reads written to `{1}`.".\
            format(nreads, idsfile))
    fw.close()

    return idsfile


def summary(args):
    """
    %prog summary cdhit.clstr

    Parse cdhit.clstr file to get distribution of cluster sizes.
    """
    from jcvi.graphics.histogram import loghistogram

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clstrfile, = args
    cf = ClstrFile(clstrfile)
    data = list(cf.iter_sizes())
    loghistogram(data, summary=True)


def deduplicate(args):
    """
    %prog deduplicate fastafile

    Wraps `cd-hit-est` to remove duplicate sequences.
    """
    p = OptionParser(deduplicate.__doc__)
    p.set_align(pctid=98)
    p.add_option("--reads", default=False, action="store_true",
                 help="Use `cd-hit-454` to deduplicate [default: %default]")
    p.add_option("--samestrand", default=False, action="store_true",
                 help="Enforce same strand alignment")
    p.set_home("cdhit")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastafile, = args
    identity = opts.pctid / 100.

    cmd = "cd-hit-454" if opts.reads else "cd-hit-est"
    cmd = op.join(opts.cdhit_home, cmd)
    cmd += " -c {0}".format(identity)
    cmd += " -d 0"  # include complete defline
    if opts.samestrand:
        cmd += " -r 0"
    cmd += " -M 0 -T {0} -i {1} -o {1}.cdhit".format(opts.cpus, fastafile)
    sh(cmd)

    dd = fastafile + ".cdhit"
    return dd


if __name__ == '__main__':
    main()
