#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Patch the sequences of one assembly using sequences from another assembly. This
is tested on merging the medicago WGS assembly with the clone-by-clone assembly.
"""

import sys
import logging

from itertools import groupby
from optparse import OptionParser

from jcvi.formats.bed import Bed
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('prepare', 'given om alignment, prepare the patchers'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def merge_ranges(beds):

    from jcvi.utils.range import range_parse
    m = [x.accn for x in beds]

    mr = [range_parse(x) for x in m]
    mc = set(x.seqid for x in mr)
    if len(mc) != 1:
        logging.error("Multiple seqid found in pocket. Aborted.")
        return

    mc = list(mc)[0]
    ms = min(x.start for x in mr)
    me = max(x.end for x in mr)

    return mc, ms, me


def prepare(args):
    """
    %prog prepare om_alignment.bed

    Given optical map alignment, prepare the patchers. Use --backbone to suggest
    which assembly is the major one, and the patchers will be extracted from
    another assembly.
    """
    p = OptionParser(prepare.__doc__)
    p.add_option("--backbone", default="Scaffold",
                 help="Prefix of the backbone assembly [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bedfile, = args
    bb = opts.backbone

    bed = Bed(bedfile)
    uniqbed = Bed()
    key = lambda x: (x.seqid, x.start, x.end)
    for k, sb in groupby(bed, key=key):
        sb = list(sb)
        backbone = [x for x in sb if x.accn.startswith(bb)]
        others = [x for x in sb if not x.accn.startswith(bb)]
        if backbone and others:
            uniqbed.extend(backbone)
        else:
            uniqbed.extend(sb)
            if others and not backbone:
                print merge_ranges(others)

    uniqbedfile = bedfile.rsplit(".", 1)[0] + ".uniq.bed"
    uniqbed.print_to_file(uniqbedfile)


if __name__ == '__main__':
    main()
