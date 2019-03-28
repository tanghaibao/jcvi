#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Functions related to processing of the pistachio genome.
"""
from __future__ import print_function

import sys


from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('agp', 'convert from the table file to agp format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def agp(args):
    """
    %prog agp Siirt_Female_pistachio_23May2017_table.txt

    The table file, as prepared by Dovetail Genomics, is not immediately useful
    to convert gene model coordinates, as assumed by formats.chain.fromagp().
    This is a quick script to do such conversion. The file structure of this
    table file is described in the .manifest file shipped in the same package::

    pistachio_b_23May2017_MeyIy.table.txt
        Tab-delimited table describing positions of input assembly scaffolds
        in the Hirise scaffolds. The table has the following format:

            1. HiRise scaffold name
            2. Input sequence name
            3. Starting base (zero-based) of the input sequence
            4. Ending base of the input sequence
            5. Strand (- or +) of the input sequence in the scaffold
            6. Starting base (zero-based) in the HiRise scaffold
            7. Ending base in the HiRise scaffold

        where '-' in the strand column indicates that the sequence is reverse
        complemented relative to the input assembly.

    CAUTION: This is NOT a proper AGP format since it does not have gaps in
    them.
    """
    p = OptionParser(agp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    tablefile, = args
    fp = open(tablefile)
    for row in fp:
        atoms = row.split()
        hr = atoms[0]
        scaf = atoms[1]
        scaf_start = int(atoms[2]) + 1
        scaf_end = int(atoms[3])
        strand = atoms[4]
        hr_start = int(atoms[5]) + 1
        hr_end = int(atoms[6])

        print("\t".join(str(x) for x in \
                (hr, hr_start, hr_end, 1, 'W',
                 scaf, scaf_start, scaf_end, strand)))


if __name__ == '__main__':
    main()
