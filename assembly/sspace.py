#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SSPACE scaffolding-related operations.
"""

import sys

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():

    actions = (
        ('anchor', 'anchor contigs to upgrade existing structure'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def anchor(args):
    """
    %prog anchor evidencefile scaffolds.fasta contigs.fasta

    Use SSPACE evidencefile to scaffold contigs into existing scaffold
    structure, as in `scaffolds.fasta`. Contigs.fasta were used by SSPACE
    directly to scaffold.
    """
    p = OptionParser(anchor.__doc__)
    p.add_option("--mingap", default=10, type="int",
                 help="Option -minGap used with gapSplit [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    evidencefile, scaffolds, contigs = args


if __name__ == '__main__':
    main()
