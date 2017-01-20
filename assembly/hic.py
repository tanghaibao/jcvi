#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Process Hi-C output into AGP for chromosomal-scale scaffolding.
"""

import sys

from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('agp', 'generate AGP file based on LACHESIS output'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def agp(args):
    """
    %prog agp main_results/ contigs.fasta

    Generate AGP file based on LACHESIS output.
    """
    p = OptionParser(agp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    odir, contigsfasta = args


if __name__ == '__main__':
    main()
