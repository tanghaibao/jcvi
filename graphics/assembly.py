"""
Report assembly statistics, using idea from:
    http://blog.malde.org/index.php/a50/

which plots the contig number versus the cumulative base pairs. useful for
assembly QC.
"""

import sys
import logging

from bisect import bisect
import numpy as np

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import ActionDispatcher, debug
debug()


def main():
    actions = (
        ('A50', 'compare A50 graphics for a set of FASTA files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_a50(fastafile):

    f = Fasta(fastafile)
    ctg_sizes = np.array([length for k, length in f.itersizes()])
    ctg_sizes = np.sort(ctg_sizes)[::-1]
    logging.debug(str(ctg_sizes))

    a50 = np.cumsum(ctg_sizes)
    logging.debug(str(a50))
    total = np.sum(ctg_sizes)
    idx = bisect(a50, total / 2)
    n50 = ctg_sizes[idx]

    return a50, n50


def A50 (args):
    """
    %prog A50 contigs_A.fasta contigs_B.fasta ...
    
    Plots the A50 graphics, see blog post (http://blog.malde.org/index.php/a50/)
    """
    from optparse import OptionParser
    
    p = OptionParser(A50.__doc__)
    opts, args = p.parse_args(args)

    if not args:
        sys.exit(p.print_help())

    a50, n50 = get_a50(args[0])
    logging.debug("`%s` N50: %d" % (args[0], n50))

    # TODO: use ggplot2
    import matplotlib.pyplot as plt

    plt.plot(np.arange(len(a50)), a50, "b-")
    plt.savefig("t.png")


if __name__ == '__main__':
    main()
