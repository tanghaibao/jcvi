"""
Report assembly statistics, using idea from:
    http://blog.malde.org/index.php/a50/

which plots the contig number versus the cumulative base pairs. useful for
assembly QC.
"""

import sys
import logging
import os.path as op

from bisect import bisect
import numpy as np

from jcvi.formats.fasta import Fasta
from jcvi.apps.base import ActionDispatcher, debug
debug()

from jcvi.apps.R import RTemplate

rplot = "A50.rplot"
rpdf = "A50.pdf"

rplot_template = """
library(ggplot2)

data <- read.table("$rplot", header=T, sep="\t")
g <- ggplot(data, aes(x=index, y=cumsize, group=fasta))
g + geom_line(aes(colour=fasta)) + 
xlab("Contigs") + ylab("Cumulative size") +
opts(title="A50 plot", legend.position="top")

ggsave(file="$rpdf")
"""


def main():
    actions = (
        ('A50', 'compare A50 graphics for a set of FASTA files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_a50(fastafile, cutoff=500):

    f = Fasta(fastafile, index=False)
    ctg_sizes = np.array([length for k, length in f.itersizes()])
    ctg_sizes = np.sort(ctg_sizes)[::-1]
    ctg_sizes = ctg_sizes[ctg_sizes >= cutoff]
    logging.debug("`%s` ctg_sizes: %s" % (fastafile, ctg_sizes))

    a50 = np.cumsum(ctg_sizes)

    total = np.sum(ctg_sizes)
    idx = bisect(a50, total / 2)
    n50 = ctg_sizes[idx]

    return a50, n50


def generate_plot(filename, rplot=rplot, rpdf=rpdf):

    rtemplate = RTemplate(rplot_template, locals())
    rtemplate.run()


def A50 (args):
    """
    %prog A50 contigs_A.fasta contigs_B.fasta ...
    
    Plots the A50 graphics, see blog post (http://blog.malde.org/index.php/a50/)
    """
    from optparse import OptionParser
    
    p = OptionParser(A50.__doc__)
    p.add_option("--overwrite", default=False, action="store_true",
            help="overwrite `%s` file if exists" % rplot)
    p.add_option("--cutoff", default=500, type="int",
            dest="cutoff",
            help="use only contigs larger than certain size [default: %default]")
    opts, args = p.parse_args(args)

    if not args:
        sys.exit(p.print_help())

    stepsize = 10 # use stepsize to speed up drawing
    if not op.exists(rplot) or opts.overwrite:
        fw = open(rplot, "w")
        header = "\t".join(("index", "cumsize", "fasta"))
        print >>fw, header
        for a in args:
            a50, n50 = get_a50(a, cutoff=opts.cutoff)
            logging.debug("`%s` N50: %d" % (a, n50))

            for i, s in zip(xrange(0, len(a50), stepsize), a50[::stepsize]):
                print >>fw, "\t".join((str(i), str(s), 
                    "%s (N50=%d)" % (op.basename(a).rsplit(".", 1)[0], n50)))
        fw.close()

    generate_plot(rplot)


if __name__ == '__main__':
    main()
