"""
Report assembly statistics, using idea from:
    http://blog.malde.org/index.php/a50/

which plots the contig number versus the cumulative base pairs. useful for
assembly QC.
"""

import sys
import logging
import os.path as op

from jcvi.formats.fasta import Fasta
from jcvi.assembly.base import calculate_A50
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
xlab("Contigs") + ylab("Cumulative size (Mb)") +
opts(title="A50 plot", legend.position="top")

ggsave(file="$rpdf")
"""


def main():
    actions = (
        ('A50', 'compare A50 graphics for a set of FASTA files'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def generate_plot(filename, rplot=rplot, rpdf=rpdf):

    rtemplate = RTemplate(rplot_template, locals())
    rtemplate.run()


def A50(args):
    """
    %prog A50 contigs_A.fasta contigs_B.fasta ...

    Plots A50 graphics, see blog post (http://blog.malde.org/index.php/a50/)
    """
    from optparse import OptionParser

    p = OptionParser(A50.__doc__)
    p.add_option("--overwrite", default=False, action="store_true",
            help="overwrite `%s` file if exists" % rplot)
    p.add_option("--cutoff", default=500, type="int", dest="cutoff",
            help="use contigs > certain size [default: %default]")
    p.add_option("--stepsize", default=10, type="int", dest="stepsize",
            help="stepsize for the distribution [default: %default]")
    opts, args = p.parse_args(args)

    if not args:
        sys.exit(p.print_help())

    stepsize = opts.stepsize  # use stepsize to speed up drawing
    if not op.exists(rplot) or opts.overwrite:
        fw = open(rplot, "w")
        header = "\t".join(("index", "cumsize", "fasta"))
        print >>fw, header
        for fastafile in args:
            f = Fasta(fastafile, index=False)
            ctgsizes = [length for k, length in f.itersizes()]

            a50, l50, n50, ctgsizes = calculate_A50(ctgsizes, cutoff=opts.cutoff)

            logging.debug("`{0}` ctg_sizes: {1}".format(fastafile, ctgsizes))

            tag = "{0} (N50={1})".format(op.basename(fastafile).rsplit(".", 1)[0], n50)
            logging.debug(tag)

            for i, s in zip(xrange(0, len(a50), stepsize), a50[::stepsize]):
                print >>fw, "\t".join((str(i), str(s / 1000000.), tag))
        fw.close()

    generate_plot(rplot)


if __name__ == '__main__':
    main()
