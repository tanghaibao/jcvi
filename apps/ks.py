#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Calculation of synonymous substitutions (Ks).
"""

import sys
import os
import os.path as op
import csv
import logging

from math import log, sqrt, pi, exp
from itertools import product, combinations
from collections import namedtuple
from optparse import OptionParser
from subprocess import Popen

import numpy as np
from Bio import SeqIO
from Bio import AlignIO
from Bio.Align.Applications import ClustalwCommandline

from jcvi.formats.base import must_open
from jcvi.apps.command import getpath, partial
from jcvi.apps.base import ActionDispatcher, debug, mkdir, set_outfile, sh
debug()


CLUSTALW_BIN = partial(getpath, name="CLUSTALW2")
PAL2NAL_BIN = partial(getpath, name="PAL2NAL")
PAML_BIN = partial(getpath, name="PAML")


class AbstractCommandline:

    def run(self):
        r = Popen(str(self), shell=True)
        return r.communicate()


class YnCommandline(AbstractCommandline):
    """Little commandline for yn00.
    """
    def __init__(self, ctl_file, command=PAML_BIN("yn00")):
        self.ctl_file = ctl_file
        self.parameters = []
        self.command = command

    def __str__(self):
        return self.command + " %s >/dev/null" % self.ctl_file


class MrTransCommandline(AbstractCommandline):
    """Simple commandline faker.
    """
    def __init__(self, prot_align_file, nuc_file, output_file,
            command=PAL2NAL_BIN("pal2nal.pl")):
        self.prot_align_file = prot_align_file
        self.nuc_file = nuc_file
        self.output_file = output_file
        self.command = command

        self.parameters = []

    def __str__(self):
        return self.command + " %s %s -output paml> %s" % (self.prot_align_file, self.nuc_file, self.output_file)


def main():

    actions = (
        ('fromgroups', 'flatten the gene families into pairs'),
        ('prepare', 'prepare pairs of sequences'),
        ('calc', 'calculate Ks between pairs of sequences'),
        ('report', 'generate a distribution of Ks values'),
        ('gc3', 'filter the Ks results to remove high GC3 genes'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def get_GC3(cdsfile):
    from jcvi.formats.fasta import Fasta

    f = Fasta(cdsfile, lazy=True)
    GC3 = {}
    for name, rec in f.iteritems_ordered():
        positions = rec.seq[2::3].upper()
        gc_counts = sum(1 for x in positions if x in "GC")
        gc_ratio = gc_counts * 1. / len(positions)
        GC3[name] = gc_ratio

    return GC3


def gc3(args):
    """
    %prog gc3 ksfile cdsfile > newksfile

    Filter the Ks results to remove high GC3 genes. High GC3 genes are
    problematic in Ks calculation - see Tang et al. 2010 PNAS. Specifically, the
    two calculation methods produce drastically different results for these
    pairs. Therefore we advise to remoeve these high GC3 genes. This is often
    the case for studying cereal genes.
    """
    import csv

    p = OptionParser(gc3.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ks_file, cdsfile = args
    header, data = read_ks_file(ks_file)
    noriginals = len(data)
    GC3 = get_GC3(cdsfile)

    writer = csv.writer(sys.stdout)
    writer.writerow(header)
    nlines = 0
    cutoff = .75
    for d in data:
        a, b = d.pair.split(";")
        aratio, bratio = GC3[a], GC3[b]
        if (aratio + bratio) / 2 > cutoff:
            continue
        writer.writerow(d)
        nlines += 1
    logging.debug("{0} records written (from {1}).".format(nlines, noriginals))


def extract_pairs(abed, bbed, groups):
    """
    Called by fromgroups(), extract pairs specific to a pair of species.
    """
    agenome = op.basename(abed.filename).split(".")[0]
    bgenome = op.basename(bbed.filename).split(".")[0]
    aorder = abed.order
    border = bbed.order
    pairsfile = "{0}.{1}.pairs".format(agenome, bgenome)
    fw = open(pairsfile, "w")

    is_self = abed.filename == bbed.filename
    npairs = 0
    for group in groups:
        iter = combinations(group, 2) if is_self \
                    else product(group, repeat=2)

        for a, b in iter:
            if a not in aorder or b not in border:
                continue

            print >> fw, "\t".join((a, b))
            npairs += 1

    logging.debug("File `{0}` written with {1} pairs.".format(pairsfile, npairs))


def fromgroups(args):
    """
    %prog fromgroups groupsfile a.bed b.bed ...

    Flatten the gene familes into pairs, the groupsfile is a file with each line
    containing the members, separated by comma. The commands also require
    several bed files in order to sort the pairs into different piles (e.g.
    pairs of species in comparison.
    """
    from jcvi.formats.bed import Bed

    p = OptionParser(fromgroups.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    groupsfile = args[0]
    bedfiles = args[1:]
    beds = [Bed(x) for x in bedfiles]
    fp = open(groupsfile)
    groups = [row.strip().split(",") for row in fp]
    for b1, b2 in product(beds, repeat=2):
        extract_pairs(b1, b2, groups)


def find_first_isoform(a, f):
    if a in f:
        return a
    for i in xrange(100):
        ia = ".".join((a, str(i)))
        if ia in f:
            return ia
    return a


def prepare(args):
    """
    %prog prepare pairsfile cdsfile > paired.cds.fasta

    Pick sequences from cdsfile to form pairs, ready to be calculated. The
    pairsfile can be generated from formats.blast.cscore(). The first two
    columns contain the pair.
    """
    from jcvi.formats.fasta import Fasta, SeqIO

    p = OptionParser(prepare.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pairsfile, cdsfile = args

    f = Fasta(cdsfile)
    fp = open(pairsfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        if row[0] == '#':
            continue
        a, b = row.split()[:2]
        if a not in f:
            a = find_first_isoform(a, f)
            assert a, a
        if b not in f:
            b = find_first_isoform(b, f)
            assert b, b

        arec = f[a]
        brec = f[b]
        SeqIO.write((arec, brec), fw, "fasta")
    fw.close()


def calc(args):
    """
    %prog calc [prot.fasta] cds.fasta > out.ks

    Protein file is optional. If only one file is given, it is assumed to
    be CDS sequences with correct frame (frame 0). Results will be written to
    stdout. Both protein file and nucleotide file are assumed to be Fasta format,
    with adjacent records as the pairs to compare.

    Author: Haibao Tang <bao@uga.edu>, Brad Chapman
    Calculate synonymous mutation rates for gene pairs

    This does the following:
        1. Fetches a protein pair.
        2. Aligns the protein pair with clustalw
        3. Convert the output to Fasta format.
        4. Use this alignment info to align gene sequences using PAL2NAL
        5. Run PAML yn00 to calculate synonymous mutation rates.
    """
    p = OptionParser(calc.__doc__)
    set_outfile(p)

    opts, args = p.parse_args(args)

    if len(args) == 1:
        protein_file, dna_file = None, args[0]
    elif len(args) == 2:
        protein_file, dna_file = args
    else:
        print >>sys.stderr, "Incorrect arguments"
        sys.exit(not p.print_help())

    output_h = must_open(opts.outfile, "w")
    output_h.write("name,dS-yn,dN-yn,dS-ng,dN-ng\n")
    work_dir = op.join(os.getcwd(), "syn_analysis")
    mkdir(work_dir)

    if not protein_file:
        protein_file = translate_dna(dna_file)

    prot_iterator = SeqIO.parse(open(protein_file), "fasta")
    dna_iterator = SeqIO.parse(open(dna_file), "fasta")
    for p_rec_1, p_rec_2, n_rec_1, n_rec_2 in \
            zip(prot_iterator, prot_iterator, dna_iterator, dna_iterator):

        print >>sys.stderr, "--------", p_rec_1.name, p_rec_2.name
        align_fasta = clustal_align_protein(p_rec_1, p_rec_2, work_dir)
        mrtrans_fasta = run_mrtrans(align_fasta, n_rec_1, n_rec_2, work_dir)
        if mrtrans_fasta:
            ds_subs_yn, dn_subs_yn, ds_subs_ng, dn_subs_ng = \
                    find_synonymous(mrtrans_fasta, work_dir)
            if ds_subs_yn is not None:
                pair_name = "%s;%s" % (p_rec_1.name, p_rec_2.name)
                output_h.write("%s\n" % (",".join(str(x) for x in (pair_name,
                        ds_subs_yn, dn_subs_yn, ds_subs_ng, dn_subs_ng))))
                output_h.flush()

    # Clean-up
    sh("rm -rf 2YN.t 2YN.dN 2YN.dS rst rub rst1 syn_analysis")


def translate_dna(dna_file):
    '''
    Translate the seqs in dna_file and produce a protein_file
    '''
    protein_file = dna_file + ".pep"
    translated = []
    for rec in SeqIO.parse(open(dna_file), "fasta"):
        rec.seq = rec.seq.translate()
        translated.append(rec)

    protein_h = open(protein_file, "w")
    SeqIO.write(translated, protein_h, "fasta")

    print >>sys.stderr, "%d records written to %s" % (len(translated),
            protein_file)
    return protein_file


def find_synonymous(input_file, work_dir):
    """Run yn00 to find the synonymous subsitution rate for the alignment.
    """
    # create the .ctl file
    ctl_file = op.join(work_dir, "yn-input.ctl")
    output_file = op.join(work_dir, "nuc-subs.yn")
    ctl_h = open(ctl_file, "w")
    ctl_h.write("seqfile = %s\noutfile = %s\nverbose = 0\n" %
                (input_file, output_file))
    ctl_h.write("icode = 0\nweighting = 0\ncommonf3x4 = 0\n")
    ctl_h.close()

    cl = YnCommandline(ctl_file)
    print >>sys.stderr, "\tyn00:", cl
    r, e = cl.run()
    ds_value_yn = None
    ds_value_ng = None
    dn_value_yn = None
    dn_value_ng = None

    # Nei-Gojobori
    output_h = open(output_file)
    row = output_h.readline()
    while row:
        if row.find("Nei & Gojobori") >=0:
            for x in xrange(5):
                row = output_h.next()
            dn_value_ng, ds_value_ng = row.split('(')[1].split(')')[0].split()
            break
        row = output_h.readline()
    output_h.close()

    # Yang
    output_h = open(output_file)
    for line in output_h:
        if line.find("+-") >= 0 and line.find("dS") == -1:
            parts = line.split(" +-")
            ds_value_yn = extract_subs_value(parts[1])
            dn_value_yn = extract_subs_value(parts[0])

    if ds_value_yn is None or ds_value_ng is None:
        h = open(output_file)
        print >>sys.stderr, "yn00 didn't work: \n%s" % h.read()

    return ds_value_yn, dn_value_yn, ds_value_ng, dn_value_ng


def extract_subs_value(text):
    """Extract a subsitution value from a line of text.

    This is just a friendly function to grab a float value for Ks and Kn
    values from the junk I get from the last line of the yn00 file.

    Line:
    2    1    52.7   193.3   2.0452  0.8979  0.0193 0.0573 +- 0.0177
    2.9732 +- 3.2002

    Parts:
        ['   2    1    52.7   193.3   2.0452  0.8979  0.0193 0.0573',
         ' 0.0177  2.9732', ' 3.2002\n']

    So we want 0.0573 for Kn and 2.9732 for Ks.
    """
    parts = text.split()
    value = float(parts[-1])

    return value


def run_mrtrans(align_fasta, rec_1, rec_2, work_dir):
    """Align two nucleotide sequences with mrtrans and the protein alignment.
    """
    align_file = op.join(work_dir, "prot-align.fasta")
    nuc_file = op.join(work_dir, "nuc.fasta")
    output_file = op.join(work_dir, "nuc-align.mrtrans")

    # make the protein alignment file
    align_h = open(align_file, "w")
    align_h.write(str(align_fasta))
    align_h.close()
    # make the nucleotide file
    SeqIO.write((rec_1, rec_2), file(nuc_file, "w"), "fasta")

    # run the program
    cl = MrTransCommandline(align_file, nuc_file, output_file)
    r, e = cl.run()
    if e is None:
        print >>sys.stderr, "\tpal2nal:", cl
        return output_file
    elif e.read().find("could not translate") >= 0:
        print >>sys.stderr, "***pal2nal could not translate"
        return None


def clustal_align_protein(rec_1, rec_2, work_dir):
    """Align the two given proteins with clustalw.
    """
    fasta_file = op.join(work_dir, "prot-start.fasta")
    align_file = op.join(work_dir, "prot.aln")
    SeqIO.write((rec_1, rec_2), file(fasta_file, "w"), "fasta")

    clustal_cl = ClustalwCommandline(CLUSTALW_BIN("clustalw2"),
            infile=fasta_file, outfile=align_file, outorder="INPUT",
            type="PROTEIN")
    stdout, stderr = clustal_cl()

    aln_file = file(clustal_cl.outfile)
    alignment = AlignIO.read(aln_file, "clustal")
    print >>sys.stderr, "\tDoing clustalw alignment: %s" % clustal_cl
    return alignment.format("fasta")


fields = "pair yn_ks yn_ka ng_ks ng_ka"
descriptions = {
        'pair': 'Gene pair',
        'yn_ks': 'Yang-Nielson method of Ks estimate',
        'yn_ka': 'Yang-Nielson method of Ka estimate',
        'ng_ks': 'Nei-Gojobori method of Ks estimate',
        'ng_ka': 'Nei-Gojobori method of Ka estimate'}

KsLine = namedtuple("KsLine", fields)


def read_ks_file(ks_file):
    reader = csv.reader(open(ks_file, "rb"))
    header = reader.next() # header
    data = []
    for row in reader:
        for i, a in enumerate(row):
            if i==0: continue
            row[i] = float(row[i])
        data.append(KsLine._make(row))

    logging.debug('File `{0}` contains a total of {1} gene pairs'.\
            format(ks_file, len(data)))

    return header, data


def my_hist(ax, l, interval, max_r, color='g'):
    if not l:
        return

    from pylab import poly_between

    n, p = [], []
    total_len = len(l)
    for i in np.arange(0, max_r, interval):
        xmin, xmax = i - .5 * interval, i + .5 * interval
        nx = [x for x in l if xmin <= x < xmax]
        n.append(i)
        p.append(len(nx) * 100. / total_len)

    xs, ys = poly_between(n, 0, p)
    return ax.fill(xs, ys, fc=color, alpha=.3)


def lognormpdf(bins, mu, sigma):
    return np.exp(-(np.log(bins) - mu) ** 2 / (2 * sigma ** 2)) / \
            (bins * sigma * sqrt(2 * pi))


def lognormpdf_mix(bins, probs, mus, sigmas, interval=.1):
    y = 0
    for prob, mu, sigma in zip(probs, mus, sigmas):
        y += prob * lognormpdf(bins, mu, sigma)
    y *= 100 * interval
    return y


def get_mixture(data, components):
    """
    probs = [.476, .509]
    mus = [.69069, -.15038]
    variances = [.468982e-1, .959052e-1]
    """
    from jcvi.apps.base import popen

    probs, mus, sigmas = [], [], []
    fw = must_open("tmp", "w")
    log_data = [log(x) for x in data if x > .05]
    data = "\n".join(["%.4f" % x for x in log_data]).replace("inf\n", "")
    fw.write(data)
    fw.close()

    cmd = "gmm-bic {0} {1} {2}".format(components, len(log_data), fw.name)
    pipe = popen(cmd)

    for row in pipe:
        if row[0] != '#':
            continue

        atoms = row.split(",")
        a, b, c = atoms[1:4]
        a = float(a)
        b = float(b)
        c = float(c)

        mus.append(a)
        sigmas.append(b)
        probs.append(c)

    os.remove(fw.name)
    return probs, mus, sigmas


def plot_ks_dist(ax, data, interval, components, ks_max, color='r'):
    from jcvi.graphics.base import _

    line, = my_hist(ax, data, interval, ks_max, color=color)
    logging.debug("Total {0} pairs after filtering.".format(len(data)))

    probs, mus, variances = get_mixture(data, components)

    bins = np.arange(0.001, ks_max, .001)
    y = lognormpdf_mix(bins, probs, mus, variances, interval)

    line_mixture, = ax.plot(bins, y, ':', color=color, lw=3)
    for i in xrange(components):
        peak_val = exp(mus[i])
        mixline = lognormpdf_mix(peak_val, probs, mus, variances, interval)
        ax.text(peak_val, mixline, _("Ks=%.2f" % peak_val), \
                color="w", size=10, bbox=dict(ec='w',fc=color, alpha=.6, boxstyle='round'))

    return line, line_mixture


def report(args):
    '''
    %prog report ksfile

    generate a report given a Ks result file (as produced by synonymous_calc.py).
    describe the median Ks, Ka values, as well as the distribution in stem-leaf plot
    '''
    from jcvi.graphics.histogram import stem_leaf_plot

    p = OptionParser(report.__doc__)
    p.add_option("--vmax", default=2., type="float",
            help="Maximum value, inclusive [default: %default]")
    p.add_option("--bins", default=20, type="int",
            help="Number of bins to plot in the histogram [default: %default]")
    p.add_option("--pdf", default=False, action="store_true",
            help="Generate graphic output for the histogram [default: %default]")
    p.add_option("--components", default=1, type="int",
            help="Number of components to decompose peaks [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) !=  1:
        sys.exit(not p.print_help())

    ks_file, = args
    header, data = read_ks_file(ks_file)
    ks_max = opts.vmax

    for f in fields.split()[1:]:
        columndata = [getattr(x, f) for x in data]
        #columndata = [x for x in columndata if 0 <= x <= ks_max]
        title = "{0}: {1:.2f}".format(descriptions[f], np.median(columndata))
        title += " ({0:.2f} +/- {1:.2f})".\
                format(np.mean(columndata), np.std(columndata))
        ks = ("ks" in f)
        if not ks:
            continue

        bins = (0, ks_max, opts.bins) if ks else (0, .6, 10)
        digit = 1 if ks else 2
        stem_leaf_plot(columndata, *bins, digit=digit, title=title)

    if not opts.pdf:
        return

    from jcvi.graphics.base import mpl, _, tex_formatter, tex_1digit_formatter
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

    fig = mpl.figure.Figure(figsize=(5, 5))

    canvas = FigureCanvas(fig)
    ax = fig.add_axes([.12, .1, .8, .8])
    components = opts.components
    data = [x.ng_ks for x in data]

    interval = ks_max / opts.bins
    line, line_mixture = plot_ks_dist(ax, data, interval, components, ks_max, color='r')
    leg = ax.legend((line, line_mixture), ("Ks", "Ks (fitted)"),
                    shadow=True, fancybox=True, prop={"size": 10})
    leg.get_frame().set_alpha(.5)

    ax.set_xlim((0, ks_max))
    ax.set_title(_('Ks distribution'), fontweight="bold")
    ax.set_xlabel(_('Synonymous substitutions per site (Ks)'))
    ax.set_ylabel(_('Percentage of gene pairs'))

    ax.xaxis.set_major_formatter(tex_1digit_formatter)
    ax.yaxis.set_major_formatter(tex_formatter)

    image_name = "Ks_plot.pdf"
    canvas.print_figure(image_name, dpi=300)
    logging.debug("Print image to `{0}`.".format(image_name))


if __name__ == '__main__':
    main()
