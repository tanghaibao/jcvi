#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Calculation of synonymous substitutions (Ks).
"""
import csv
import logging
import os
import os.path as op
import sys
from functools import partial
from itertools import combinations, product
from math import exp, log, pi, sqrt

import numpy as np
from Bio import AlignIO, SeqIO
from Bio.Align.Applications import ClustalwCommandline, MuscleCommandline

from jcvi.apps.base import (
    ActionDispatcher,
    OptionParser,
    Popen,
    getpath,
    iglob,
    mkdir,
    sh,
)
from jcvi.formats.base import LineFile, must_open
from jcvi.graphics.base import AbstractLayout, adjust_spines, markup, plt, savefig
from jcvi.utils.cbook import gene_name
from jcvi.utils.table import write_csv

CLUSTALW_BIN = partial(getpath, name="CLUSTALW2", warn="warn")
MUSCLE_BIN = partial(getpath, name="MUSCLE", warn="warn")
PAL2NAL_BIN = partial(getpath, name="PAL2NAL", warn="warn")
PAML_BIN = partial(getpath, name="PAML", warn="warn")


class AbstractCommandline:
    def run(self):
        r = Popen(str(self))
        return r.communicate()


class YnCommandline(AbstractCommandline):
    """Little commandline for yn00."""

    def __init__(self, ctl_file, command=PAML_BIN("yn00")):
        self.ctl_file = ctl_file
        self.parameters = []
        self.command = command

    def __str__(self):
        return self.command + " %s >/dev/null" % self.ctl_file


class MrTransCommandline(AbstractCommandline):
    """Simple commandline faker."""

    def __init__(
        self,
        prot_align_file,
        nuc_file,
        output_file,
        outfmt="paml",
        command=PAL2NAL_BIN("pal2nal.pl"),
    ):
        self.prot_align_file = prot_align_file
        self.nuc_file = nuc_file
        self.output_file = output_file
        self.outfmt = outfmt
        self.command = command

        self.parameters = []

    def __str__(self):
        return self.command + " %s %s -output %s > %s" % (
            self.prot_align_file,
            self.nuc_file,
            self.outfmt,
            self.output_file,
        )


def main():

    actions = (
        ("batch", "compute ks for a set of anchors file"),
        ("fromgroups", "flatten the gene families into pairs"),
        ("prepare", "prepare pairs of sequences"),
        ("calc", "calculate Ks between pairs of sequences"),
        ("subset", "subset pre-calculated Ks according to pairs file"),
        ("gc3", "filter the Ks results to remove high GC3 genes"),
        ("report", "generate a distribution of Ks values"),
        ("multireport", "generate several Ks value distributions in same figure"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def batch(args):
    """
    %prog batch all.cds *.anchors

    Compute Ks values for a set of anchors file. This will generate a bunch of
    work directories for each comparisons. The anchorsfile should be in the form
    of specie1.species2.anchors.
    """
    from jcvi.apps.grid import MakeManager

    p = OptionParser(batch.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    cdsfile = args[0]
    anchors = args[1:]
    workdirs = [".".join(op.basename(x).split(".")[:2]) for x in anchors]
    for wd in workdirs:
        mkdir(wd)

    mm = MakeManager()
    for wd, ac in zip(workdirs, anchors):
        pairscdsfile = wd + ".cds.fasta"
        cmd = "python -m jcvi.apps.ks prepare {} {} -o {}".format(
            ac, cdsfile, pairscdsfile
        )
        mm.add((ac, cdsfile), pairscdsfile, cmd)
        ksfile = wd + ".ks"
        cmd = "python -m jcvi.apps.ks calc {} -o {} --workdir {}".format(
            pairscdsfile, ksfile, wd
        )
        mm.add(pairscdsfile, ksfile, cmd)
    mm.write()


class LayoutLine(object):
    def __init__(self, row, delimiter=","):
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.ksfile = args[0]
        self.components = int(args[1])
        self.label = args[2]
        self.color = args[3]
        self.marker = args[4]

    def __str__(self):
        return ", ".join(
            str(x)
            for x in (self.ksfile, self.components, self.label, self.color, self.marker)
        )


class Layout(AbstractLayout):
    def __init__(self, filename, delimiter=","):
        super(Layout, self).__init__(filename)
        if not op.exists(filename):
            ksfiles = iglob(".", "*.ks")
            header = "Ks file|ncomponents|label|color|marker".split("|")
            contents = []
            for ksfile in ksfiles:
                leg = op.basename(ksfile).rsplit(".", 1)[0]
                if leg.count(".") == 1:
                    leg = leg.replace(".", " *vs.* ")
                contents.append((ksfile, "1", leg, "", ""))
            write_csv(header, contents, comment=True, filename=filename)

        fp = open(filename)
        for row in fp:
            if row[0] == "#":
                continue
            self.append(LayoutLine(row, delimiter=delimiter))

        self.assign_colors()
        self.assign_markers()


class KsPlot(object):
    def __init__(self, ax, ks_max, bins, legendp="upper left"):

        self.ax = ax
        self.ks_max = ks_max
        self.interval = ks_max / bins
        self.legendp = legendp
        self.lines = []
        self.labels = []

    def add_data(
        self,
        data,
        components=1,
        label="Ks",
        color="r",
        marker=".",
        fill=False,
        fitted=True,
        kde=False,
    ):

        ax = self.ax
        ks_max = self.ks_max
        interval = self.interval
        if kde:
            marker = None

        line, line_mixture = plot_ks_dist(
            ax,
            data,
            interval,
            components,
            ks_max,
            color=color,
            marker=marker,
            fill=fill,
            fitted=fitted,
            kde=kde,
        )
        self.lines.append(line)
        self.labels.append(label)

        if fitted:
            self.lines.append(line_mixture)
            self.labels.append(label + " (fitted)")

    def draw(self, title="*Ks* distribution", filename="Ks_plot.pdf"):

        ax = self.ax
        ks_max = self.ks_max
        lines = self.lines
        labels = [markup(x) for x in self.labels]
        legendp = self.legendp
        if len(lines) > 1:
            leg = ax.legend(
                lines,
                labels,
                loc=legendp,
                shadow=True,
                fancybox=True,
                prop={"size": 10},
            )
            leg.get_frame().set_alpha(0.5)

        ax.set_xlim((0, ks_max - self.interval))
        ylim = ax.get_ylim()[-1]
        ax.set_ylim(0, ylim)
        ax.set_title(markup(title), fontweight="bold")
        ax.set_xlabel(markup("Synonymous substitutions per site (*Ks*)"))
        ax.set_ylabel("Percentage of gene pairs (bin={})".format(self.interval))

        ax.set_xticklabels(ax.get_xticks(), family="Helvetica")
        ax.set_yticklabels(ax.get_yticks(), family="Helvetica")

        adjust_spines(ax, ["left", "bottom"], outward=True)

        if filename:
            savefig(filename, dpi=300)


def multireport(args):
    """
    %prog multireport layoutfile

    Generate several Ks value distributions in the same figure. If the layout
    file is missing then a template file listing all ks files will be written.

    The layout file contains the Ks file, number of components, colors, and labels:

    # Ks file, ncomponents, label, color, marker
    LAP.sorghum.ks, 1, LAP-sorghum, r, o
    SES.sorghum.ks, 1, SES-sorghum, g, +
    MOL.sorghum.ks, 1, MOL-sorghum, m, ^

    If color or marker is missing, then a random one will be assigned.
    """
    p = OptionParser(multireport.__doc__)
    p.set_outfile(outfile="Ks_plot.pdf")
    add_plot_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="8x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (layoutfile,) = args
    ks_min = opts.vmin
    ks_max = opts.vmax
    bins = opts.bins
    fill = opts.fill
    layout = Layout(layoutfile)
    print(layout, file=sys.stderr)

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax = fig.add_axes([0.12, 0.13, 0.8, 0.8])

    kp = KsPlot(ax, ks_max, bins, legendp=opts.legendp)
    for lo in layout:
        data = KsFile(lo.ksfile)
        data = [x.ng_ks for x in data]
        data = [x for x in data if ks_min <= x <= ks_max]
        kp.add_data(
            data,
            lo.components,
            label=lo.label,
            color=lo.color,
            marker=lo.marker,
            fill=fill,
            fitted=opts.fit,
            kde=opts.kde,
        )

    kp.draw(title=opts.title, filename=opts.outfile)


def get_GC3(cdsfile):
    from jcvi.formats.fasta import Fasta

    f = Fasta(cdsfile, lazy=True)
    GC3 = {}
    for name, rec in f.iteritems_ordered():
        positions = rec.seq[2::3].upper()
        gc_counts = sum(1 for x in positions if x in "GC")
        gc_ratio = gc_counts * 1.0 / len(positions)
        GC3[name] = gc_ratio

    return GC3


def plot_GC3(GC3, cdsfile, fill="white"):
    from jcvi.graphics.histogram import histogram

    numberfile = "{0}.gc3".format(cdsfile)
    fw = must_open(numberfile, "w")
    fw.write("\n".join(map(str, GC3.values())))
    fw.close()
    histogram(
        numberfile,
        vmin=0,
        vmax=1,
        xlabel="GC3",
        title=cdsfile,
        bins=50,
        skip=0,
        ascii=False,
        fill=fill,
    )

    logging.debug("{0} GC3 values plotted to {1}.pdf".format(len(GC3), numberfile))


def gc3(args):
    """
    %prog gc3 ksfile cdsfile [cdsfile2] -o newksfile

    Filter the Ks results to remove high GC3 genes. High GC3 genes are
    problematic in Ks calculation - see Tang et al. 2010 PNAS. Specifically, the
    two calculation methods produce drastically different results for these
    pairs. Therefore we advise to remoeve these high GC3 genes. This is often
    the case for studying cereal genes.

    If 2 genomes are involved, the cdsfile of the 2nd genome can be provided
    concatenated or separated.
    """
    p = OptionParser(gc3.__doc__)
    p.add_option(
        "--plot", default=False, action="store_true", help="Also plot the GC3 histogram"
    )
    p.set_outfile()

    opts, args = p.parse_args(args)

    outfile = opts.outfile
    plot = opts.plot

    if not 1 < len(args) < 4:
        sys.exit(not p.print_help())

    ks_file, cdsfile = args[:2]
    GC3 = get_GC3(cdsfile)
    if plot:
        plot_GC3(GC3, cdsfile, fill="green")

    if len(args) == 3:
        cdsfile2 = args[2]
        GC3_2 = get_GC3(cdsfile2)
        GC3.update(GC3_2)
        if plot:
            plot_GC3(GC3_2, cdsfile2, fill="lightgreen")

    data = KsFile(ks_file)
    noriginals = len(data)

    fw = must_open(outfile, "w")
    writer = csv.writer(fw)
    writer.writerow(fields.split(","))
    nlines = 0
    cutoff = 0.75
    for d in data:
        a, b = d.name.split(";")
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
        iter = combinations(group, 2) if is_self else product(group, repeat=2)

        for a, b in iter:
            if a not in aorder or b not in border:
                continue

            print("\t".join((a, b)), file=fw)
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
    for i in range(100):
        ia = ".".join((a, str(i)))
        if ia in f:
            return ia
    return a


def prepare(args):
    """
    %prog prepare pairsfile cdsfile [pepfile] -o paired.cds.fasta

    Pick sequences from cdsfile to form pairs, ready to be calculated. The
    pairsfile can be generated from formats.blast.cscore(). The first two
    columns contain the pair.
    """
    from jcvi.formats.fasta import Fasta

    p = OptionParser(prepare.__doc__)
    p.set_outfile()

    opts, args = p.parse_args(args)
    outfile = opts.outfile

    if len(args) == 2:
        pairsfile, cdsfile = args
        pepfile = None
    elif len(args) == 3:
        pairsfile, cdsfile, pepfile = args
    else:
        sys.exit(not p.print_help())

    f = Fasta(cdsfile)
    fp = open(pairsfile)
    fw = must_open(outfile, "w")
    if pepfile:
        assert outfile != "stdout", "Please specify outfile name."
        f2 = Fasta(pepfile)
        fw2 = must_open(outfile + ".pep", "w")
    for row in fp:
        if row[0] == "#":
            continue
        a, b = row.split()[:2]
        if a == b:
            logging.debug("Self pairs found: {0} - {1}. Ignored".format(a, b))
            continue

        if a not in f:
            a = find_first_isoform(a, f)
            assert a, a
        if b not in f:
            b = find_first_isoform(b, f)
            assert b, b

        acds = f[a]
        bcds = f[b]
        SeqIO.write((acds, bcds), fw, "fasta")
        if pepfile:
            apep = f2[a]
            bpep = f2[b]
            SeqIO.write((apep, bpep), fw2, "fasta")
    fw.close()
    if pepfile:
        fw2.close()


def calc(args):
    """
    %prog calc [prot.fasta] cds.fasta > out.ks

    Protein file is optional. If only one file is given, it is assumed to
    be CDS sequences with correct frame (frame 0). Results will be written to
    stdout. Both protein file and nucleotide file are assumed to be Fasta format,
    with adjacent records as the pairs to compare.

    Author: Haibao Tang <bao@uga.edu>, Brad Chapman, Jingping Li
    Calculate synonymous mutation rates for gene pairs

    This does the following:
        1. Fetches a protein pair.
        2. Aligns the protein pair with clustalw (default) or muscle.
        3. Convert the output to Fasta format.
        4. Use this alignment info to align gene sequences using PAL2NAL
        5. Run PAML yn00 to calculate synonymous mutation rates.
    """
    from jcvi.formats.fasta import translate

    p = OptionParser(calc.__doc__)
    p.add_option(
        "--longest",
        action="store_true",
        help="Get longest ORF, only works if no pep file, e.g. ESTs",
    )
    p.add_option(
        "--msa",
        default="clustalw",
        choices=("clustalw", "muscle"),
        help="software used to align the proteins",
    )
    p.add_option("--workdir", default=os.getcwd(), help="Work directory")
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) == 1:
        protein_file, dna_file = None, args[0]
    elif len(args) == 2:
        protein_file, dna_file = args
    else:
        print("Incorrect arguments", file=sys.stderr)
        sys.exit(not p.print_help())

    output_h = must_open(opts.outfile, "w")
    print(fields, file=output_h)
    work_dir = op.join(opts.workdir, "syn_analysis")
    mkdir(work_dir)

    if not protein_file:
        protein_file = dna_file + ".pep"
        translate_args = [dna_file, "--outfile=" + protein_file]
        if opts.longest:
            translate_args += ["--longest"]
        dna_file, protein_file = translate(translate_args)

    prot_iterator = SeqIO.parse(open(protein_file), "fasta")
    dna_iterator = SeqIO.parse(open(dna_file), "fasta")
    for p_rec_1, p_rec_2, n_rec_1, n_rec_2 in zip(
        prot_iterator, prot_iterator, dna_iterator, dna_iterator
    ):

        print("--------", p_rec_1.name, p_rec_2.name, file=sys.stderr)
        if opts.msa == "clustalw":
            align_fasta = clustal_align_protein((p_rec_1, p_rec_2), work_dir)
        elif opts.msa == "muscle":
            align_fasta = muscle_align_protein((p_rec_1, p_rec_2), work_dir)
        mrtrans_fasta = run_mrtrans(align_fasta, (n_rec_1, n_rec_2), work_dir)
        if mrtrans_fasta:
            ds_subs_yn, dn_subs_yn, ds_subs_ng, dn_subs_ng = find_synonymous(
                mrtrans_fasta, work_dir
            )
            if ds_subs_yn is not None:
                pair_name = "%s;%s" % (p_rec_1.name, p_rec_2.name)
                output_h.write(
                    "%s\n"
                    % (
                        ",".join(
                            str(x)
                            for x in (
                                pair_name,
                                ds_subs_yn,
                                dn_subs_yn,
                                ds_subs_ng,
                                dn_subs_ng,
                            )
                        )
                    )
                )
                output_h.flush()

    # Clean-up
    sh("rm -rf 2YN.t 2YN.dN 2YN.dS rst rub rst1 syn_analysis")


def find_synonymous(input_file, work_dir):
    """Run yn00 to find the synonymous subsitution rate for the alignment."""
    cwd = os.getcwd()
    os.chdir(work_dir)
    # create the .ctl file
    ctl_file = "yn-input.ctl"
    output_file = "nuc-subs.yn"
    ctl_h = open(ctl_file, "w")
    ctl_h.write(
        "seqfile = %s\noutfile = %s\nverbose = 0\n"
        % (op.basename(input_file), output_file)
    )
    ctl_h.write("icode = 0\nweighting = 0\ncommonf3x4 = 0\n")
    ctl_h.close()

    cl = YnCommandline(ctl_file)
    print("\tyn00:", cl, file=sys.stderr)
    r, e = cl.run()
    ds_value_yn = None
    ds_value_ng = None
    dn_value_yn = None
    dn_value_ng = None

    # Nei-Gojobori
    output_h = open(output_file)
    row = output_h.readline()
    while row:
        if row.find("Nei & Gojobori") >= 0:
            for x in range(5):
                row = next(output_h)
            dn_value_ng, ds_value_ng = row.split("(")[1].split(")")[0].split()
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
        print("yn00 didn't work: \n%s" % h.read(), file=sys.stderr)

    os.chdir(cwd)
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


def run_mrtrans(align_fasta, recs, work_dir, outfmt="paml"):
    """Align nucleotide sequences with mrtrans and the protein alignment."""
    align_file = op.join(work_dir, "prot-align.fasta")
    nuc_file = op.join(work_dir, "nuc.fasta")
    output_file = op.join(work_dir, "nuc-align.mrtrans")

    # make the prot_align file and nucleotide file
    align_h0 = open(align_file + "0", "w")
    align_h0.write(str(align_fasta))
    align_h0.close()
    prot_seqs = {}
    i = 0
    for rec in SeqIO.parse(align_h0.name, "fasta"):
        prot_seqs[i] = rec.seq
        i += 1
    align_h = open(align_file, "w")
    for i, rec in enumerate(recs):
        if len(rec.id) > 30:
            rec.id = rec.id[:28] + "_" + str(i)
            rec.description = ""
        print(">{0}\n{1}".format(rec.id, prot_seqs[i]), file=align_h)
    align_h.close()
    SeqIO.write(recs, open(nuc_file, "w"), "fasta")

    # run the program
    cl = MrTransCommandline(align_file, nuc_file, output_file, outfmt=outfmt)
    r, e = cl.run()
    if e is None:
        print("\tpal2nal:", cl, file=sys.stderr)
        return output_file
    elif e.read().find("could not translate") >= 0:
        print("***pal2nal could not translate", file=sys.stderr)
        return None


def clustal_align_protein(recs, work_dir, outfmt="fasta"):
    """
    Align given proteins with clustalw.
    recs are iterable of Biopython SeqIO objects
    """
    fasta_file = op.join(work_dir, "prot-start.fasta")
    align_file = op.join(work_dir, "prot.aln")
    SeqIO.write(recs, open(fasta_file, "w"), "fasta")

    clustal_cl = ClustalwCommandline(
        cmd=CLUSTALW_BIN("clustalw2"),
        infile=fasta_file,
        outfile=align_file,
        outorder="INPUT",
        type="PROTEIN",
    )
    stdout, stderr = clustal_cl()

    aln_file = open(clustal_cl.outfile)
    alignment = AlignIO.read(aln_file, "clustal")
    print("\tDoing clustalw alignment: %s" % clustal_cl, file=sys.stderr)
    if outfmt == "fasta":
        return alignment.format("fasta")
    if outfmt == "clustal":
        return alignment


def muscle_align_protein(recs, work_dir, outfmt="fasta", inputorder=True):
    """
    Align given proteins with muscle.
    recs are iterable of Biopython SeqIO objects
    """
    fasta_file = op.join(work_dir, "prot-start.fasta")
    align_file = op.join(work_dir, "prot.aln")
    SeqIO.write(recs, open(fasta_file, "w"), "fasta")

    muscle_cl = MuscleCommandline(
        cmd=MUSCLE_BIN("muscle"),
        input=fasta_file,
        out=align_file,
        seqtype="protein",
        clwstrict=True,
    )
    stdout, stderr = muscle_cl()
    alignment = AlignIO.read(muscle_cl.out, "clustal")

    if inputorder:
        try:
            muscle_inputorder(muscle_cl.input, muscle_cl.out)
        except ValueError:
            return ""
        alignment = AlignIO.read(muscle_cl.out, "fasta")

    print("\tDoing muscle alignment: %s" % muscle_cl, file=sys.stderr)
    if outfmt == "fasta":
        return alignment.format("fasta")
    if outfmt == "clustal":
        return alignment.format("clustal")


def muscle_inputorder(inputfastafile, alnfile, trunc_name=True):
    """
    Fix for muscle -stable option according to here:
    http://drive5.com/muscle/stable.html
    """
    sh("cp {0} {0}.old".format(alnfile), log=False)
    maxi = 30 if trunc_name else 1000

    aa = AlignIO.read(alnfile, "clustal")
    alignment = dict((a.id[:maxi], a) for a in aa)
    if trunc_name and len(alignment) < len(aa):
        raise ValueError("ERROR: The first 30 chars of your seq names are not unique")

    fw = must_open(alnfile, "w")
    for rec in SeqIO.parse(inputfastafile, "fasta"):
        a = alignment[rec.id[:maxi]]
        fw.write(">{0}\n{1}\n".format(a.id[:maxi], a.seq))

    fw.close()
    sh("rm {0}.old".format(alnfile), log=False)


def subset(args):
    """
    %prog subset pairsfile ksfile1 ksfile2 ... -o pairs.ks

    Subset some pre-calculated ks ka values (in ksfile) according to pairs
    in tab delimited pairsfile/anchorfile.
    """
    p = OptionParser(subset.__doc__)
    p.add_option(
        "--noheader", action="store_true", help="don't write ksfile header line"
    )
    p.add_option(
        "--block", action="store_true", help="preserve block structure in input"
    )
    p.set_stripnames()
    p.set_outfile()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pairsfile, ksfiles = args[0], args[1:]
    noheader = opts.noheader
    block = opts.block
    if block:
        noheader = True
    outfile = opts.outfile

    ksvals = {}
    for ksfile in ksfiles:
        ksvals.update(
            dict(
                (line.name, line)
                for line in KsFile(ksfile, strip_names=opts.strip_names)
            )
        )

    fp = open(pairsfile)
    fw = must_open(outfile, "w")

    if not noheader:
        print(fields, file=fw)

    i = j = 0
    for row in fp:
        if row[0] == "#":
            if block:
                print(row.strip(), file=fw)
            continue
        a, b = row.split()[:2]
        name = ";".join((a, b))
        if name not in ksvals:
            name = ";".join((b, a))
            if name not in ksvals:
                j += 1
                print("\t".join((a, b, ".", ".")), file=fw)
                continue
        ksline = ksvals[name]
        if block:
            print("\t".join(str(x) for x in (a, b, ksline.ks)), file=fw)
        else:
            ksline.name = ";".join((a, b))
            print(ksline, file=fw)
        i += 1
    fw.close()

    logging.debug("{0} pairs not found in ksfiles".format(j))
    logging.debug("{0} ks records written to `{1}`".format(i, outfile))
    return outfile


fields = "name,yn_ks,yn_ka,ng_ks,ng_ka"
descriptions = {
    "name": "Gene pair",
    "yn_ks": "Yang-Nielson Ks estimate",
    "yn_ka": "Yang-Nielson Ka estimate",
    "ng_ks": "Nei-Gojobori Ks estimate",
    "ng_ka": "Nei-Gojobori Ka estimate",
}


class KsLine:
    def __init__(self, row, strip_names=False):
        args = row.strip().split(",")
        self.name = args[0]
        self.yn_ks = self.get_float(args[1])
        self.yn_ka = self.get_float(args[2])
        self.ng_ks = self.get_float(args[3])
        self.ng_ka = self.get_float(args[4])
        self.ks = self.ng_ks
        if ";" in self.name:
            self.gene_a, self.gene_b = self.name.split(";")
            if strip_names:
                self.gene_a = gene_name(self.gene_a)
                self.gene_b = gene_name(self.gene_b)

    def get_float(self, x):
        try:
            x = float(x)
        except:
            x = -1
        return x

    def __str__(self):
        return ",".join(
            str(x) for x in (self.name, self.yn_ks, self.yn_ka, self.ng_ks, self.ng_ka)
        )

    @property
    def anchorline(self):
        return "\t".join(
            (gene_name(self.gene_a), gene_name(self.gene_b), "{:.3f}".format(self.ks))
        )


class KsFile(LineFile):
    def __init__(self, filename, strip_names=False):
        super(KsFile, self).__init__(filename)

        fp = open(filename)
        for row in fp:
            ksline = KsLine(row, strip_names=strip_names)
            if ksline.name == "name":  # header
                continue
            self.append(ksline)

        logging.debug(
            "File `{0}` contains a total of {1} gene pairs".format(filename, len(self))
        )

    def print_to_anchors(self, outfile):
        fw = must_open(outfile, "w")
        for row in self:
            print(row.anchorline, file=fw)
        fw.close()


def my_hist(ax, l, interval, max_r, color="g", marker=".", fill=False, kde=False):
    if not l:
        return

    n, p = [], []
    total_len = len(l)
    for i in np.arange(0, max_r, interval):
        xmin, xmax = i - 0.5 * interval, i + 0.5 * interval
        nx = [x for x in l if xmin <= x < xmax]
        n.append(i)
        p.append(len(nx) * 100.0 / total_len)

    if kde:
        from scipy import stats

        kernel = stats.gaussian_kde(l)
        n = np.arange(0, max_r, interval)
        kn = kernel(n)
        p = kn / sum(kn) * 100

    if fill:
        from pylab import poly_between

        xs, ys = poly_between(n, 0, p)
        line = ax.fill(xs, ys, fc=color, alpha=0.5)

    else:
        line = ax.plot(
            n, p, color=color, lw=2, ms=3, marker=marker, mfc="w", mec=color, mew=2
        )

    return line


def lognormpdf(bins, mu, sigma):
    return np.exp(-((np.log(bins) - mu) ** 2) / (2 * sigma ** 2)) / (
        bins * sigma * sqrt(2 * pi)
    )


def lognormpdf_mix(bins, probs, mus, sigmas, interval=0.1):
    y = 0
    for prob, mu, sigma in zip(probs, mus, sigmas):
        y += prob * lognormpdf(bins, mu, sigma)
    y *= 100 * interval  # Percentage

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
    log_data = [log(x) for x in data if x > 0.05]
    data = "\n".join(["%.4f" % x for x in log_data]).replace("inf\n", "")
    fw.write(data)
    fw.close()

    cmd = "gmm-bic {0} {1} {2}".format(components, len(log_data), fw.name)
    pipe = popen(cmd)

    for row in pipe:
        if row[0] != "#":
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


def plot_ks_dist(
    ax,
    data,
    interval,
    components,
    ks_max,
    color="r",
    marker=".",
    fill=False,
    fitted=True,
    kde=False,
):

    (line,) = my_hist(
        ax, data, interval, ks_max, color=color, marker=marker, fill=fill, kde=kde
    )
    logging.debug("Total {0} pairs after filtering.".format(len(data)))

    line_mixture = None
    if fitted:
        probs, mus, variances = get_mixture(data, components)

        iv = 0.001
        bins = np.arange(iv, ks_max, iv)
        y = lognormpdf_mix(bins, probs, mus, variances, interval)

        (line_mixture,) = ax.plot(bins, y, ":", color=color, lw=3)

        for i in range(components):
            peak_val = exp(mus[i])
            mixline = lognormpdf_mix(peak_val, probs, mus, variances, interval)
            ax.text(
                peak_val,
                mixline,
                "Ks=%.2f" % peak_val,
                color="w",
                size=10,
                bbox=dict(ec="w", fc=color, alpha=0.6, boxstyle="round"),
            )

    return line, line_mixture


def add_plot_options(p):
    p.add_option("--fit", default=False, action="store_true", help="Plot fitted lines")
    p.add_option("--kde", default=False, action="store_true", help="Use KDE smoothing")
    p.add_option("--vmin", default=0.0, type="float", help="Minimum value, inclusive")
    p.add_option("--vmax", default=3.0, type="float", help="Maximum value, inclusive")
    p.add_option(
        "--bins", default=60, type="int", help="Number of bins to plot in the histogram"
    )
    p.add_option("--legendp", default="upper right", help="Place of the legend")
    p.add_option(
        "--fill",
        default=False,
        action="store_true",
        help="Do not fill the histogram area",
    )
    p.add_option("--title", default="*Ks* distribution", help="Title of the plot")


def report(args):
    """
    %prog report ksfile

    generate a report given a Ks result file (as produced by synonymous_calc.py).
    describe the median Ks, Ka values, as well as the distribution in stem-leaf plot
    """
    from jcvi.utils.cbook import SummaryStats
    from jcvi.graphics.histogram import stem_leaf_plot

    p = OptionParser(report.__doc__)
    p.add_option(
        "--pdf",
        default=False,
        action="store_true",
        help="Generate graphic output for the histogram",
    )
    p.add_option(
        "--components",
        default=1,
        type="int",
        help="Number of components to decompose peaks",
    )
    add_plot_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="5x5")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (ks_file,) = args
    data = KsFile(ks_file)
    ks_min = opts.vmin
    ks_max = opts.vmax
    bins = opts.bins

    for f in fields.split(",")[1:]:
        columndata = [getattr(x, f) for x in data]
        ks = "ks" in f
        if not ks:
            continue

        columndata = [x for x in columndata if ks_min <= x <= ks_max]

        st = SummaryStats(columndata)
        title = "{0} ({1}): ".format(descriptions[f], ks_file)
        title += "Median:{0:.3f} (1Q:{1:.3f}|3Q:{2:.3f}||".format(
            st.median, st.firstq, st.thirdq
        )
        title += "Mean:{0:.3f}|Std:{1:.3f}||N:{2})".format(st.mean, st.sd, st.size)

        tbins = (0, ks_max, bins) if ks else (0, 0.6, 10)
        digit = 2 if (ks_max * 1.0 / bins) < 0.1 else 1
        stem_leaf_plot(columndata, *tbins, digit=digit, title=title)

    if not opts.pdf:
        return

    components = opts.components
    data = [x.ng_ks for x in data]
    data = [x for x in data if ks_min <= x <= ks_max]

    fig = plt.figure(1, (iopts.w, iopts.h))
    ax = fig.add_axes([0.12, 0.1, 0.8, 0.8])
    kp = KsPlot(ax, ks_max, opts.bins, legendp=opts.legendp)
    kp.add_data(data, components, fill=opts.fill, fitted=opts.fit, kde=opts.kde)
    kp.draw(title=opts.title)


if __name__ == "__main__":
    main()
