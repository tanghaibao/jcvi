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

from collections import namedtuple
from optparse import OptionParser
from subprocess import Popen

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
        ('prepare', 'prepare pairs of sequences'),
        ('calc', 'calculate Ks between pairs of sequences'),
        ('report', 'generate a distribution of Ks values'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def prepare(args):
    """
    %prog prepare pairsfile cdsfile > paired.cds.fasta

    Pick sequences from cdsfile to form pairs, ready to be calculated. The
    pairsfile can be generated from formats.blast.cscore(). The first two
    columns contain the pair.
    """
    from jcvi.formats.fasta import Fasta, SeqIO

    p = OptionParser(prepare.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pairsfile, cdsfile = args

    f = Fasta(cdsfile)
    fp = open(pairsfile)
    fw = sys.stdout
    for row in fp:
        a, b = row.split()[:2]
        arec = f[a]
        brec = f[b]
        SeqIO.write((arec, brec), fw, "fasta")


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
    reader.next() # header
    data = []
    for row in reader:
        for i, a in enumerate(row):
            if i==0: continue
            row[i] = float(row[i])
        data.append(KsLine._make(row))
    return data


def report(args):
    '''
    %prog report ksfile

    generate a report given a Ks result file (as produced by synonymous_calc.py).
    describe the median Ks, Ka values, as well as the distribution in stem-leaf plot
    '''
    import numpy as np

    from jcvi.graphics.histogram import stem_leaf_plot

    p = OptionParser(report.__doc__)
    opts, args = p.parse_args(args)

    if len(args) !=  1:
        sys.exit(not p.print_help())

    ks_file, = args

    data = read_ks_file(ks_file)
    logging.debug('File `{0}` contains a total of {1} gene pairs'.\
            format(ks_file, len(data)))

    for f in fields.split()[1:]:
        columndata = [getattr(x, f) for x in data]
        title = "{0}: {1:.2f}".format(descriptions[f], np.median(columndata))
        title += " ({0:.2f} +/- {1:.2f})".\
                format(np.mean(columndata), np.std(columndata))
        ks = ("ks" in f)
        if not ks:
            continue

        bins = (0, 2., 20) if ks else (0, .6, 10)
        digit = 1 if ks else 2
        stem_leaf_plot(columndata, *bins, digit=digit, title=title)


if __name__ == '__main__':
    main()
