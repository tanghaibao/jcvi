#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utility to run Automated Human Readable Description (AHRD) pipeline.

<https://github.com/groupschoof/AHRD>
"""

import os.path as op
import sys
import re

from glob import glob
from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug, mkdir
debug()


Template = """
proteins_fasta: {2}
blast_dbs:
  swissprot:
    weight: 100
    file: swissprot/{1}.swissprot.pairwise
    blacklist: {0}/blacklist_descline.txt
    filter: {0}/filter_descline_sprot.txt
    token_blacklist: {0}/blacklist_token.txt
    description_score_bit_score_weight: 0.2

  tair:
    weight: 50
    file: tair/{1}.tair.pairwise
    blacklist: {0}/blacklist_descline.txt
    filter: {0}/filter_descline_tair.txt
    token_blacklist: {0}/blacklist_token.txt
    description_score_bit_score_weight: 0.4

  trembl:
    weight: 10
    file: trembl/{1}.trembl.pairwise
    blacklist: {0}/blacklist_descline.txt
    filter: {0}/filter_descline_trembl.txt
    token_blacklist: {0}/blacklist_token.txt
    description_score_bit_score_weight: 0.4

token_score_bit_score_weight: 0.5
token_score_database_score_weight: 0.3
token_score_overlap_score_weight: 0.2
description_score_relative_description_frequency_weight: 0.6
output: {3}
"""

def main():

    actions = (
        ('batch', 'batch run AHRD'),
        ('merge', 'merge AHRD run results'),
        ('fix', 'fix AHRD names'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


Unknown = "Unknown protein"
Hypothetical = "hypothetical protein"


def fix_text(s):

    s = s.split(";")[0]

    # Cellular locations
    loc_pat = re.compile(r",\s*(chloroplastic|cytoplasmic|mitochondrial).*?\s$", re.I)
    # Any word that matches e.g. AT5G54690
    atg_pat = re.compile(r"\bAT[1-5M]G\d{5}.*?\s", re.I)
    # Any word that matches e.g. Os02g0234800
    osg_pat = re.compile(r"\bOs\d{2}g\d{7}.*?\s", re.I)
    # (fragment)
    frag_pat = re.compile(r"\(fragment[s]?\)", re.I)
    # Trailing protein numeric copy (e.g. Myb 1)
    trail_pat = re.compile(r"(?<!type)\s\d+\s*$", re.I)
    # UPF
    upf_pat = re.compile(r"^UPF.*$")

    for pat in (loc_pat, atg_pat, osg_pat, frag_pat, trail_pat, upf_pat):
        # below is a hack since word boundaries don't work on /
        s = s.strip() + " "
        s = re.sub(pat, "", s)

    # &apos;? => '
    apos_pat = re.compile(r"&apos;?")
    s = re.sub(apos_pat, "'", s)
    # &gt => none
    gt_pat = re.compile(r"&gt")
    s = re.sub(gt_pat, "", s)
    # reduce runs such as -- '''
    s = re.sub(r"[-]+", "-", s)
    s = re.sub(r"[']+", "'", s)

    s = s.strip()

    # -like to -like protein
    like_pat = re.compile(r"[-]like$", re.I)
    s = re.sub(like_pat, "-like protein", s)

    sl = s.lower()

    # Any mention of `clone` or `contig` is not informative
    if "clone" in sl or "contig" in sl:
        s = Unknown

    # All that's left is `protein` is not informative
    if sl in ("protein", ""):
        s = Unknown

    if Unknown.lower() in sl:
        s = Unknown

    if "FUNCTIONS IN" in sl and "unknown" in sl:
        s = Unknown

    if s == Unknown:
        s = Hypothetical

    # Compact all spaces
    s = ' '.join(s.split())

    return s


def fix(args):
    """
    %prog fix ahrd.csv > ahrd.fixed.csv

    Fix ugly names from Uniprot.
    """
    p = OptionParser(fix.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    csvfile, = args
    fp = open(csvfile)
    for row in fp:
        atoms = row.rstrip("\r\n").split("\t")
        name = atoms[3]
        newname = fix_text(name)
        print "\t".join(atoms[:4] + [newname])


def merge(args):
    """
    %prog merge output/*.csv > ahrd.csv

    Merge AHRD results, remove redundant headers, empty lines, etc.
    """
    p = OptionParser(merge.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    csvfiles = args
    cf = csvfiles[0]
    fp = open(cf)
    for row in fp:
        if row.startswith("Protein"):
            break
    header = row.rstrip()
    print header

    for cf in csvfiles:
        fp = open(cf)
        for row in fp:
            if row[0] == '#':
                continue
            if row.strip() == "":
                continue
            if row.strip() == header:
                continue
            print row.strip()


def batch(args):
    """
    %prog batch splits output

    The arguments are two folders.
    Input FASTA sequences are in splits/.
    Output csv files are in output/.

    Must have folders swissprot/, tair/, trembl/ that contains the respective
    BLAST output. Once finished, you can run, for example:

    $ parallel java -Xmx2g -jar ~/code/AHRD/dist/ahrd.jar {} ::: output/*.yml
    """
    p = OptionParser(batch.__doc__)
    p.add_option("--path", default="~/code/AHRD/",
                 help="Path where AHRD is installed [default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    splits, output = args
    mkdir(output)

    for f in glob("{0}/*.fasta".format(splits)):
        fb = op.basename(f).split(".")[0]
        fw = open(op.join(output, fb + ".yml"), "w")

        path = op.expanduser(opts.path)
        dir = op.join(path, "test/resources")
        outfile = op.join(output, fb + ".csv")
        print >> fw, Template.format(dir, fb, f, outfile)


if __name__ == '__main__':
    main()
