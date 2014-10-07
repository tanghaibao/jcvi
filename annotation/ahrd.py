#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utility to run Automated Human Readable Description (AHRD) pipeline.

<https://github.com/groupschoof/AHRD>
"""

import os.path as op
from os import symlink
import sys
import re
import logging

from jcvi.formats.base import must_open
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, glob


##### Compiled RegExps #####
# Cellular locations
loc_pat = re.compile(r",\s*(chloroplastic|cytoplasmic|mitochondrial).*?\s$", re.I)
# Any word that matches e.g. Os02g0234800
osg_pat = re.compile(r"\bOs\d{2}g\d{7}.*?\s", re.I)
# (fragment)
frag_pat = re.compile(r"\(fragment[s]?\)", re.I)
# Trailing protein numeric copy (e.g. Myb 1)
trail_pat = re.compile(r"(?<!type)\s\d+\s*$", re.I)
# UPF
upf_pat = re.compile(r"^UPF.*$")
# Remove 'DDB_G\d+' ID
ddb_pat = re.compile(r"\s+DDB_G\d+", re.I)

# Any AHRD that matches e.g. "AT5G54690-like protein"
atg_pat = re.compile(r"\bAT[1-5M]G\d{5}-like protein", re.I)

# remove 'arabidopsis thaliana'
atg_id_pat = re.compile(r"[_]*AT\d{1}G\d+[/]*", re.I)
athila_pat1 = re.compile(r"Belongs to|^Encodes|^Expression|^highly", re.I)
athila_pat2 = re.compile(r"^Arabidopsis thaliana ", re.I)
athila_pat3 = re.compile(r"^Arabidopsis ", re.I)
athila_pat4 = re.compile(r"BEST Arabidopsis thaliana protein match is: ", re.I)

# &apos;? => '
apos_pat = re.compile(r"&apos;?")

# &gt => none
gt_pat = re.compile(r"&gt")

# -like to -like protein
like_pat = re.compile(r"[-]like$", re.I)

# 'repeat$' to 'repeat protein'
repeat_pat = re.compile(r"repeat$", re.I)

# re used by the following 3 cases
Protein_pat = re.compile(r"Protein\s+", re.I)

# 'binding$' to 'binding protein'
binding_pat = re.compile(r"binding$", re.I)

# 'domain$' to 'domain-containing protein'
domain_pat = re.compile(r"domain$", re.I)

# 'related$' to '-like protein'
related_pat = re.compile(r"[,\s+]*[\s+|-]*related$", re.I)

# '[0-9]+ homolog' to '-like protein'
homolog_pat1 = re.compile(r"(?<!type)\s+\d+\s+homolog.*$", re.I)

# 'Protein\s+(.*)\s+homolog' to '$1-like protein'
homolog_pat2 = re.compile(r"^Protein([\s+\S+]+)\s+homolog.*", re.I)

# 'homolog protein' to '-like protein'
homolog_pat3 = re.compile(r"\s+homolog\s+protein.*", re.I)
# 'homolog \S+' to '-like protein'
homolog_pat4 = re.compile(r"\s+homolog\s+\S+$", re.I)
# 'homologue$' to '-like protein'
homolog_pat5 = re.compile(r"\s+homologue[\s+\S+]$", re.I)
# 'homolog$' to '-like protein'
homolog_pat6 = re.compile(r"\s+homolog$", re.I)

# 'Agenet domain-containing protein / bromo-adjacent homology (BAH) domain-containing protein'
# to 'Agenet and bromo-adjacent homology (BAH) domain-containing protein'
agenet_pat = re.compile(r"Agenet domain-containing protein \/ ", re.I)

# plural to singular
plural_pat = re.compile(r"[deinr]s$", re.I)

# 'like_TBP' or 'likeTBP' to 'like TBP'
tbp_pat = re.compile(r"like[_]*TBP", re.I)

# 'protein protein' to 'protein'
prot_pat = re.compile(r" protein protein", re.I)

# 'Candidate|Hypothetical|Novel|Predicted|Possible' to 'Putative'
put_pat = re.compile(r"Candidate|Hypothetical|Novel|Predicted|Possible", re.I)

# 'dimerisation' to 'dimerization'
dimer_pat = re.compile(r"dimerisation", re.I)

# '\s+LENGTH=\d+' to ''
length_pat = re.compile(r"\s+LENGTH\=\d+", re.I)

# disallowed words
disallow = ("genome", "annotation", "project")
disallow_pat = re.compile("|".join(str(x) for x in disallow))

# disallowed organism names
organism = ("thaliana", "rickettsia", "rice", "yeast")
organism_pat = re.compile("|".join("^.*{0}".format(str(x)) for x in organism))

# consolidate glycosidic links
glycosidic_link_pat = re.compile("\d+,\d+")

# Kevin Silverstein suggested names (exclude list)
spada = ("LCR", "RALF", "SCR")
spada_pat = re.compile("|".join("^{0}$".format(str(x)) for x in spada))

# names with all capital letters (maybe followed by numbers)
sym_pat = re.compile(r"^[A-Z]+[A-Z0-9\-]{0,}$")
lc_sym_pat = re.compile(r"^[A-z]{1}[a-z]+[0-9]{1,}$")
eol_sym_pat = re.compile(r"\([A-Z]+[A-Z0-9\-]{0,}\)$")

# sulfer -> sulfur
sulfer_pat = re.compile(r"sulfer")

# assessory -> accessory
assessory_pat = re.compile(r"assessory")

# british to american spelling conversion
# -ise -> -ize
# -isation -> ization
ise_pat = re.compile(r"\b([A-z]+)ise\b")
isation_pat = re.compile(r"\b([A-z]+)isation\b")

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
{7}
token_score_bit_score_weight: {4}
token_score_database_score_weight: {5}
token_score_overlap_score_weight: {6}
description_score_relative_description_frequency_weight: 0.6
output: {3}
"""

iprscanTemplate = """
interpro_database: ./interpro.xml
interpro_result: {0}
"""

# Necessary for the script to know the location of `interpro.xml` and `interpro.dtd`
iprscan_datadir = "/usr/local/devel/ANNOTATION/iprscan/iprscan_v4.7/data"


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

    # Fix descriptions like D7TDB1 (
    s = re.sub("([A-Z0-9]){6} \(", "", s)
    s = s.translate(None, "[]")
    s = s.replace("(-)", "[-]")
    s = s.replace("(+)", "[+]")
    s = s.replace("(Uncharacterized protein)", "")
    s = s.translate(None, "()")

    # before trimming off at the first ";", check if name has glycosidic
    # linkage information (e.g 1,3 or 1,4). If so, also check if multiple
    # linkages are separated by ";". If so, replace ";" by "-"
    m = re.findall(glycosidic_link_pat, s)
    if m and ";" in s:
        s = re.sub(";\s*", "-", s)

    s = s.split(";")[0]

    # Cellular locations
    # Any word that matches e.g. AT5G54690
    # Any word that matches e.g. Os02g0234800
    # (fragment)
    # Trailing protein numeric copy (e.g. Myb 1)
    # UPF
    # Remove 'DDB_G\d+' ID
    # '_At[0-9]+g[0-9]+' to ''
    for pat in (loc_pat, osg_pat, frag_pat, trail_pat, upf_pat, ddb_pat):
        # below is a hack since word boundaries don't work on /
        s = s.strip() + " "
        s = re.sub(pat, "", s)

    # &apos;? => '
    s = re.sub(apos_pat, "'", s)
    # &gt => none
    s = re.sub(gt_pat, "", s)
    # reduce runs such as -- '''
    s = re.sub(r"[-]+", "-", s)
    s = re.sub(r"[']+", "'", s)

    s = s.strip()

    # -like to -like protein
    s = re.sub(like_pat, "-like protein", s)

    # 'repeat$' to 'repeat protein'
    if re.search(repeat_pat, s):
        s += "-containing protein"

    # 'binding$' to 'binding protein'
    if re.search(binding_pat, s):
        s += " protein"
        if re.match(Protein_pat, s):
            s = re.sub(Protein_pat, "", s)

    # 'domain$' to 'domain-containing protein'
    if re.search(domain_pat, s):
        s += "-containing protein"
        if re.search(r"-domain", s):
            s = re.sub(r"-domain", " domain", s)
        if re.match(Protein_pat, s):
            s = re.sub(Protein_pat, "", s)

    # 'related$' to '-like protein'
    if re.search(related_pat, s):
        s = re.sub(related_pat, "-like protein", s)
        if re.match(Protein_pat, s) and not re.match(r"Protein kinase", s):
            s = re.sub(Protein_pat, "", s)

    # '[0-9]+ homolog' to '-like protein'
    if re.search(homolog_pat1, s):
        s = re.sub(homolog_pat1, "-like protein", s)
        if re.match(Protein_pat, s):
            s = re.sub(Protein_pat, "", s)

    # 'Protein\s+(.*)\s+homolog' to '$1-like protein'
    match = re.search(homolog_pat2, s)
    if match and not re.match(r"Protein kinase", s):
        ret = match.group(1)
        s = re.sub(homolog_pat2, ret + "-like protein", s)
        s = re.sub(r"^\s+", "", s)
        s = s.capitalize()

    # 'homolog protein' to '-like protein'
    # 'homolog \S+' to '-like protein'
    # 'homologue$' to '-like protein'
    # 'homolog$' to '-like protein'
    for pat in (homolog_pat3, homolog_pat4, homolog_pat5, homolog_pat6):
        if re.search(pat, s):
            s = re.sub(pat, "-like protein", s)

    # 'Agenet domain-containing protein / bromo-adjacent homology (BAH) domain-containing protein'
    # to 'Agenet and bromo-adjacent homology (BAH) domain-containing protein'
    if re.search(agenet_pat, s):
        s = re.sub(agenet_pat, "Agenet and ", s)

    # plural to singular
    if re.search(plural_pat, s):
        #if s.find('biogenesis') == -1 and s.find('Topors') == -1 and s.find('allergens') == -1:
        if s.find('biogenesis') == -1 and s.find('Topors') == -1:
            s = re.sub(r"s$", "", s)

    # 'like_TBP' or 'likeTBP' to 'like TBP'
    if re.search(tbp_pat, s):
        s = re.sub(tbp_pat, "like TBP", s)

    # 'protein protein' to 'protein'
    if re.search(prot_pat, s):
        s = re.sub(prot_pat, " protein", s)

    # 'Candidate|Hypothetical|Novel|Predicted|Possible' to 'Putative'
    if re.search(put_pat, s):
        s = re.sub(put_pat, "Putative", s)

    # 'dimerisation' to 'dimerization'
    if re.search(dimer_pat, s):
        s = re.sub(dimer_pat, "dimerization", s)

    # Any AHRD that matches e.g. "AT5G54690-like protein"
    # Any AHRD that contains the words '^Belongs|^Encoded|^Expression|^highly'
    for pat in (atg_pat, athila_pat1):
        if re.search(pat, s):
            s = Unknown

    # remove 'arabidopsis[ thaliana]' and/or embedded Atg IDs
    for pat in (atg_id_pat, athila_pat2, athila_pat3, athila_pat4):
        # below is a hack since word boundaries don't work on /
        s = s.strip() + " "
        s = re.sub(pat, "", s)

    # remove "\s+LENGTH=\d+" from TAIR deflines
    if re.search(length_pat, s):
        s = re.sub(length_pat, "", s)

    # if name has a dot followed by a space (". ") in it and contains multiple
    # parts separated by a comma, strip name starting from first occurrence of ","
    if re.search(r"\. ", s):
        if re.search(r",", s):
            s = s.split(",")[0]

    # if name contains any of the disallowed words,
    # remove word occurrence from name
    # if name contains references to any other organism, trim name upto
    # that occurrence
    for pat in (disallow_pat, organism_pat):
        if re.search(pat, s):
            s = re.sub(pat, "", s)

    s = s.strip()

    # if name is entirely a gene symbol-like (all capital letters, maybe followed by numbers)
    # add a "-like protein" at the end
    if (re.search(sym_pat, s) or re.search(lc_sym_pat, s)) \
            and not re.search(spada_pat, s):
        s = s + "-like protein"

    # if gene symbol in parantheses at EOL, remove symbol
    if re.search(eol_sym_pat, s):
        s = re.sub(eol_sym_pat, "", s)

    # if name terminates at a symbol([^A-Za-z0-9_]), trim it off
    if re.search(r"\W{1,}$", s) and not re.search(r"\)$", s):
        s = re.sub("\W{1,}$", "", s)

    # change sulfer to sulfur
    if re.search(sulfer_pat, s):
        s = re.sub(sulfer_pat, "sulfur", s)

    # change assessory to accessory
    if re.search(assessory_pat, s):
        s = re.sub(assessory_pat, "accessory", s)

    # change -ise/-isation to -ize/-ization
    match = re.search(ise_pat, s)
    if match:
        ret = match.group(1)
        s = re.sub(ise_pat, "{0}ize".format(ret), s)

    match = re.search(isation_pat, s)
    if match:
        ret = match.group(1)
        s = re.sub(isation_pat, "{0}ization".format(ret), s)

    """
    case (qr/^Histone-lysine/) { $ahrd =~ s/,\s+H\d{1}\s+lysine\-\d+//gs; }
    """
    sl = s.lower()

    # Any mention of `clone` or `contig` is not informative
    if "clone" in sl or "contig" in sl:
        s = Unknown

    # All that's left is `protein` is not informative
    if sl in ("protein", "protein, putative", ""):
        s = Unknown

    if Unknown.lower() in sl:
        s = Unknown

    if "FUNCTIONS IN".lower() in sl and "unknown" in sl:
        s = Unknown

    if "uncharacterized" in sl:
        s = "uncharacterized protein"

    s = s.replace(", putative", "")
    s = s.replace("Putative ", "")

    if s == Unknown or s.strip() == "protein":
        s = Hypothetical

    # Compact all spaces
    s = ' '.join(s.split())

    assert s.strip()

    return s


def fix(args):
    """
    %prog fix ahrd.csv > ahrd.fixed.csv

    Fix ugly names from Uniprot.
    """
    p = OptionParser(fix.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    csvfile, = args
    fp = open(csvfile)
    fw = must_open(opts.outfile, "w")
    for row in fp:
        if row[0] == '#':
            continue
        if row.strip() == "":
            continue
        atoms = row.rstrip("\r\n").split("\t")
        name, hit, ahrd_code, desc = atoms[:4] \
                if len(atoms) > 2 else \
                atoms[0], None, None, atoms[-1]
        newdesc = fix_text(desc)
        if hit and hit.strip() != "" and newdesc == Hypothetical:
            newdesc = "conserved " + newdesc
        print >> fw, "\t".join(atoms[:4] + [newdesc] + atoms[4:])


def merge(args):
    """
    %prog merge output/*.csv > ahrd.csv

    Merge AHRD results, remove redundant headers, empty lines, etc. If there are
    multiple lines containing the same ID (first column). Then whatever comes
    the first will get retained.
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

    seen = set()
    for cf in csvfiles:
        fp = open(cf)
        for row in fp:
            if row[0] == '#':
                continue
            if row.strip() == "":
                continue
            if row.strip() == header:
                continue

            atoms = row.rstrip().split("\t")
            id = atoms[0]
            if id in seen:
                logging.error("ID `{0}` ignored.".format(id))
                continue

            seen.add(id)
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

    ahrd_weights = { "blastp": [0.5, 0.3, 0.2],
                     "blastx": [0.6, 0.4, 0.0]
                   }
    blast_progs = tuple(ahrd_weights.keys())

    p.add_option("--path", default="~/code/AHRD/",
                 help="Path where AHRD is installed [default: %default]")
    p.add_option("--blastprog", default="blastp", choices=blast_progs,
                help="Specify the blast program being run. Based on this option," \
                   + " the AHRD parameters (score_weights) will be modified." \
                   + " [default: %default]")
    p.add_option("--iprscan", default=None,
                help="Specify path to InterProScan results file if available." \
                   + " If specified, the yml conf file will be modified" \
                   + " appropriately. [default: %default]")

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    splits, output = args
    mkdir(output)

    bit_score, db_score, ovl_score = ahrd_weights[opts.blastprog]

    for f in glob("{0}/*.fasta".format(splits)):
        fb = op.basename(f).rsplit(".", 1)[0]
        fw = open(op.join(output, fb + ".yml"), "w")

        path = op.expanduser(opts.path)
        dir = op.join(path, "test/resources")
        outfile = op.join(output, fb + ".csv")
        interpro = iprscanTemplate.format(opts.iprscan) if opts.iprscan else ""

        print >> fw, Template.format(dir, fb, f, outfile, bit_score, db_score, ovl_score, interpro)

    if opts.iprscan:
        if not op.lexists("interpro.xml"):
            symlink(op.join(iprscan_datadir, "interpro.xml"), "interpro.xml")

        if not op.lexists("interpro.dtd"):
            symlink(op.join(iprscan_datadir, "interpro.dtd"), "interpro.dtd")

if __name__ == '__main__':
    main()
