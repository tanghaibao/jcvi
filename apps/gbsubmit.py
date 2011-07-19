#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Prepare the data for Genbank submission
"""

import os.path as op
import sys
import string
import logging

from glob import glob
from urlparse import parse_qs
from optparse import OptionParser
from collections import defaultdict

from Bio import SeqIO

from jcvi.utils.cbook import memoized, fill
from jcvi.formats.base import DictFile
from jcvi.apps.console import print_red
from jcvi.apps.base import ActionDispatcher, sh, mkdir, debug
debug()

"""
GSS submission template files

<http://www.ncbi.nlm.nih.gov/dbGSS/how_to_submit.html>
"""

# Modify the following if a different submission
# TODO: make this generic and exist outside source code
Title = """Comparative Genomics of Sisymbrium irio"""

Authors = """Town,C.D., Tang,H., Paterson,A.H. and Pires,J.C."""

Libname = "Sisymbrium irio BAC library SIL"
Contact = "Chris D. Town"

PublicationTemplate = """TYPE: Pub
MEDUID:
TITLE:
{Title}
AUTHORS:
{Authors}
JOURNAL:
VOLUME:
ISSUE:
PAGES:
YEAR: 2011
STATUS: 1
||"""

LibraryTemplate = """TYPE: Lib
NAME: {Libname}
ORGANISM: Sisymbrium irio
STRAIN: Gomez-Campo 1146-67
SEX:
STAGE:
TISSUE:
CELL_TYPE:
VECTOR: pCC1BAC
RE_1: HindIII
DESCR:
Constructed by Amplicon Express;
Transformed into Invitrogen DH10b phage resistant E. coli.
||"""

ContactTemplate = """TYPE: Cont
NAME: {Contact}
FAX: 301-795-7070
TEL: 301-795-7523
EMAIL: cdtown@jcvi.org
LAB: Plant Genomics
INST: J. Craig Venter Institute
ADDR: 9704 Medical Center Dr., Rockville, MD 20850, USA
||"""

Directions = {"forward": "TR",
              "reverse": "TV"}

Primers = {"TR": "M13 Universal For 18bp Primer (TGTAAAACGACGGCCAGT)",
           "TV": "T7 Rev 20bp Primer (TAATACGACTCACTATAGGG)"}

GSSTemplate = """TYPE: GSS
STATUS: New
CONT_NAME: {Contact}
GSS#: {gssID}
CLONE: {cloneID}
SOURCE: JCVI
OTHER_GSS: {othergss}
CITATION:
{Title}
INSERT: 120000
PLATE: {plate}
ROW: {row}
COLUMN: {column}
SEQ_PRIMER: {primer}
DNA_TYPE: Genomic
CLASS: BAC ends
LIBRARY: {Libname}
PUBLIC:
PUT_ID:
COMMENT:
SEQUENCE:
{seq}
||"""

Nrows, Ncols = 16, 24
vars = globals()


def main():

    actions = (
        ('gss', 'prepare package for genbank gss submission'),
        ('htg', 'prepare sqn for genbank htg submission'),
        ('t384', 'print out a table converting between 96 well to 384 well'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def htg(args):
    """
    %prog htg fastafile phasefile template.sbt

    Prepare sqnfiles for Genbank HTG submission.

    `fastafile` contains the records to update, multiple records are allowed
    (with each one generating separate sqn file in the sqn/ folder). The record
    defline has the accession ID. For example,
    >AC148290.3

    `phasefile` is a 2-column file containing the mappings. For example:
    AC148290.3  2

    which means this is a Phase-2 BAC. Record with only a single contig will be
    labeled as Phase-3 regardless of the info in the `phasefile`. Template file
    is the Genbank sbt template. See jcvi.formats.sbt for generation of such
    files.
    """
    from jcvi.formats.fasta import Fasta, sequin

    p = OptionParser(htg.__doc__)
    p.add_option("--comment",
            help="Comments for this update [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    fastafile, phasefile, sbtfile = args
    comment = opts.comment

    fastadir = "fasta"
    sqndir = "sqn"
    mkdir(fastadir, overwrite=True)
    mkdir(sqndir, overwrite=True)

    cmd = "faSplit byname {0} {1}/".format(fastafile, fastadir)
    sh(cmd)

    acmd = 'tbl2asn -a z -p fasta -r {sqndir}'
    acmd += ' -i {splitfile} -t {sbtfile} -C tigr'
    acmd += ' -j "[tech=htgs {phase}] [organism=Medicago truncatula] [strain=A17]"'
    acmd += ' -A {accession_nv} -o {sqndir}/{accession_nv}.sqn -V Vbr'
    acmd += ' -y "{comment}" -W T -T T'
    for fafile in glob(fastadir + "/*.fa"):
        splitfile, gaps = sequin([fafile])
        splitfile = op.basename(splitfile)
        accession = op.basename(fafile).rsplit(".", 1)[0]
        accession_nv = accession.split(".", 1)[0]
        phase = 2

        if gaps == 0 and phase != 3:
            newphase = 3
            logging.debug("Promote {0} from Phase {1} => {2}".\
                    format(accession, phase, newphase))
            phase = newphase

        if gaps != 0 and phase == 3:
            newphase = 2
            logging.debug("Demote {0} from Phase {1} => {2}".\
                    format(accession, phase, newphase))
            phase = newphase

        cmd = acmd.format(accession=accession, accession_nv=accession_nv,
                sqndir=sqndir, sbtfile=sbtfile, splitfile=splitfile,
                phase=phase, comment=comment)
        sh(cmd)

        valfile = "{0}/{1}.val".format(sqndir, accession)
        contents = open(valfile).read().strip()
        assert not contents, "Validation error:\n{0}".format(contents)


@memoized
def get_rows_cols(nrows=Nrows, ncols=Ncols):
    rows, cols = string.ascii_uppercase[:nrows], range(1, ncols + 1)
    return rows, cols


@memoized
def get_plate(nrows=Nrows, ncols=Ncols):

    rows, cols = get_rows_cols(nrows, ncols)
    plate = [[""] * ncols for x in xrange(nrows)]
    n = 0
    # 384 to (96+quadrant)
    for i in xrange(0, nrows, 2):
        for j in xrange(0, ncols, 2):
            n += 1
            prefix = "{0:02d}".format(n)
            plate[i][j] = prefix + 'A'
            plate[i][j + 1] = prefix + 'B'
            plate[i + 1][j] = prefix + 'C'
            plate[i + 1][j + 1] = prefix + 'D'

    # (96+quadrant) to 384
    splate = {}
    for i in xrange(nrows):
        for j in xrange(ncols):
            c = plate[i][j]
            splate[c] = "{0}{1}".format(rows[i], j + 1)

    return plate, splate


def convert_96_to_384(c96, quad, nrows=Nrows, ncols=Ncols):
    """
    Convert the 96-well number and quad number to 384-well number

    >>> convert_96_to_384("B02", 1)
    'C3'
    >>> convert_96_to_384("H09", 4)
    'P18'
    """
    rows, cols = get_rows_cols()
    plate, splate = get_plate()

    n96 = rows.index(c96[0]) * ncols / 2 + int(c96[1:])
    q = "{0:02d}{1}".format(n96, "ABCD"[quad - 1])
    return splate[q]


def t384(args):
    """
    %prog t384

    Print out a table converting between 96 well to 384 well
    """
    p = OptionParser(t384.__doc__)
    opts, args = p.parse_args(args)

    plate, splate = get_plate()

    fw = sys.stdout
    for i in plate:
        for j, p in enumerate(i):
            if j != 0:
                fw.write('|')
            fw.write(p)
        fw.write('\n')


def parse_description(s):
    """
    Returns a dictionary based on the FASTA header, assuming JCVI data
    """
    s = "".join(s.split()[1:]).replace("/", ";")
    a = parse_qs(s)
    return a


def gss(args):
    """
    %prog gss fastafile plateMapping

    Generate sequence files and metadata templates suited for gss submission.
    The FASTA file is assumed to be exported from the JCVI data delivery folder
    which looks like:

    >1127963806024 /library_name=SIL1T054-B-01-120KB /clear_start=0
    /clear_end=839 /primer_id=1049000104196 /trace_id=1064147620169
    /trace_file_id=1127963805941 /clone_insert_id=1061064364776
    /direction=reverse /sequencer_run_id=1064147620155
    /sequencer_plate_barcode=B906423 /sequencer_plate_well_coordinates=C3
    /sequencer_plate_96well_quadrant=1 /sequencer_plate_96well_coordinates=B02
    /template_plate_barcode=CC0251602AB /growth_plate_barcode=BB0273005AB
    AGCTTTAGTTTCAAGGATACCTTCATTGTCATTCCCGGTTATGATGATATCATCAAGATAAACAAGAATG
    ACAATGATACCTGTTTGGTTCTGAAGTGTAAAGAGGGTATGTTCAGCTTCAGATCTTCTAAACCCTTTGT
    CTAGTAAGCTGGCACTTAGCTTCCTATACCAAACCCTTTGTGATTGCTTCAGTCCATAAATTGCCTTTTT

    Plate mapping file maps the JTC `sequencer_plate_barcode` to external IDs.
    For example:
    B906423 SIL-001
    """
    p = OptionParser(gss.__doc__)

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(p.print_help())

    fastafile, mappingfile = args
    seen = defaultdict(int)
    clone = defaultdict(set)

    plateMapping = DictFile(mappingfile)

    fw = open("MetaData.txt", "w")
    print >> fw, PublicationTemplate.format(**vars)
    print >> fw, LibraryTemplate.format(**vars)
    print >> fw, ContactTemplate.format(**vars)
    logging.debug("Meta data written to `{0}`".format(fw.name))

    fw = open("GSS.txt", "w")
    fw_log = open("GSS.log", "w")
    for rec in SeqIO.parse(fastafile, "fasta"):
        # First pass just check well number matchings and populate sequences in
        # the same clone
        description = rec.description
        a = parse_description(description)
        direction = a["direction"][0]
        sequencer_plate_barcode = a["sequencer_plate_barcode"][0]
        sequencer_plate_well_coordinates = \
            a["sequencer_plate_well_coordinates"][0]
        sequencer_plate_96well_quadrant = \
            a["sequencer_plate_96well_quadrant"][0]
        sequencer_plate_96well_coordinates = \
            a["sequencer_plate_96well_coordinates"][0]

        # Check the 96-well ID is correctly converted to 384-well ID
        w96 = sequencer_plate_96well_coordinates
        w96quad = int(sequencer_plate_96well_quadrant)
        w384 = sequencer_plate_well_coordinates
        assert convert_96_to_384(w96, w96quad) == w384

        plate = sequencer_plate_barcode
        assert plate in plateMapping, \
            "{0} not found in `{1}` !".format(plate, mappingfile)

        plate = plateMapping[plate]
        d = Directions[direction]

        cloneID = "{0}{1}".format(plate, w384)
        gssID = "{0}{1}".format(cloneID, d)
        seen[gssID] += 1

        if seen[gssID] > 1:
            gssID = "{0}{1}".format(gssID, seen[gssID])

        seen[gssID] += 1
        clone[cloneID].add(gssID)

    seen = defaultdict(int)
    for rec in SeqIO.parse(fastafile, "fasta"):
        # need to populate gssID, mateID, cloneID, seq, plate, row, column
        description = rec.description
        a = parse_description(description)
        direction = a["direction"][0]
        sequencer_plate_barcode = a["sequencer_plate_barcode"][0]
        sequencer_plate_well_coordinates = \
            a["sequencer_plate_well_coordinates"][0]
        w384 = sequencer_plate_well_coordinates

        plate = sequencer_plate_barcode
        plate = plateMapping[plate]
        d = Directions[direction]

        row = w384[0]
        column = int(w384[1:])
        seq = fill(str(rec.seq), width=70)

        cloneID = "{0}{1}".format(plate, w384)
        gssID = "{0}{1}".format(cloneID, d)
        primer = Primers[d]
        seen[gssID] += 1

        if seen[gssID] > 1:
            logging.error("duplicate key {0} found".format(gssID))
            gssID = "{0}{1}".format(gssID, seen[gssID])

        othergss = clone[cloneID] - set([gssID])
        othergss = ", ".join(sorted(othergss))
        vars.update(locals())

        print >> fw, GSSTemplate.format(**vars)

        # Write conversion logs to log file
        print >> fw_log, "{0}\t{1}".format(gssID, description)
        print >> fw_log, "=" * 60

    logging.debug("A total of {0} seqs written to `{1}`".\
            format(len(seen), fw.name))
    fw.close()
    fw_log.close()


if __name__ == '__main__':
    import doctest
    doctest.testmod()
    main()
