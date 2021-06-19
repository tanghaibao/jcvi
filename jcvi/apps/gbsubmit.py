#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Prepare the data for Genbank submission
"""
import os.path as op
import sys
import string
import logging

from collections import defaultdict
from functools import lru_cache

from Bio import SeqIO

from jcvi.utils.orderedcollections import parse_qs
from jcvi.formats.base import DictFile
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, mkdir, glob


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

Directions = {"forward": "TR", "reverse": "TV"}

Primers = {
    "TR": "M13 Universal For 18bp Primer (TGTAAAACGACGGCCAGT)",
    "TV": "T7 Rev 20bp Primer (TAATACGACTCACTATAGGG)",
}

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
        ("fcs", "process the results from Genbank contaminant screen"),
        ("gss", "prepare package for genbank gss submission"),
        ("htg", "prepare sqn to update existing genbank htg records"),
        ("htgnew", "prepare sqn to submit new genbank htg records"),
        ("asn", "get the name tags from a bunch of asn.1 files"),
        ("t384", "print out a table converting between 96 well to 384 well"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def fcs(args):
    """
    %prog fcs fcsfile

    Process the results from Genbank contaminant screen. An example of the file
    looks like:

    contig name, length, span(s), apparent source
    contig0746      11760   1..141  vector
    contig0751      14226   13476..14226    vector
    contig0800      124133  30512..30559    primer/adapter
    """
    p = OptionParser(fcs.__doc__)
    p.add_option(
        "--cutoff",
        default=200,
        help="Skip small components less than",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fcsfile,) = args
    cutoff = opts.cutoff
    fp = open(fcsfile)
    for row in fp:
        if row[0] == "#":
            continue

        sep = "\t" if "\t" in row else None
        atoms = row.rstrip().split(sep, 3)
        contig, length = atoms[:2]
        length = int(length)
        label = atoms[-1]
        label = label.replace(" ", "_")

        if len(atoms) == 3:
            ranges = "{0}..{1}".format(1, length)
        else:
            assert len(atoms) == 4
            ranges = atoms[2]

        for ab in ranges.split(","):
            a, b = ab.split("..")
            a, b = int(a), int(b)
            assert a <= b
            ahang = a - 1
            bhang = length - b
            if ahang < cutoff:
                a = 1
            if bhang < cutoff:
                b = length
            print("\t".join(str(x) for x in (contig, a - 1, b, label)))


def asn(args):
    """
    %prog asn asnfiles

    Mainly to get this block, and extract `str` field:

        general {
          db "TIGR" ,
          tag
            str "mtg2_12952" } ,
        genbank {
          accession "AC148996" ,
    """
    from jcvi.formats.base import must_open

    p = OptionParser(asn.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fw = must_open(opts.outfile, "w")
    for asnfile in args:
        fp = open(asnfile)
        ingeneralblock = False
        ingenbankblock = False
        gb, name = None, None
        for row in fp:
            if row.strip() == "":
                continue

            tag = row.split()[0]

            if tag == "general":
                ingeneralblock = True
            if ingeneralblock and tag == "str":
                if name is None:  # Only allow first assignment
                    name = row.split('"')[1]
                ingeneralblock = False

            if tag == "genbank":
                ingenbankblock = True
            if ingenbankblock and tag == "accession":
                if gb is None:
                    gb = row.split('"')[1]
                ingenbankblock = False

        assert gb and name
        print("{0}\t{1}".format(gb, name), file=fw)


def verify_sqn(sqndir, accession):
    valfile = "{0}/{1}.val".format(sqndir, accession)
    contents = open(valfile).read().strip()
    assert not contents, "Validation error:\n{0}".format(contents)

    cmd = "gb2fasta {0}/{1}.gbf".format(sqndir, accession)
    outfile = "{0}/{1}.fasta".format(sqndir, accession)
    sh(cmd, outfile=outfile)


def htgnew(args):
    """
    %prog htgnew fastafile phasefile template.sbt

    Prepare sqnfiles for submitting new Genbank HTG records.

    `fastafile` contains the sequences.
    `phasefile` contains the phase information, it is a two column file:

    mth2-45h12    3

    `template.sbt` is the Genbank submission template.

    This function is simpler than htg, since the record names have not be
    assigned yet (so less bookkeeping).
    """
    from jcvi.formats.fasta import sequin

    p = OptionParser(htgnew.__doc__)
    p.add_option("--comment", default="", help="Comments for this submission")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    fastafile, phasefile, sbtfile = args
    comment = opts.comment

    fastadir = "fasta"
    sqndir = "sqn"
    mkdir(fastadir)
    mkdir(sqndir)

    cmd = "faSplit byname {0} {1}/".format(fastafile, fastadir)
    sh(cmd, outfile="/dev/null", errfile="/dev/null")

    acmd = "tbl2asn -a z -p fasta -r {sqndir}"
    acmd += " -i {splitfile} -t {sbtfile} -C tigr"
    acmd += ' -j "[tech=htgs {phase}] [organism=Medicago truncatula] [strain=A17]"'
    acmd += " -o {sqndir}/{accession_nv}.sqn -V Vbr"
    acmd += ' -y "{comment}" -W T -T T'

    nupdated = 0
    for row in open(phasefile):
        name, phase = row.split()[:2]
        fafile = op.join(fastadir, name + ".fa")
        cloneopt = "--clone={0}".format(name)
        splitfile, gaps = sequin([fafile, cloneopt])
        splitfile = op.basename(splitfile)
        accession = accession_nv = name

        phase = int(phase)
        assert phase in (1, 2, 3)

        cmd = acmd.format(
            accession_nv=accession_nv,
            sqndir=sqndir,
            sbtfile=sbtfile,
            splitfile=splitfile,
            phase=phase,
            comment=comment,
        )
        sh(cmd)

        verify_sqn(sqndir, accession)
        nupdated += 1

    print("A total of {0} records updated.".format(nupdated), file=sys.stderr)


def htg(args):
    """
    %prog htg fastafile template.sbt

    Prepare sqnfiles for Genbank HTG submission to update existing records.

    `fastafile` contains the records to update, multiple records are allowed
    (with each one generating separate sqn file in the sqn/ folder). The record
    defline has the accession ID. For example,
    >AC148290.3

    Internally, this generates two additional files (phasefile and namesfile)
    and download records from Genbank. Below is implementation details:

    `phasefile` contains, for each accession, phase information. For example:
    AC148290.3      3       HTG     2       mth2-45h12

    which means this is a Phase-3 BAC. Record with only a single contig will be
    labeled as Phase-3 regardless of the info in the `phasefile`. Template file
    is the Genbank sbt template. See jcvi.formats.sbt for generation of such
    files.

    Another problem is that Genbank requires the name of the sequence to stay
    the same when updating and will kick back with a table of name conflicts.
    For example:

    We are unable to process the updates for these entries
    for the following reason:

    Seqname has changed

    Accession Old seq_name New seq_name
    --------- ------------ ------------
    AC239792 mtg2_29457 AC239792.1

    To prepare a submission, this script downloads genbank and asn.1 format,
    and generate the phase file and the names file (use formats.agp.phase() and
    apps.gbsubmit.asn(), respectively). These get automatically run.

    However, use --phases if the genbank files contain outdated information.
    For example, the clone name changes or phase upgrades. In this case, run
    formats.agp.phase() manually, modify the phasefile and use --phases to override.
    """
    from jcvi.formats.fasta import sequin, ids
    from jcvi.formats.agp import phase
    from jcvi.apps.fetch import entrez

    p = OptionParser(htg.__doc__)
    p.add_option(
        "--phases",
        default=None,
        help="Use another phasefile to override",
    )
    p.add_option("--comment", default="", help="Comments for this update")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, sbtfile = args
    pf = fastafile.rsplit(".", 1)[0]

    idsfile = pf + ".ids"
    phasefile = pf + ".phases"
    namesfile = pf + ".names"

    ids([fastafile, "--outfile={0}".format(idsfile)])

    asndir = "asn.1"
    mkdir(asndir)
    entrez([idsfile, "--format=asn.1", "--outdir={0}".format(asndir)])
    asn(glob("{0}/*".format(asndir)) + ["--outfile={0}".format(namesfile)])

    if opts.phases is None:
        gbdir = "gb"
        mkdir(gbdir)
        entrez([idsfile, "--format=gb", "--outdir={0}".format(gbdir)])
        phase(glob("{0}/*".format(gbdir)) + ["--outfile={0}".format(phasefile)])
    else:
        phasefile = opts.phases

    assert op.exists(namesfile) and op.exists(phasefile)

    newphasefile = phasefile + ".new"
    newphasefw = open(newphasefile, "w")
    comment = opts.comment

    fastadir = "fasta"
    sqndir = "sqn"
    mkdir(fastadir)
    mkdir(sqndir)

    from jcvi.graphics.histogram import stem_leaf_plot

    names = DictFile(namesfile)
    assert len(set(names.keys())) == len(set(names.values()))

    phases = DictFile(phasefile)
    ph = [int(x) for x in phases.values()]
    # vmin 1, vmax 4, bins 3
    stem_leaf_plot(ph, 1, 4, 3, title="Counts of phases before updates")
    logging.debug("Information loaded for {0} records.".format(len(phases)))
    assert len(names) == len(phases)

    newph = []

    cmd = "faSplit byname {0} {1}/".format(fastafile, fastadir)
    sh(cmd, outfile="/dev/null", errfile="/dev/null")

    acmd = "tbl2asn -a z -p fasta -r {sqndir}"
    acmd += " -i {splitfile} -t {sbtfile} -C tigr"
    acmd += ' -j "{qualifiers}"'
    acmd += " -A {accession_nv} -o {sqndir}/{accession_nv}.sqn -V Vbr"
    acmd += ' -y "{comment}" -W T -T T'

    qq = "[tech=htgs {phase}] [organism=Medicago truncatula] [strain=A17]"

    nupdated = 0
    for row in open(phasefile):
        atoms = row.rstrip().split("\t")
        # see formats.agp.phase() for column contents
        accession, phase, clone = atoms[0], atoms[1], atoms[-1]
        fafile = op.join(fastadir, accession + ".fa")
        accession_nv = accession.split(".", 1)[0]

        newid = names[accession_nv]
        newidopt = "--newid={0}".format(newid)
        cloneopt = "--clone={0}".format(clone)
        splitfile, gaps = sequin([fafile, newidopt, cloneopt])
        splitfile = op.basename(splitfile)
        phase = int(phase)
        assert phase in (1, 2, 3)

        oldphase = phase
        if gaps == 0 and phase != 3:
            phase = 3

        if gaps != 0 and phase == 3:
            phase = 2

        print("{0}\t{1}\t{2}".format(accession_nv, oldphase, phase), file=newphasefw)
        newph.append(phase)

        qualifiers = qq.format(phase=phase)
        if ";" in clone:
            qualifiers += " [keyword=HTGS_POOLED_MULTICLONE]"

        cmd = acmd.format(
            accession=accession,
            accession_nv=accession_nv,
            sqndir=sqndir,
            sbtfile=sbtfile,
            splitfile=splitfile,
            qualifiers=qualifiers,
            comment=comment,
        )
        sh(cmd)

        verify_sqn(sqndir, accession)
        nupdated += 1

    stem_leaf_plot(newph, 1, 4, 3, title="Counts of phases after updates")
    print("A total of {0} records updated.".format(nupdated), file=sys.stderr)


@lru_cache(maxsize=None)
def get_rows_cols(nrows=Nrows, ncols=Ncols):
    rows, cols = string.ascii_uppercase[:nrows], range(1, ncols + 1)
    return rows, cols


@lru_cache(maxsize=None)
def get_plate(nrows=Nrows, ncols=Ncols):

    rows, cols = get_rows_cols(nrows, ncols)
    plate = [[""] * ncols for _ in range(nrows)]
    n = 0
    # 384 to (96+quadrant)
    for i in range(0, nrows, 2):
        for j in range(0, ncols, 2):
            n += 1
            prefix = "{0:02d}".format(n)
            plate[i][j] = prefix + "A"
            plate[i][j + 1] = prefix + "B"
            plate[i + 1][j] = prefix + "C"
            plate[i + 1][j + 1] = prefix + "D"

    # (96+quadrant) to 384
    splate = {}
    for i in range(nrows):
        for j in range(ncols):
            c = plate[i][j]
            splate[c] = "{0}{1}".format(rows[i], j + 1)

    return plate, splate


def convert_96_to_384(c96, quad, ncols=Ncols):
    """
    Convert the 96-well number and quad number to 384-well number

    >>> convert_96_to_384("B02", 1)
    'C3'
    >>> convert_96_to_384("H09", 4)
    'P18'
    """
    rows, cols = get_rows_cols()
    plate, splate = get_plate()

    n96 = rows.index(c96[0]) * ncols // 2 + int(c96[1:])
    q = "{0:02d}{1}".format(n96, "ABCD"[quad - 1])
    return splate[q]


def t384(args):
    """
    %prog t384

    Print out a table converting between 96 well to 384 well
    """
    p = OptionParser(t384.__doc__)
    p.parse_args(args)

    plate, splate = get_plate()

    fw = sys.stdout
    for i in plate:
        for j, p in enumerate(i):
            if j != 0:
                fw.write("|")
            fw.write(p)
        fw.write("\n")


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
    print(PublicationTemplate.format(**vars), file=fw)
    print(LibraryTemplate.format(**vars), file=fw)
    print(ContactTemplate.format(**vars), file=fw)
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
        sequencer_plate_well_coordinates = a["sequencer_plate_well_coordinates"][0]
        sequencer_plate_96well_quadrant = a["sequencer_plate_96well_quadrant"][0]
        sequencer_plate_96well_coordinates = a["sequencer_plate_96well_coordinates"][0]

        # Check the 96-well ID is correctly converted to 384-well ID
        w96 = sequencer_plate_96well_coordinates
        w96quad = int(sequencer_plate_96well_quadrant)
        w384 = sequencer_plate_well_coordinates
        assert convert_96_to_384(w96, w96quad) == w384

        plate = sequencer_plate_barcode
        assert plate in plateMapping, "{0} not found in `{1}` !".format(
            plate, mappingfile
        )

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
        sequencer_plate_well_coordinates = a["sequencer_plate_well_coordinates"][0]
        w384 = sequencer_plate_well_coordinates

        plate = sequencer_plate_barcode
        plate = plateMapping[plate]
        d = Directions[direction]

        cloneID = "{0}{1}".format(plate, w384)
        gssID = "{0}{1}".format(cloneID, d)
        seen[gssID] += 1

        if seen[gssID] > 1:
            logging.error("duplicate key {0} found".format(gssID))
            gssID = "{0}{1}".format(gssID, seen[gssID])

        othergss = clone[cloneID] - {gssID}
        othergss = ", ".join(sorted(othergss))
        vars.update(locals())

        print(GSSTemplate.format(**vars), file=fw)

        # Write conversion logs to log file
        print("{0}\t{1}".format(gssID, description), file=fw_log)
        print("=" * 60, file=fw_log)

    logging.debug("A total of {0} seqs written to `{1}`".format(len(seen), fw.name))
    fw.close()
    fw_log.close()


if __name__ == "__main__":
    import doctest

    doctest.testmod()
    main()
