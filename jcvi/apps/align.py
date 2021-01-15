#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Perform DNA-DNA alignment using BLAST, NUCMER and BLAT. Keep the interface the
same and does parallelization both in core and on grid.
"""
import os.path as op
import sys
import shutil
import logging

from jcvi.utils.cbook import depends
from jcvi.apps.base import (
    OptionParser,
    ActionDispatcher,
    sh,
    get_abs_path,
    which,
    mkdir,
)


@depends
def run_formatdb(infile=None, outfile=None, dbtype="nucl"):
    cmd = "makeblastdb"
    cmd += " -dbtype {0} -in {1}".format(dbtype, infile)
    sh(cmd)


@depends
def run_blat(
    infile=None,
    outfile=None,
    db="UniVec_Core",
    pctid=95,
    hitlen=50,
    cpus=16,
    overwrite=True,
):

    cmd = "pblat -threads={0}".format(cpus) if which("pblat") else "blat"
    cmd += " {0} {1} -out=blast8 {2}".format(db, infile, outfile)
    sh(cmd)

    blatfile = outfile
    filtered_blatfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
    run_blast_filter(
        infile=blatfile, outfile=filtered_blatfile, pctid=pctid, hitlen=hitlen
    )
    if overwrite:
        shutil.move(filtered_blatfile, blatfile)


@depends
def run_vecscreen(infile=None, outfile=None, db="UniVec_Core", pctid=None, hitlen=None):
    """
    BLASTN parameters reference:
    http://www.ncbi.nlm.nih.gov/VecScreen/VecScreen_docs.html
    """
    db = get_abs_path(db)
    nin = db + ".nin"
    run_formatdb(infile=db, outfile=nin)

    cmd = "blastn"
    cmd += " -task blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -penalty -5 -gapopen 4 -gapextend 4 -dust yes -soft_masking true"
    cmd += " -searchsp 1750000000000 -evalue 0.01 -outfmt 6 -num_threads 8"
    sh(cmd)


@depends
def run_megablast(
    infile=None,
    outfile=None,
    db=None,
    wordsize=None,
    pctid=98,
    hitlen=100,
    best=None,
    evalue=0.01,
    task="megablast",
    cpus=16,
):

    assert db, "Need to specify database fasta file."

    db = get_abs_path(db)
    nin = db + ".nin"
    nin00 = db + ".00.nin"
    nin = nin00 if op.exists(nin00) else (db + ".nin")
    run_formatdb(infile=db, outfile=nin)

    cmd = "blastn"
    cmd += " -query {0} -db {1} -out {2}".format(infile, db, outfile)
    cmd += " -evalue {0} -outfmt 6 -num_threads {1}".format(evalue, cpus)
    cmd += " -task {0}".format(task)
    if wordsize:
        cmd += " -word_size {0}".format(wordsize)
    if pctid:
        cmd += " -perc_identity {0}".format(pctid)
    if best:
        cmd += " -max_target_seqs {0}".format(best)
    sh(cmd)

    if pctid and hitlen:
        blastfile = outfile
        filtered_blastfile = outfile + ".P{0}L{1}".format(pctid, hitlen)
        run_blast_filter(
            infile=blastfile, outfile=filtered_blastfile, pctid=pctid, hitlen=hitlen
        )
        shutil.move(filtered_blastfile, blastfile)


def run_blast_filter(infile=None, outfile=None, pctid=95, hitlen=50):
    from jcvi.formats.blast import filter

    logging.debug("Filter BLAST result (pctid={0}, hitlen={1})".format(pctid, hitlen))
    pctidopt = "--pctid={0}".format(pctid)
    hitlenopt = "--hitlen={0}".format(hitlen)
    filter([infile, pctidopt, hitlenopt])


def main():

    actions = (
        ("blast", "run blastn using query against reference"),
        ("blat", "run blat using query against reference"),
        ("blasr", "run blasr on a set of pacbio reads"),
        ("nucmer", "run nucmer using query against reference"),
        ("last", "run last using query against reference"),
        ("lastgenome", "run whole genome LAST"),
        ("lastgenomeuniq", "run whole genome LAST and screen for 1-to-1 matches"),
        ("minimap", "run minimap2 aligner"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def minimap(args):
    """
    %prog minimap ref.fasta query.fasta

    Wrap minimap2 aligner using query against sequences. When query and ref
    is the same, we are in "self-scan" mode (e.g. useful for finding internal
    duplications resulted from mis-assemblies).
    """
    from jcvi.apps.grid import MakeManager
    from jcvi.formats.fasta import Fasta

    p = OptionParser(minimap.__doc__)
    p.add_option(
        "--chunks",
        type="int",
        default=2000000,
        help="Split ref.fasta into chunks of size in self-scan mode",
    )
    p.set_outdir(outdir="outdir")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ref, query = args
    chunks = opts.chunks
    outdir = opts.outdir
    if ref != query:
        raise NotImplementedError

    # "self-scan" mode
    # build faidx (otherwise, parallel make may complain)
    sh("samtools faidx {}".format(ref))
    f = Fasta(ref)
    mkdir(outdir)
    mm = MakeManager()
    for name, size in f.itersizes():
        start = 0
        for end in range(chunks, size, chunks):
            fafile = op.join(outdir, "{}_{}_{}.fa".format(name, start + 1, end))
            cmd = "samtools faidx {} {}:{}-{} -o {}".format(
                ref, name, start + 1, end, fafile
            )
            mm.add(ref, fafile, cmd)

            paffile = fafile.rsplit(".", 1)[0] + ".paf"
            cmd = "minimap2 -P {} {} > {}".format(fafile, fafile, paffile)
            mm.add(fafile, paffile, cmd)

            epsfile = fafile.rsplit(".", 1)[0] + ".eps"
            cmd = "minidot {} > {}".format(paffile, epsfile)
            mm.add(paffile, epsfile, cmd)
            start += chunks

    mm.write()


def nucmer(args):
    """
    %prog nucmer ref.fasta query.fasta

    Run NUCMER using query against reference. Parallel implementation derived
    from: <https://github.com/fritzsedlazeck/sge_mummer>
    """
    from itertools import product

    from jcvi.apps.grid import MakeManager
    from jcvi.formats.base import split

    p = OptionParser(nucmer.__doc__)
    p.add_option(
        "--chunks", type="int", help="Split both query and subject into chunks"
    )
    p.set_params(prog="nucmer", params="-l 100 -c 500")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ref, query = args
    cpus = opts.cpus
    nrefs = nqueries = opts.chunks or int(cpus ** 0.5)
    refdir = ref.split(".")[0] + "-outdir"
    querydir = query.split(".")[0] + "-outdir"
    reflist = split([ref, refdir, str(nrefs)]).names
    querylist = split([query, querydir, str(nqueries)]).names

    mm = MakeManager()
    for i, (r, q) in enumerate(product(reflist, querylist)):
        pf = "{0:04d}".format(i)
        cmd = "nucmer -maxmatch"
        cmd += " {0}".format(opts.extra)
        cmd += " {0} {1} -p {2}".format(r, q, pf)
        deltafile = pf + ".delta"
        mm.add((r, q), deltafile, cmd)
        print(cmd)

    mm.write()


def blasr(args):
    """
    %prog blasr ref.fasta fofn

    Run blasr on a set of PacBio reads. This is based on a divide-and-conquer
    strategy described below.
    """
    from more_itertools import grouper
    from jcvi.apps.grid import MakeManager

    p = OptionParser(blasr.__doc__)
    p.set_cpus(cpus=8)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    reffasta, fofn = args
    flist = sorted([x.strip() for x in open(fofn)])
    h5list = []
    mm = MakeManager()
    for i, fl in enumerate(grouper(flist, 3)):
        chunkname = "chunk{0:03d}".format(i)
        fn = chunkname + ".fofn"
        h5 = chunkname + ".cmp.h5"
        fw = open(fn, "w")
        print("\n".join(fl), file=fw)
        fw.close()

        cmd = "pbalign {0} {1} {2}".format(fn, reffasta, h5)
        cmd += " --nproc {0} --forQuiver --tmpDir .".format(opts.cpus)
        mm.add((fn, reffasta), h5, cmd)
        h5list.append(h5)

    # Merge h5, sort and repack
    allh5 = "all.cmp.h5"
    tmph5 = "tmp.cmp.h5"
    cmd_merge = "cmph5tools.py merge --outFile {0}".format(allh5)
    cmd_merge += " " + " ".join(h5list)
    cmd_sort = "cmph5tools.py sort --deep {0} --tmpDir .".format(allh5)
    cmd_repack = "h5repack -f GZIP=1 {0} {1}".format(allh5, tmph5)
    cmd_repack += " && mv {0} {1}".format(tmph5, allh5)
    mm.add(h5list, allh5, [cmd_merge, cmd_sort, cmd_repack])

    # Quiver
    pf = reffasta.rsplit(".", 1)[0]
    variantsgff = pf + ".variants.gff"
    consensusfasta = pf + ".consensus.fasta"
    cmd_faidx = "samtools faidx {0}".format(reffasta)
    cmd = "quiver -j 32 {0}".format(allh5)
    cmd += " -r {0} -o {1} -o {2}".format(reffasta, variantsgff, consensusfasta)
    mm.add(allh5, consensusfasta, [cmd_faidx, cmd])

    mm.write()


def get_outfile(reffasta, queryfasta, suffix="blast", outdir=None):
    q = op.basename(queryfasta).split(".")[0]
    r = op.basename(reffasta).split(".")[0]
    outfile = ".".join((q, r, suffix))
    if outdir:
        outfile = op.join(outdir, outfile)

    return outfile


def blat(args):
    """
    %prog blat ref.fasta query.fasta

    Calls blat and filters BLAST hits.
    """
    p = OptionParser(blat.__doc__)
    p.set_align(pctid=95, hitlen=30)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    reffasta, queryfasta = args
    blastfile = get_outfile(reffasta, queryfasta, suffix="blat")

    run_blat(
        infile=queryfasta,
        outfile=blastfile,
        db=reffasta,
        pctid=opts.pctid,
        hitlen=opts.hitlen,
        cpus=opts.cpus,
        overwrite=False,
    )

    return blastfile


def blast(args):
    """
    %prog blast ref.fasta query.fasta

    Calls blast and then filter the BLAST hits. Default is megablast.
    """
    task_choices = ("blastn", "blastn-short", "dc-megablast", "megablast", "vecscreen")
    p = OptionParser(blast.__doc__)
    p.set_align(pctid=0, evalue=0.01)
    p.add_option("--wordsize", type="int", help="Word size")
    p.add_option("--best", default=1, type="int", help="Only look for best N hits")
    p.add_option(
        "--task", default="megablast", choices=task_choices, help="Task of the blastn"
    )
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    reffasta, queryfasta = args
    blastfile = get_outfile(reffasta, queryfasta)

    run_megablast(
        infile=queryfasta,
        outfile=blastfile,
        db=reffasta,
        wordsize=opts.wordsize,
        pctid=opts.pctid,
        evalue=opts.evalue,
        hitlen=None,
        best=opts.best,
        task=opts.task,
        cpus=opts.cpus,
    )

    return blastfile


def lastgenome(args):
    """
    %prog genome_A.fasta genome_B.fasta

    Run LAST by calling LASTDB, LASTAL. The script runs the following steps:
    $ lastdb -P0 -uNEAR -R01 Chr10A-NEAR Chr10A.fa
    $ lastal -E0.05 -C2 Chr10A-NEAR Chr10A.fa -fTAB > Chr10A.Chr10A.tab
    $ last-dotplot Chr10A.Chr10A.tab
    """
    from jcvi.apps.grid import MakeManager

    p = OptionParser(lastgenome.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gA, gB = args
    mm = MakeManager()
    bb = lambda x: op.basename(x).rsplit(".", 1)[0]
    gA_pf, gB_pf = bb(gA), bb(gB)

    # Build LASTDB
    dbname = "-".join((gA_pf, "NEAR"))
    dbfile = dbname + ".suf"
    build_db_cmd = "lastdb -P0 -uNEAR -R01 {} {}".format(dbfile, gA)
    mm.add(gA, dbfile, build_db_cmd)

    # Run LASTAL
    tabfile = "{}.{}.tab".format(gA_pf, gB_pf)
    lastal_cmd = "lastal -E0.05 -C2 {} {}".format(dbname, gB)
    lastal_cmd += " -fTAB > {}".format(tabfile)
    mm.add([dbfile, gB], tabfile, lastal_cmd)

    mm.write()


def lastgenomeuniq(args):
    """
    %prog genome_A.fasta genome_B.fasta

    Run LAST by calling LASTDB, LASTAL and LAST-SPLIT. The recipe is based on
    tutorial here:

    <https://github.com/mcfrith/last-genome-alignments>

    The script runs the following steps:
    $ lastdb -P0 -uNEAR -R01 Chr10A-NEAR Chr10A.fa
    $ lastal -E0.05 -C2 Chr10A-NEAR Chr10B.fa | last-split -m1 | maf-swap | last-split -m1 -fMAF > Chr10A.Chr10B.1-1.maf
    $ maf-convert -n blasttab Chr10A.Chr10B.1-1.maf > Chr10A.Chr10B.1-1.blast

    Works with LAST v959.
    """
    from jcvi.apps.grid import MakeManager

    p = OptionParser(lastgenome.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    gA, gB = args
    mm = MakeManager()
    bb = lambda x: op.basename(x).rsplit(".", 1)[0]
    gA_pf, gB_pf = bb(gA), bb(gB)

    # Build LASTDB
    dbname = "-".join((gA_pf, "NEAR"))
    dbfile = dbname + ".suf"
    build_db_cmd = "lastdb -P0 -uNEAR -R01 {} {}".format(dbfile, gA)
    mm.add(gA, dbfile, build_db_cmd)

    # Run LASTAL
    maffile = "{}.{}.1-1.maf".format(gA_pf, gB_pf)
    lastal_cmd = "lastal -E0.05 -C2 {} {}".format(dbname, gB)
    lastal_cmd += " | last-split -m1"
    lastal_cmd += " | maf-swap"
    lastal_cmd += " | last-split -m1 -fMAF > {}".format(maffile)
    mm.add([dbfile, gB], maffile, lastal_cmd)

    # Convert to BLAST format
    blastfile = maffile.replace(".maf", ".blast")
    convert_cmd = "maf-convert -n blasttab {} > {}".format(maffile, blastfile)
    mm.add(maffile, blastfile, convert_cmd)

    mm.write()


@depends
def run_lastdb(
    infile=None, outfile=None, mask=False, lastdb_bin="lastdb", dbtype="nucl"
):
    outfilebase = outfile.rsplit(".", 1)[0]
    db = "-p " if dbtype == "prot" else ""
    mask = "-c " if mask else ""
    cmd = "{0} {1}{2}{3} {4}".format(lastdb_bin, db, mask, outfilebase, infile)
    sh(cmd)


def last(args, dbtype=None):
    """
    %prog database.fasta query.fasta

    Run LAST by calling LASTDB and LASTAL. LAST program available:
    <http://last.cbrc.jp>

    Works with LAST-719.
    """
    p = OptionParser(last.__doc__)
    p.add_option(
        "--dbtype",
        default="nucl",
        choices=("nucl", "prot"),
        help="Molecule type of subject database",
    )
    p.add_option("--path", help="Specify LAST path")
    p.add_option(
        "--mask", default=False, action="store_true", help="Invoke -c in lastdb"
    )
    p.add_option(
        "--format",
        default="BlastTab",
        choices=("TAB", "MAF", "BlastTab", "BlastTab+"),
        help="Output format",
    )
    p.add_option(
        "--minlen",
        default=0,
        type="int",
        help="Filter alignments by how many bases match",
    )
    p.add_option("--minid", default=0, type="int", help="Minimum sequence identity")
    p.set_cpus()
    p.set_outdir()
    p.set_params()

    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    subject, query = args
    path = opts.path
    cpus = opts.cpus
    if not dbtype:
        dbtype = opts.dbtype
    getpath = lambda x: op.join(path, x) if path else x
    lastdb_bin = getpath("lastdb")
    lastal_bin = getpath("lastal")

    subjectdb = subject.rsplit(".", 1)[0]
    run_lastdb(
        infile=subject,
        outfile=subjectdb + ".prj",
        mask=opts.mask,
        lastdb_bin=lastdb_bin,
        dbtype=dbtype,
    )

    u = 2 if opts.mask else 0
    cmd = "{0} -u {1}".format(lastal_bin, u)
    cmd += " -P {0} -i3G".format(cpus)
    cmd += " -f {0}".format(opts.format)
    cmd += " {0} {1}".format(subjectdb, query)

    minlen = opts.minlen
    minid = opts.minid
    extra = opts.extra
    assert minid != 100, "Perfect match not yet supported"
    mm = minid / (100 - minid)

    if minlen:
        extra += " -e{0}".format(minlen)
    if minid:
        extra += " -r1 -q{0} -a{0} -b{0}".format(mm)
    if extra:
        cmd += " " + extra.strip()

    lastfile = get_outfile(subject, query, suffix="last", outdir=opts.outdir)
    sh(cmd, outfile=lastfile)
    return lastfile


if __name__ == "__main__":
    main()
