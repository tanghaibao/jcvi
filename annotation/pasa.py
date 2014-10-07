#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utilities for submitting PASA jobs and processing PASA results.
"""

import os
import os.path as op
import sys
import logging

from jcvi.formats.base import write_file, must_open, FileMerger
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, \
        which, mkdir


alignAssembly_conf = """
# MySQL settings
MYSQLDB={0}

#######################################################
# Parameters to specify to specific scripts in pipeline
# create a key = "script_name" + ":" + "parameter"
# assign a value as done above.

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED={1}
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID={2}
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY={3}

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
"""

tdn, tfasta = "tdn.accs", "transcripts.fasta"
aaconf, acconf = "alignAssembly.conf", "annotCompare.conf"
ALLOWED_ALIGNERS = ("blat", "gmap")


def main():

    actions = (
        ('assemble', 'run pasa alignment assembly pipeline'),
        ('longest', 'label longest transcript per gene as full-length'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def assemble(args):
    """
    %prog assemble pasa_db_name genome.fasta transcripts-dn.fasta [transcript-gg.fasta]

    Run the PASA alignment assembly pipeline

    If two transcript fasta files (Trinity denovo and genome guided) are provided,
    the PASA Comprehensive Transcriptome protocol is followed
    <http://pasa.sourceforge.net/#A_ComprehensiveTranscriptome>

    Using the `--prepare` option creates a shell script with the run commands without
    executing the pipeline
    """
    p = OptionParser(assemble.__doc__)
    p.set_home("pasa")
    p.set_align(pctid=95, pctcov=90, intron=2000, bpsplice=3, compreh_pctcov=30)
    p.add_option("--aligners", default="blat,gmap",
            help="Specify splice aligners to use for mapping [default: %default]")
    p.add_option("--clean", default=False, action="store_true",
            help="Clean transcripts using tgi seqclean [default: %default]")
    p.set_cpus()
    p.set_grid()
    p.set_grid_opts()
    p.add_option("--prepare", default=False, action="store_true",
            help="Prepare PASA run script with commands [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) not in (3, 4):
        sys.exit(not p.print_help())

    pasa_db, genome, dnfasta, = args[:3]
    ggfasta = args[3] if len(args) == 4 else None

    PASA_HOME = opts.pasa_home
    if not op.isdir(PASA_HOME):
        logging.error("PASA_HOME={0} directory does not exist".format(PASA_HOME))
        sys.exit()

    aligners = opts.aligners.split(",")
    for aligner in aligners:
        if aligner not in ALLOWED_ALIGNERS:
            logging.error("Error: Unknown aligner `{0}`".format(aligner))
            logging.error("Can be any of {0}, ".format("|".join(ALLOWED_ALIGNERS)) + \
                    "combine multiple aligners in list separated by comma")
            sys.exit()

    clean = opts.clean
    seqclean = which("seqclean")
    if clean and not seqclean:
        logging.error("Cannot find tgi seqclean in PATH")
        sys.exit()

    accn_extract = which(op.join(PASA_HOME, "misc_utilities", "accession_extractor.pl"))
    launch_pasa = which(op.join(PASA_HOME, "scripts", "Launch_PASA_pipeline.pl"))
    build_compreh_trans = which(op.join(PASA_HOME, "scripts", "build_comprehensive_transcriptome.dbi"))

    cpus = opts.cpus
    grid = opts.grid
    prepare, runfile, meta = opts.prepare, "run.sh", "run script"
    pctcov, pctid = opts.pctcov, opts.pctid
    compreh_pctcov, bpsplice = opts.compreh_pctcov, opts.bpsplice

    mkdir(pasa_db)
    os.chdir(pasa_db)

    if prepare:
        write_file(runfile, "", meta=meta)  # initialize run script

    if ggfasta:
        transcripts = FileMerger([dnfasta, ggfasta], tfasta).merge()
        accn_extract_cmd = "cat {0} | {1} > {2}".format(dnfasta, accn_extract, tdn)
        write_file(runfile, accn_extract_cmd, meta=meta, append=True) \
                if prepare else sh(accn_extract_cmd)
    else:
        transcripts = dnfasta

    if opts.grid and not opts.threaded:
        opts.threaded = opts.cpus

    prjobid = None
    if clean:
        cleancmd = "{0} {1} -c {2} -l 60".format(seqclean, transcripts, cpus)
        if prepare:
            write_file(runfile, cleancmd, meta=meta, append=True)
        else:
            prjobid = sh(cleancmd, grid=grid, grid_opts=opts)

    aafw = must_open(aaconf, "w")
    print >> aafw, alignAssembly_conf.format("{0}_pasa".format(pasa_db), pctcov, pctid, bpsplice)
    aafw.close()

    aacmd = "{0} -c {1} -C -R -g {2}".format(launch_pasa, aaconf, genome)
    aacmd += " -t {0}.clean -T -u {0} ".format(transcripts) if clean else \
             " -t {0} ".format(transcripts)
    if ggfasta:
        aacmd += " --TDN {0} ".format(tdn)
    aacmd += " --ALIGNERS {0} -I {1}".format(",".join(aligners), opts.intron)

    if prepare:
        write_file(runfile, aacmd, meta=meta, append=True)
    else:
        opts.hold_jid = prjobid
        prjobid = sh(aacmd, grid=grid, grid_opts=opts)

    if ggfasta:
        comprehcmd = "{0} -c {1} -t {2}".format(build_compreh_trans, aaconf, transcripts)
        comprehcmd += "--min_per_ID {0} --min_per_aligned {1}".format(pctid, pctcov)

        if prepare:
            write_file(runfile, comprehcmd, meta=meta, append=True)
        else:
            opts.hold_jid = prjobid
            prjobid = sh(comprehcmd, grid=grid, grid_opts=opts)


def longest(args):
    """
    %prog longest pasa.fasta output.subclusters.out

    Find the longest PASA assembly and label it as full-length. Also removes
    transcripts shorter than half the length of the longest, or shorter than
    200bp. The assemblies for the same locus is found in
    `output.subclusters.out`. In particular the lines that look like:

    sub-cluster: asmbl_25 asmbl_26 asmbl_27
    """
    from jcvi.formats.fasta import Fasta, SeqIO
    from jcvi.formats.sizes import Sizes

    p = OptionParser(longest.__doc__)
    p.add_option("--prefix", default="pasa",
                 help="Replace asmbl_ with prefix [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, subclusters = args
    prefix = fastafile.rsplit(".", 1)[0]

    idsfile = prefix + ".fl.ids"
    fw = open(idsfile, "w")
    sizes = Sizes(fastafile).mapping

    name_convert = lambda x: x.replace("asmbl", opts.prefix)

    keep = set()  # List of IDs to write
    fp = open(subclusters)
    nrecs = 0
    for row in fp:
        if not row.startswith("sub-cluster:"):
            continue
        asmbls = row.split()[1:]
        longest_asmbl = max(asmbls, key=lambda x: sizes[x])
        longest_size = sizes[longest_asmbl]
        print >> fw, name_convert(longest_asmbl)
        nrecs += 1
        cutoff = max(longest_size / 2, 200)
        keep.update(set(x for x in asmbls if sizes[x] >= cutoff))

    fw.close()
    logging.debug("{0} fl-cDNA records written to `{1}`.".format(nrecs, idsfile))

    f = Fasta(fastafile, lazy=True)
    newfastafile = prefix + ".clean.fasta"
    fw = open(newfastafile, "w")
    nrecs = 0
    for name, rec in f.iteritems_ordered():
        if name not in keep:
            continue

        rec.id = name_convert(name)
        rec.description = ""
        SeqIO.write([rec], fw, "fasta")
        nrecs += 1

    fw.close()
    logging.debug("{0} valid records written to `{1}`.".format(nrecs, newfastafile))


if __name__ == '__main__':
    main()
