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

#script validate_alignments_in_db.dbi
validate_alignments_in_db.dbi:--MIN_PERCENT_ALIGNED={1}
validate_alignments_in_db.dbi:--MIN_AVG_PER_ID={2}
validate_alignments_in_db.dbi:--NUM_BP_PERFECT_SPLICE_BOUNDARY={3}

#script subcluster_builder.dbi
subcluster_builder.dbi:-m=50
"""

annotCompare_conf = """
# MySQL settings
MYSQLDB={0}

#script cDNA_annotation_comparer.dbi
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP={1}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_PROT_CODING={2}
cDNA_annotation_comparer.dbi:--MIN_PERID_PROT_COMPARE={3}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_FL_COMPARE={4}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_LENGTH_NONFL_COMPARE={5}
cDNA_annotation_comparer.dbi:--MIN_FL_ORF_SIZE={6}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_ALIGN_LENGTH={7}
cDNA_annotation_comparer.dbi:--MIN_PERCENT_OVERLAP_GENE_REPLACE={8}
cDNA_annotation_comparer.dbi:--STOMP_HIGH_PERCENTAGE_OVERLAPPING_GENE={9}
cDNA_annotation_comparer.dbi:--TRUST_FL_STATUS={10}
cDNA_annotation_comparer.dbi:--MAX_UTR_EXONS={11}
"""

tdn, tfasta = "tdn.accs", "transcripts.fasta"
aaconf, acconf = "alignAssembly.conf", "annotCompare.conf"
ALLOWED_ALIGNERS = ("blat", "gmap")


def main():

    actions = (
        ('assemble', 'run pasa alignment assembly pipeline'),
        ('compare', 'run pasa annotation comparison pipeline'),
        ('longest', 'label longest transcript per gene as full-length'),
        ('consolidate', 'generate consolidated annotation set from 2 or more annot compare results')
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def assemble(args):
    """
    %prog assemble pasa_db_name genome.fasta transcripts-dn.fasta [transcript-gg.fasta]

    Run the PASA alignment assembly pipeline

    If two transcript fasta files (Trinity denovo and genome guided) are provided
    and the `--compreh` param is enabled, the PASA Comprehensive Transcriptome DB
    protocol is followed <http://pasa.sourceforge.net/#A_ComprehensiveTranscriptome>

    Using the `--prepare` option creates a shell script with the run commands without
    executing the pipeline
    """
    p = OptionParser(assemble.__doc__)
    p.set_pasa_opts()
    p.add_option("--prepare", default=False, action="store_true",
            help="Prepare PASA run script with commands [default: %default]")
    p.set_grid()
    p.set_grid_opts()
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
    seqclean = op.join(opts.tgi_home, "seqclean")

    accn_extract = which(op.join(PASA_HOME, "misc_utilities", \
            "accession_extractor.pl"))
    launch_pasa = which(op.join(PASA_HOME, "scripts", \
            "Launch_PASA_pipeline.pl"))
    build_compreh_trans = which(op.join(PASA_HOME, "scripts", \
            "build_comprehensive_transcriptome.dbi"))

    cpus = opts.cpus
    grid = opts.grid
    prepare, runfile = opts.prepare, "run.sh"
    pctcov, pctid = opts.pctcov, opts.pctid
    compreh_pctcov, bpsplice = opts.compreh_pctcov, opts.bpsplice

    mkdir(pasa_db)
    os.chdir(pasa_db)

    if prepare:
        write_file(runfile, "")  # initialize run script

    if ggfasta:
        transcripts = FileMerger([dnfasta, ggfasta], tfasta).merge()
        accn_extract_cmd = "cat {0} | {1} > {2}".format(dnfasta, accn_extract, tdn)
        write_file(runfile, accn_extract_cmd, append=True) \
                if prepare else sh(accn_extract_cmd)
    else:
        transcripts = dnfasta

    if opts.grid and not opts.threaded:
        opts.threaded = opts.cpus

    prjobid = None
    if clean:
        cleancmd = "{0} {1} -c {2} -l 60".format(seqclean, transcripts, cpus)
        if prepare:
            write_file(runfile, cleancmd, append=True)
        else:
            prjobid = sh(cleancmd, grid=grid, grid_opts=opts)

    aafw = must_open(aaconf, "w")
    print >> aafw, alignAssembly_conf.format("{0}_pasa".format(pasa_db), \
            pctcov, pctid, bpsplice)
    aafw.close()

    aacmd = "{0} -c {1} -C -R -g {2}".format(launch_pasa, aaconf, genome)
    aacmd += " -t {0}.clean -T -u {0} ".format(transcripts) if clean else \
             " -t {0} ".format(transcripts)
    if ggfasta:
        aacmd += " --TDN {0} ".format(tdn)
    aacmd += " --ALIGNERS {0} -I {1} --CPU {2}".format(",".join(aligners), \
            opts.intron, cpus)

    if prepare:
        write_file(runfile, aacmd, append=True)
    else:
        opts.hold_jid = prjobid
        prjobid = sh(aacmd, grid=grid, grid_opts=opts)

    if opts.compreh and ggfasta:
        comprehcmd = "{0} -c {1} -t {2}".format(build_compreh_trans, aaconf, transcripts)
        comprehcmd += " --min_per_ID {0} --min_per_aligned {1}".format(compreh_pctid, compreh_pctcov)

        if prepare:
            write_file(runfile, comprehcmd, append=True)
        else:
            opts.hold_jid = prjobid
            prjobid = sh(comprehcmd, grid=grid, grid_opts=opts)


def compare(args):
    """
    %prog compare pasa_db_name genome.fasta transcripts.fasta [annotation.gff]

    Run the PASA annotation comparison pipeline

    If annotation.gff file is provided, the PASA database is loaded with the annotations
    first before starting annotation comparison. Otherwise, it uses previously
    loaded annotation data.

    Using the `--prepare` option creates a shell script with the run commands without
    executing the pipeline
    """
    p = OptionParser(compare.__doc__)
    p.set_pasa_opts(action="compare")
    p.add_option("--prepare", default=False, action="store_true",
            help="Prepare PASA run script with commands [default: %default]")
    p.set_grid()
    p.set_grid_opts()
    opts, args = p.parse_args(args)

    if len(args) not in (3, 4):
        sys.exit(not p.print_help())

    pasa_db, genome, transcripts, = args[:3]
    annotation = args[3] if len(args) == 4 else None

    PASA_HOME = opts.pasa_home
    if not op.isdir(PASA_HOME):
        logging.error("PASA_HOME={0} directory does not exist".format(PASA_HOME))
        sys.exit()

    launch_pasa = which(op.join(PASA_HOME, "scripts", \
            "Launch_PASA_pipeline.pl"))

    grid = opts.grid
    prepare, runfile = opts.prepare, "run.sh"

    os.chdir(pasa_db)

    if prepare:
        write_file(runfile, "")  # initialize run script

    if opts.grid and not opts.threaded:
        opts.threaded = opts.cpus

    acfw = must_open(acconf, "w")
    print >> acfw, annotCompare_conf.format("{0}_pasa".format(pasa_db), \
            opts.pctovl, opts.pct_coding, opts.pctid_prot, opts.pctlen_FL, \
            opts.pctlen_nonFL, opts.orf_size, opts.pct_aln, opts.pctovl_gene, \
            opts.stompovl, opts.trust_FL, opts.utr_exons)
    acfw.close()

    if op.exists("{0}.clean".format(transcripts)):
        transcripts = "{0}.clean".format(transcripts)

    accmd = "{0} -c {1} -A -g {2} -t {3} --GENETIC_CODE {4}".format(launch_pasa, \
            acconf, genome, transcripts, opts.genetic_code)
    if annotation:
        accmd += " -L --annots_gff3 {0}".format(annotation)
    if prepare:
        write_file(runfile, accmd, append=True)
    else:
        sh(accmd, grid=grid, grid_opts=opts)


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


def consolidate(args):
    """
    %prog consolidate gffile1 gffile2 ... > consolidated.out

    Given 2 or more gff files generated by pasa annotation comparison,
    iterate through every gene locus and identify all cases of same and
    different isoforms across the different input datasets.
    """
    from jcvi.formats.gff import make_index
    from jcvi.utils.cbook import AutoVivification
    from jcvi.utils.grouper import Grouper
    from itertools import combinations, product

    p = OptionParser(consolidate.__doc__)
    p.add_option("--slop", default=False, action="store_true",
            help="allow minor variation in terminal 5'/3' UTR" + \
                 " start/stop position [default: %default]")
    p.set_outfile()

    opts, args = p.parse_args(args)
    slop = opts.slop

    if len(args) < 2:
        sys.exit(not p.print_help())

    gffdbx = {}
    mrna = AutoVivification()
    for gffile in args:
        dbn = gffile.rsplit(".", 1)[0]
        gffdbx[dbn] = make_index(gffile)
        for gene in gffdbx[dbn].features_of_type('gene', order_by=('seqid', 'start')):
            c = list(gffdbx[dbn].children(gene, featuretype='mRNA', order_by='start'))
            if len(c) > 0:
                mrna[gene][dbn] = c

    fw = must_open(opts.outfile, "w")
    outln = ["id"]
    outln.extend(gffdbx.keys())
    print >> fw, "\t".join(str(x) for x in outln)
    for gene in mrna:
        g = Grouper()
        dbns = list(combinations(mrna[gene], 2))
        if len(dbns) > 0:
            for dbn1, dbn2 in dbns:
                for mrna1, mrna2 in product(mrna[gene][dbn1], mrna[gene][dbn2]):
                    g.join((dbn1, mrna1.id))
                    g.join((dbn2, mrna2.id))

                    fUTR, tUTR = None, None
                    if match_subfeats(mrna1, mrna2, gffdbx[dbn1], gffdbx[dbn2]):
                        fUTR = match_subfeats(mrna1, mrna2, gffdbx[dbn1], gffdbx[dbn2], \
                                featuretype='five_prime_UTR', slop=slop)
                        tUTR = match_subfeats(mrna1, mrna2, gffdbx[dbn1], gffdbx[dbn2], \
                                featuretype='three_prime_UTR', slop=slop)

                    if fUTR and tUTR:
                        g.join((dbn1, mrna1.id), (dbn2, mrna2.id))
        else:
            for dbn1 in mrna[gene]:
                for mrna1 in mrna[gene][dbn1]:
                    g.join((dbn1, mrna1.id))

        print >> sys.stderr, list(g)
        for group in g:
            outln = []
            outln.append("_".join(str(m) for m in set([elem[1] for elem in group])))
            dbn = set([elem[0] for elem in group])
            for db in gffdbx:
                d = 'Y' if db in dbn else 'N'
                outln.append(d)
            print >> fw, "\t".join(str(x) for x in outln)

    fw.close()


def match_subfeats(f1, f2, dbx1, dbx2, featuretype='CDS', slop=False):
    from jcvi.formats.gff import match_span, match_nchildren, \
            match_Nth_child

    f1c, f2c = list(dbx1.children(f1, featuretype=featuretype)), \
            list(dbx2.children(f2, featuretype=featuretype))

    if len(f1c) > 0 and len(f2c) > 0:
        if match_nchildren(f1c, f2c):
            if featuretype.endswith('UTR'):
                n = 1 if featuretype.startswith('five_prime') else len(f1c)
                if match_Nth_child(f1c, f2c, N=n, slop=slop):
                    del f1c[n-1], f2c[n-1]
                else:
                    return False

            for cf1, cf2 in zip(f1c, f2c):
                if not match_span(cf1, cf2):
                    return False

        else:
            return False

    return True


if __name__ == '__main__':
    main()
