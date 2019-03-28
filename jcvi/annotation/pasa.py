#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Utilities for submitting PASA jobs and processing PASA results.
"""
from __future__ import print_function

import os
import os.path as op
import sys
import logging

from jcvi.formats.base import write_file, must_open, FileMerger
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, symlink, which


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

annotation = "annotation.gff3"
tdn, flaccs = "tdn.accs", "FL_accs.txt"
tfasta, gfasta = "transcripts.fasta", "genome.fasta"
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

    fl_accs = opts.fl_accs
    cpus = opts.cpus
    grid = opts.grid
    prepare, runfile = opts.prepare, "run.sh"
    pctcov, pctid = opts.pctcov, opts.pctid
    compreh_pctid = opts.compreh_pctid
    compreh_pctcov, bpsplice = opts.compreh_pctcov, opts.bpsplice

    cmds = []

    # set PASAHOME env variable if preparing shell script
    if prepare:
        env_cmd = 'export PASAHOME="{0}"'.format(PASA_HOME)
        cmds.append(env_cmd)

    if ggfasta:
        transcripts = FileMerger([dnfasta, ggfasta], tfasta).merge()
        accn_extract_cmd = "cat {0} | {1} > {2}".format(dnfasta, accn_extract, tdn)
        cmds.append(accn_extract_cmd)
        if not prepare:
            sh(accn_extract_cmd)
    else:
        symlink(dnfasta, tfasta)
        transcripts = tfasta

    if opts.grid and not opts.threaded:
        opts.threaded = opts.cpus

    prjobid = None
    if clean:
        ccpus = 16 if cpus >= 16 else cpus
        cleancmd = "{0} {1} -c {2} -l 60".format(seqclean, transcripts, ccpus)
        if prepare:
            cmds.append(cleancmd)
        else:
            prjobid = sh(cleancmd, grid=grid, grid_opts=opts)

    aafw = must_open(aaconf, "w")
    print(alignAssembly_conf.format("{0}_pasa".format(pasa_db), \
            pctcov, pctid, bpsplice), file=aafw)
    aafw.close()

    symlink(genome, gfasta)

    aacmd = "{0} -c {1} -C -R -g {2}".format(launch_pasa, aaconf, gfasta)
    aacmd += " -t {0}.clean -T -u {0}".format(transcripts) if clean else \
             " -t {0}".format(transcripts)
    if fl_accs:
        symlink(fl_accs, flaccs)
        aacmd += " -f {0}".format(flaccs)
    if ggfasta:
        aacmd += " --TDN {0}".format(tdn)
    aacmd += " --ALIGNERS {0} -I {1} --CPU {2}".format(",".join(aligners), \
            opts.intron, cpus)

    if prepare:
        cmds.append(aacmd)
    else:
        opts.hold_jid = prjobid
        prjobid = sh(aacmd, grid=grid, grid_opts=opts)

    if opts.compreh and ggfasta:
        comprehcmd = "{0} -c {1} -t {2}".format(build_compreh_trans, aaconf, transcripts)
        comprehcmd += " --min_per_ID {0} --min_per_aligned {1}".format(compreh_pctid, compreh_pctcov)

        if prepare:
            cmds.append(comprehcmd)
        else:
            opts.hold_jid = prjobid
            prjobid = sh(comprehcmd, grid=grid, grid_opts=opts)

    if prepare:
        write_file(runfile, "\n".join(cmds))  # initialize run script


def compare(args):
    """
    %prog compare pasa_db_name [--annots_gff3=annotation.gff3]

    Run the PASA annotation comparison pipeline

    This assumes that PASA alignment assembly has alredy been completed and
    run directory contains `genome.fasta` and `transcript.fasta` files.

    If `--annots_gff3` is specified, the PASA database is loaded with the annotations
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

    if len(args) < 1:
        sys.exit(not p.print_help())

    pasa_db, = args

    PASA_HOME = opts.pasa_home
    if not op.isdir(PASA_HOME):
        logging.error("PASA_HOME={0} directory does not exist".format(PASA_HOME))
        sys.exit()

    launch_pasa = which(op.join(PASA_HOME, "scripts", \
            "Launch_PASA_pipeline.pl"))

    annots_gff3 = opts.annots_gff3
    grid = opts.grid
    prepare, runfile = opts.prepare, "run.sh"

    os.chdir(pasa_db)

    if prepare:
        write_file(runfile, "", append=True, skipcheck=True)  # initialize run script

    acfw = must_open(acconf, "w")
    print(annotCompare_conf.format("{0}_pasa".format(pasa_db), \
            opts.pctovl, opts.pct_coding, opts.pctid_prot, opts.pctlen_FL, \
            opts.pctlen_nonFL, opts.orf_size, opts.pct_aln, opts.pctovl_gene, \
            opts.stompovl, opts.trust_FL, opts.utr_exons), file=acfw)
    acfw.close()

    if not op.exists(gfasta):
        sys.exit("Genome fasta file `{0}` does not exist".format(gfasta))

    transcripts = tfasta
    if not op.exists(transcripts):
        sys.exit("Transcript fasta file `{0}` does not exist".format(transcripts))

    if op.exists("{0}.clean".format(transcripts)):
        transcripts = "{0}.clean".format(transcripts)

    accmd = "{0} -c {1} -A -g {2} -t {3} --GENETIC_CODE {4}".format(launch_pasa, \
            acconf, gfasta, transcripts, opts.genetic_code)

    if annots_gff3:
        if not op.exists(annots_gff3):
            sys.exit("Annotation gff3 file `{0}` does not exist".format(annots_gff3))
        symlink(annots_gff3, annotation)
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
        print(name_convert(longest_asmbl), file=fw)
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
    iterate through each locus (shared locus name or overlapping CDS)
    and identify same/different isoforms (shared splicing structure)
    across the input datasets.

    If `slop` is enabled, consolidation will collapse any variation
    in terminal UTR lengths, keeping the longest as representative.
    """
    from jcvi.formats.base import longest_unique_prefix
    from jcvi.formats.gff import make_index, match_subfeats
    from jcvi.utils.cbook import AutoVivification
    from jcvi.utils.grouper import Grouper
    from itertools import combinations, product

    supported_modes = ["name", "coords"]
    p = OptionParser(consolidate.__doc__)
    p.add_option("--slop", default=False, action="store_true",
            help="allow minor variation in terminal 5'/3' UTR" + \
                 " start/stop position [default: %default]")
    p.add_option("--inferUTR", default=False, action="store_true",
            help="infer presence of UTRs from exon coordinates")
    p.add_option("--mode", default="name", choices=supported_modes,
            help="method used to determine overlapping loci")
    p.add_option("--summary", default=False, action="store_true",
            help="Generate summary table of consolidation process")
    p.add_option("--clusters", default=False, action="store_true",
            help="Generate table of cluster members after consolidation")
    p.set_outfile()

    opts, args = p.parse_args(args)
    slop = opts.slop
    inferUTR = opts.inferUTR
    mode = opts.mode

    if len(args) < 2:
        sys.exit(not p.print_help())

    gffdbx = {}
    for gffile in args:
        dbn = longest_unique_prefix(gffile, args)
        gffdbx[dbn] = make_index(gffile)

    loci = Grouper()
    for dbn in gffdbx:
        odbns = [odbn for odbn in gffdbx if dbn != odbn]
        for gene in gffdbx[dbn].features_of_type('gene', order_by=('seqid', 'start')):
            if mode == "name":
                loci.join(gene.id, (gene.id, dbn))
            else:
                if (gene.id, dbn) not in loci:
                    loci.join((gene.id, dbn))
                    gene_cds = list(gffdbx[dbn].children(gene, \
                        featuretype='CDS', order_by=('start')))
                    gene_cds_start, gene_cds_stop = gene_cds[0].start, \
                        gene_cds[-1].stop
                    for odbn in odbns:
                        for ogene_cds in gffdbx[odbn].region(seqid=gene.seqid, \
                                start=gene_cds_start, end=gene_cds_stop, \
                                strand=gene.strand, featuretype='CDS'):
                            for ogene in gffdbx[odbn].parents(ogene_cds, featuretype='gene'):
                                loci.join((gene.id, dbn), (ogene.id, odbn))

    gfeats = {}
    mrna = AutoVivification()
    for i, locus in enumerate(loci):
        gene = "gene_{0:0{pad}}".format(i, pad=6) \
                if mode == "coords" else None

        for elem in locus:
            if type(elem) == tuple:
                _gene, dbn = elem
                if gene is None: gene = _gene

                g = gffdbx[dbn][_gene]
                if gene not in gfeats:
                    gfeats[gene] = g
                    gfeats[gene].attributes['ID'] = [gene]
                else:
                    if g.start < gfeats[gene].start:
                        gfeats[gene].start = g.start
                    if g.stop > gfeats[gene].stop:
                        gfeats[gene].stop = g.stop

                c = list(gffdbx[dbn].children(_gene, featuretype='mRNA', order_by='start'))
                if len(c) > 0:
                    mrna[gene][dbn] = c

    fw = must_open(opts.outfile, "w")
    print("##gff-version	3", file=fw)
    seen = {}
    if opts.summary:
        summaryfile = "{0}.summary.txt".format(opts.outfile.rsplit(".")[0])
        sfw = must_open(summaryfile, "w")
        summary = ["id"]
        summary.extend(gffdbx.keys())
        print("\t".join(str(x) for x in summary), file=sfw)
    if opts.clusters:
        clustersfile = "{0}.clusters.txt".format(opts.outfile.rsplit(".")[0])
        cfw = must_open(clustersfile, "w")
        clusters = ["id", "dbns", "members", "trlens"]
        print("\t".join(str(x) for x in clusters), file=cfw)
    for gene in mrna:
        g = Grouper()
        dbns = list(combinations(mrna[gene], 2))
        if len(dbns) > 0:
            for dbn1, dbn2 in dbns:
                dbx1, dbx2 = gffdbx[dbn1], gffdbx[dbn2]
                for mrna1, mrna2 in product(mrna[gene][dbn1], mrna[gene][dbn2]):
                    mrna1s, mrna2s = mrna1.stop - mrna1.start + 1, \
                            mrna2.stop - mrna2.start + 1
                    g.join((dbn1, mrna1.id, mrna1s))
                    g.join((dbn2, mrna2.id, mrna2s))

                    if match_subfeats(mrna1, mrna2, dbx1, dbx2, featuretype='CDS'):
                        res = []
                        ftypes = ['exon'] if inferUTR else ['five_prime_UTR', 'three_prime_UTR']
                        for ftype in ftypes:
                            res.append(match_subfeats(mrna1, mrna2, dbx1, dbx2, featuretype=ftype, slop=slop))

                        if all(r == True for r in res):
                            g.join((dbn1, mrna1.id, mrna1s), (dbn2, mrna2.id, mrna2s))
        else:
            for dbn1 in mrna[gene]:
                for mrna1 in mrna[gene][dbn1]:
                    g.join((dbn1, mrna1.id, mrna1.stop - mrna1.start + 1))

        print(gfeats[gene], file=fw)

        for group in g:
            group.sort(key=lambda x: x[2], reverse=True)
            dbs, mrnas = [el[0] for el in group], [el[1] for el in group]
            d, m = dbs[0], mrnas[0]

            dbid, _mrnaid = "|".join(str(x) for x in set(dbs)), []
            for x in mrnas:
                if x not in _mrnaid: _mrnaid.append(x)
            mrnaid = "{0}|{1}".format(dbid, "-".join(_mrnaid))
            if mrnaid not in seen:
                seen[mrnaid] = 0
            else:
                seen[mrnaid] += 1
                mrnaid = "{0}-{1}".format(mrnaid, seen[mrnaid])

            _mrna = gffdbx[d][m]
            _mrna.attributes['ID'] = [mrnaid]
            _mrna.attributes['Parent'] = [gene]
            children = gffdbx[d].children(m, order_by='start')
            print(_mrna, file=fw)
            for child in children:
                child.attributes['ID'] = ["{0}|{1}".format(dbid, child.id)]
                child.attributes['Parent'] = [mrnaid]
                print(child, file=fw)

            if opts.summary:
                summary = [mrnaid]
                summary.extend(['Y' if db in set(dbs) else 'N' for db in gffdbx])
                print("\t".join(str(x) for x in summary), file=sfw)

            if opts.clusters:
                clusters = [mrnaid]
                clusters.append(",".join(str(el[0]) for el in group))
                clusters.append(",".join(str(el[1]) for el in group))
                clusters.append(",".join(str(el[2]) for el in group))
                print("\t".join(str(x) for x in clusters), file=cfw)

    fw.close()
    if opts.summary: sfw.close()
    if opts.clusters: cfw.close()


if __name__ == '__main__':
    main()
