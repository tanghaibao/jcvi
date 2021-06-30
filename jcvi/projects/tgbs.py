#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Reference-free tGBS related functions.
"""

import logging
import os
import os.path as op
import sys

from collections import Counter
from pickle import dump, load

from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.fastq import iter_fastq
from jcvi.formats.base import must_open, write_file
from jcvi.formats.sam import get_prefix
from jcvi.apps.cdhit import deduplicate
from jcvi.apps.gmap import check_index
from jcvi.apps.grid import MakeManager
from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, iglob, mkdir


speedupsh = r"""
cd {0}

find *.native | sed 's/\..*//' | sort -u | \
    awk '{{ printf("split_by_chromosome.pl -s %s -o splitted_%s -native %s.*.native -x 5\n", \
                  $0, $0, $0); }}' > split.sh
parallel -j {1} < split.sh

find splitted_* -name "*.native" | \
    awk '{{ printf("SNP_Discovery-short.pl -native %s -o %s.SNPs_Het.txt -a 2 -ac 0.3 -c 0.8\n", \
                  $0, $0); }}' > snps.sh
parallel -j {1} < snps.sh

find splitted_*.log | \
    awk '{{ gsub("splitted_|.log", "", $0); \
           printf("combine_snps_single_file.pl -d splitted_%s -p \"*.txt\" -o %s.SNPs_Het.txt\n", \
                  $0, $0); }}' > combine.sh
parallel -j {1} < combine.sh

cd ..
"""


def main():

    actions = (
        ("snpflow", "run SNP calling pipeline from reads to allele_counts"),
        ("count", "count the number of reads in all clusters"),
        ("snpplot", "illustrate the SNP sites in CDT"),
        ("weblogo", "extract base composition for reads"),
        ("novo", "reference-free tGBS pipeline v1"),
        ("novo2", "reference-free tGBS pipeline v2"),
        ("mstmap", "convert LMDs to MSTMAP input"),
        ("query", "random access to loci file"),
        ("synteny", "plot mst map against reference genome"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def build_index(locifile):
    idxfile = locifile + ".idx"
    if need_update(locifile, idxfile):
        fp = open(locifile)
        fw = open(idxfile, "w")
        idx = {}
        while True:
            pos = fp.tell()
            line = fp.readline()
            if not line:
                break
            if not line.startswith("//"):
                continue
            tag, contig = line.split()[:2]
            idx[contig] = pos
        dump(idx, fw)
        fw.close()
        return idx

    idx = load(open(idxfile))
    return idx


def query(args):
    """
    %prog query out.loci contig

    Random access to loci file. This script helps speeding up debugging.
    """
    p = OptionParser(query.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    locifile, contig = args
    idx = build_index(locifile)
    pos = idx[contig]
    logging.debug("Contig {0} found at pos {1}".format(contig, pos))
    fp = open(locifile)
    fp.seek(pos)
    section = []
    while True:
        row = fp.readline()
        if row.startswith("//") and row.split()[1] != contig:
            break
        section.append(row)
    print("".join(section))


def synteny(args):
    """
    %prog synteny mstmap.out novo.final.fasta reference.fasta

    Plot MSTmap against reference genome.
    """
    from jcvi.assembly.geneticmap import bed as geneticmap_bed
    from jcvi.apps.align import blat
    from jcvi.formats.blast import bed as blast_bed, best

    p = OptionParser(synteny.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    mstmapout, novo, ref = args
    pf = mstmapout.split(".")[0]
    rf = ref.split(".")[0]
    mstmapbed = geneticmap_bed([mstmapout])
    cmd = "cut -d. -f1 {0}".format(mstmapbed)
    tmpbed = mstmapbed + ".tmp"
    sh(cmd, outfile=tmpbed)
    os.rename(tmpbed, pf + ".bed")

    cmd = "cut -f4 {0} | cut -d. -f1 | sort -u".format(mstmapbed)
    idsfile = pf + ".ids"
    sh(cmd, outfile=idsfile)
    fastafile = pf + ".fasta"
    cmd = "faSomeRecords {0} {1} {2}".format(novo, idsfile, fastafile)
    sh(cmd)
    blastfile = blat([ref, fastafile])
    bestblastfile = best([blastfile])
    blastbed = blast_bed([bestblastfile])
    os.rename(blastbed, rf + ".bed")

    anchorsfile = "{0}.{1}.anchors".format(pf, rf)
    cmd = "paste {0} {0}".format(idsfile)
    sh(cmd, outfile=anchorsfile)


def mstmap(args):
    """
    %prog mstmap LMD50.snps.genotype.txt

    Convert LMDs to MSTMAP input.
    """
    from jcvi.assembly.geneticmap import MSTMatrix

    p = OptionParser(mstmap.__doc__)
    p.add_option(
        "--population_type",
        default="RIL6",
        help="Type of population, possible values are DH and RILd",
    )
    p.add_option(
        "--missing_threshold",
        default=0.5,
        help="Missing threshold, .25 excludes any marker with >25% missing",
    )
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (lmd,) = args
    fp = open(lmd)
    next(fp)  # Header
    table = {"0": "-", "1": "A", "2": "B", "3": "X"}
    mh = ["locus_name"] + next(fp).split()[4:]
    genotypes = []
    for row in fp:
        atoms = row.split()
        chr, pos, ref, alt = atoms[:4]
        locus_name = ".".join((chr, pos))
        codes = [table[x] for x in atoms[4:]]
        genotypes.append([locus_name] + codes)

    mm = MSTMatrix(genotypes, mh, opts.population_type, opts.missing_threshold)
    mm.write(opts.outfile, header=True)


def weblogo(args):
    """
    %prog weblogo [fastafile|fastqfile]

    Extract base composition for reads
    """
    import numpy as np
    from rich.progress import Progress

    p = OptionParser(weblogo.__doc__)
    p.add_option("-N", default=10, type="int", help="Count the first and last N bases")
    p.add_option("--nreads", default=1000000, type="int", help="Parse first N reads")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastqfile,) = args
    N = opts.N
    nreads = opts.nreads

    pat = "ATCG"
    L = np.zeros((4, N), dtype="int32")
    R = np.zeros((4, N), dtype="int32")
    p = dict((a, i) for (i, a) in enumerate(pat))
    L4, R3 = Counter(), Counter()

    k = 0
    fw_L = open("L.fasta", "w")
    fw_R = open("R.fasta", "w")
    fastq = fastqfile.endswith(".fastq")
    it = iter_fastq(fastqfile) if fastq else SeqIO.parse(must_open(fastqfile), "fasta")

    with Progress() as progress:
        progress.add_task("[green] Processing ...", start=False, total=nreads)
        for rec in it:
            k += 1
            if k > nreads:
                break
            if rec is None:
                break
            s = str(rec.seq)
            for i, a in enumerate(s[:N]):
                if a in p:
                    a = p[a]
                    L[a][i] += 1
            for j, a in enumerate(s[-N:][::-1]):
                if a in p:
                    a = p[a]
                    R[a][N - 1 - j] += 1
            l4, r3 = s[:4], s[-3:]
            L4[l4] += 1
            R3[r3] += 1
            print(">{0}\n{1}".format(k, s[:N]), file=fw_L)
            print(">{0}\n{1}".format(k, s[-N:]), file=fw_R)

    fw_L.close()
    fw_R.close()

    cmd = "weblogo -F png -s large -f {0}.fasta -o {0}.png"
    cmd += " --color-scheme classic --composition none -U probability"
    cmd += " --title {1}"
    sh(cmd.format("L", "First_10_bases"))
    sh(cmd.format("R", "Last_10_bases"))

    np.savetxt("L.{0}.csv".format(pat), L, delimiter=",", fmt="%d")
    np.savetxt("R.{0}.csv".format(pat), R, delimiter=",", fmt="%d")

    fw = open("L4.common", "w")
    for p, c in L4.most_common(N):
        print("\t".join((p, str(c))), file=fw)
    fw.close()

    fw = open("R3.common", "w")
    for p, c in R3.most_common(N):
        print("\t".join((p, str(c))), file=fw)
    fw.close()


def count(args):
    """
    %prog count cdhit.consensus.fasta

    Scan the headers for the consensus clusters and count the number of reads.
    """
    from jcvi.graphics.histogram import stem_leaf_plot
    from jcvi.utils.cbook import SummaryStats

    p = OptionParser(count.__doc__)
    p.add_option("--csv", help="Write depth per contig to file")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastafile,) = args
    csv = open(opts.csv, "w") if opts.csv else None

    f = Fasta(fastafile, lazy=True)
    sizes = []
    for desc, rec in f.iterdescriptions_ordered():
        if desc.startswith("singleton"):
            sizes.append(1)
            continue

        # consensus_for_cluster_0 with 63 sequences
        if "with" in desc:
            name, w, size, seqs = desc.split()
            if csv:
                print("\t".join(str(x) for x in (name, size, len(rec))), file=csv)
            assert w == "with"
            sizes.append(int(size))
        # MRD85:00603:02472;size=167;
        else:
            name, size, tail = desc.split(";")
            sizes.append(int(size.replace("size=", "")))

    if csv:
        csv.close()
        logging.debug("File written to `%s`.", opts.csv)

    s = SummaryStats(sizes)
    print(s, file=sys.stderr)
    stem_leaf_plot(s.data, 0, 100, 20, title="Cluster size")


def novo(args):
    """
    %prog novo reads.fastq

    Reference-free tGBS pipeline v1.
    """
    from jcvi.assembly.kmer import jellyfish, histogram
    from jcvi.assembly.preprocess import diginorm
    from jcvi.formats.fasta import filter as fasta_filter, format
    from jcvi.apps.cdhit import filter as cdhit_filter

    p = OptionParser(novo.__doc__)
    p.add_option(
        "--technology",
        choices=("illumina", "454", "iontorrent"),
        default="iontorrent",
        help="Sequencing platform",
    )
    p.set_depth(depth=50)
    p.set_align(pctid=96)
    p.set_home("cdhit", default="/usr/local/bin/")
    p.set_home("fiona", default="/usr/local/bin/")
    p.set_home("jellyfish", default="/usr/local/bin/")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (fastqfile,) = args
    cpus = opts.cpus
    depth = opts.depth
    pf, sf = fastqfile.rsplit(".", 1)

    diginormfile = pf + ".diginorm." + sf
    if need_update(fastqfile, diginormfile):
        diginorm([fastqfile, "--single", "--depth={0}".format(depth)])
        keepabund = fastqfile + ".keep.abundfilt"
        sh("cp -s {0} {1}".format(keepabund, diginormfile))

    jf = pf + "-K23.histogram"
    if need_update(diginormfile, jf):
        jellyfish(
            [
                diginormfile,
                "--prefix={0}".format(pf),
                "--cpus={0}".format(cpus),
                "--jellyfish_home={0}".format(opts.jellyfish_home),
            ]
        )

    genomesize = histogram([jf, pf, "23"])
    fiona = pf + ".fiona.fa"
    if need_update(diginormfile, fiona):
        cmd = op.join(opts.fiona_home, "fiona")
        cmd += " -g {0} -nt {1} --sequencing-technology {2}".format(
            genomesize, cpus, opts.technology
        )
        cmd += " -vv {0} {1}".format(diginormfile, fiona)
        logfile = pf + ".fiona.log"
        sh(cmd, outfile=logfile, errfile=logfile)

    dedup = "cdhit"
    pctid = opts.pctid
    cons = fiona + ".P{0}.{1}.consensus.fasta".format(pctid, dedup)
    if need_update(fiona, cons):
        deduplicate(
            [
                fiona,
                "--consensus",
                "--reads",
                "--pctid={0}".format(pctid),
                "--cdhit_home={0}".format(opts.cdhit_home),
            ]
        )

    filteredfile = pf + ".filtered.fasta"
    if need_update(cons, filteredfile):
        covfile = pf + ".cov.fasta"
        cdhit_filter(
            [cons, "--outfile={0}".format(covfile), "--minsize={0}".format(depth / 5)]
        )
        fasta_filter([covfile, "50", "--outfile={0}".format(filteredfile)])

    finalfile = pf + ".final.fasta"
    if need_update(filteredfile, finalfile):
        format(
            [
                filteredfile,
                finalfile,
                "--sequential=replace",
                "--prefix={0}_".format(pf),
            ]
        )


def scan_read_files(trimmed, patterns):
    reads = iglob(trimmed, patterns)
    samples = sorted(set(op.basename(x).split(".")[0] for x in reads))
    logging.debug(
        "Total {0} read files from {1} samples".format(len(reads), len(samples))
    )
    return reads, samples


def novo2(args):
    """
    %prog novo2 trimmed projectname

    Reference-free tGBS pipeline v2.
    """
    p = OptionParser(novo2.__doc__)
    p.set_fastq_names()
    p.set_align(pctid=95)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    trimmed, pf = args
    pctid = opts.pctid
    reads, samples = scan_read_files(trimmed, opts.names)

    # Set up directory structure
    clustdir = "uclust"
    acdir = "allele_counts"
    for d in (clustdir, acdir):
        mkdir(d)

    mm = MakeManager()
    clustfiles = []
    # Step 0 - clustering within sample
    for s in samples:
        flist = [x for x in reads if op.basename(x).split(".")[0] == s]
        outfile = s + ".P{0}.clustS".format(pctid)
        outfile = op.join(clustdir, outfile)
        cmd = "python -m jcvi.apps.uclust cluster --cpus=8"
        cmd += " {0} {1}".format(s, " ".join(flist))
        cmd += " --outdir={0}".format(clustdir)
        cmd += " --pctid={0}".format(pctid)
        mm.add(flist, outfile, cmd)
        clustfiles.append(outfile)

    # Step 1 - make consensus within sample
    allcons = []
    for s, clustfile in zip(samples, clustfiles):
        outfile = s + ".P{0}.consensus".format(pctid)
        outfile = op.join(clustdir, outfile)
        cmd = "python -m jcvi.apps.uclust consensus"
        cmd += " {0}".format(clustfile)
        mm.add(clustfile, outfile, cmd)
        allcons.append(outfile)

    # Step 2 - clustering across samples
    clustSfile = pf + ".P{0}.clustS".format(pctid)
    cmd = "python -m jcvi.apps.uclust mcluster {0}".format(" ".join(allcons))
    cmd += " --prefix={0}".format(pf)
    mm.add(allcons, clustSfile, cmd)

    # Step 3 - make consensus across samples
    locifile = pf + ".P{0}.loci".format(pctid)
    cmd = "python -m jcvi.apps.uclust mconsensus {0}".format(" ".join(allcons))
    cmd += " --prefix={0}".format(pf)
    mm.add(allcons + [clustSfile], locifile, cmd)

    mm.write()


def snpflow(args):
    """
    %prog snpflow trimmed reference.fasta

    Run SNP calling pipeline until allele_counts are generated. This includes
    generation of native files, SNP_Het file. Speedup for fragmented genomes
    are also supported.
    """
    p = OptionParser(snpflow.__doc__)
    p.set_fastq_names()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    trimmed, ref = args
    nseqs = len(Fasta(ref))
    supercat = nseqs >= 1000
    if supercat:
        logging.debug("Total seqs in ref: {0} (supercat={1})".format(nseqs, supercat))

    reads, samples = scan_read_files(trimmed, opts.names)

    # Set up directory structure
    nativedir, countsdir = "native", "allele_counts"
    for d in (nativedir, countsdir):
        mkdir(d)

    mm = MakeManager()
    # Step 0 - index database
    db = op.join(*check_index(ref, supercat=supercat, go=False))
    cmd = "python -m jcvi.apps.gmap index {0}".format(ref)
    if supercat:
        cmd += " --supercat"
        coordsfile = db + ".coords"
        supercatfile = ref.rsplit(".", 1)[0] + ".supercat.fasta"
        mm.add(ref, (db, coordsfile), cmd)
    else:
        mm.add(ref, db, cmd)

    # Step 1 - GSNAP alignment and conversion to native file
    allnatives = []
    allsamstats = []
    gmapdb = supercatfile if supercat else ref
    for f in reads:
        prefix = get_prefix(f, ref)
        gsnapfile = op.join(nativedir, prefix + ".gsnap")
        nativefile = op.join(nativedir, prefix + ".unique.native")
        samstatsfile = op.join(nativedir, prefix + ".unique.sam.stats")
        cmd = "python -m jcvi.apps.gmap align {0} {1}".format(gmapdb, f)
        cmd += " --outdir={0} --native --cpus=1".format(nativedir)
        mm.add((f, db), nativefile, cmd)

        cmd = "python -m jcvi.apps.gmap bam {0} {1} --cpus=1".format(gsnapfile, gmapdb)
        mm.add(nativefile, samstatsfile, cmd)
        allnatives.append(nativefile)
        allsamstats.append(samstatsfile)

    # Step 2 - call SNP discovery
    if supercat:
        nativeconverted = nativedir + "-converted"
        mkdir(nativeconverted)
        allnativesc = [op.join(nativeconverted, op.basename(x)) for x in allnatives]
        cmd = "tGBS-Convert_Pseudo_Genome_NATIVE_Coordinates.pl"
        cmd += " -i {0}/*.native -o {1}".format(nativedir, nativeconverted)
        cmd += " -c {0}".format(coordsfile)
        cmds = ["rm -rf {0}".format(nativeconverted), cmd]
        mm.add(allnatives + [coordsfile], allnativesc, cmds)

        runfile = "speedup.sh"
        write_file(runfile, speedupsh.format(nativeconverted, opts.cpus))
        nativedir = nativeconverted
        allsnps = [op.join(nativedir, "{0}.SNPs_Het.txt".format(x)) for x in samples]
        mm.add(allnativesc, allsnps, "./{0}".format(runfile))
    else:
        for s in samples:
            snpfile = op.join(nativedir, "{0}.SNPs_Het.txt".format(s))
            cmd = "SNP_Discovery-short.pl"
            cmd += " -native {0}/{1}.*unique.native".format(nativedir, s)
            cmd += " -o {0} -a 2 -ac 0.3 -c 0.8".format(snpfile)
            flist = [x for x in allnatives if op.basename(x).split(".")[0] == s]
            mm.add(flist, snpfile, cmd)

    # Step 3 - generate equal file
    allsnps = [op.join(nativedir, "{0}.SNPs_Het.txt".format(x)) for x in samples]
    for s in samples:
        equalfile = op.join(nativedir, "{0}.equal".format(s))
        cmd = "extract_reference_alleles.pl"
        cmd += " --native {0}/{1}.*unique.native".format(nativedir, s)
        cmd += " --genotype {0}/{1}.SNPs_Het.txt".format(nativedir, s)
        cmd += " --allgenotypes {0}/*.SNPs_Het.txt".format(nativedir)
        cmd += " --fasta {0} --output {1}".format(ref, equalfile)
        mm.add(allsnps, equalfile, cmd)

    # Step 4 - generate snp matrix
    allequals = [op.join(nativedir, "{0}.equal".format(x)) for x in samples]
    matrix = "snps.matrix.txt"
    cmd = "generate_matrix.pl"
    cmd += " --tables {0}/*SNPs_Het.txt --equal {0}/*equal".format(nativedir)
    cmd += " --fasta {0} --output {1}".format(ref, matrix)
    mm.add(allsnps + allequals, matrix, cmd)

    # Step 5 - generate allele counts
    allcounts = []
    for s in samples:
        allele_counts = op.join(countsdir, "{0}.SNPs_Het.allele_counts".format(s))
        cmd = "count_reads_per_allele.pl -m snps.matrix.txt"
        cmd += " -s {0} --native {1}/{0}.*unique.native".format(s, nativedir)
        cmd += " -o {0}".format(allele_counts)
        mm.add(matrix, allele_counts, cmd)
        allcounts.append(allele_counts)

    # Step 6 - generate raw snps
    rawsnps = "Genotyping.H3.txt"
    cmd = "/home/shared/scripts/delin/SamplesGenotyping.pl --homo 3"
    cmd += " -pf allele_counts -f {0} --outfile {1}".format(countsdir, rawsnps)
    cmds = ["rm -f {0}".format(rawsnps), cmd]
    mm.add(allcounts, rawsnps, cmds)

    # Step 7 - generate alignment report
    sam_summary = "sam.summary"
    cmd = "/home/shared/scripts/eddyyeh/alignment_stats.pl"
    cmd += " -f {0} -o {1}".format(" ".join(allsamstats), sam_summary)
    mm.add(allsamstats, sam_summary, cmd)

    native_summary = "native.summary"
    cmd = "/home/shared/scripts/eddyyeh/alignment_stats.pl"
    cmd += " -n {0} -o {1}".format(" ".join(allnatives), native_summary)
    mm.add(allnatives, native_summary, cmd)

    mm.write()


def snpplot(args):
    """
    %prog counts.cdt

    Illustrate the histogram per SNP site.
    """
    p = OptionParser(snpplot.__doc__)
    opts, args, iopts = p.set_image_options(args, format="png")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (datafile,) = args
    # Read in CDT file
    fp = open(datafile)
    next(fp)
    next(fp)
    data = []
    for row in fp:
        atoms = row.split()[4:]
        nval = len(atoms)
        values = [float(x) for x in atoms]
        # normalize
        values = [x * 1.0 / sum(values) for x in values]
        data.append(values)

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    xmin, xmax = 0.1, 0.9
    ymin, ymax = 0.1, 0.9
    yinterval = (ymax - ymin) / len(data)
    colors = "rbg" if nval == 3 else ["lightgray"] + list("rbg")
    ystart = ymax
    for d in data:
        xstart = xmin
        for dd, c in zip(d, colors):
            xend = xstart + (xmax - xmin) * dd
            root.plot((xstart, xend), (ystart, ystart), "-", color=c)
            xstart = xend
        ystart -= yinterval

    root.text(
        0.05,
        0.5,
        "{0} LMD50 SNPs".format(len(data)),
        ha="center",
        va="center",
        rotation=90,
        color="lightslategray",
    )

    for x, t, c in zip((0.3, 0.5, 0.7), ("REF", "ALT", "HET"), "rbg"):
        root.text(x, 0.95, t, color=c, ha="center", va="center")
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
