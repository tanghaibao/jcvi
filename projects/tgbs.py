#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Reference-free tGBS related functions.
"""

import os.path as op
import logging
import sys

from collections import defaultdict
from itertools import combinations

from jcvi.formats.fasta import Fasta, SeqIO
from jcvi.formats.fastq import iter_fastq
from jcvi.formats.base import must_open, write_file
from jcvi.formats.bed import Bed, mergeBed
from jcvi.formats.sam import get_prefix
from jcvi.utils.counter import Counter
from jcvi.apps.cdhit import uclust, deduplicate
from jcvi.apps.gmap import check_index
from jcvi.apps.grid import MakeManager
from jcvi.graphics.base import plt, savefig, normalize_axes
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh, \
            iglob, mkdir


class HaplotypeResolver (object):

    def __init__(self, haplotype_set, maf=.1):
        self.haplotype_set = haplotype_set
        self.nind = len(haplotype_set)
        self.notmissing = sum(1 for x in haplotype_set if x)
        counter = Counter()
        for haplotypes in haplotype_set:
            counter.update(Counter(haplotypes))
        self.counter = {}
        for h, c in counter.items():
            if c >= self.notmissing * maf:
                self.counter[h] = c

    def __str__(self):
        return "N={0} M={1} C={2}".format(len(self.counter), \
                                self.notmissing, self.counter)

    def solve(self, fw):
        haplotype_counts = self.counter.items()
        for (a, ai), (b, bi) in combinations(haplotype_counts, 2):
            abi = sum(1 for haplotypes in self.haplotype_set \
                if a in haplotypes and b in haplotypes)
            pct = max(abi * 100 / ai, abi * 100 / bi)
            print >> fw, a, b, "A={0}".format(ai), "B={0}".format(bi), \
                               "AB={0}".format(abi), "{0}%".format(pct), \
                               "compatible" if pct < 50 else ""
            fw.flush()


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
        ('snpflow', 'run SNP calling pipeline from reads to allele_counts'),
        ('novo', 'reference-free tGBS pipeline'),
        ('resolve', 'separate repeats on collapsed contigs'),
        ('count', 'count the number of reads in all clusters'),
        ('snpplot', 'illustrate the SNP sites in CDT'),
        ('track', 'track and contrast read mapping in two bam files'),
        ('weblogo', 'extract base composition for reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def weblogo(args):
    """
    %prog weblogo [fastafile|fastqfile]

    Extract base composition for reads
    """
    import numpy as np
    from jcvi.utils.progressbar import ProgressBar, Percentage, Bar, ETA

    p = OptionParser(weblogo.__doc__)
    p.add_option("-N", default=10, type="int",
                 help="Count the first and last N bases")
    p.add_option("--nreads", default=1000000, type="int",
                 help="Parse first N reads")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
    N = opts.N
    nreads = opts.nreads

    pat = "ATCG"
    L = np.zeros((4, N), dtype="int32")
    R = np.zeros((4, N), dtype="int32")
    p = dict((a, i) for (i, a) in enumerate(pat))
    L4, R3 = Counter(), Counter()
    widgets = ['Parse reads: ', Percentage(), ' ',
               Bar(marker='>', left='[', right=']'), ' ', ETA()]
    pr = ProgressBar(maxval=nreads, term_width=60, widgets=widgets).start()

    k = 0
    fw_L = open("L.fasta", "w")
    fw_R = open("R.fasta", "w")
    fastq = fastqfile.endswith(".fastq")
    it = iter_fastq(fastqfile) if fastq else \
           SeqIO.parse(must_open(fastqfile), "fasta")
    for rec in it:
        k += 1
        if k % 1000 == 0:
            pr.update(k)
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
        print >> fw_L, ">{0}\n{1}".format(k, s[:N])
        print >> fw_R, ">{0}\n{1}".format(k, s[-N:])

    fw_L.close()
    fw_R.close()

    cmd = "weblogo -F png -s large -f {0}.fasta -o {0}.png"
    cmd += " --color-scheme classic --composition none -U probability"
    cmd += " --title {1}"
    sh(cmd.format('L', "First_10_bases"))
    sh(cmd.format('R', "Last_10_bases"))

    np.savetxt("L.{0}.csv".format(pat), L, delimiter=',', fmt="%d")
    np.savetxt("R.{0}.csv".format(pat), R, delimiter=',', fmt="%d")

    fw = open("L4.common", "w")
    for p, c in L4.most_common(N):
        print >> fw, "\t".join((p, str(c)))
    fw.close()

    fw = open("R3.common", "w")
    for p, c in R3.most_common(N):
        print >> fw, "\t".join((p, str(c)))
    fw.close()


def bed_store(bedfile, sorted=False):
    bedfile = mergeBed(bedfile, s=True, nms=True, sorted=sorted)
    bed = Bed(bedfile)
    reads, reads_r = {}, defaultdict(list)
    for b in bed:
        target = "{0}:{1}".format(b.seqid, b.start)
        for accn in b.accn.split(","):
            reads[accn] = target
            reads_r[target].append(accn)
    return reads, reads_r


def contrast_stores(bed1_store_r, bed2_store, minreads=10, minpct=.1, prefix="AB"):
    for target, reads in bed1_store_r.iteritems():
        nreads = len(reads)
        if nreads < minreads:
            continue
        good_mapping = max(minreads / 2, minpct * nreads)
        bed2_targets = Counter(bed2_store.get(r) for r in reads)
        c = dict((k, v) for (k, v) in bed2_targets.items() if v >= good_mapping)
        ctag = "|".join("{0}({1})".format(k, v) for (k, v) in c.items())
        print prefix, target, nreads, ctag, len(set(c.keys()) - set([None]))


def track(args):
    """
    %prog track bed1 bed2

    Track and contrast read mapping in two bam files.
    """
    p = OptionParser(track.__doc__)
    p.add_option("--sorted", default=False, action="store_true",
                 help="BED already sorted")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bed1, bed2 = args
    sorted = opts.sorted
    bed1_store, bed1_store_r = bed_store(bed1, sorted=sorted)
    bed2_store, bed2_store_r = bed_store(bed2, sorted=sorted)
    contrast_stores(bed1_store_r, bed2_store)
    contrast_stores(bed2_store_r, bed1_store, prefix="BA")


def resolve(args):
    """
    %prog resolve matrixfile fastafile bamfolder

    Separate repeats along collapsed contigs. First scan the matrixfile for
    largely heterozygous sites. For each heterozygous site, we scan each bam to
    retrieve distinct haplotypes. The frequency of each haplotype is then
    computed, the haplotype with the highest frequency, assumed to be
    paralogous, is removed.
    """
    import pysam
    from collections import defaultdict
    from itertools import groupby

    p = OptionParser(resolve.__doc__)
    p.add_option("--missing", default=.5, type="float",
                 help="Max level of missing data")
    p.add_option("--het", default=.5, type="float",
                 help="Min level of heterozygous calls")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    matrixfile, fastafile, bamfolder = args
    #f = Fasta(fastafile)
    fp = open(matrixfile)
    for row in fp:
        if row[0] != '#':
            break
    header = row.split()
    ngenotypes = len(header) - 4
    nmissing = int(round(opts.missing * ngenotypes))
    logging.debug("A total of {0} individuals scanned".format(ngenotypes))
    logging.debug("Look for markers with < {0} missing and > {1} het".\
                    format(opts.missing, opts.het))
    bamfiles = iglob(bamfolder, "*.bam")
    logging.debug("Folder `{0}` contained {1} bam files".\
                    format(bamfolder, len(bamfiles)))

    data = []
    for row in fp:
        if row[0] == '#':
            continue
        atoms = row.split()
        seqid, pos, ref, alt = atoms[:4]
        genotypes = atoms[4:]
        c = Counter(genotypes)
        c0 = c.get('0', 0)
        c3 = c.get('3', 0)
        if c0 >= nmissing:
            continue
        hetratio = c3 * 1. / (ngenotypes - c0)
        if hetratio <= opts.het:
            continue
        pos = int(pos)
        data.append((seqid, pos, ref, alt, c, hetratio))

    data.sort()
    logging.debug("A total of {0} target markers in {1} contigs.".\
                    format(len(data), len(set(x[0] for x in data))))
    samfiles = [pysam.AlignmentFile(x, "rb") for x in bamfiles]
    samfiles = [(op.basename(x.filename).split(".")[0], x) for x in samfiles]
    samfiles.sort()
    logging.debug("BAM files grouped to {0} individuals".\
                    format(len(set(x[0] for x in samfiles))))

    fw = must_open(opts.outfile, "w")
    for seqid, d in groupby(data, lambda x: x[0]):
        d = list(d)
        nmarkers = len(d)
        logging.debug("Process contig {0} ({1} markers)".format(seqid, nmarkers))
        haplotype_set = []
        for pf, sf in groupby(samfiles, key=lambda x: x[0]):
            haplotypes = []
            for pfi, samfile in sf:
                reads = defaultdict(list)
                positions = []
                for s, pos, ref, alt, c, hetratio in d:
                    for c in samfile.pileup(seqid):
                        if c.reference_pos != pos - 1:
                            continue
                        for r in c.pileups:
                            rname = r.alignment.query_name
                            rbase = r.alignment.query_sequence[r.query_position]
                            reads[rname].append((pos, rbase))
                    positions.append(pos)
                for read in reads.values():
                    hap = ['-'] * nmarkers
                    for p, rbase in read:
                        hap[positions.index(p)] = rbase
                    hap = "".join(hap)
                    if "-" in hap:
                        continue
                    haplotypes.append(hap)
            haplotypes = set(haplotypes)
            haplotype_set.append(haplotypes)
        hr = HaplotypeResolver(haplotype_set)
        print >> fw, seqid, hr
        hr.solve(fw)


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

    fastafile, = args
    csv = open(opts.csv, "w") if opts.csv else None

    f = Fasta(fastafile, lazy=True)
    sizes = []
    for desc, rec in f.iterdescriptions_ordered():
        if desc.startswith("singleton"):
            sizes.append(1)
            continue
        # consensus_for_cluster_0 with 63 sequences
        name, w, size, seqs = desc.split()
        if csv:
            print >> csv, "\t".join(str(x) for x in (name, size, len(rec)))
        assert w == "with"
        sizes.append(int(size))

    if csv:
        csv.close()
        logging.debug("File written to `{0}`".format(opts.csv))

    s = SummaryStats(sizes)
    print >> sys.stderr, s
    stem_leaf_plot(s.data, 0, 100, 20, title="Cluster size")


def novo(args):
    """
    %prog novo reads.fastq

    Reference-free tGBS pipeline.
    """
    from jcvi.assembly.kmer import jellyfish, histogram
    from jcvi.assembly.preprocess import diginorm
    from jcvi.formats.fasta import filter as fasta_filter, format
    from jcvi.apps.cdhit import filter as cdhit_filter

    p = OptionParser(novo.__doc__)
    p.add_option("--technology", choices=("illumina", "454", "iontorrent"),
                 default="iontorrent", help="Sequencing platform")
    p.add_option("--dedup", choices=("uclust", "cdhit"),
                 default="cdhit", help="Dedup algorithm")
    p.set_depth(depth=50)
    p.set_align(pctid=96)
    p.set_home("cdhit", default="/usr/local/bin/")
    p.set_home("fiona", default="/usr/local/bin/")
    p.set_home("jellyfish", default="/usr/local/bin/")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fastqfile, = args
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
        jellyfish([diginormfile, "--prefix={0}".format(pf),
                    "--cpus={0}".format(cpus),
                    "--jellyfish_home={0}".format(opts.jellyfish_home)])

    genomesize = histogram([jf, pf, "23"])
    fiona = pf + ".fiona.fa"
    if need_update(diginormfile, fiona):
        cmd = op.join(opts.fiona_home, "fiona")
        cmd += " -g {0} -nt {1} --sequencing-technology {2}".\
                    format(genomesize, cpus, opts.technology)
        cmd += " -vv {0} {1}".format(diginormfile, fiona)
        logfile = pf + ".fiona.log"
        sh(cmd, outfile=logfile, errfile=logfile)

    dedup = opts.dedup
    pctid = opts.pctid
    cons = fiona + ".P{0}.{1}.consensus.fasta".format(pctid, dedup)
    if need_update(fiona, cons):
        if dedup == "cdhit":
            deduplicate([fiona, "--consensus", "--reads",
                         "--pctid={0}".format(pctid),
                         "--cdhit_home={0}".format(opts.cdhit_home)])
        else:
            uclust([fiona, "--pctid={0}".format(pctid)])

    filteredfile = pf + ".filtered.fasta"
    if need_update(cons, filteredfile):
        covfile = pf + ".cov.fasta"
        cdhit_filter([cons, "--outfile={0}".format(covfile),
                      "--minsize={0}".format(depth / 5)])
        fasta_filter([covfile, "50", "--outfile={0}".format(filteredfile)])

    finalfile = pf + ".final.fasta"
    if need_update(filteredfile, finalfile):
        format([filteredfile, finalfile, "--sequential=replace",
                    "--prefix={0}_".format(pf)])


def snpflow(args):
    """
    %prog snpflow trimmed reference.fasta

    Run SNP calling pipeline until allele_counts are generated. This includes
    generation of native files, SNP_Het file. Speedup for fragmented genomes
    are also supported.
    """
    p = OptionParser(snpflow.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    trimmed, ref = args
    nseqs = len(Fasta(ref))
    supercat = nseqs >= 1000
    if supercat:
        logging.debug("Total seqs in ref: {0} (supercat={1})".\
                      format(nseqs, supercat))

    flist = iglob(trimmed, "*.fq", "*.fq.gz")
    nreads = len(flist)
    samples = sorted(set(op.basename(x).split(".")[0] for x in flist))
    nsamples = len(samples)
    logging.debug("Total {0} read files from {1} samples".\
                    format(nreads, nsamples))

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
    for f in flist:
        prefix = get_prefix(f, ref)
        gsnapfile = op.join(nativedir, prefix + ".gsnap")
        nativefile = op.join(nativedir, prefix + ".unique.native")
        samstatsfile = op.join(nativedir, prefix + ".unique.sam.stats")
        cmd = "python -m jcvi.apps.gmap align {0} {1}".format(gmapdb, f)
        cmd += " --outdir={0} --native --cpus=1".format(nativedir)
        mm.add((f, db), nativefile, cmd)

        cmd = "python -m jcvi.apps.gmap bam {0} {1} --cpus=1".\
                format(gsnapfile, gmapdb)
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
        mm.add(flist + allsnps, equalfile, cmd)

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
        mm.add([matrix] + flist, allele_counts, cmd)
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

    datafile, = args
    # Read in CDT file
    fp = open(datafile)
    fp.next(); fp.next()
    data = []
    for row in fp:
        atoms = row.split()[4:]
        nval = len(atoms)
        values = [float(x) for x in atoms]
        # normalize
        values = [x * 1. / sum(values) for x in values]
        data.append(values)

    pf = datafile.rsplit(".", 1)[0]
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    xmin, xmax = .1, .9
    ymin, ymax = .1, .9
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

    root.text(.05, .5, "{0} LMD50 SNPs".format(len(data)),
              ha="center", va="center", rotation=90, color="lightslategray")

    for x, t, c in zip((.3, .5, .7), ("REF", "ALT", "HET"), "rbg"):
        root.text(x, .95, t, color=c, ha="center", va="center")
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
