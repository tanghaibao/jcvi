#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Identify repeat numbers in STR repeats.
"""

import re
import os
import os.path as op
import sys
import vcf
import logging
import pyfasta
import numpy as np
import pandas as pd

from math import log, ceil
from collections import defaultdict
from multiprocessing import Pool

from jcvi.utils.cbook import percentage, uniqify
from jcvi.formats.bed import natsorted
from jcvi.apps.grid import MakeManager
from jcvi.formats.base import LineFile, must_open
from jcvi.utils.aws import push_to_s3, pull_from_s3, check_exists_s3, ls_s3
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, sh


READLEN = 150
MINSCORE = 30
YSEARCH_HAPLOTYPE = """
DYS393  DYS390 DYS19/DYS394  DYS19b        DYS391        DYS385a       DYS385b DYS426  DYS388  DYS439
DYS389I DYS392 DYS389B       DYS458        DYS459a/b     DYS459a/b     DYS455  DYS454  DYS447  DYS437
DYS448  DYS449 DYS464a/b/c/d DYS464a/b/c/d DYS464a/b/c/d DYS464a/b/c/d DYS464e DYS464f DYS464g DYS460
GATA-H4 YCAIIa YCAIIb        DYS456        DYS607        DYS576        DYS570  CDYa    CDYb    DYS442
DYS438  DYS531 DYS578        DYS395S1a/b   DYS395S1a/b   DYS590        DYS537  DYS641  DYS472  DYS406S1
DYS511  DYS425 DYS413a       DYS413b       DYS557        DYS594        DYS436  DYS490  DYS534  DYS450
DYS444  DYS481 DYS520        DYS446        DYS617        DYS568        DYS487  DYS572  DYS640  DYS492
DYS565  DYS461 DYS462        GATA-A10      DYS635        GAAT1B07      DYS441  DYS445  DYS452  DYS463
DYS434  DYS435 DYS485        DYS494        DYS495        DYS505        DYS522  DYS533  DYS549  DYS556
DYS575  DYS589 DYS636        DYS638        DYS643        DYS714        DYS716  DYS717  DYS726  DXYS156-Y
""".split()
YSEARCH_LL = """
L1  L2  L3  L4  L5  L6  L7  L8  L9  L10
L11 L12 L13 L14 L15 L16 L17 L18 L19 L20
L21 L22 L23 L24 L25 L26 L27 L28 L29 L30
L31 L32 L33 L34 L35 L36 L37 L38 L39 L40
L41 L54 L55 L56 L57 L58 L59 L60 L61 L62
L63 L42 L64 L65 L66 L67 L68 L69 L70 L71
L49 L72 L73 L51 L74 L75 L76 L77 L78 L79
L80 L43 L44 L45 L46 L47 L48 L50 L52 L53
L81 L82 L83 L84 L85 L86 L87 L88 L89 L90
L91 L92 L93 L94 L95 L96 L97 L98 L99 L100
""".split()
YHRD_YFILER = """
DYS456 DYS389I DYS390 DYS389B DYS458 DYS19/DYS394 DYS385
DYS393 DYS391 DYS439 DYS635 DYS392 GATA-H4 DYS437 DYS438 DYS448
""".split()
YHRD_YFILERPLUS = """
DYS576 DYS389I DYS635 DYS389B DYS627 DYS460 DYS458 DYS19/DYS394 GATA-H4 DYS448 DYS391
DYS456 DYS390 DYS438 DYS392 DYS518 DYS570 DYS437 DYS385a DYS449
DYS393 DYS439 DYS481 DYF387S1 DYS533
""".split()
USYSTR_ALL = """
DYF387S1 DYS19/DYS394 DYS385 DYS389I
DYS389B DYS390 DYS391 DYS392
DYS393 DYS437 DYS438 DYS439
DYS448 DYS449 DYS456 DYS458
DYS460 DYS481 DYS518 DYS533
DYS549 DYS570 DYS576 DYS627
DYS635 DYS643 GATA-H4
""".split()


class STRLine(object):

    def __init__(self, line, named=False):
        args = line.split()
        self.seqid = args[0]
        self.start = int(args[1])
        self.end = int(args[2])
        self.period = int(args[3])
        self.copynum = args[4]
        self.consensusSize = int(args[5])
        self.pctmatch = int(args[6])
        self.pctindel = int(args[7])
        self.score = args[8]
        self.A = args[9]
        self.C = args[10]
        self.G = args[11]
        self.T = args[12]
        self.entropy = float(args[13])
        self.motif = args[14]
        assert self.period == len(self.motif)
        if named:
            self.name = args[15] if len(args) == 16 else None
        else:
            self.name = None

    def __str__(self):
        fields = [self.seqid, self.start, self.end,
                  self.period, self.copynum, self.consensusSize,
                  self.pctmatch, self.pctindel, self.score,
                  self.A, self.C, self.G, self.T,
                  "{0:.2f}".format(self.entropy), self.motif]
        if self.name is not None:
            fields += [self.name]
        return "\t".join(str(x) for x in fields)

    @property
    def longname(self):
        name = "_".join(str(x) for x in \
                    (self.seqid, self.start, self.motif,
                    int(float(self.copynum))))
        if self.name is not None:
            name += "_" + self.name
        return name

    def is_valid(self, maxperiod=6, maxlength=READLEN, minscore=MINSCORE):
        return 1 <= self.period <= maxperiod and \
               (self.end - self.start + 1) <= maxlength and \
               self.score >= minscore

    def calc_entropy(self):
        total = self.A + self.C + self.G + self.T
        if total == 0:  # Perhaps they are all Ns - might crash in lobstrindex()
            return 0
        fractions = [x * 1.0 / total for x in [self.A, self.C, self.G, self.T]]
        entropy = sum([-1.0 * x * log(x, 2) for x in fractions if x != 0])
        return entropy

    def iter_exact_str(self, genome):
        pat = re.compile("(({0}){{2,}})".format(self.motif))
        start = self.start
        s = genome[self.seqid][self.start - 1: self.end].upper()
        for m in re.finditer(pat, s):
            self.start = start + m.start()
            length = m.end() - m.start()
            subseq = m.group(0)
            assert length % self.period == 0
            assert subseq.startswith(self.motif)

            self.end = self.start - 1 + length
            self.copynum = length / self.period
            self.pctmatch = 100
            self.pctindel = 0
            self.score = 2 * length
            self.fix_counts(subseq)
            yield self

    def fix_counts(self, subseq):
        length = int(ceil(self.period * float(self.copynum)))
        # Sanity check for length, otherwise lobSTR misses units
        self.end = max(self.end, self.start + length - 1)
        self.A = subseq.count('A')
        self.C = subseq.count('C')
        self.G = subseq.count('G')
        self.T = subseq.count('T')
        self.entropy = self.calc_entropy()


class STRFile (LineFile):

    def __init__(self, lobstr_home, db="hg38"):
        filename = op.join(lobstr_home, "{0}/index.info".format(db))
        super(STRFile, self).__init__(filename)
        fp = open(filename)
        for row in fp:
            self.append(STRLine(row))

    @property
    def ids(self):
        return [s.longname for s in self]

    @property
    def register(self):
        return dict(((s.seqid, s.start), s.name) for s in self)


class LobSTRvcf(dict):

    def __init__(self, columnidsfile="STR.ids"):
        self.samplekey = None
        fp = open(columnidsfile)
        self.columns = [x.strip() for x in fp]
        self.evidence = {}  # name: (supporting reads, stutter reads)
        logging.debug("A total of {} markers imported".format(len(self.columns)))

    def parse(self, filename, cleanup=False):
        self.samplekey = op.basename(filename).split(".")[0]
        fp = must_open(filename)
        reader = vcf.Reader(fp)
        for record in reader:
            info = record.INFO
            ref = float(info["REF"])
            rpa = info.get("RPA", ref)
            motif = info["MOTIF"]
            name = "_".join(str(x) for x in (record.CHROM, record.POS,
                                             motif, int(ref)))
            for sample in record.samples:
                gt = sample["GT"]
                if gt == "0/0":
                    alleles = (ref, ref)
                elif gt in ("0/1", "1/0"):
                    alleles = (ref, rpa[0])
                elif gt == "1/1":
                    alleles = (rpa[0], rpa[0])
                elif gt == "1/2":
                    alleles = rpa
                self[name] = ",".join(str(int(x)) for x in sorted(alleles))

                # Collect supporting read evidence
                motif_length = len(motif)
                adjusted_alleles = [(x - ref) * motif_length for x in alleles]
                support = stutter = 0
                allreads = sample["ALLREADS"]
                for r in allreads.split(";"):
                    k, v = r.split("|")
                    k, v = int(k), int(v)
                    min_dist = min([abs(k - x) for x in adjusted_alleles])
                    if motif_length * .5 < min_dist < motif_length * 1.5:
                        stutter += v
                    support += v
                #print name, ref, allreads, alleles, adjusted_alleles, \
                #      "stutter={} support={}".format(stutter, support)
                self.evidence[name] = "{},{}".format(stutter, support)

        if cleanup:
            sh("rm -f {}".format(op.basename(filename)))

    @property
    def csvline(self):
        return ",".join([self.get(c, "-1,-1") for c in self.columns])

    @property
    def evline(self):
        return ",".join([self.evidence.get(c, "-1,-1") for c in self.columns])


def main():

    actions = (
        ('compile', "compile vcf results into master spreadsheet"),
        ('mergecsv', "combine csv into binary array"),
        ('count', 'count alleles in master spreadsheet'),
        ('filterloci', 'select subset of loci based on allele counts'),
        ('filterdata', 'select subset of data based on filtered loci'),
        ('batchlobstr', "run batch lobSTR"),
        ('batchhtt', "run batch HTT caller"),
        ('htt', 'extract HTT region and run lobSTR'),
        ('liftover', 'liftOver CODIS/Y-STR markers'),
        ('lobstr', 'run lobSTR on a big BAM'),
        ('lobstrindex', 'make lobSTR index'),
        ('trf', 'run TRF on FASTA files'),
        ('ystr', 'print out Y-STR info given VCF'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def counts_to_percentile(counts):
    percentile = {}
    s = 0
    for k, v in sorted(counts.items(), reverse=True):
        s += v
        percentile[k] = s
    for k, v in percentile.items():
        v = "{:.6f}".format(v * 1. / s)
        percentile[k] = v
    return percentile


def convert_to_percentile(arg):
    i, a, percentile = arg
    xb = a % 1000
    pp = np.array([percentile.get(x, "1.000000") for x in xb], dtype="S8")
    if i % 1000 == 0:
        print >> sys.stderr, i
        print >> sys.stderr, xb
        print >> sys.stderr, pp
    return i, pp


def write_csv(csvfile, m, index, columns, sep="\t", index_label="SampleKey"):
    fw = open(csvfile, "w")
    print >> fw, "\t".join([index_label] + columns)
    for i, a in enumerate(m):
        print >> fw, index[i] + "\t" + "\t".join(str(x) for x in a)
    fw.close()


def filterdata(args):
    """
    %prog filterdata data.bin samples.ids STR.ids allele_freq remove.ids final.ids

    Filter subset of data after dropping remove.ids.
    """
    p = OptionParser(filterdata.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 6:
        sys.exit(not p.print_help())

    binfile, sampleids, strids, af, remove, final = args
    df, m, samples, loci = read_binfile(binfile, sampleids, strids)
    remove = [x.strip() for x in open(remove)]
    removes = set(remove)
    final = [x.strip() for x in open(final)]
    assert len(loci) == len(remove) + len(final)

    fp = open(af)
    percentiles = {}
    for row in fp:
        sname, counts = row.split()
        countst = [x for x in counts.strip("{}").split(",") if x]
        countsd = {}
        for x in countst:
            a, b = x.split(":")
            countsd[int(a)] = int(b)
        percentile = counts_to_percentile(countsd)
        percentiles[sname] = percentile

    p = Pool(processes=opts.cpus)
    run_args = []
    for i, sname in enumerate(loci):
        if sname in removes:
            continue
        a = m[:, i]
        percentile = percentiles[sname]
        run_args.append((i, a, percentile))

    res = []
    for r in p.map_async(convert_to_percentile, run_args).get():
        res.append(r)
    res.sort()

    # Write mask (P-value) matrix
    ii, pvalues = zip(*res)
    m = np.vstack(pvalues).T
    write_csv("final.mask.tsv", m, samples, final)

    df.drop(remove, inplace=True, axis=1)
    df.columns = final

    # Save a copy of the raw numpy array
    filtered_bin = "filtered.bin"
    m = df.as_matrix()
    m[m < 0] = -1
    m.tofile(filtered_bin)
    logging.debug("Binary matrix written to `{}`".format(filtered_bin))

    # Write data output
    df.to_csv("final.data.tsv", sep="\t", index_label="SampleKey")


def filterloci(args):
    """
    %prog filterloci allele_freq STR-exons.wo.bed SAMPLES

    Filter subset of loci that satisfy:
    1. no redundancy (unique chr:pos)
    2. variable (n_alleles > 1)
    3. low level of missing data (>= 50% autosomal + X, > 25% for Y)

    Write meta file with the following infor:
    1. id
    2. title
    3. gene_name
    4. variant_type
    5. motif
    6. allele_frequency

    `STR-exons.wo.bed` can be generated like this:
    $ tail -n 854476 /mnt/software/lobSTR-4.0.0/hg38/index.tab | cut -f1-3 > all-STR.bed
    $ intersectBed -a all-STR.bed -b all-exons.bed -wo > STR-exons.wo.bed
    """
    p = OptionParser(filterloci.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    af, wobed, samples = args
    nsamples = len([x.strip() for x in open(samples)])
    nalleles = nsamples * 2
    logging.debug("Load {} samples ({} alleles) from `{}`".\
                    format(nsamples, nalleles, samples))

    logging.debug("Load gene intersections from `{}`".format(wobed))
    fp = open(wobed)
    gene_map = defaultdict(set)
    for row in fp:
        chr1, start1, end1, chr2, start2, end2, name, ov = row.split()
        gene_map[(chr1, start1)] |= set(name.split(","))
    for k, v in gene_map.items():
        non_enst = sorted(x for x in v if not x.startswith("ENST"))
        enst = sorted(x.rsplit(".", 1)[0] for x in v if x.startswith("ENST"))
        gene_map[k] = ",".join(non_enst + enst)

    logging.debug("Filtering loci from `{}`".format(af))
    fp = open(af)
    treds = """chr19_45770205_CAG
    chr6_170561926_CAG
    chrX_67545318_CAG
    chr9_69037287_GAA
    chrX_147912051_CGG
    chr4_3074877_CAG
    chrX_148500639_CCG
    chr12_6936729_CAG
    chr13_70139384_CTG
    chr6_16327636_CTG
    chr14_92071011_CTG
    chr12_111598951_CTG
    chr3_63912686_CAG
    chr19_13207859_CTG""".split()
    seen = set(treds)
    remove = []
    fw = open("meta.tsv", "w")
    header = "id title gene_name variant_type motif allele_frequency".\
                replace(" ",  "\t")
    print >> fw, header
    variant_type = "short tandem repeats"
    title = "Short tandem repeats ({})n"
    for row in fp:
        sname, counts = row.split()
        name = sname.rsplit("_", 1)[0]
        seqid, pos, motif = name.split("_")
        countst = [x for x in counts.strip("{}").split(",") if x]
        countsd = {}
        for x in countst:
            a, b = x.split(":")
            countsd[int(a)] = int(b)

        if counts_filter(countsd, nalleles, seqid):
            remove.append(sname)
            continue

        if name in seen:
            remove.append(sname)
            continue
        seen.add(name)

        gene_name = gene_map.get((seqid, pos), "")
        print >> fw, "\t".join((name, title.format(motif),
                                gene_name, variant_type, motif, counts))
    fw.close()

    removeidsfile = "remove.ids"
    fw = open(removeidsfile, "w")
    print >> fw, "\n".join(remove)
    fw.close()
    logging.debug("A total of {} filtered loci written to `{}`".\
                    format(len(remove), removeidsfile))


def counts_filter(countsd, nalleles, seqid):
    # Check for variability
    if len(countsd) < 2:
        return True

    # Check for missingness
    observed = sum(countsd.values())
    observed_pct = observed * 100  / nalleles
    if observed_pct < 50:
        if not (seqid == "chrY" and observed_pct >= 25):
            return True


def read_binfile(binfile, sampleids, strids):
    m = np.fromfile(binfile, dtype=np.int32)
    samples = [x.strip() for x in open(sampleids)]
    loci = [x.strip() for x in open(strids)]
    nsamples, nloci = len(samples), len(loci)
    print >> sys.stderr, "{} x {} entries imported".format(nsamples, nloci)

    m.resize(nsamples, nloci)
    df = pd.DataFrame(m, index=samples, columns=loci)
    return df, m, samples, loci


def count(args):
    """
    %prog count data.bin samples.ids STR.ids

    Count alleles in master spreadsheet per locus.
    """
    from collections import Counter

    p = OptionParser(count.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    binfile, sampleids, strids = args
    df, m, samples, loci = read_binfile(binfile, sampleids, strids)
    fw = open("allele_freq", "w")
    for i, columnname in enumerate(loci):
        a = m[:, i]
        counts = Counter()
        xa = a / 1000
        xb = a % 1000
        counts.update(xa)
        counts.update(xb)
        del counts[-1]
        af = "{" + ",".join("{}:{}".format(k, v) for k, v in sorted(counts.items())) + "}"
        print >> fw, "\t".join((columnname, af))
    fw.close()


def mergecsv(args):
    """
    %prog mergecsv *.csv

    Combine CSV into binary array.
    """
    p = OptionParser(mergecsv.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    csvfiles = args
    arrays = []
    samplekeys = []
    for csvfile in csvfiles:
        samplekey = op.basename(csvfile).split(".")[0]
        a = np.fromfile(csvfile, sep=",", dtype=np.int32)
        x1 = a[::2]
        x2 = a[1::2]
        a = x1 * 1000 + x2
        a[a < 0] = -1
        arrays.append(a)
        samplekeys.append(samplekey)
        print >> sys.stderr, samplekey, a
    print >> sys.stderr, "Merging"
    b = np.concatenate(arrays)
    b.tofile("data.bin")

    fw = open("samples", "w")
    print >> fw, "\n".join(samplekeys)
    fw.close()


def write_csv_ev(filename, store, cleanup):
    lv = LobSTRvcf()
    lv.parse(filename, cleanup=cleanup)
    lv.parse(filename.replace(".hg38.", ".hg38-named."), cleanup=cleanup)

    csvfile = op.basename(filename) + ".csv"
    fw = open(csvfile, "w")
    print >> fw, lv.csvline
    fw.close()

    evfile = op.basename(filename) + ".ev"
    fw = open(evfile, "w")
    print >> fw, lv.evline
    fw.close()

    # Save to s3
    push_to_s3(store, csvfile)
    push_to_s3(store, evfile)


def run(arg):
    filename, store, cleanup = arg
    csvfile = filename + ".csv"
    #evfile = filename + ".ev"
    try:
        if check_exists_s3(csvfile):
            logging.debug("{} exists. Skipped.".format(csvfile))
        else:
            write_csv_ev(filename, store, cleanup)
            logging.debug("{} written and uploaded.".format(csvfile))
    except:
        logging.debug("Thread failed!")


def compile(args):
    """
    %prog compile samples.csv

    Compile vcf results into master spreadsheet.
    """
    p = OptionParser(compile.__doc__)
    p.add_option("--db", default="hg38,hg38-named",
                 help="Use these lobSTR db")
    p.set_home("lobstr")
    p.set_cpus()
    p.set_aws_opts(store="hli-mv-data-science/htang/str")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    samples, = args
    workdir = opts.workdir
    dbs = opts.db.split(",")
    cwd = os.getcwd()
    mkdir(workdir)
    os.chdir(workdir)
    samples = op.join(cwd, samples)

    stridsfile = "STR.ids"
    vcffiles = [x.strip() for x in must_open(samples)]
    if not op.exists(stridsfile):
        ids = []
        for db in dbs:
            ids.extend(STRFile(opts.lobstr_home, db=db).ids)
        uids = uniqify(ids)
        logging.debug("Combined: {} Unique: {}".format(len(ids), len(uids)))

        fw = open(stridsfile, "w")
        print >> fw, "\n".join(uids)
        fw.close()

        # Generate two alleles
        dipuids = []
        for uid in uids:
            dipuids.extend([uid + ".1", uid + ".2"])
        fw = open("header.ids", "w")
        print >> fw, ",".join(dipuids)
        fw.close()

    p = Pool(processes=opts.cpus)
    run_args = [(x, opts.store, opts.cleanup) for x in vcffiles]
    for res in p.map_async(run, run_args).get():
        continue


def build_ysearch_link(r, ban=["DYS520", "DYS413a", "DYS413b"]):
    template = \
    "http://www.ysearch.org/search_search.asp?fail=2&uid=&freeentry=true&"
    markers = []
    for i, marker in zip(YSEARCH_LL, YSEARCH_HAPLOTYPE):
        z = r.get(marker, "null")
        if "a/b" in marker or marker in ban:
            z = "null"
        m = "{0}={1}".format(i, z)
        markers.append(m)
    print template + "&".join(markers)


def build_yhrd_link(r, panel, ban=["DYS385"]):
    L = []
    for marker in panel:
        z = r.get(marker, "--")
        if marker in ban:
            z = "--"
        L.append(z)
    print " ".join(str(x) for x in L)


def ystr(args):
    """
    %prog ystr chrY.vcf

    Print out Y-STR info given VCF. Marker name extracted from tabfile.
    """
    from jcvi.utils.table import write_csv

    p = OptionParser(ystr.__doc__)
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    vcffile, = args
    si = STRFile(opts.lobstr_home, db="hg38-named")
    register = si.register

    header = "Marker|Reads|Ref|Genotype|Motif".split("|")
    contents = []
    fp = must_open(vcffile)
    reader = vcf.Reader(fp)
    simple_register = {}
    for record in reader:
        name = register[(record.CHROM, record.POS)]
        info = record.INFO
        ref = int(float(info["REF"]))
        rpa = info.get("RPA", ref)
        if isinstance(rpa, list):
            rpa = "|".join(str(int(float(x))) for x in rpa)
        ru = info["RU"]
        simple_register[name] = rpa
        for sample in record.samples:
            contents.append((name, sample["ALLREADS"], ref, rpa, ru))

    # Multi-part markers
    a, b, c = "DYS389I", "DYS389B.1", "DYS389B"
    if a in simple_register and b in simple_register:
        simple_register[c] = int(simple_register[a]) + int(simple_register[b])

    # Multi-copy markers
    mm = ["DYS385", "DYS413", "YCAII"]
    for m in mm:
        ma, mb = m + 'a', m + 'b'
        if ma not in simple_register or mb not in simple_register:
            simple_register[ma] = simple_register[mb] = None
            del simple_register[ma]
            del simple_register[mb]
            continue
        if simple_register[ma] > simple_register[mb]:
            simple_register[ma], simple_register[mb] = \
                    simple_register[mb], simple_register[ma]

    write_csv(header, contents, sep=" ")
    print "[YSEARCH]"
    build_ysearch_link(simple_register)
    print "[YFILER]"
    build_yhrd_link(simple_register, panel=YHRD_YFILER)
    print "[YFILERPLUS]"
    build_yhrd_link(simple_register, panel=YHRD_YFILERPLUS)
    print "[YSTR-ALL]"
    build_yhrd_link(simple_register, panel=USYSTR_ALL)


def get_motif(s, motif_length):
    sl = len(s)
    kmers = set()
    # Get all kmers
    for i in xrange(sl - motif_length):
        ss = s[i: i + motif_length]
        kmers.add(ss)

    kmer_counts = []
    for kmer in kmers:
        kmer_counts.append((s.count(kmer), -s.index(kmer), kmer))

    return sorted(kmer_counts, reverse=True)[0][-1]


def liftover(args):
    """
    %prog liftover lobstr_v3.0.2_hg38_ref.bed hg38.upper.fa

    LiftOver CODIS/Y-STR markers.
    """
    p = OptionParser(liftover.__doc__)
    p.add_option("--checkvalid", default=False, action="store_true",
                help="Check minscore, period and length")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    refbed, fastafile = args
    genome = pyfasta.Fasta(fastafile)
    edits = []
    fp = open(refbed)
    for i, row in enumerate(fp):
        s = STRLine(row)
        seq = genome[s.seqid][s.start - 1: s.end].upper()
        s.motif = get_motif(seq, len(s.motif))
        s.fix_counts(seq)
        if opts.checkvalid and not s.is_valid():
            continue
        edits.append(s)
        if i % 10000 == 0:
            print >> sys.stderr, i, "lines read"

    edits = natsorted(edits, key=lambda x: (x.seqid, x.start))
    for e in edits:
        print str(e)


def trf(args):
    """
    %prog trf outdir

    Run TRF on FASTA files.
    """
    from jcvi.apps.base import iglob

    p = OptionParser(trf.__doc__)
    p.add_option("--mismatch", default=31, type="int",
                 help="Mismatch and gap penalty")
    p.add_option("--minscore", default=MINSCORE, type="int",
                 help="Minimum score to report")
    p.add_option("--period", default=6, type="int",
                 help="Maximum period to report")
    p.add_option("--telomeres", default=False, action="store_true",
                 help="Run telomere search: minscore=140 period=7")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    outdir, = args
    mm = MakeManager()
    if opts.telomeres:
        opts.minscore, opts.period = 140, 7

    params = "2 {0} {0} 80 10 {1} {2}".\
            format(opts.mismatch, opts.minscore, opts.period).split()
    bedfiles = []
    for fastafile in natsorted(iglob(outdir, "*.fa,*.fasta")):
        pf = op.basename(fastafile).split(".")[0]
        cmd1 = "trf {0} {1} -d -h".format(fastafile, " ".join(params))
        datfile = op.basename(fastafile) + "." + ".".join(params) + ".dat"
        bedfile = "{0}.trf.bed".format(pf)
        cmd2 = "cat {0} | awk '($9 >0)' | sed 's/ /\\t/g'".format(datfile)
        cmd2 += " | awk '{{print \"{0}\\t\" $0}}' > {1}".format(pf, bedfile)
        mm.add(fastafile, datfile, cmd1)
        mm.add(datfile, bedfile, cmd2)
        bedfiles.append(bedfile)

    bedfile = "trf.bed"
    cmd = "cat {0} > {1}".format(" ".join(natsorted(bedfiles)), bedfile)
    mm.add(bedfiles, bedfile, cmd)

    mm.write()


def batchlobstr(args):
    """
    %prog batchlobstr samples.csv

    Run lobSTR sequentially on list of samples. Each line contains:
    sample-name,s3-location
    """
    p = OptionParser(batchlobstr.__doc__)
    p.add_option("--sep", default=",", help="Separator for building commandline")
    p.set_aws_opts(store="hli-mv-data-science/htang/str")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    samplesfile, = args
    store = opts.store
    computed = ls_s3(store)
    fp = open(samplesfile)
    skipped = total = 0
    for row in fp:
        total += 1
        sample, s3file = row.strip().split(",")[:2]
        bamfile = s3file.replace(".gz", "").replace(".vcf", ".bam")

        gzfile = sample + ".{0}.vcf.gz".format("hg38-named")
        if gzfile in computed:
            skipped += 1
            continue

        print opts.sep.join("python -m jcvi.variation.str lobstr".split() + \
                            [bamfile, "hg38", "hg38-named",
                            "--prefix", sample,
                            "--workdir", opts.workdir, "--cleanup"])
    fp.close()
    logging.debug("Total skipped: {0}".format(percentage(skipped, total)))


def lobstr(args):
    """
    %prog lobstr bamfile lobstr_index1 lobstr_index2 ...

    Run lobSTR on a big BAM file. There can be multiple lobSTR indices.
    """
    p = OptionParser(lobstr.__doc__)
    p.add_option("--chr", help="Run only this chromosome")
    p.add_option("--prefix", help="Use prefix file name")
    p.set_home("lobstr")
    p.set_cpus()
    p.set_aws_opts(store="hli-mv-data-science/htang/str")
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    bamfile = args[0]
    lbindices = args[1:]
    s3mode = bamfile.startswith("s3")
    store = opts.store
    workdir = opts.workdir
    mkdir(workdir)
    os.chdir(workdir)

    pf = opts.prefix or bamfile.split("/")[-1].split(".")[0]
    if s3mode:
        gzfile = pf + ".{0}.vcf.gz".format(lbindices[-1])
        remotegzfile = "s3://{0}/{1}".format(store, gzfile)
        if check_exists_s3(remotegzfile):
            logging.debug("Object `{0}` exists. Computation skipped."\
                            .format(remotegzfile))
            return
        localbamfile = pf + ".bam"
        localbaifile = localbamfile + ".bai"
        if op.exists(localbamfile):
            logging.debug("BAM file already downloaded.")
        else:
            pull_from_s3(bamfile, localbamfile)
        if op.exists(localbaifile):
            logging.debug("BAM index file already downloaded.")
        else:
            remotebaifile = bamfile + ".bai"
            if check_exists_s3(remotebaifile):
                pull_from_s3(remotebaifile, localbaifile)
            else:
                remotebaifile = bamfile.rsplit(".")[0] + ".bai"
                if check_exists_s3(remotebaifile):
                    pull_from_s3(remotebaifile, localbaifile)
                else:
                    logging.debug("BAM index cannot be found in S3!")
                    sh("samtools index {0}".format(localbamfile))
        bamfile = localbamfile

    lhome = opts.lobstr_home
    chrs = [opts.chr] if opts.chr else (range(1, 23) + ["X", "Y"])
    for lbidx in lbindices:
        mm = MakeManager(filename="makefile.{0}".format(lbidx))
        vcffiles = []
        for chr in chrs:
            cmd, vcffile = allelotype_on_chr(bamfile, chr, lhome, lbidx)
            mm.add(bamfile, vcffile, cmd)
            vcffiles.append(vcffile)

        gzfile = bamfile.split(".")[0] + ".{0}.vcf.gz".format(lbidx)
        cmd = "vcf-concat {0} | vcf-sort".format(" ".join(vcffiles))
        cmd += " | bgzip -c > {0}".format(gzfile)
        mm.add(vcffiles, gzfile, cmd)
        mm.run(cpus=opts.cpus)

        if s3mode:
            push_to_s3(store, gzfile)

    if opts.cleanup:
        sh("rm -f *")


def allelotype_on_chr(bamfile, chr, lhome, lbidx):
    outfile = "{0}.chr{1}".format(bamfile.split(".")[0], chr)
    cmd = "allelotype --command classify --bam {1}".format(lhome, bamfile)
    cmd += " --noise_model {0}/models/illumina_v3.pcrfree".format(lhome)
    cmd += " --strinfo {0}/{1}/index.tab".format(lhome, lbidx)
    cmd += " --index-prefix {0}/{1}/lobSTR_".format(lhome, lbidx)
    cmd += " --chrom chr{0} --out {1}.{2}".format(chr, outfile, lbidx)
    cmd += " --max-diff-ref {0}".format(READLEN)
    cmd += " --haploid chrY"
    return cmd, ".".join((outfile, lbidx, "vcf"))


def batchhtt(args):
    """
    %prog htt samples.csv

    Batch HTT caller. The samples are:
    SampleId,ExecID,Norm VCF,gVCF,BAM
    """
    p = OptionParser(batchhtt.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    samplesfile, = args
    fp = open(samplesfile)
    fp.next()  # header
    for row in fp:
        sample, ex, vcf, gvcf, bam = row.strip().split(",")
        htt([bam])
    fp.close()


def htt(args):
    """
    %prog htt bamfile

    Extract HTT region and run lobSTR.
    """
    p = OptionParser(htt.__doc__)
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bamfile, = args
    lhome = opts.lobstr_home

    minibamfile = bamfile.split("/")[-1]
    cmd = "samtools view {0} chr4:3070000-3080000 -b".format(bamfile)
    cmd += " -o {0}".format(minibamfile)
    sh(cmd)

    sh("rm {0}.bai".format(minibamfile))
    sh("samtools index {0}".format(minibamfile))

    cmd = allelotype_on_chr(minibamfile, 4, lhome, "hg38-named")
    sh(cmd)


def lobstrindex(args):
    """
    %prog lobstrindex hg38.trf.bed hg38.upper.fa hg38

    Make lobSTR index. Make sure the FASTA contain only upper case (so use
    fasta.format --upper to convert from UCSC fasta). The bed file is generated
    by str().
    """
    p = OptionParser(lobstrindex.__doc__)
    p.add_option("--fixseq", action="store_true", default=False,
                 help="Scan sequences to extract perfect STRs")
    p.set_home("lobstr")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    trfbed, fastafile, pf = args
    lhome = opts.lobstr_home
    mkdir(pf)

    if opts.fixseq:
        genome = pyfasta.Fasta(fastafile)
        newbedfile = trfbed + ".new"
        newbed = open(newbedfile, "w")
        fp = open(trfbed)
        retained = total = 0
        for row in fp:
            s = STRLine(row)
            total += 1
            for ns in s.iter_exact_str(genome):
                if not ns.is_valid():
                    continue
                print >> newbed, ns
                retained += 1
        newbed.close()
        logging.debug("Retained: {0}".format(percentage(retained, total)))
    else:
        newbedfile = trfbed

    mm = MakeManager()
    cmd = "python {0}/scripts/lobstr_index.py".format(lhome)
    cmd += " --str {0} --ref {1} --out {2}".format(newbedfile, fastafile, pf)
    mm.add((newbedfile, fastafile), op.join(pf, "lobSTR_ref.fasta.rsa"), cmd)

    tabfile = "{0}/index.tab".format(pf)
    cmd = "python {0}/scripts/GetSTRInfo.py".format(lhome)
    cmd += " {0} {1} > {2}".format(newbedfile, fastafile, tabfile)
    mm.add((newbedfile, fastafile), tabfile, cmd)

    infofile = "{0}/index.info".format(pf)
    cmd = "cp {0} {1}".format(trfbed, infofile)
    mm.add(trfbed, infofile, cmd)
    mm.write()


if __name__ == '__main__':
    main()
