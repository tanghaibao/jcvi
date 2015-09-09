#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using VCLUST to derep, cluster, and make consensus from duplicate reads.
The VCLUST implementation borrows ideas and code from PyRAD. PyRAD link:

<https://github.com/dereneaton/pyrad>
"""

import os.path as op
import shutil
import sys
import logging
import numpy as np
import scipy
import scipy.stats
import scipy.optimize

from collections import defaultdict
from copy import deepcopy
from functools import partial
from itertools import groupby
from subprocess import Popen, PIPE, STDOUT
from tempfile import mkdtemp

from jcvi.formats.base import BaseFile, FileMerger, must_open, split
from jcvi.formats.fasta import parse_fasta
from jcvi.formats.fastq import fasta
from jcvi.utils.cbook import memoized
from jcvi.utils.orderedcollections import DefaultOrderedDict
from jcvi.utils.iter import grouper
from jcvi.apps.base import OptionParser, ActionDispatcher, datadir, listify, mkdir, \
            need_update, sh


SEP = "//"
BASES = "ACTGN_-"  # CAUTION: change of this line must also change AlleleCounts
NBASES = len(BASES)
AMB = "RKSYWM"
AMBL = "rksywm"
EXTENDEDBASES = BASES + AMB
ACHEADER = """
TAXON     CHR   POS     REF_NT  REF_ALLELE      ALT_ALLELE      REF_COUNT
ALT_COUNT       OTHER_COUNT     TOTAL_READS     A       G       C       T
READ_INS        READ_DEL        TOTAL_READS
""".split()
ACHEADER_NO_TAXON = ACHEADER[1:]


alleles = lambda x: (",".join(x).replace("-", "*") if x else "N")
getsize = lambda name: (0 if ";" not in name else \
            int(name.split(";")[1].replace("size=", "")))


class ClustFile (BaseFile):

    def __init__(self, filename):
        super(ClustFile, self).__init__(filename)

    def __iter__(self):
        nstacks = 0
        fp = must_open(self.filename)
        for tag, contents in groupby(fp, lambda row: row[0] == '/'):
            if tag:
                continue
            data = Clust()
            for name, seq in grouper(contents, 2):
                name, seq = name.strip(), seq.strip()
                nrep = getsize(name)
                data.append((name, seq, nrep))
            yield data
            nstacks += 1
            if nstacks % 5000 == 0:
                logging.debug("{0} stacks parsed".format(nstacks))


class Clust (list):

    def __init__(self):
        super(Clust, self).__init__(self)

    def __str__(self):
        s = []
        for d in self:
            s.append("\n".join(d[:2]))
        return "\n".join(s) + "\n" + SEP


class ClustStore (BaseFile):
    def __init__(self, consensfile):
        super(ClustStore, self).__init__(consensfile)
        binfile = consensfile + ".bin"
        idxfile = consensfile + ".idx"
        self.bin = np.fromfile(binfile, dtype=np.uint16)
        assert self.bin.size % NBASES == 0

        self.bin = self.bin.reshape((self.bin.size / NBASES, NBASES))
        self.index = {}
        fp = open(idxfile)
        for row in fp:
            name, start, end = row.split()
            start, end = int(start), int(end)
            self.index[name.strip(">")] = (start, end)

    def __getitem__(self, name):
        start, end = self.index[name]
        return self.bin[start:end, :]


class AlleleCount (object):
    """
    Each record represents a line in the .allele_count file

    Fields are:
    # CHR   POS     REF_NT  REF_ALLELE      ALT_ALLELE      REF_COUNT
    # ALT_COUNT       OTHER_COUNT     TOTAL_READS     A       G       C       T
    # READ_INS        READ_DEL        TOTAL_READS
    """
    def __init__(self, taxon, chr, pos, ref_allele, alt_allele, profile):
        self.taxon = taxon
        self.chr = chr
        self.pos = pos
        self.ref_nt = listify(ref_allele)
        self.ref_allele = listify(ref_allele)
        alts = []
        for a in listify(alt_allele):
            alts.extend(unhetero(a))
        self.alt_allele = sorted(set(alts))
        self.update(profile)

    def tostring(self, taxon=False):
        ref_allele = alleles(self.ref_allele)
        ar = [self.chr, self.pos,
                ref_allele, ref_allele, alleles(self.alt_allele),
                self.ref_count, self.alt_count, self.other_count, self.total_count,
                self.A, self.G, self.C, self.T,
                self.read_ins, self.read_del, self.total_count]
        if taxon:
            ar = [self.taxon] + ar
        return "\t".join(str(x) for x in ar)

    def update(self, profile):
        self.ref_count = sum(profile[BASES.index(x)] for x in self.ref_allele)
        self.alt_count = sum(profile[BASES.index(x)] for x in self.alt_allele)
        self.A, self.C, self.T, self.G, N, tgaps, gaps = profile
        self.total_count = sum(profile) - tgaps
        others = set(BASES) - set(self.ref_allele) - set(self.alt_allele)
        self.other_count = sum(profile[BASES.index(x)] for x in others) - tgaps
        self.read_ins = self.total_count if '-' in self.ref_allele else 0
        self.read_del = gaps

    def clear(self):
        self.update([0] * NBASES)


class ClustStores (dict):
    """
    ClustStores provides random access to any consensus read
    """
    def __init__(self, consensfiles):
        super(ClustStores, self).__init__(self)
        for cs in consensfiles:
            name = op.basename(cs).split(".")[0]
            self[name] = ClustStore(cs)


def main():

    actions = (
        ('align', 'align clustfile to clustSfile'),
        ('estimateHE', 'estimate heterozygosity and error rate for stacks'),
        ('cluster', 'cluster within samples'),
        ('consensus', 'call consensus bases within samples'),
        ('mcluster', 'cluster across samples'),
        ('mconsensus', 'call consensus bases across samples'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add_consensus_options(p):
    p.add_option("--prefix", default="mcluster", help="Output prefix")
    p.add_option("--minlength", default=30, type="int", help="Min contig length")
    p.add_option("--mindepth", default=3, type="int", help="Min depth for each stack")
    p.add_option("--minsamp", default=3, type="int", help="Min number of samples")


@memoized
def binom_consens(n1, n2, E, H):
    """
    Given two bases are observed at a site n1 and n2, and the error rate E, the
    probability the site is aa, bb, ab is calculated using binomial distribution
    as in Li_et al 2009, 2011, and 500 reads were randomly if high coverage.
    """
    prior_homo = (1 - H) / 2.
    prior_het = H
    ab = scipy.misc.comb(n1 + n2, n1) / (2. ** (n1 + n2))
    aa = scipy.stats.binom.pmf(n1, n1 + n2, E)
    bb = scipy.stats.binom.pmf(n2, n1 + n2, E)
    Q = [prior_homo * aa, prior_homo * bb, prior_het * ab]
    Qn = ['aa', 'bb', 'ab']
    P = max(Q) / sum(Q)
    return P, Qn[Q.index(max(Q))]


def hetero(n1,n2):
    """
    Returns IUPAC symbol for ambiguity bases, used for polymorphic sites.
    """
    D = {('G', 'A'): 'R',
         ('G', 'T'): 'K',
         ('G', 'C'): 'S',
         ('T', 'C'): 'Y',
         ('T', 'A'): 'W',
         ('C', 'A'): 'M'}
    a = D.get((n1, n2))
    b = D.get((n2, n1))
    return a or b


def unhetero(amb):
    """
    Returns bases from ambiguity code.
    """
    amb = amb.upper()
    D = {'R': ('G', 'A'),
         'K': ('G', 'T'),
         'S': ('G', 'C'),
         'Y': ('T', 'C'),
         'W': ('T', 'A'),
         'M': ('C', 'A')}
    return D.get(amb, amb)


def uplow(b):
    """
    Precedence G > T > C > A
    """
    D = {('G', 'A'): 'G',
         ('A', 'G'): 'G',
         ('G', 'T'): 'G',
         ('T', 'G'): 'G',
         ('G', 'C'): 'G',
         ('C', 'G'): 'G',
         ('T', 'C'): 'T',
         ('C', 'T'): 'T',
         ('T', 'A'): 'T',
         ('A', 'T'): 'T',
         ('C', 'A'): 'C',
         ('A', 'C'): 'C'}
    r = D.get(b)
    if not r:
        r = b[0]
    return r


def findalleles(consensus, alleles, AL):
    cons = list(consensus)
    bigbase = uplow(tuple([i.split("_")[0] for i in AL]))
    bigallele = AL.index([i for i in AL if i.split("_")[0] == bigbase][0])
    for k in range(1, len(alleles)):
        c = uplow(tuple([i.split("_")[k] for i in AL]))
        which = AL.index([i for i in AL if i.split("_")[k] == c][0])
        if AL[bigallele] != AL[which]:
            cons[alleles[k]] = cons[alleles[k]].lower()

    return "".join(cons)


def breakalleles(consensus):
    """
    Break ambiguity code consensus seqs into two alleles.
    """
    a1 = ""
    a2 = ""
    for base in consensus:
        if base in AMB:
            a, b = unhetero(base)
        elif base in AMBL:
            b, a = unhetero(base)
        else:
            a = b = base
        a1 += a
        a2 += b
    return a1, a2


def find_pctid(consensusfiles):
    pctid = min([int(op.basename(x).split(".")[-2].replace("P", "")) \
            for x in consensusfiles])
    logging.debug("Set pctid={0}".format(pctid))
    return pctid


def mcluster(args):
    """
    %prog mcluster *.consensus

    Cluster across samples using consensus sequences.
    """
    p = OptionParser(mcluster.__doc__)
    add_consensus_options(p)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    consensusfiles = args
    minlength = opts.minlength
    cpus = opts.cpus
    pf = opts.prefix
    pctid = find_pctid(consensusfiles)

    pf += ".P{0}".format(pctid)
    consensusfile = pf + ".consensus.fasta"
    haplotypefile = pf + ".haplotype.fasta"
    if need_update(consensusfiles, (consensusfile, haplotypefile)):
        fw_cons = must_open(consensusfile, "w")
        fw_haps = must_open(haplotypefile, "w")
        totalseqs = 0
        for cf in consensusfiles:
            nseqs = 0
            s = op.basename(cf).split(".")[0]
            for name, seq in parse_fasta(cf):
                name = '.'.join((s, name))
                a1, a2 = breakalleles(seq)
                print >> fw_cons, ">{0}\n{1}".format(name, seq)
                print >> fw_haps, ">{0}\n{1}".format(name, a1)
                nseqs += 1
            logging.debug("Read `{0}`: {1} seqs".format(cf, nseqs))
            totalseqs += nseqs
        logging.debug("Total: {0} seqs".format(totalseqs))
        fw_cons.close()
        fw_haps.close()

    userfile = pf + ".u"
    notmatchedfile = pf + ".notmatched"
    if need_update(haplotypefile, userfile):
        cluster_smallmem(haplotypefile, userfile, notmatchedfile,
                         minlength, pctid, cpus)

    clustfile = pf + ".clust"
    if need_update((consensusfile, userfile, notmatchedfile), clustfile):
        makeclust(consensusfile, userfile, notmatchedfile, clustfile)

    clustSfile = pf + ".clustS"
    if need_update(clustfile, clustSfile):
        parallel_musclewrap(clustfile, cpus, minsamp=opts.minsamp)


def makeloci(clustSfile, store, prefix):
    C = ClustFile(clustSfile)
    pf = clustSfile.rsplit(".", 1)[0]
    locifile = pf + ".loci"
    finalfastafile = pf + ".final.fasta"
    fw = open(locifile, "w")
    fw_finalfasta = open(finalfastafile, "w")
    locid = 0
    AC = []
    for data in C:
        names, seqs, nreps = zip(*data)
        # Strip off cut site
        seqs = [x.upper() for x in seqs]
        fname = "{0}_{1}".format(prefix, locid)

        # TODO: apply number of shared heteros paralog filter

        # Record variable sites
        seed_seq = seqs[0]
        ncols = len(seed_seq)
        snpsite = [' '] * ncols
        seed_ungapped_pos = []
        ref_alleles = []
        alt_alleles = []
        ungapped_i = 0
        for i in xrange(ncols):
            r = seed_seq[i]
            ref_allele = unhetero(r)[0]
            ref_alleles.append(ref_allele)
            seed_ungapped_pos.append(ungapped_i)
            if r in '_-':                # Skip if reference is a deletion
                alt_alleles.append([])
                continue
            else:
                ungapped_i += 1

            site = [s[i] for s in seqs]   # Column slice in MSA
            reals = []
            for x in site:
                if x not in "N_-":
                    reals.extend(unhetero(x))
            realcounts = sorted([(reals.count(x), x) for x in BASES], reverse=True)
            altcount = realcounts[1][0]
            if altcount > 1:
                snpsite[i] = '*'
            elif altcount == 1:
                snpsite[i] = '-'
            nonzeros = [x for c, x in realcounts if (c and x != ref_allele)]
            alt_alleles.append(nonzeros)

        assert len(seed_ungapped_pos) == ncols
        assert len(ref_alleles) == ncols
        assert len(alt_alleles) == ncols

        for name, seq in zip(names, seqs):
            name = name.strip(">")
            taxon, readname = name.split(".", 1)
            profile = store[taxon][readname]
            assert len(seq) == ncols

            ungapped_i = 0
            gap_p = [0, 0, 0, 0, 0, 0, sum(profile[0])]
            for pos, ref_allele, alt_allele, r, ispoly in \
                zip(seed_ungapped_pos, ref_alleles, alt_alleles, seq, snpsite):
                if r in '_-':     # insertion in ref, deletion in read
                    p = gap_p
                else:
                    p = profile[ungapped_i]
                    ungapped_i += 1

                if ispoly != '*':
                    continue
                ac = AlleleCount(taxon, fname, pos + 1,  # 1-based coordinate
                                 ref_allele, alt_allele, p)
                AC.append(ac)

        longname = max(len(x) for x in names)
        longname = max(len(fname) + 3, longname) + 1
        for name, seq, nrep in data:
            print >> fw, name.ljust(longname) + seq
        print >> fw, "// {0}".format(fname).ljust(longname) + "".join(snpsite) + "|"

        seed_seq = seed_seq.strip("_N").replace("-", "")
        print >> fw_finalfasta, ">{0} with {1} sequences\n{2}".\
                    format(fname, sum(nreps), seed_seq)
        locid += 1

    logging.debug("Stacks written to `{0}`".format(locifile))
    logging.debug("Final consensus sequences written to `{0}` (n={1})".\
                    format(finalfastafile, locid))
    fw.close()
    fw_finalfasta.close()

    return AC


def mconsensus(args):
    """
    %prog mconsensus *.consensus

    Call consensus along the stacks from cross-sample clustering.
    """
    p = OptionParser(mconsensus.__doc__)
    p.add_option("--allele_counts", default="allele_counts",
                 help="Directory to generate allele counts")
    add_consensus_options(p)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    consensusfiles = args
    prefix = opts.prefix
    acdir = opts.allele_counts
    store = ClustStores(consensusfiles)
    pctid = find_pctid(consensusfiles)
    pf = prefix + ".P{0}".format(pctid)

    clustSfile = pf + ".clustS"
    AC = makeloci(clustSfile, store, prefix)

    mkdir(acdir)
    acfile = pf + ".allele_counts"
    fw = open(acfile, "w")
    seen = DefaultOrderedDict(list)        # chr, pos => taxa
    print >> fw, "# " + "\t".join(ACHEADER)
    # Sort allele counts into separate files
    for ac in AC:
        chrpos = ac.chr, ac.pos
        seen[chrpos].append(ac)
        print >> fw, ac.tostring(taxon=True)
    fw.close()

    logging.debug("Populate all taxa and instantiate empty vector if missing")
    all_taxa = set([op.basename(x).split(".")[0] for x in consensusfiles])
    taxon_to_ac = defaultdict(list)
    for chrpos, aclist in seen.items():
        included_taxa = set([x.taxon for x in aclist])
        missing_taxa = all_taxa - included_taxa
        template = deepcopy(aclist[0])
        template.clear()
        for ac in aclist:
            taxon_to_ac[ac.taxon].append(ac)
        for tx in missing_taxa:
            taxon_to_ac[tx].append(template)

    logging.debug("Write allele counts for all taxa")
    for tx, aclist in sorted(taxon_to_ac.items()):
        tx_acfile = op.join(acdir, tx + ".allele_counts")
        fw = open(tx_acfile, "w")
        print >> fw, "# " + "\t".join(ACHEADER_NO_TAXON)
        for ac in aclist:
            print >> fw, ac.tostring()
        fw.close()
        logging.debug("Written {0} sites in `{1}`".\
                format(len(aclist), tx_acfile))


def consensus(args):
    """
    %prog consensus clustSfile

    Call consensus along the stacks. Tabulate bases at each site, tests for
    errors according to error rate, calls consensus. HEfile contains the
    heterozygosity and error rate as calculated by estimateHE().
    """
    p = OptionParser(consensus.__doc__)
    p.add_option("--ploidy", default=2, type="int",
                 help="Number of haplotypes per locus")
    p.add_option("--estimate_errors", default=False, action="store_true",
                 help="Estimate H and E from data")
    add_consensus_options(p)
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clustSfile, = args
    pf = clustSfile.rsplit(".", 1)[0]
    HEfile = pf + ".HE"
    mindepth = opts.mindepth

    H, E = .01, .001
    if opts.estimate_errors:
        try:
            HEfile = estimateHE([clustSfile])
            H, E = open(HEfile).readline().split()
            H, E = float(H), float(E)
        except:
            pass
    logging.debug("H={0} E={1}".format(H, E))

    bases = BASES[:4]
    C = ClustFile(clustSfile)
    output = []
    bins = []
    indices = []
    start = end = 0  # Index into base count array
    for data in C:
        names, seqs, nreps = zip(*data)
        total_nreps = sum(nreps)
        # Depth filter
        if total_nreps < mindepth:
            continue

        name, seq, nrep = data[0]
        fname = name.split(";")[0] + ";size={0};".format(total_nreps)
        cons_name, cons_seq, cons_nrep = data[-1]
        if not cons_name.startswith("CONSENS"):
            cons_name, cons_seq, cons_nrep = name, seq, nrep

        # List for sequence data
        S = [(seq, nrep) for name, seq, nrep in data if nrep]
        # Pileups for base counting
        RAD = stack(S)

        if len(data) == 1:   # No computation needed
            output.append((fname, seq))
            bins.extend(RAD)
            start = end
            end += len(seq)
            indices.append((fname, start, end))
            continue

        # Strip N's from either end and gaps
        leftjust, rightjust = get_left_right(cons_seq)
        terminalNs = set(range(leftjust) + range(rightjust + 1, len(cons_seq)))
        gaps = set(i for i, j in enumerate(cons_seq) if j == '-')

        filtered = terminalNs | gaps
        shortcon = "".join(j for (i, j) in enumerate(cons_seq) \
                            if i not in filtered)
        shortRAD = [j for (i, j) in enumerate(RAD) if i not in filtered]
        conlen = len(shortcon)
        assert conlen == len(shortRAD)

        # Correct consensus by converting to top voting bases
        corr_shortcon = ""
        for base, site in zip(shortcon, shortRAD):
            site = site[:4]
            n1 = max(site)  # Base with highest count
            # Check count for original base for possible ties
            n0 = site[BASES.index(base)]
            corr_base = bases[site.index(n1)] if n1 > n0 else base
            corr_shortcon += corr_base

        if opts.verbose:
            print_list = lambda L: ",".join(str(x) for x in sorted(L))
            print fname
            print "\n".join(["{0} {1}".format(*x) for x in S])
            print cons_seq
            print "Terminal Ns:", print_list(terminalNs)
            print "Gaps:", print_list(gaps)
            print shortcon
            diff = "".join((" " if a == b else "*") \
                        for a, b in zip(shortcon, corr_shortcon))
            if diff.strip():
                print corr_shortcon
                print diff
            print "|".join(["{0}{1}:{2}".\
                        format(i, corr_shortcon[i], " ".join(str(x) for x in j)) \
                        for i, j in enumerate(shortRAD)])
            print "-" * 60

        output.append((fname, corr_shortcon))
        bins.extend(shortRAD)

        start = end
        end += len(shortcon)
        indices.append((fname, start, end))

    consensfile = pf + ".consensus"
    consens = open(consensfile, 'w')
    for k, v in output:
        print >> consens, "\n".join((k, v))
    consens.close()
    logging.debug("Consensus sequences written to `{0}`".format(consensfile))

    binfile = consensfile + ".bin"
    bins = np.array(bins, dtype=np.uint32)
    ulimit = 65535
    bins[bins > ulimit] = ulimit
    bins = np.array(bins, dtype=np.uint16)  # Compact size
    bins.tofile(binfile)
    logging.debug("Allele counts written to `{0}`".format(binfile))

    idxfile = consensfile + ".idx"
    fw = open(idxfile, "w")
    for fname, start, end in indices:
        print >> fw, "\t".join(str(x) for x in (fname, start, end))
    fw.close()
    logging.debug("Serializing indices to `{0}`".format(idxfile))

    return consensfile, binfile, idxfile


def stack(S):
    """
    From list of bases at a site D,  make counts of bases
    """
    S, nreps = zip(*S)
    S = np.array([list(x) for x in S])
    rows, cols = S.shape
    counts = []
    for c in xrange(cols):
        freq = [0] * NBASES
        for b, nrep in zip(S[:, c], nreps):
            freq[BASES.index(b)] += nrep
        counts.append(freq)
    return counts


def get_left_right(seq):
    """
    Find position of the first and last base
    """
    cseq = seq.strip("_-N")
    leftjust = seq.index(cseq[0])
    rightjust = seq.rindex(cseq[-1])

    return leftjust, rightjust


def cons(f, mindepth):
    """
    Makes a list of lists of reads at each site
    """
    C = ClustFile(f)
    for data in C:
        names, seqs, nreps = zip(*data)
        total_nreps = sum(nreps)
        # Depth filter
        if total_nreps < mindepth:
            continue

        S = []
        for name, seq, nrep in data:
            # Append sequence * number of dereps
            S.append([seq, nrep])

        # Make list for each site in sequences
        res = stack(S)
        yield [x[:4] for x in res if sum(x[:4]) >= mindepth]


def makeP(N):
    # Make list of freq. for BASES
    sump = float(sum([sum(i) for i in N]))
    if sump:
        p1 = sum([i[0] for i in N]) / sump
        p2 = sum([i[1] for i in N]) / sump
        p3 = sum([i[2] for i in N]) / sump
        p4 = sum([i[3] for i in N]) / sump
    else:
        p1 = p2 = p3 = p4 = 0.0
    return [p1, p2, p3, p4]


def makeC(N):
    """
    Makes a dictionary with counts of base counts [x,x,x,x]:x,
    speeds up Likelihood calculation
    """
    C = defaultdict(int)
    for d in N:
        C[tuple(d)] += 1

    return [i for i in C.items() if (0, 0, 0, 0) not in i]


def L1(E, P, N):
    # Probability of homozygous
    h = []
    s = sum(N)
    for i, l in enumerate(N):
        p = P[i]
        b = scipy.stats.binom.pmf(s - l, s, E)
        h.append(p * b)
    return sum(h)


def L2(E, P, N):
    # Probability of heterozygous
    h = []
    s = sum(N)
    for l, i in enumerate(N):
        for j, k in enumerate(N):
            if j > l:
                one = 2. * P[l] * P[j]
                two = scipy.stats.binom.pmf(s - i - k, s, (2. * E) / 3.)
                three = scipy.stats.binom.pmf(i, k + i, 0.5)
                four = 1. - (sum([q ** 2. for q in P]))
                h.append(one * two * (three / four))
    return sum(h)


def totlik(E, P, H, N):
    # Total probability
    lik = ((1 - H) * L1(E, P, N)) + (H * L2(E, P, N))
    return lik


def LL(x0, P, C):
    # Log likelihood score given values [H, E]
    H, E = x0
    L = []
    if H <= 0. or E <= 0.:
        r = np.exp(100)
    else:
        for i in C:
            ll = totlik(E, P, H, i[0])
            if ll > 0:
                L.append(i[1] * np.log(ll))
        r = -sum(L)
    return r


def estimateHE(args):
    """
    %prog estimateHE clustSfile

    Estimate heterozygosity (H) and error rate (E). Idea borrowed heavily from
    the PyRad paper.
    """
    p = OptionParser(estimateHE.__doc__)
    add_consensus_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clustSfile, = args
    HEfile = clustSfile.rsplit(".", 1)[0] + ".HE"
    if not need_update(clustSfile, HEfile):
        logging.debug("File `{0}` found. Computation skipped.".format(HEfile))
        return HEfile

    D = []
    for d in cons(clustSfile, opts.mindepth):
        D.extend(d)

    logging.debug("Computing base frequencies ...")
    P = makeP(D)
    C = makeC(D)
    logging.debug("Solving log-likelihood function ...")
    x0 = [.01, .001]  # initital values
    H, E = scipy.optimize.fmin(LL, x0, args=(P, C))

    fw = must_open(HEfile, "w")
    print >> fw, H, E
    fw.close()

    return HEfile


def alignfast(names, seqs):
    """
    Performs MUSCLE alignments on cluster and returns output as string
    """
    matfile = op.join(datadir, "blosum80.mat")
    cmd = "poa -read_fasta - -pir stdout {0} -tolower -silent -hb".format(matfile)
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    s = ""
    for i, j in zip(names, seqs):
        s += "\n".join((i, j)) + "\n"
    return p.communicate(s)[0]


def replace_terminal(seq):
    leftjust, rightjust = get_left_right(seq)
    seq = "_" * leftjust + seq[leftjust: rightjust + 1] \
                + "_" * (len(seq) - rightjust - 1)
    return seq


def sortalign(stringnames):
    G = stringnames.split("\n>")
    aligned = [('>' + i.split("\n")[0].strip('>'),
               replace_terminal("".join(i.split("\n")[1:]).upper())) for i in G]
    return aligned


def parallel_musclewrap(clustfile, cpus, minsamp=0):
    musclewrap_minsamp = partial(musclewrap, minsamp=minsamp)
    if cpus == 1:
        return musclewrap_minsamp(clustfile)

    from jcvi.apps.grid import Jobs

    outdir = mkdtemp(dir=".")
    fs = split([clustfile, outdir, str(cpus), "--format=clust"])
    g = Jobs(musclewrap_minsamp, fs.names)
    g.run()

    clustnames = [x.replace(".clust", ".clustS") for x in fs.names]
    clustSfile = clustfile.replace(".clust", ".clustS")
    FileMerger(clustnames, outfile=clustSfile).merge()
    shutil.rmtree(outdir)


def filter_samples(names, seqs, sep='.'):
    """
    When there are uncollapsed contigs within the same sample, only retain the
    first seq, or the seq that is most abundant (with cluster_size).
    """
    seen = set()
    filtered_names, filtered_seqs = [], []
    for name, seq in zip(names, seqs):
        samp = name.split(sep, 1)[0]
        if samp in seen:
            continue
        seen.add(samp)
        filtered_names.append(name)
        filtered_seqs.append(seq)

    nfiltered, nnames = len(filtered_names), len(names)
    assert nfiltered == len(seen)

    return filtered_names, filtered_seqs, seen


def musclewrap(clustfile, minsamp=0):
    cnts = 0
    C = ClustFile(clustfile)
    clustSfile = clustfile.replace(".clust", ".clustS")
    fw = open(clustSfile, 'w')
    for data in C:
        STACK = Clust()
        names = []
        seqs = []
        names, seqs, nreps = zip(*data)
        if minsamp:  # Filter based on samples, applicable in mcluster()
            names, seqs, samples = filter_samples(names, seqs)
            if len(samples) < minsamp:
                continue
        else:
            names, seqs = names[:255], seqs[:255]  # Reduce high coverage data

        if len(names) == 1:
            STACK.append((names[0], seqs[0]))
        else:
            stringnames = alignfast(names, seqs)
            aligned = sortalign(stringnames)
            # Reorder keys by derep number
            D1 = [(getsize(name), name, seq) for name, seq in aligned]
            D1.sort(key=lambda x: (-x[0], x[1]))
            for size, name, seq in D1:
                STACK.append((name, seq))

        if STACK:
            print >> fw, STACK
        cnts += 1

    fw.close()


def stats(clustSfile, statsfile, mindepth):
    C = ClustFile(clustSfile)
    depth = []
    for data in C:
        d = 0
        for name, seq, nrep in data:
            d += nrep
        depth.append(d)
    keep = [i for i in depth if i >= mindepth]
    namecheck = op.basename(clustSfile).split(".")[0]
    if depth:
        me = round(np.mean(depth), 3)
        std = round(np.std(depth), 3)
    else:
        me = std = 0.0
    if keep:
        mek = round(np.mean(keep), 3)
        stdk = round(np.std(keep), 3)
    else:
        mek = stdk = 0.0
    out = dict(label=namecheck, cnts=len(depth), mean=me, std=std,
               keep=len(keep), meank=mek, stdk=stdk)
    header = "label cnts mean std keep meank stdk".split()

    bins = [0, 5, 10, 15, 20, 25, 30, 35, 40, 50, 100, 250, 500, 99999]
    ohist, edges = np.histogram(depth, bins)
    hist = [float(i) / sum(ohist) for i in ohist]
    hist = [int(round(i * 30)) for i in hist]

    logging.debug("Sample {0} finished, {1} loci".\
                    format(clustSfile, len(depth)))

    fw = open(statsfile,'w')
    print >> fw, "# Params: mindepth={0}".format(mindepth)
    print >> fw, " ".join("{0}={1}".format(k, out[k]) for k in header)
    print >> fw, "\nbins\tdepth_histogram\tcnts"
    print >> fw, "   :\t0------------50-------------100%"

    for i, j, k in zip(edges, hist, ohist):
        firststar = " "
        if k > 0:
            firststar = "*"
        print >> fw, i,'\t', firststar + "*" * j + " " * (34 - j), k
    fw.close()


def makeclust(derepfile, userfile, notmatchedfile, clustfile,
              mindepth=3, include=.95):
    D = dict(parse_fasta(derepfile))
    U = defaultdict(list)  # Clusters
    fp = open(userfile)
    for row in fp:
        query, target, id, qcov, tcov = row.rstrip().split("\t")
        U[target].append((query, getsize(query),
                          float(id) * float(qcov) * float(tcov)))

    fw = open(clustfile, "w")
    for key, members in U.items():
        keysize = getsize(key)
        members.sort(key=lambda x: (-x[1], -x[2]))
        totalsize = keysize + sum(x[1] for x in members)
        cutoff = int(round(totalsize * include))
        cumulative = keysize

        # Remove outliers within each cluster
        seqs = [('>' + key, D[key])]
        for name, size, id in members:
            seqs.append(('>' + name, D[name]))
            cumulative += size
            if cumulative >= cutoff:
                break
        if cumulative < mindepth:
            continue
        seq = "\n".join("\n".join(x) for x in seqs)
        print >> fw, "\n".join((seq, SEP))

    I = dict(parse_fasta(notmatchedfile))
    singletons = set(I.keys()) - set(U.keys())
    for key in singletons:
        if getsize(key) < mindepth:
            continue
        print >> fw, "\n".join(('>' + key, I[key], SEP))
    fw.close()


def derep(fastafile, derepfile, minlength, cpus, usearch="vsearch"):
    cmd = usearch + " -minseqlength {0}".format(minlength)
    cmd += " -derep_fulllength {0}".format(fastafile)
    cmd += " -output {0} -sizeout".format(derepfile)
    cmd += " -threads {0}".format(cpus)
    sh(cmd)


def cluster_smallmem(derepfile, userfile, notmatchedfile, minlength, pctid,
                     cpus, usearch="vsearch"):
    identity = pctid / 100.
    cmd = usearch + " -minseqlength {0}".format(minlength)
    cmd += " -leftjust"
    cmd += " -cluster_size {0}".format(derepfile)
    cmd += " -id {0}".format(identity)
    cmd += " -query_cov {0}".format(identity)
    cmd += " -userout {0}".format(userfile)
    cmd += " -userfields query+target+id+qcov+tcov"
    cmd += " -maxaccepts 1 -maxrejects 8"  # Decrease maxrejects for speed
    cmd += " -usersort -sizein"
    cmd += " -notmatched {0}".format(notmatchedfile)
    cmd += " -threads {0}".format(cpus)
    sh(cmd)


def cluster(args):
    """
    %prog cluster prefix fastqfiles

    Use `vsearch` to remove duplicate reads. This routine is heavily influenced
    by PyRAD: <https://github.com/dereneaton/pyrad>.
    """
    p = OptionParser(cluster.__doc__)
    add_consensus_options(p)
    p.set_align(pctid=94)
    p.set_outdir()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    prefix = args[0]
    fastqfiles = args[1:]
    cpus = opts.cpus
    pctid = opts.pctid
    mindepth = opts.mindepth
    minlength = opts.minlength
    fastafile, qualfile = fasta(fastqfiles + ["--seqtk",
                                "--outdir={0}".format(opts.outdir),
                                "--outfile={0}".format(prefix + ".fasta")])

    prefix = op.join(opts.outdir, prefix)
    pf = prefix + ".P{0}".format(pctid)
    derepfile = prefix + ".derep"
    if need_update(fastafile, derepfile):
        derep(fastafile, derepfile, minlength, cpus)

    userfile = pf + ".u"
    notmatchedfile = pf + ".notmatched"
    if need_update(derepfile, userfile):
        cluster_smallmem(derepfile, userfile, notmatchedfile,
                         minlength, pctid, cpus)

    clustfile = pf + ".clust"
    if need_update((derepfile, userfile, notmatchedfile), clustfile):
        makeclust(derepfile, userfile, notmatchedfile, clustfile,
                  mindepth=mindepth)

    clustSfile = pf + ".clustS"
    if need_update(clustfile, clustSfile):
        parallel_musclewrap(clustfile, cpus)

    statsfile = pf + ".stats"
    if need_update(clustSfile, statsfile):
        stats(clustSfile, statsfile, mindepth=mindepth)


def align(args):
    """
    %prog align clustfile

    Align clustfile to clustSfile. Useful for benchmarking aligners.
    """
    p = OptionParser(align.__doc__)
    p.set_cpus()
    add_consensus_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clustfile, = args
    parallel_musclewrap(clustfile, opts.cpus)


if __name__ == '__main__':
    main()
