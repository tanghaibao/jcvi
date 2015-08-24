#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Using VCLUST to derep, cluster, and make consensus from duplicate reads.
The VCLUST implementation borrows ideas and code from PyRAD. PyRAD link:

<https://github.com/dereneaton/pyrad>
"""

import os.path as op
import sys
import logging
import numpy as np
import scipy
import scipy.stats
import scipy.optimize

from collections import defaultdict
from itertools import groupby
from random import sample
from subprocess import Popen, PIPE, STDOUT

from jcvi.formats.base import BaseFile, must_open
from jcvi.formats.fasta import parse_fasta
from jcvi.formats.fastq import fasta
from jcvi.utils.cbook import memoized
from jcvi.utils.iter import grouper
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


SEP = "//"


class ClustFile (BaseFile):

    def __init__(self, filename):
        super(ClustFile, self).__init__(filename)

    def iter_seqs(self):
        nstacks = 0
        fp = must_open(self.filename)
        for tag, contents in groupby(fp, lambda row: row[0] == '/'):
            if tag:
                continue
            data = []
            for name, seq in grouper(contents, 2):
                name, seq = name.strip(), seq.strip()
                nrep = int(name.split(";")[1].replace("size=", ""))
                data.append((name, seq, nrep))
            yield data
            nstacks += 1
            if nstacks % 1000 == 0:
                logging.debug("{0} stacks parsed".format(nstacks))


def main():

    actions = (
        # UCLUST/VCLUST related
        ('estimateHE', 'estimate heterozygosity and error rate for stacks'),
        ('cluster', 'cluster within samples'),
        ('consensus', 'call consensus bases within samples'),
        ('mcluster', 'cluster across samples'),
        ('mconsensus', 'call consensus bases across samples'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def add_consensus_options(p):
    p.add_option("--minlength", default=30, type="int", help="Min contig length")
    p.add_option("--mindepth", default=3, type="int", help="Min depth for each stack")
    p.add_option("--cut", default=4, type="int",
                 help="Size of restriction site to trim on the left, e.g. CATG is 4")


@memoized
def binom_consens(n1, n2, E, H):
    """
    Given two bases are observed at a site n1 and n2, and the error rate E, the
    probability the site is aa, bb, ab is calculated using binomial distribution
    as in Li_et al 2009, 2011, and if coverage > 500, 500 reads were randomly
    sampled.
    """
    maf = n1 / (n1 + n2)
    prior_homo = (1 - H) / 2.
    prior_het = H
    ab = scipy.misc.comb(n1 + n2, n1) / (2. ** (n1 + n2))
    aa = scipy.stats.binom.pmf(n1, n1 + n2, E)
    bb = scipy.stats.binom.pmf(n2, n1 + n2, E)
    Q = [prior_homo * aa, prior_homo * bb, prior_het * ab]
    Qn = ['aa', 'bb', 'ab']
    P = max(Q) / sum(Q)
    return P, maf, Qn[Q.index(max(Q))]


def naive_consens(n1, n2):
    """
    Majority consensus calling for sites with too low of coverage for
    statistical calling. Only used with 'lowcounts' option.
    """
    maf = n1 * 1. / (n1 + n2)
    return [1.0, maf, 'aa']


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
    D = {'R': ('G', 'A'),
         'K': ('G', 'T'),
         'S': ('G', 'C'),
         'Y': ('T', 'C'),
         'W': ('T', 'A'),
         'M': ('C', 'A')}
    return D.get(amb.upper())


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
        if base in tuple("RKSYWM"):
            a, b = unhetero(base)
        elif base in tuple("rksywm"):
            b, a = unhetero(base)
        else:
            a = b = base
        a1 += a
        a2 += b
    return a1, a2


def removerepeat_Ns(shortcon):
    """
    Checks for interior Ns in consensus seqs remove those that arise next to
    *single repeats* of at least 3 bases on either side, which may be
    sequencing errors on deep coverage repeats
    """
    Nlocs = [i for i, j in enumerate(shortcon) if j == 'N']
    repeats = set()
    for n in Nlocs:
        r1 = len(set(list(shortcon)[n - 3: n]))
        if r1 < 2:
            repeats.add(n)
        r2 = len(set(list(shortcon)[n + 1: n + 4]))
        if r2 < 2:
            repeats.add(n)
    return "".join([j for (i, j) in enumerate(shortcon) if i not in repeats])


def mcluster(args):
    """
    %prog mcluster *.consensus

    Cluster across samples using consensus sequences.
    """
    p = OptionParser(mcluster.__doc__)
    add_consensus_options(p)
    p.set_align(pctid=96)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    consensusfiles = args
    minlength = opts.minlength
    identity = opts.pctid / 100.
    cpus = opts.cpus

    fastafile = "mcluster.fasta"
    if need_update(consensusfiles, fastafile):
        fw = must_open(fastafile, "w")
        totalseqs = 0
        for consensusfile in consensusfiles:
            nseqs = 0
            s = op.basename(consensusfile).split(".")[0]
            for name, seq in parse_fasta(consensusfile, upper=False):
                a1, a2 = breakalleles(seq)
                print >> fw, ">{0}-{1}\n{2}".format(s, name, a1)
                nseqs += 1
            logging.debug("Read `{0}`: {1} seqs".format(consensusfile, nseqs))
            totalseqs += nseqs
        logging.debug("Total: {0} seqs".format(totalseqs))
        fw.close()

    pf = fastafile
    userfile = pf + ".u"
    notmatchedfile = pf + ".notmatched"
    if need_update(fastafile, userfile):
        cluster_smallmem(fastafile, userfile, notmatchedfile,
                         minlength, identity, cpus)

    clustfile = pf + ".clust"
    if need_update((fastafile, userfile, notmatchedfile), clustfile):
        makeclust(fastafile, userfile, notmatchedfile, clustfile)

    clustSfile = pf + ".clustS"
    if need_update(clustfile, clustSfile):
        musclewrap(clustfile)


def makealign(clustSfile, locifile, CUT1):
    C = ClustFile(clustSfile)
    fw = open(locifile, "w")
    for data in C.iter_seqs():
        names, seqs, nreps = zip(*data)
        # Strip off cut site
        seqs = [x[CUT1:] for x in seqs]
        longname = max(len(x) for x in names) + 2

        # Apply number of shared heteros paralog filter

        # Record variable sites
        ncols = len(seqs[0])
        snpsite = [' '] * ncols

        for i in xrange(ncols):
            site = [s[i] for s in seqs]
            reals = [x for x in site if x not in "N-"]
            if len(set(reals)) <= 1:
                continue

            realcounts = sorted([reals.count(x) for x in set(reals)],
                                 reverse=True)
            snpsite[i] = '*' if realcounts[1] > 1 else '-'

        for name, seq, nrep in data:
            print >> fw, name.ljust(longname) + seq
        print >> fw, "//".ljust(longname + CUT1) + "".join(snpsite) + "|"

    logging.debug("Stacks written to `{0}`".format(locifile))
    fw.close()


def mconsensus(args):
    """
    %prog mconsensus clustSfile

    Call consensus along the stacks from cross-sample clustering.
    """
    p = OptionParser(mconsensus.__doc__)
    add_consensus_options(p)
    p.add_option("--minsamp", default=3, type="int",
                 help="Min number of samples")
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clustSfile, = args
    locifile = clustSfile.rsplit(".", 1)[0] + ".loci"
    if need_update(clustSfile, locifile):
        makealign(clustSfile, locifile, opts.cut)


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
    p.add_option("--maxH", default=10, type="int",
                 help="Max number of heterozygous sites allowed")
    p.add_option("--maxN", default=10, type="int",
                 help="Max number of Ns allowed")
    add_consensus_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    clustSfile, = args
    HEfile = clustSfile.rsplit(".", 1)[0] + ".HE"
    mindepth = opts.mindepth
    haplos = opts.ploidy
    maxH = opts.maxH
    maxN = opts.maxN

    HEfile = estimateHE([clustSfile])

    H, E = open(HEfile).readline().split()
    try:
        H, E = float(H), float(E)
    except:
        H, E = .01, .001
    logging.debug("H={0} E={1}".format(H, E))

    bases = "ACTG"
    locus = minsamplocus = npoly = P = 0
    C = ClustFile(clustSfile)
    output = []
    for data in C.iter_seqs():
        names, seqs, nreps = zip(*data)
        name, seq, nrep = data[0]
        fname = name.split(";")[0] + ";size={0};".format(sum(nreps))
        if len(data) == 1:
            if nrep >= mindepth:   # Same thing
                output.append((fname, seq))
            continue

        S = []               # List for sequence data
        alleles = []         # Measures # alleles, detect paralogs
        locus += 1           # Measures # of loci
        nHs = 0              # Measures heterozygous sites in this locus
        cons_seq = ""        # Consensus sequence
        basenumber = 0       # Tracks het locations
        rights = []
        leftjust = rightjust = None
        for name, seq, nrep in data:
            # Append sequence * number of dereps
            for i in xrange(nrep):
                S.append(tuple(seq.strip()))

            # Record left and right most index of seed and hits (for GBS)
            # leftjust is seed's left, rightjust is the shortest reverse hit
            if name[-1] == ";":
                leftjust = seq.index([i for i in seq if i not in list("-N")][0])

            if name[-1] == "-":
                rights.append(max(-1, [seq.rindex(i) for i in seq if i in bases]))

            # Trim off overhang edges of gbs reads
            if rights:
                # Record in name that there was a reverse hit
                fname = "_".join(fname.split("_")[0:-1]) + "_c1"
                try:
                    rightjust = min([min(i) for i in rights])
                except ValueError:
                    S = ""

            for s in xrange(len(S)):
                S[s] = S[s][leftjust:]
                if rightjust:
                    S[s] = S[s][:rightjust + 1]

        # Apply paralog filters, depth filter disabled
        if len(S) < mindepth:
            continue

        minsamplocus += 1
        RAD = stack(S)
        paralog = False
        for site in RAD:
            site, Ns, gaps = site

            # Minimum depth of coverage for base calling
            depthofcoverage = sum(site)
            if depthofcoverage < gaps:
                cons = '-'
            elif depthofcoverage < mindepth:
                cons = 'N'
                n1 = depthofcoverage - 1
                n2 = 0   # Prevents zero division error
            else:
                n1, n2, n3, n4 = sorted(site, reverse=True)

                # Speed hack = if diploid exclude if a third base present at > 20%
                quickthirdbasetest = False
                if haplos == 2:
                    if float(n3) / (n1 + n2 + n3 + n4) > .2:
                        quickthirdbasetest = True

                if quickthirdbasetest:
                    cons = "@"   # Paralog
                else:
                    m1, m2 = n1, n2
                    # For high cov data, reduce for base calling
                    if n1 + n2 >= 500:
                        s = sample('A' * n1 + 'B' * n2, 500)
                        m1, m2 = s.count('A'), s.count('B')

                    # Make base calls, two different methods available:
                    # binom_consens and naive_consens
                    if n1 + n2 >= mindepth:
                        P, maf, who = binom_consens(m1, m2, E, H)

                    # High conf if base could be called with 95% post. prob.
                    if P < .95:
                        cons = 'N'
                    else:
                        if who in 'ab':
                            a = [i for i, l in enumerate(site) if l == n1]
                            if len(a) == 2:       # alleles came up equal freq
                                cons = hetero(bases[a[0]], bases[a[1]])
                            else:                 # alleles came up diff freq
                                b = [i for i, l in enumerate(site) if l == n2]

                                # If three alleles came up equal, only need if diploid paralog filter off
                                if a == b:
                                    cons = hetero(bases[a[0]], bases[a[1]])
                                else:
                                    cons = hetero(bases[a[0]], bases[b[0]])
                            alleles.append(basenumber)
                            nHs += 1
                        else:
                            cons = bases[site.index(n1)]

            cons_seq += cons
            basenumber += 1

            # Only allow maxH polymorphic sites in a locus
            if cons == '@' or nHs > maxH:
                paralog = True
                break

        if paralog:
            continue

        # Filter to limit to N haplotypes
        al = []
        if len(alleles) > 1:
            for s in S:
                d = []
                for z in alleles:
                    if s[z] in unhetero(cons_seq[z]):
                        d.append(s[z])
                if "N" not in d and len(d) == len(alleles):
                    al.append('_'.join(d))

            AL = sorted(set(al), key=al.count)

            # Set correct alleles relative to first polymorphic base
            if AL:  # and len(AL) <= haplos   - TODO: check ploidy level
                cons_seq = findalleles(cons_seq, alleles, AL)

        # strip N's from either end
        shortcon = cons_seq.lstrip("N").rstrip("N").replace("-", "")
        shortcon = removerepeat_Ns(shortcon)

        # Only allow maxN internal "N"s in a locus
        if shortcon.count("N") <= maxN and len(shortcon) >= opts.minlength:
            npoly += nHs
            output.append((fname, shortcon))

    consens = open(clustSfile.replace(".clustS", ".consensus"), 'w+')
    for k, v in output:
        print >> consens, "\n".join((k, v))
    consens.close()

    nsites = sum([len(v) - len(opts.cut) for k, v in output])
    ldic = len(output)
    NP = 0 if not nsites else npoly / float(nsites)

    return [clustSfile.split('/')[-1], locus, minsamplocus, ldic, \
            nsites, npoly, round(NP, 7)]


def stack(D):
    """
    from list of bases at a site D,
    returns an ordered list of counts of bases
    """
    L = len(D)
    counts = []
    for i in range(len(D[0])):
        A = C = T = G = N = S = 0
        for nseq in range(L):
            s = D[nseq][i]
            A += s.count("A")
            C += s.count("C")
            T += s.count("T")
            G += s.count("G")
            N += s.count("N")
            S += s.count("-")
        counts.append([[A, C, T, G], N, S])
    return counts


def cons(f, minsamp, CUT1):
    """ makes a list of lists of reads at each site """
    C = ClustFile(f)
    for data in C.iter_seqs():
        S = []
        rights = []
        lefts = []
        leftjust = rightjust = None
        for name, seq, nrep in data:
            # Record left and right most for cutting
            if name.split(";")[-1] == "":
                leftjust = seq.index([i for i in seq if i not in list("-N")][0])
                rightjust = seq.rindex([i for i in seq if i not in list("-N")][0])
            lefts.append(seq.index([i for i in seq if i not in list("-N")][0]))
            rights.append(seq.rindex([i for i in seq if i not in list("-N")][0]))

            # Append sequence * number of dereps
            for i in range(nrep):
                S.append(tuple(seq))

        # Trim off overhang edges of gbs reads
        if any([i < leftjust for i in lefts]):
            rightjust = min(rights)
        if any([i < rightjust for i in rights]):
            leftjust = max(lefts)

        for s in range(len(S)):
            if leftjust or rightjust:
                S[s] = S[s][leftjust: rightjust + 1]

        # Trim off restriction sites from ends
        for s in range(len(S)):
            S[s] = S[s][CUT1:]

        if len(S) >= minsamp:
            # Make list for each site in sequences
            res = stack(S)
            # Exclude sites with indels
            yield [i[0] for i in res if i[2] == 0]


def makeP(N):
    """ returns a list of freq. for ATGC"""
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
    """ Makes a dictionary with counts of base counts [x,x,x,x]:x,
    speeds up Likelihood calculation"""
    C = defaultdict(int)
    k = iter(N)
    while True:
        try:
            d = k.next()
        except StopIteration:
            break
        C[tuple(d)] += 1

    L = [(i, j) for i, j in C.items()]
    return [i for i in L if (0,0,0,0) not in i]


def L1(E, P, N):
    """probability homozygous"""
    h = []
    s = sum(N)
    for i, l in enumerate(N):
        p = P[i]
        b = scipy.stats.binom.pmf(s - l, s, E)
        h.append(p*b)
    return sum(h)


def L2(E, P, N):
    """probability of heterozygous"""
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
    """ total probability """
    lik = ((1 - H) * L1(E, P, N)) + (H * L2(E, P, N))
    return lik


def LL(x0, P, C):
    """ Log likelihood score given values [H, E] """
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
    for d in cons(clustSfile, opts.mindepth, opts.cut):
        D.extend(d)

    logging.debug("Computing base frequencies ...")
    P = makeP(D)
    logging.debug("Computing base vector counts ...")
    C = makeC(D)
    logging.debug("Solving log-likelihood function ...")
    x0 = [.01, .001]  # initital values
    H, E = scipy.optimize.fmin(LL, x0, args=(P, C))

    fw = must_open(HEfile, "w")
    print >> fw, H, E
    fw.close()

    return HEfile


def alignfast(names, seqs, bigfile=10000):
    """
    Performs MUSCLE alignments on cluster and returns output as string
    """
    cmd = "muscle -quiet -in -"
    p = Popen(cmd, shell=True, stdin=PIPE, stdout=PIPE, stderr=STDOUT, close_fds=True)
    for i, j in zip(names, seqs):
        s = "\n".join((i, j)) + "\n"
        p.stdin.write(s)
    p.stdin.close()
    p.wait()
    return p.stdout.read()


def sortalign(stringnames):
    G = stringnames.split("\n>")
    aligned = [('>' + i.split("\n")[0].strip('>'),
               "".join(i.split("\n")[1:])) for i in G]
    return aligned


def musclewrap(clustfile):
    clustSfile = clustfile.replace(".clust", ".clustS")
    cnts = 0
    C = ClustFile(clustfile)
    fw = open(clustSfile, 'w')
    for data in C.iter_seqs():
        STACK = []
        names = []
        seqs = []
        names, seqs, nreps = zip(*data)
        if len(names) == 1:
            STACK = [names[0] + "\n" + seqs[0]]
        else:
            # Keep only the 200 most common dereps, aligning more is surely junk
            stringnames = alignfast(names[0:200], seqs[0:200])
            aligned = sortalign(stringnames)
            D1 = {}
            leftlimit = 0
            for name, seq in aligned:
                D1[name] = seq

                # Do not allow seqeuence to the left of the seed (may include adapter/barcodes)
                if not name.split(";")[-1]:
                    leftlimit = min([seq.index(j) for j in seq if j != "-"])

            # Reorder keys by derep number
            keys = D1.keys()
            keys.sort(key=lambda x: int(x.split(";")[1].replace("size=", "")),
                      reverse=True)
            for key in keys:
                STACK.append(key + "\n" + D1[key][leftlimit:])

        if STACK:
            print >> fw, "\n".join(STACK)
            print >> fw, SEP
        cnts += 1

    fw.close()


def stats(clustSfile, statsfile, mindepth=0):
    C = ClustFile(clustSfile)
    depth = []
    for data in C.iter_seqs():
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


def makeclust(derepfile, userfile, notmatchedfile, clustfile):
    D = dict(parse_fasta(derepfile))
    U = defaultdict(list)  # Clusters
    fp = open(userfile)
    for row in fp:
        query, target, id, gaps, qstrand, qcov = row.rstrip().split("\t")
        U[target].append([query, qstrand, qcov])

    fw = open(clustfile, "w")
    for key, values in U.items():
        seqs = [('>' + key, D[key])]
        for name, strand, cov in values:
            cov = float(cov)
            # Only match forward reads if high Cov
            if strand == "+" and cov >= 90:
                seqs.append(('>' + name, D[name]))
        seq = "\n".join("\n".join(x) for x in seqs)
        print >> fw, "\n".join((seq, SEP))

    I = dict(parse_fasta(notmatchedfile))
    singletons = set(I.keys()) - set(U.keys())
    for key in singletons:
        print >> fw, "\n".join(('>' + key, I[key], SEP))
    fw.close()


def derep(fastafile, derepfile, minlength, cpus, usearch="vsearch"):
    cmd = usearch + " -minseqlength {0}".format(minlength)
    cmd += " -derep_fulllength {0}".format(fastafile)
    cmd += " -output {0} -sizeout".format(derepfile)
    cmd += " -threads {0}".format(cpus)
    sh(cmd)


def cluster_smallmem(derepfile, userfile, notmatchedfile, minlength, identity,
                     cpus, usearch="vsearch"):
    cmd = usearch + " -minseqlength {0}".format(minlength)
    cmd += " -leftjust"
    cmd += " -cluster_smallmem {0}".format(derepfile)
    cmd += " -id {0}".format(identity)
    cmd += " -userout {0}".format(userfile)
    cmd += " -userfields query+target+id+gaps+qstrand+qcov"
    cmd += " -maxaccepts 1 -maxrejects 0"
    cmd += " -minsl .5 -fulldp"
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
    p.set_align(pctid=96)
    p.set_outdir()
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    prefix = args[0]
    fastqfiles = args[1:]
    cpus = opts.cpus
    identity = opts.pctid / 100.
    minlength = opts.minlength
    fastafile, qualfile = fasta(fastqfiles + ["--seqtk",
                                "--outdir={0}".format(opts.outdir),
                                "--outfile={0}".format(prefix + ".fasta")])

    pf = prefix + ".P{0}".format(opts.pctid)
    pf = op.join(opts.outdir, pf)
    derepfile = pf + ".derep"
    if need_update(fastafile, derepfile):
        derep(fastafile, derepfile, minlength, cpus)

    userfile = pf + ".u"
    notmatchedfile = pf + ".notmatched"
    if need_update(derepfile, userfile):
        cluster_smallmem(derepfile, userfile, notmatchedfile,
                         minlength, identity, cpus)

    clustfile = pf + ".clust"
    if need_update((derepfile, userfile, notmatchedfile), clustfile):
        makeclust(derepfile, userfile, notmatchedfile, clustfile)

    clustSfile = pf + ".clustS"
    if need_update(clustfile, clustSfile):
        musclewrap(clustfile)

    statsfile = pf + ".stats"
    if need_update(clustSfile, statsfile):
        stats(clustSfile, statsfile, mindepth=opts.mindepth)


if __name__ == '__main__':
    main()
