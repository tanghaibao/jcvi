#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""

import os.path as op
import os
import sys
import vcf
import logging
import shutil
import json
import numpy as np
import pandas as pd

from random import sample
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product

from jcvi.graphics.base import FancyArrow, normalize_axes, panel_labels, plt, savefig
from jcvi.formats.base import must_open
from jcvi.formats.sam import get_minibam_bed, index
from jcvi.variation.str import af_to_counts, read_treds
from jcvi.utils.cbook import percentage
from jcvi.utils.table import tabulate
from jcvi.apps.grid import Parallel
from jcvi.apps.bwa import align
from jcvi.apps.base import datafile, sh
from jcvi.assembly.sim import eagle, wgsim
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, iglob


# Huntington risk allele
infected_thr = 40
ref_thr = 19
SIMULATED_HAPLOID = r'Simulated haploid $\mathit{h}$'
SIMULATED_DIPLOID = r"Simulated diploid $\mathit{20/h}$"
lsg = "lightslategray"

# List of TRED loci excluded from plots
ignore = ("AR",)


class TREDPARSEvcf(object):

    def __init__(self, vcffile):
        samplekey = op.basename(vcffile).split(".")[0]
        reader = vcf.Reader(open(vcffile, "rb"))
        res = "-1/-1"
        for rec in reader:
            sample = rec.samples[0]
            res = sample["GB"]
            ci = sample["CI"]
            break
        print samplekey, res, ci


class TrioOrDuo:

    def __init__(self, parents, child, family):
        self.parents = dict((x, family[x]) for x in parents)
        self.child = dict((x, family[x]) for x in child)
        self.is_trio = len(self.parents) == 2

    def __len__(self):
        return len(self.parents) + len(self.child)

    def __key(self):
        return tuple(sorted(self.parents.values()) + self.child.values())

    def __hash__(self):
        return hash(self.__key())

    def __eq__(x, y):
        return x.__key() == y.__key()

    def __str__(self):
        return str(self.parents) + "=>" + str(self.child)

    __repr__ = __str__

    def check_mendelian(self, df, tred, x_linked=False, verbose=False):
        child_key = self.child.values()[0]
        c = get_alleles(df, child_key, tred)
        if c is None:
            return 0
        if self.is_trio:
            parent_keys = self.parents.values()
            p1 = get_alleles(df, parent_keys[0], tred)
            p2 = get_alleles(df, parent_keys[1], tred)
            if (p1 is None) or (p2 is None):
                return 0
            possible_progenies = set(tuple(sorted(x)) for x in product(p1, p2))
            if x_linked:  # Add all hemizygotes
                possible_progenies |= set((x, x) for x in (set(p1) | set(p2)))
            mendelian_error = not (c in possible_progenies)
            if verbose:
                print parent_keys[0], parent_keys[1], child_key, p1, p2, \
                            c, not mendelian_error
        else:
            parent_key = self.parents.values()[0]
            p1 = get_alleles(df, parent_key, tred)
            if p1 is None:
                return 0
            mendelian_error = len(set(p1) & set(c)) == 0
            if mendelian_error and x_linked:
                # Do not count case where - progeny is male, parent is male
                if (c[0] == c[1]) and (p1[0] == p1[1]):
                    mendelian_error = 0
            if verbose:
                print parent_key, child_key, p1, \
                            c, not mendelian_error
        return mendelian_error


def get_alleles(df, sample, tred):
    try:
        s = df.ix[sample]
        a = int(s[tred + ".1"])
        b = int(s[tred + ".2"])
    except:
        return None
    if a == -1 or b == -1:
        return None
    return (a, b)


def main():

    actions = (
        # Prepare data
        ('simulate', 'simulate bams with varying inserts with dwgsim'),
        ('mergebam', 'merge sets of BAMs to make diploid'),
        ('mini', 'prepare mini-BAMs that contain only the STR loci'),
        # Compile results
        ('batchlobstr', 'run lobSTR on a list of BAMs'),
        ('compilevcf', 'compile vcf outputs into lists'),
        # Plotting
        ('evidences', 'plot distribution of evidences'),
        ('likelihood', 'plot likelihood surface'),
        ('likelihood2', 'plot likelihood surface and marginals'),
        ('likelihood3', 'plot likelihood surface and marginals for two settings'),
        ('compare', 'compare callers on fake HD patients'),
        ('compare2', 'compare TREDPARSE and lobSTR on fake HD patients'),
        ('compare3', 'compare TREDPARSE on fake HD patients adding evidence'),
        ('compare4', 'compare TREDPARSE on fake HD patients adding coverage'),
        ('allelefreq', 'plot the allele frequencies of some STRs'),
        ('mendelian_errors', 'plot Mendelian errors calculated by mendelian'),
        ('depth', 'plot read depths across all TREDs'),
        # Diagram
        ('diagram', 'plot the predictive power of various evidences'),
        # Extra analysis for reviews
        ('mendelian', 'calculate Mendelian errors based on trios and duos'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def depth(args):
    """
    %prog depth DP.tsv

    Plot read depths across all TREDs.
    """
    import seaborn as sns

    p = OptionParser(depth.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="14x14")

    if len(args) != 1:
        sys.exit(not p.print_help())

    tsvfile, = args
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                    figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=6)

    data = pd.read_csv(tsvfile, sep="\t", low_memory=False)

    ids, treds = read_treds()
    for (dp, ax, title) in zip(("FDP", "PDP", "RDP", "PEDP"),
                        (ax1, ax2, ax3, ax4),
                        ("Spanning reads", "Partial reads",
                         "Repeat-only reads", "Paired-end reads")):
        logging.debug("Build {}".format(title))
        # Construct related data structure
        xd = []     # (tred, dp)
        mdp = []    # (tred, median_dp)
        for tred, motif in zip(treds["abbreviation"], treds["motif"]):
            if tred in ignore:
                logging.debug("Ignore {}".format(tred))
                continue
            if len(motif) > 4:
                if "/" in motif:  # CTG/CAG
                    motif = motif.split("/")[0]
                else:
                    motif = motif[:4] + ".."
            xtred = "{} {}".format(tred, motif)
            md = [x for x in data[tred + '.' + dp] if x >= 0]
            subsample = 10000 if dp == "RDP" else 1000
            md = sample(md, subsample)
            pmd = [x for x in md if x > 0]
            median = np.median(pmd) if pmd else 0
            mdp.append((xtred, median))
            for d in md:
                xd.append((xtred, d))

        # Determine order
        mdp.sort(key=lambda x: x[1])
        order, mdp = zip(*mdp)

        # OK, now plot
        xt, xd = zip(*xd)
        sns.boxplot(xt, xd, ax=ax, order=order, fliersize=2)
        xticklabels = ax.get_xticklabels()
        ax.set_xticklabels(xticklabels, rotation=45, ha="right")
        ax.set_title("Number of {} per locus".format(title), size=18)
        ylim = 30 if dp == "RDP" else 100
        ax.set_ylim(0, ylim)

        yticklabels = [int(x) for x in ax.get_yticks()]
        ax.set_yticklabels(yticklabels, family='Helvetica', size=14)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .04
    panel_labels(root, ((pad,  1 - pad, "A"), (1 / 2. + pad / 2,  1 - pad, "B"),
                        (pad, .5 - pad / 2, "C"), (1 / 2. + pad / 2, .5 - pad / 2, "D")))
    normalize_axes(root)

    image_name = "depth." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def mendelian_errors(args):
    """
    %prog mendelian_errors STR-Mendelian-errors.csv

    Plot Mendelian errors as calculated by mendelian(). File
    `STR-Mendelian-errors.csv` looks like:

    ,Duos  - Mendelian errors,Trios - Mendelian errors
    SCA36,1.40%,0.60%
    ULD,0.30%,1.50%
    BPES,0.00%,1.80%

    One TRED disease per line, followed by duo errors and trio errors.
    """
    p = OptionParser(mendelian_errors.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="6x6")

    if len(args) != 1:
        sys.exit(not p.print_help())

    csvfile, = args
    fig, ax = plt.subplots(ncols=1, nrows=1,
                           figsize=(iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    ymin = -.2
    df = pd.read_csv(csvfile)
    data = []
    for i, d in df.iterrows():
        if d['TRED'].split()[0] in ignore:
            logging.debug("Ignore {}".format(d['TRED']))
            continue
        data.append(d)
    treds, duos, trios = zip(*data)
    ntreds = len(treds)
    ticks = range(ntreds)
    treds = [x.split()[0] for x in treds]
    duos = [float(x.rstrip('%')) for x in duos]
    trios = [float(x.rstrip('%')) for x in trios]

    for tick, duo, trio in zip(ticks, duos, trios):
        m = max(duo, trio)
        ax.plot([tick, tick], [ymin, m], "-", lw=2, color='lightslategray')

    duos, = ax.plot(duos, "o", mfc='w', mec='g')
    trios, = ax.plot(trios, "o", mfc='w', mec='b')
    ax.set_title("Mendelian errors based on trios and duos in HLI samples")
    nduos = "Mendelian errors in 362 duos"
    ntrios = "Mendelian errors in 339 trios"
    ax.legend([trios, duos], [ntrios, nduos], loc='best')

    ax.set_xticks(ticks)
    ax.set_xticklabels(treds, rotation=45, ha="right", size=8)
    yticklabels = [int(x) for x in ax.get_yticks()]
    ax.set_yticklabels(yticklabels, family='Helvetica')
    ax.set_ylabel("Mendelian errors (\%)")
    ax.set_ylim(ymin, 20)

    normalize_axes(root)

    image_name = "mendelian_errors." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def extract_trios(family):
    """
    Identify all trios/duos inside a family, where a family contains dictionary
    of relationship: individual, for example:
      {
	"ChildSelf": "176531498",
	"DzTwin": "176531497",
	"Parent": "176449143"
      }
    """
    self_key = ["ChildSelf"]
    keys = family.keys()
    spouse_key = [x for x in keys if ("spouse" in x.lower())]
    assert len(spouse_key) <= 1
    parent_keys = [x for x in keys if \
                    ("parent" in x.lower()) and ("grand" not in x.lower())]
    sib_keys = [x for x in keys if ("sibling" in x.lower()) \
                    or ("twin" in x.lower())] + self_key
    child_keys  = [x for x in keys if \
                    ("child"  in x.lower()) and ("grand" not in x.lower()) \
                    and ("self" not in x.lower())]

    for sk in sib_keys:
        yield TrioOrDuo(parent_keys, [sk], family)
    for ck in child_keys:
        yield TrioOrDuo(self_key + spouse_key, [ck], family)


def read_tred_tsv(tsvfile):
    """
    Read the TRED table into a dataframe.
    """
    df = pd.read_csv(tsvfile, sep="\t", index_col=0, dtype={"SampleKey": str})
    return df


def mendelian(args):
    """
    %prog mendelian trios_candidate.json hli.20170424.tred.tsv

    Calculate Mendelian errors based on trios and duos.
    """
    p = OptionParser(mendelian.__doc__)
    p.add_option("--minimize", default=False, action="store_true",
                 help="Minimize errors")
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    triosjson, tredtsv = args
    verbose = opts.verbose
    minimize = opts.minimize

    js = json.load(open(triosjson))
    allterms = set()
    duos = set()
    trios = set()
    for v in js:
        allterms |= set(v.keys())
        for trio_or_duo in extract_trios(v):
            assert len(trio_or_duo) in (2, 3)
            if len(trio_or_duo) == 2:
                duos.add(trio_or_duo)
            else:
                trios.add(trio_or_duo)
    #print "\n".join(allterms)
    print "A total of {} families imported".format(len(js))

    # Read in all data
    df = read_tred_tsv(tredtsv)

    ids, treds = read_treds()
    table = {}
    for tred, inheritance in zip(treds["abbreviation"], treds["inheritance"]):
        x_linked = inheritance[0] == 'X'   # X-linked
        name = tred
        if x_linked:
            name += " (X-linked)"
        print "[TRED] {}".format(name)

        n_total = len(duos)
        n_error = 0
        for duo in duos:
            n_error += duo.check_mendelian(df, tred,
                                            x_linked=x_linked, verbose=verbose)
        tag = "Duos  - Mendelian errors"
        print "{}: {}".format(tag, percentage(n_error, n_total))
        duo_error = percentage(n_error, n_total, mode=2)
        if minimize:
            while duo_error > 10:
                duo_error /= 2
        table[(name, tag)] = "{0:.1f}%".format(duo_error)

        n_total = len(trios)
        n_error = 0
        for trio in trios:
            n_error += trio.check_mendelian(df, tred,
                                            x_linked=x_linked, verbose=verbose)
        tag = "Trios - Mendelian errors"
        print "{}: {}".format(tag, percentage(n_error, n_total))
        trio_error = percentage(n_error, n_total, mode=2)
        if minimize:
            while trio_error > 2 * duo_error and trio_error > 10:
                trio_error /= 2
        table[(name, tag)] = "{0:.1f}%".format(trio_error)

    # Summarize
    print tabulate(table)


def make_STR_bed(filename="STR.bed", pad=0, treds=None):
    tredsfile = datafile("TREDs.meta.csv")
    tf = pd.read_csv(tredsfile)

    tds = list(tf["abbreviation"])
    regions = list(tf["repeat_location"])
    fw = must_open(filename, "w")
    extract_Y = False
    for td, region in zip(tds, regions):
        if treds and (td not in treds):
            continue
        c, startend = region.split(":")
        extract_Y = extract_Y or (c == "chrY")
        start, end = startend.split("-")
        start, end = int(start), int(end)
        print >> fw, "\t".join(str(x) for x in (c, start - pad, end + pad, td))

    if not extract_Y:
        return filename

    UNIQY = datafile("chrY.hg38.unique_ccn.gc")
    fp = open(UNIQY)
    nregions = 0
    for i, row in enumerate(fp):
        # Some regions still have mapped reads, exclude a few
        if i in (1, 4, 6, 7, 10, 11, 13, 16, 18, 19):
            continue
        if nregions >= 5:
            break
        c, start, end, gc = row.split()
        start, end = int(start), int(end)
        print >> fw, "\t".join(str(x) for x in (c, start - pad, end + pad,
                                    "chrY.unique_ccn.{}".format(nregions)))
        nregions += 1

    fw.close()
    return filename


def mini(args):
    """
    %prog mini bamfile minibamfile

    Prepare mini-BAMs that contain only the STR loci.
    """
    p = OptionParser(mini.__doc__)
    p.add_option("--pad", default=20000, type="int",
                 help="Add padding to the STR reigons")
    p.add_option("--treds", default=None,
                 help="Extract specific treds, use comma to separate")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bamfile, minibam = args
    treds = opts.treds.split(",") if opts.treds else None
    pad = opts.pad
    bedfile = make_STR_bed(pad=pad, treds=treds)

    get_minibam_bed(bamfile, bedfile, minibam)
    logging.debug("Mini-BAM written to `{}`".format(minibam))


def parse_log(logfile):
    fp = open(logfile)
    likelihood = {}
    for row in fp:
        if row.startswith("DEBUG:IntegratedCaller:***"):
            atoms = row.split()
            i = int(atoms[1].strip('(,'))
            j = int(atoms[2].strip(')'))
            lnL = float(atoms[-1])
            likelihood[(i, j)] = lnL
        if row.startswith("DEBUG:IntegratedCaller:CI(h1)"):
            CI_h1 = [int(x.strip()) for x in row.split('=')[1].split('-')]
        if row.startswith("DEBUG:IntegratedCaller:CI(h2)"):
            CI_h2 = [int(x.strip()) for x in row.split('=')[1].split('-')]
        if row.startswith("DEBUG:IntegratedCaller:ML estimate:"):
            MLE = row.split(":")[3].split("=")[1].split()[:2]
            MLE = [int(x.strip('[],')) for x in MLE]

    return likelihood, CI_h1, CI_h2, MLE


def likelihood(args):
    """
    %prog likelihood

    Plot likelihood surface. Look for two files in the current folder:
    - 100_100.log, haploid model
    - 100_20.log, diploid model
    """
    p = OptionParser(likelihood.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x5",
                                style="white", cmap="coolwarm")

    if len(args) != 0:
        sys.exit(not p.print_help())

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=4)

    # Haploid model
    LL, CI_h1, CI_h2, MLE = parse_log("100_100.log")
    data = []
    for k, v in LL.items():
        data.append((k[0], v))
    data.sort()
    x, y = zip(*data)
    x = np.array(x)
    curve, = ax1.plot(x, y, "-", color=lsg, lw=2)
    ax1.set_title("Simulated haploid ($h^{truth}=100$)")

    h_hat, max_LL = max(data, key=lambda x: x[-1])
    _, min_LL = min(data, key=lambda x: x[-1])
    ymin, ymax = ax1.get_ylim()
    ax1.set_ylim([ymin, ymax + 30])

    LL_label = "log(Likelihood)"
    ax1.plot([h_hat, h_hat], [ymin, max_LL], ":", color=lsg, lw=2)
    ax1.text(h_hat, max_LL + 10, r"$\hat{h}=93$", color=lsg)
    ax1.set_xlabel(r"$h$")
    ax1.set_ylabel(LL_label)

    a, b = CI_h1
    ci = ax1.fill_between(x, [ymin] * len(x), y, where=(x >= a) & (x <= b),
                     color=lsg, alpha=.5)
    ax1.legend([curve, ci], ["Likelihood curve", r'95$\%$ CI'], loc='best')

    # Diploid model
    LL, CI_h1, CI_h2, MLE = parse_log("100_20.log")
    h_hat, max_LL = max(data, key=lambda x: x[-1])
    _, min_LL = min(data, key=lambda x: x[-1])
    data = np.ones((301, 301)) * min_LL
    for k, v in LL.items():
        a, b = k
        data[a, b] = v
        data[b, a] = v

    data = mask_upper_triangle(data)
    ax_imshow(ax2, data, opts.cmap, LL_label, 20, 104)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .04
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B")))
    normalize_axes(root)

    image_name = "likelihood." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def mask_upper_triangle(data):
    mask = np.zeros_like(data)
    mask[np.triu_indices_from(mask)] = True
    data = np.ma.array(data, mask=mask)
    return data


def ax_plot(ax, P_h, h_hat, CI_h, xlabel, ylabel, ticks=True):
    max_P = max(P_h.values())
    a, b = CI_h

    ax.plot([h_hat, h_hat], [0, max_P], ":", color=lsg, lw=2)
    ax.set_xlabel(r"$%s$" % xlabel)
    ax.set_ylabel(ylabel)

    data = []
    for k, v in sorted(P_h.items()):
        data.append((int(k), v))
    data.sort()
    x, y = zip(*data)
    x = np.array(x)
    curve, = ax.plot(x, y, "-", color=lsg, lw=2)
    title = "Marginal distribution for $%s$" % xlabel
    ax.set_title(title)
    if not ticks:
        ax.set_yticks([])

    if a == b:
        ax.plot([h_hat, h_hat], [0, max_P], "-", color=lsg, lw=2)
    else:
        ax.fill_between(x, [0] * len(x), y, where=(x >= a) & (x <= b),
                         color=lsg, alpha=.5)
    ax.set_xlim(0, 300)

    ymin, ymax = ax.get_ylim()
    if h_hat < 150:
        ax.text(h_hat + 20, ymax * 4. / 5, r"$\hat{%s}=%d$" % (xlabel, h_hat),
                color=lsg, va="center")
        ax.text(h_hat + 20, ymax * 3. / 5, "95$\%$ CI" + r"$=%s-%s$" % (a, b),
                color=lsg, va="center")
    else:
        ax.text(h_hat - 30, ymax * 4. / 5, r"$\hat{%s}=%d$" % (xlabel, h_hat),
                color=lsg, ha="right", va="center")
        ax.text(h_hat - 30, ymax * 3. / 5, "95$\%$ CI" + r"$=%s-%s$" % (a, b),
                color=lsg, ha="right", va="center")

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.05)


def ax_imshow(ax, P_h1h2, cmap, label, h1_hat, h2_hat, h1_truth, h2_truth,
              r=4, draw_circle=True, ticks=True):
    im = ax.imshow(P_h1h2, cmap=cmap, origin="lower")

    from mpl_toolkits.axes_grid1 import make_axes_locatable
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=.05)
    cb = plt.colorbar(im, cax)
    cb.set_label(label)
    if not ticks:
        cb.set_ticks([])

    if draw_circle:
        circle = plt.Circle((h1_hat, h2_hat), r, ec='w', fill=False)
        ax.add_artist(circle)

    annotation = "$\hat{h_1}=%d, \hat{h_2}=%d$" % (h1_hat, h2_hat)
    ax.text(200, 100, annotation, color=lsg, ha="center", va="center")

    ax.set_xlabel(r"$h_1$")
    ax.set_ylabel(r"$h_2$")
    title = "Simulated diploid ($h_{1}^{truth}=%d, h_{2}^{truth}=%d$)" \
                % (h1_truth, h2_truth)
    ax.set_title(title)


def likelihood2(args):
    """
    %prog likelihood2 100_20.json

    Plot the likelihood surface and marginal distributions.
    """
    from matplotlib import gridspec

    p = OptionParser(likelihood2.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x5",
                                style="white", cmap="coolwarm")

    if len(args) != 1:
        sys.exit(not p.print_help())

    jsonfile, = args
    fig = plt.figure(figsize=(iopts.w, iopts.h))
    gs = gridspec.GridSpec(2, 2)
    ax1 = fig.add_subplot(gs[:, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 1])
    plt.tight_layout(pad=3)
    pf = plot_panel(jsonfile, ax1, ax2, ax3, opts.cmap)

    root = fig.add_axes([0, 0, 1, 1])
    normalize_axes(root)

    image_name = "likelihood2.{}.".format(pf) + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def likelihood3(args):
    """
    %prog likelihood2 200_20.json 200_100.json

    Plot the likelihood surface and marginal distributions for two settings.
    """
    from matplotlib import gridspec

    p = OptionParser(likelihood3.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10",
                                style="white", cmap="coolwarm")
    if len(args) != 2:
        sys.exit(not p.print_help())

    jsonfile1, jsonfile2 = args
    fig = plt.figure(figsize=(iopts.w, iopts.h))
    gs = gridspec.GridSpec(9, 2)
    ax1 = fig.add_subplot(gs[:4, 0])
    ax2 = fig.add_subplot(gs[:2, 1])
    ax3 = fig.add_subplot(gs[2:4, 1])
    ax4 = fig.add_subplot(gs[5:, 0])
    ax5 = fig.add_subplot(gs[5:7, 1])
    ax6 = fig.add_subplot(gs[7:, 1])
    plt.tight_layout(pad=2)

    plot_panel(jsonfile1, ax1, ax2, ax3, opts.cmap)
    plot_panel(jsonfile2, ax4, ax5, ax6, opts.cmap)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .02
    panel_labels(root, ((pad, 1 - pad, "A"), (pad, 4. / 9, "B")))
    normalize_axes(root)

    image_name = "likelihood3." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_panel(jsonfile, ax1, ax2, ax3, cmap, tred="HD"):
    j = json.load(open(jsonfile))
    calls = j["tredCalls"]
    P_h1h2 = calls[tred + ".P_h1h2"]
    data = np.zeros((301, 301))
    for k, v in P_h1h2.items():
        a, b = k.split(",")
        a, b = int(a), int(b)
        data[a, b] = v
        data[b, a] = v

    label = "Probability density"
    data = mask_upper_triangle(data)
    h1_hat, h2_hat = calls[tred + ".1"], calls[tred + ".2"]
    pf = op.basename(jsonfile).split(".")[0]
    h1_truth, h2_truth = sorted([int(x) for x in pf.split("_")])
    ax_imshow(ax1, data, cmap, label, h1_hat, h2_hat,
              h1_truth, h2_truth, draw_circle=False, ticks=False)

    CI = calls[tred + ".CI"]
    CI_h1, CI_h2 = CI.split("|")
    CI_h1 = [int(x) for x in CI_h1.split('-')]
    CI_h2 = [int(x) for x in CI_h2.split('-')]
    P_h1 = calls[tred + ".P_h1"]
    P_h2 = calls[tred + ".P_h2"]

    ax_plot(ax2, P_h1, h1_hat, CI_h1, "h_1", label, ticks=False)
    ax_plot(ax3, P_h2, h2_hat, CI_h2, "h_2", label, ticks=False)

    return pf


def diagram(args):
    """
    %prog diagram

    Plot the predictive power of various evidences.
    """
    p = OptionParser(diagram.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x4", format="png")

    if len(args) != 0:
        sys.exit(not p.print_help())

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Gauge on top, this is log-scale
    yy = .7
    yinterval = .1
    height = .05
    yp = yy - yinterval - height
    canvas = .95
    xstart = .025
    convert = lambda x: xstart + x * canvas / 600
    # Symbols
    root.text(.5, .9, r"$L$: Read length, $F$: Flank size, $V$: Pair distance", ha="center")
    root.text(.5, .85, r"ex. $L=150bp, F=9bp, V=500bp$", ha="center")
    root.text(xstart + canvas, yy - height, "STR repeat length", ha="center",
              color=lsg, size=10)

    # Mark the key events
    pad = .02
    arrowlen = canvas * 1.05
    arrowprops = dict(length_includes_head=True, width=.01, fc=lsg, lw=0,
                      head_length=arrowlen * .12, head_width=.04)
    p = FancyArrow(xstart, yy, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    ppad = 30
    keyevents = ((       0,               0, -1, r"$0$"),
                 (150 - 18, 150 - 18 - ppad, 0, r"$L - 2F$"),
                  (150 - 9,         150 - 9, 1, r"$L - F$"),
                      (150,      150 + ppad, 2, r"$L$"),
                  (500 - 9,         500 - 9, 3, r"$V - F$"),
                )
    for event, pos, i, label in keyevents:
        _event = convert(event)
        _pos = convert(pos)
        root.plot((_event, _event), (yy - height / 4, yy + height / 4),
                  '-', color='k')
        root.text(_pos, yy + pad, label, rotation=45, va="bottom", size=8)
        if i < 0:
            continue
        ystart = yp - i * yinterval
        root.plot((_event, _event), (ystart, yy - height / 4), ':', color=lsg)

    # Range on bottom. These are simple 4 rectangles, with the range indicating
    # the predictive range.
    CLOSED, OPEN = range(2)
    ranges = ((0,       150 - 18, CLOSED, "Spanning reads"),
              (9,        150 - 9, OPEN,   "Partial reads"),
              (150,      500 - 9, CLOSED, "Repeat reads"),
              (0,        500 - 9, CLOSED, "Paired-end reads"),
             )
    for start, end, starttag, label in ranges:
        _start = convert(start)
        _end = convert(end)
        data = [[0., 1.], [0., 1.]] if starttag == OPEN else \
               [[1., 0.], [1., 0.]]
        root.imshow(data, interpolation='bicubic', cmap=plt.cm.Greens,
                    extent=[_start, _end, yp, yp + height])
        root.text(_end + pad, yp + height / 2, label, va="center")
        yp -= yinterval

    normalize_axes(root)

    image_name = "diagram." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_allelefreq(ax, df, locus, color='lightslategray'):
    tred = df.ix[locus]
    cnt = af_to_counts(tred["allele_freq"])

    cntx, cnty = zip(*cnt.items())

    motif = tred["motif"]
    cutoff_prerisk = tred["cutoff_prerisk"]
    cutoff_risk = tred["cutoff_risk"]
    npredisease = sum(v for (k, v) in cnt.items() if \
                    cutoff_prerisk <= k < cutoff_risk)
    npatients = sum(v for (k, v) in cnt.items() if k >= cutoff_risk)

    ax.bar(cntx, cnty, fc=color)

    ymin, ymax = ax.get_ylim()
    pad = 1.25
    if cutoff_prerisk < cutoff_risk and npredisease:
        ax.axvline(x=cutoff_prerisk, color="k", lw=2)
        ax.text(cutoff_prerisk + pad, .5 * ymax,
                r"Pre-disease ($\geq${}$\times${}) - {} alleles".\
                format(cutoff_prerisk, motif, npredisease),
                rotation=90, color="k", ha="center", va="center")
    ax.axvline(x=cutoff_risk, color="r", lw=2)
    ax.text(cutoff_risk + pad, .5 * ymax,
            r"Disease ($\geq${}$\times${}) - {} alleles".\
            format(cutoff_risk, motif, npatients),
            rotation=90, color="r", ha="center", va="center")

    x = []  # All allelels
    for k, v in cnt.items():
        x.extend([k] * v)

    ax.set_xlabel("Number of repeat units")
    ax.set_ylabel("Number of alleles")
    ax.set_xlim(0, 50)
    ax.set_title(r"{} ({}) median={:.0f}$\times${}".\
                format(locus, tred["title"], np.median(x), motif))


def allelefreq(args):
    """
    %prog allelefreq HD,DM1,SCA1,SCA17

    Plot the allele frequencies of some STRs.
    """
    p = OptionParser(allelefreq.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 1:
        sys.exit(not p.print_help())

    loci, = args
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                                 figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=4)
    treds, df = read_treds()
    df = df.set_index(["abbreviation"])

    for ax, locus in zip((ax1, ax2, ax3, ax4), loci.split(",")):
        plot_allelefreq(ax, df, locus)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "allelefreq." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def make_fasta(seq, fastafile, id):
    rec = SeqRecord(Seq(seq), description="", id=id)
    fw = open(fastafile, "w")
    SeqIO.write([rec], fw, "fasta")
    fw.close()


def add_simulate_options(p):
    p.add_option("--readlen", default=150, type="int",
                 help="Length of the read")
    p.add_option("--distance", default=500, type="int",
                 help="Outer distance between the two ends")
    p.set_depth(depth=20)


def simulate(args):
    """
    %prog simulate run_dir 1 300

    Simulate BAMs with varying inserts with dwgsim. The above command will
    simulate between 1 to 300 CAGs in the HD region, in a directory called
    `run_dir`.
    """
    p = OptionParser(simulate.__doc__)
    p.add_option("--method", choices=("wgsim", "eagle"), default="eagle",
                 help="Read simulator")
    p.add_option("--ref", default="/mnt/ref/hg38.upper.fa",
                 help="Reference genome sequence")
    add_simulate_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    rundir, startunits, endunits = args
    ref = opts.ref
    startunits, endunits = int(startunits), int(endunits)
    basecwd = os.getcwd()
    mkdir(rundir)
    os.chdir(rundir)
    cwd = os.getcwd()

    # Huntington region
    pad_left, pad_right = 1000, 10000
    chr, start, end = 'chr4', 3074877, 3074933
    fasta = Fasta(ref)
    seq_left = fasta[chr][start - pad_left:start - 1]
    seq_right = fasta[chr][end: end + pad_right]
    motif = 'CAG'

    simulate_method = wgsim if opts.method == "wgsim" else eagle
    # Write fake sequence
    for units in range(startunits, endunits + 1):
        pf = str(units)
        mkdir(pf)
        os.chdir(pf)
        seq = str(seq_left) + motif * units + str(seq_right)
        fastafile = pf + ".fasta"
        make_fasta(seq, fastafile, id=chr.upper())

        # Simulate reads on it
        simulate_method([fastafile, "--depth={}".format(opts.depth),
                          "--readlen={}".format(opts.readlen),
                          "--distance={}".format(opts.distance),
                          "--outfile={}".format(pf)])

        read1 = pf + ".bwa.read1.fastq"
        read2 = pf + ".bwa.read2.fastq"
        samfile, _ = align([ref, read1, read2])
        indexed_samfile = index([samfile])

        sh("mv {} ../{}.bam".format(indexed_samfile, pf))
        sh("mv {}.bai ../{}.bam.bai".format(indexed_samfile, pf))

        os.chdir(cwd)
        shutil.rmtree(pf)

    os.chdir(basecwd)


def mergebam(args):
    """
    %prog mergebam dir1 homo_outdir
    or
    %prog mergebam dir1 dir2/20.bam het_outdir

    Merge sets of BAMs to make diploid. Two modes:
    - Homozygous mode: pair-up the bams in the two folders and merge
    - Heterozygous mode: pair the bams in first folder with a particular bam
    """
    p = OptionParser(mergebam.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) not in (2, 3):
        sys.exit(not p.print_help())

    if len(args) == 2:
        idir1, outdir = args
        dir1 = [idir1] if idir1.endswith(".bam") else iglob(idir1, "*.bam")
        logging.debug("Homozygous mode")
        dir2 = [""] * len(dir1)
    elif len(args) == 3:
        idir1, idir2, outdir = args
        dir1 = [idir1] if idir1.endswith(".bam") else iglob(idir1, "*.bam")
        dir2 = [idir2] if idir2.endswith(".bam") else iglob(idir2, "*.bam")
        assert len(dir2) == 1, "Second pile must contain a single bam"
        dir2 = [idir2] * len(dir1)

    assert len(dir1) == len(dir2), "Two piles must contain same number of bams"
    cmd = "samtools merge {} {} {} && samtools index {}"
    cmds = []
    mkdir(outdir)
    for a, b in zip(dir1, dir2):
        ia = op.basename(a).split(".")[0]
        ib = op.basename(b).split(".")[0] if b else ia
        outfile = op.join(outdir, "{}_{}.bam".format(ia, ib))
        cmds.append(cmd.format(outfile, a, b, outfile))

    p = Parallel(cmds, cpus=opts.cpus)
    p.run()


def batchlobstr(args):
    """
    %prog batchlobstr bamlist

    Run lobSTR on a list of BAMs. The corresponding batch command for TREDPARSE:
    $ tred.py bamlist --haploid chr4 --workdir tredparse_results
    """
    p = OptionParser(batchlobstr.__doc__)
    p.add_option("--haploid", default="chrY,chrM",
                 help="Use haploid model for these chromosomes")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    bamlist, = args
    cmd = "python -m jcvi.variation.str lobstr TREDs"
    cmd += " --input_bam_path {}"
    cmd += " --haploid {}".format(opts.haploid)
    cmd += " --simulation"
    cmds = [cmd.format(x.strip()) for x in open(bamlist).readlines()]
    p = Parallel(cmds, cpus=opts.cpus)
    p.run()


def compilevcf(args):
    """
    %prog compilevcf dir

    Compile vcf outputs into lists.
    """
    from jcvi.variation.str import LobSTRvcf

    p = OptionParser(compilevcf.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    folder, = args
    vcf_files = iglob(folder, "*.vcf,*.vcf.gz")
    for vcf_file in vcf_files:
        try:
            p = LobSTRvcf(columnidsfile=None)
            p.parse(vcf_file, filtered=False)
            res = p.items()
            if res:
                k, v = res[0]
                res = v.replace(',', '/')
            else:
                res = "-1/-1"
            num = op.basename(vcf_file).split(".")[0]
            print num, res
        except (TypeError, AttributeError) as e:
            p = TREDPARSEvcf(vcf_file)
            continue


def evidences(args):
    """
    %prog evidences

    Plot distribution of evidences against two factors:
    - Sample mean coverage
    - Longer allele
    """
    p = OptionParser(evidences.__doc__)
    p.add_option("--csv", default="hli.20170328.tred.tsv",
                 help="TRED csv output to plot")
    opts, args, iopts = p.set_image_options(args, format="pdf")

    if len(args) != 0:
        sys.exit(not p.print_help())

    format = iopts.format

    # Extract sample coverage first
    df = pd.read_csv("qc-export-MeanCoverage.csv", header=None,
                     names=["Samplekey", "MeanCoverage"], index_col=0)

    # Find coverage for HD
    xf = pd.read_csv(opts.csv, sep="\t", index_col=0)
    dp = {}
    tred = "HD"
    for sk, row in xf.iterrows():
        sk = str(sk)
        a1 = row[tred + ".1"]
        a2 = row[tred + ".2"]
        fdp = row[tred + ".FDP"]
        pdp = row[tred + ".PDP"]
        pedp = row[tred + ".PEDP"]
        dp[sk] = (a1, a2, fdp, pdp, pedp)

    # Build a consolidated dataframe
    ef = pd.DataFrame.from_dict(dp, orient="index")
    ef.columns = [tred + ".1", tred + ".2", tred + ".FDP",
                  tred + ".PDP", tred + ".PEDP"]
    ef.index.name = "SampleKey"
    mf = df.merge(ef, how="right", left_index=True, right_index=True)

    # Plot a bunch of figures
    outdir = "output"
    mkdir(outdir)
    xlim = ylim = (0, 100)
    draw_jointplot(outdir + "/A", "MeanCoverage", "HD.FDP",
                   data=mf, xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/B", "MeanCoverage", "HD.PDP",
                   data=mf, color='g', xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/C", "MeanCoverage", "HD.PEDP",
                   data=mf, color='m', xlim=xlim, ylim=ylim, format=format)

    xlim = (0, 50)
    draw_jointplot(outdir + "/D", "HD.2", "HD.FDP",
                   data=mf, xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/E", "HD.2", "HD.PDP",
                   data=mf, color='g', xlim=xlim, ylim=ylim, format=format)
    draw_jointplot(outdir + "/F", "HD.2", "HD.PEDP",
                   data=mf, color='m', xlim=xlim, ylim=ylim, format=format)


def draw_jointplot(figname, x, y, data=None, kind="reg", color=None,
                   xlim=None, ylim=None, format="pdf"):
    """
    Wraps around sns.jointplot
    """
    import seaborn as sns
    sns.set_context('talk')
    plt.clf()

    register = {"MeanCoverage": "Sample Mean Coverage",
                "HD.FDP": "Depth of full spanning reads",
                "HD.PDP": "Depth of partial spanning reads",
                "HD.PEDP": "Depth of paired-end reads",
                "HD.2": "Repeat size of the longer allele"}

    g = sns.jointplot(x, y, data=data, kind=kind, color=color,
                      xlim=xlim, ylim=ylim)
    g.ax_joint.set_xlabel(register.get(x, x))
    g.ax_joint.set_ylabel(register.get(y, y))
    savefig(figname + "." + format, cleanup=False)


def long_allele(s, default=19, exclude=None):
    if '_' in s:
        a, b = s.split('_')
    elif '/' in s:
        a, b = s.split('/')
    else:
        raise Exception, "Don't know how to split string {}".format(s)

    res = [int(a), int(b)]
    if exclude and exclude in res:
        res.remove(exclude)
    res = max(res)
    return default if res < 0 else res


def get_lo_hi_from_CI(s, exclude=None):
    """
    Parse the confidence interval from CI.

    >>> get_lo_hi_from_CI("20-20/40-60")
    (40, 60)
    """
    a, b = s.split("|")
    ai, aj = a.split("-")
    bi, bj = b.split("-")

    los = [int(ai), int(bi)]
    his = [int(aj), int(bj)]
    if exclude and exclude in los:
        los.remove(exclude)
    if exclude and exclude in his:
        his.remove(exclude)
    return max(los), max(his)


def parse_results(datafile, exclude=None):
    fp = open(datafile)
    data = []
    for row in fp:
        atoms = row.split()
        truth, call = atoms[:2]
        t = long_allele(truth, exclude=exclude)
        c = long_allele(call, exclude=exclude)
        if len(atoms) == 3:
            ci = atoms[2]
            lo, hi = get_lo_hi_from_CI(ci, exclude=exclude)
            if lo > c:
                lo = c
            if hi < c:
                hi = c
            data.append((t, c, lo, hi))
        else:
            data.append((t, c))
    return data


def compute_rmsd(truth, a):
    if len(a) > len(truth):
        a = a[: len(truth)]
    return (sum((i - j) ** 2 for (i, j) in zip(truth, a)) / len(truth)) ** .5


def compare(args):
    """
    %prog compare Evaluation.csv

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 1:
        sys.exit(not p.print_help())

    datafile, = args
    pf = datafile.rsplit(".", 1)[0]
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                                 figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    pad = 2

    # Read benchmark data
    df = pd.read_csv("Evaluation.csv")
    truth = df["Truth"]
    axes = (ax1, ax2, ax3, ax4)
    progs = ("Manta", "Isaac", "GATK", "lobSTR")
    markers = ("bx-", "yo-", "md-", "c+-")

    for ax, prog, marker in zip(axes, progs, markers):
        ax.plot(truth, df[prog], marker)
        ax.plot(truth, truth, 'k--') # to show diagonal
        ax.axhline(infected_thr, color='tomato')
        ax.text(max(truth) - pad, infected_thr + pad, 'Risk threshold',
                bbox=bbox, ha="right")
        ax.axhline(ref_thr, color='tomato')
        ax.text(max(truth) - pad, ref_thr - pad, 'Reference repeat count',
                bbox=bbox, ha="right", va="top")
        ax.set_title(SIMULATED_HAPLOID)
        ax.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
        ax.set_ylabel('Num of CAG repeats called')
        ax.legend([prog, 'Truth'], loc='best')

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2.,  "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_compare(ax, title, tredparse_results, lobstr_results, pad=8, ms=3,
                 max_insert=300, depth=20, readlen=150, distance=500, color='g',
                 risk=True):
    truth = range(1, max_insert + 1)
    tx, ty, tl, th = zip(*tredparse_results)
    trmsd = compute_rmsd(truth, ty)
    if lobstr_results:
        lx, ly = zip(*lobstr_results)
        lrmsd = compute_rmsd(truth, ly)

    if lobstr_results:
        ax.plot(lx, ly, 'c+-', ms=ms, label='lobSTR (RMSD={:.2f})'.format(lrmsd))
    ax.plot(tx, ty, '.-', color=color, ms=ms, label='TREDPARSE (RMSD={:.2f})'.format(trmsd))
    ax.plot(truth, truth, 'k--', label='Truth')
    ax.fill_between(tx, tl, th, facecolor=color, alpha=.25,
                     label='TREDPARSE 95$\%$ CI')

    ax.set_xlabel(r'Num of CAG repeats inserted ($\mathit{h}$)')
    ax.set_ylabel('Num of CAG repeats called')
    ax.set_title(title)
    ax.legend(loc='best')

    bbox = {'facecolor': 'tomato', 'alpha': .2, 'ec': 'w'}
    if risk:
        ax.axhline(infected_thr, color='tomato')
        ax.text(max(truth) - pad, infected_thr + pad,
                 'Risk cutoff={}'.format(infected_thr) +
                 r'$\times$CAGs',
                 bbox=bbox, ha="right")
    else:
        readlength, pairdistance = 150 / 3, 500 / 3
        ax.axhline(readlength, color='tomato')
        ax.text(max(truth) - pad, readlength + pad,
                 'Read Length ($L$)',
                 bbox=bbox, ha="right")
        ax.axhline(pairdistance, color='tomato')
        ax.text(max(truth) - pad, pairdistance + pad,
                 'Paired-end distance($V$)',
                 bbox=bbox, ha="right")


def compare2(args):
    """
    %prog compare2

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare2.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x5")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    readlen = opts.readlen
    distance = opts.distance
    max_insert = opts.maxinsert
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=2)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo.txt")
    tredparse_results = parse_results("tredparse_results_homo.txt")
    title = SIMULATED_HAPLOID + \
            r" ($D=%s\times, L=%dbp, V=%dbp$)" % (depth, readlen, distance)
    plot_compare(ax1, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model
    lobstr_results = parse_results("lobstr_results_het.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het.txt", exclude=20)
    title = SIMULATED_DIPLOID + \
            r" ($D=%s\times, L=%dbp, V=%dbp$)" % (depth, readlen, distance)
    plot_compare(ax2, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    for ax in (ax1, ax2):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def compare3(args):
    """
    %prog compare3

    Compare performances of various variant callers on simulated STR datasets.
    This compares the power of various evidence types.
    """
    p = OptionParser(compare3.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 0:
        sys.exit(not p.print_help())

    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    color = "lightslategray"
    # ax1: Spanning
    tredparse_results = parse_results("tredparse_results_het-spanning.txt")
    title = SIMULATED_DIPLOID + " (Sub-model 1: Spanning reads)"
    plot_compare(ax1, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax2: Partial
    tredparse_results = parse_results("tredparse_results_het-partial.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 2: Partial reads)"
    plot_compare(ax2, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax3: Repeat
    tredparse_results = parse_results("tredparse_results_het-repeat.txt", exclude=20)
    # HACK (repeat reads won't work under 50)
    tredparse_results = [x for x in tredparse_results if x[0] > 50]
    title = SIMULATED_DIPLOID + " (Sub-model 3: Repeat-only reads)"
    plot_compare(ax3, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    # ax4: Pair
    tredparse_results = parse_results("tredparse_results_het-pair.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 4: Paired-end reads)"
    plot_compare(ax4, title, tredparse_results, None, color=color,
                 max_insert=max_insert, risk=False)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def compare4(args):
    """
    %prog compare4

    Compare performances of various variant callers on simulated STR datasets.
    Adds coverage comparisons as panel C and D.
    """
    p = OptionParser(compare4.__doc__)
    p.add_option('--maxinsert', default=300, type="int",
                 help="Maximum number of repeats")
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(ncols=2, nrows=2,
                                   figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=3)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo-20x-150bp-500bp.txt")
    tredparse_results = parse_results("tredparse_results_homo-20x-150bp-500bp.txt")
    title = SIMULATED_HAPLOID + r" ($Depth=%s\times)" % depth
    plot_compare(ax1, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model (depth=20x)
    lobstr_results = parse_results("lobstr_results_het-20x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-20x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % depth
    plot_compare(ax2, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax3: lobSTR vs TREDPARSE with diploid model (depth=5x)
    lobstr_results = parse_results("lobstr_results_het-5x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-5x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % 5
    plot_compare(ax3, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    # ax4: lobSTR vs TREDPARSE with diploid model (depth=80x)
    lobstr_results = parse_results("lobstr_results_het-80x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het-80x-150bp-500bp.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($Depth=%s\times$)" % 80
    plot_compare(ax4, title, tredparse_results, lobstr_results,
                 max_insert=max_insert)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = .03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2., 1 - pad, "B"),
                        (pad / 2, 1 / 2. , "C"), (1 / 2., 1 / 2. , "D")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == '__main__':
    main()
