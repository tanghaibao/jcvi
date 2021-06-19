#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Related scripts for the HLI-STR (TREDPARSE) paper.
"""
import os.path as op
import os
import csv
import sys
import vcf
import logging
import shutil
import json
import numpy as np
import pandas as pd

from collections import defaultdict
from random import sample
from pyfaidx import Fasta
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from itertools import product
from natsort import natsorted

from jcvi.graphics.base import (
    FancyArrow,
    normalize_axes,
    panel_labels,
    plt,
    savefig,
    set_helvetica_axis,
)
from jcvi.formats.base import is_number, must_open
from jcvi.formats.sam import get_minibam_bed, index
from jcvi.variation.str import TREDsRepo, af_to_counts, read_treds
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
SIMULATED_HAPLOID = r"Simulated haploid $\mathit{h}$"
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
        print(samplekey, res, ci)


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

    def __eq__(self, other):
        return self.__key() == other.__key()

    def __str__(self):
        return str(self.parents) + "=>" + str(self.child)

    __repr__ = __str__

    def check_mendelian(self, df, tred, tolerance=0, x_linked=False, verbose=False):
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
            possible_progenies = get_progenies(
                p1, p2, x_linked=x_linked, tolerance=tolerance
            )
            mendelian_error = not (c in possible_progenies)
            if verbose:
                print(
                    parent_keys[0],
                    parent_keys[1],
                    child_key,
                    p1,
                    p2,
                    c,
                    not mendelian_error,
                )
        else:
            parent_key = self.parents.values()[0]
            p1 = get_alleles(df, parent_key, tred)
            if p1 is None:
                return 0
            _p1 = expand_alleles(p1, tolerance=tolerance)
            mendelian_error = len(set(_p1) & set(c)) == 0
            if mendelian_error and x_linked:
                # Do not count case where - progeny is male, parent is male
                if (c[0] == c[1]) and (p1[0] == p1[1]):
                    mendelian_error = 0
            if verbose:
                print(parent_key, child_key, p1, c, not mendelian_error)
        return mendelian_error


def expand_alleles(p, tolerance=0):
    """
    Returns expanded allele set given the tolerance.
    """
    _p = set()
    for x in p:
        _p |= set(range(x - tolerance, x + tolerance + 1))
    return _p


def get_progenies(p1, p2, x_linked=False, tolerance=0):
    """
    Returns possible progenies in a trio.
    """
    _p1 = expand_alleles(p1, tolerance=tolerance)
    _p2 = expand_alleles(p2, tolerance=tolerance)
    possible_progenies = set(tuple(sorted(x)) for x in product(_p1, _p2))
    if x_linked:  # Add all hemizygotes
        possible_progenies |= set((x, x) for x in (set(_p1) | set(_p2)))
    return possible_progenies


def get_alleles(df, sample, tred):
    try:
        s = df.ix[sample]
        a = int(s[tred + ".1"])
        b = int(s[tred + ".2"])
    except:
        return None
    if a == -1 or b == -1:
        return None
    return a, b


def main():

    actions = (
        # Prepare data
        ("simulate", "simulate bams with varying inserts with dwgsim"),
        ("mergebam", "merge sets of BAMs to make diploid"),
        ("mini", "prepare mini-BAMs that contain only the STR loci"),
        ("alts", "build alternative loci based on simulation data"),
        # Compile results
        ("batchlobstr", "run lobSTR on a list of BAMs"),
        ("compilevcf", "compile vcf outputs into lists"),
        # Plotting
        ("evidences", "plot distribution of evidences"),
        ("likelihood", "plot likelihood surface"),
        ("likelihood2", "plot likelihood surface and marginals"),
        ("likelihood3", "plot likelihood surface and marginals for two settings"),
        ("compare", "compare callers on fake HD patients"),
        ("compare2", "compare TREDPARSE and lobSTR on fake HD patients"),
        ("power", "compare TREDPARSE on fake HD patients adding evidence"),
        ("tredparse", "compare TREDPARSE on fake HD patients adding coverage"),
        ("allelefreq", "plot the allele frequencies of some STRs"),
        ("allelefreqall", "plot all 30 STR allele frequencies"),
        ("depth", "plot read depths across all TREDs"),
        # Diagram
        ("diagram", "plot the predictive power of various evidences"),
        # Extra analysis for reviews
        ("mendelian", "calculate Mendelian errors based on trios and duos"),
        ("mendelian2", "second iteration of Mendelian error calculation"),
        ("mendelian_errors", "plot Mendelian errors calculated by mendelian"),
        ("mendelian_errors2", "plot Mendelian errors calculated by mendelian2"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def mendelian_errors2(args):
    """
    %prog mendelian_errors2 Trios.summary.csv

    Plot Mendelian errors as calculated by mendelian(). File
    `Trios.summary.csv` looks like:

    Name,Motif,Inheritance,N_Correct,N_Error,N_missing,ErrorRate [N_Error / (N_Correct + N_Error))]
    DM1,CTG,AD,790,12,0,1.5%
    DM2,CCTG,AD,757,45,0,5.6%
    DRPLA,CAG,AD,791,11,0,1.4%
    """
    p = OptionParser(mendelian_errors2.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="7x7", format="png")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (csvfile,) = args
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    ymin = -0.2
    df = pd.read_csv(csvfile)
    data = []
    for i, d in df.iterrows():
        tred = d["Name"]
        motif = d["Motif"]
        if tred in ignore:
            logging.debug("Ignore {}".format(d["TRED"]))
            continue

        if len(motif) > 6:
            if "/" in motif:  # CTG/CAG
                motif = motif.split("/")[0]
            else:
                motif = motif[:6] + ".."
        xtred = "{} {}".format(tred, motif)
        accuracy = d[-1]
        data.append((xtred, accuracy))

    key = lambda x: float(x.rstrip("%"))
    data.sort(key=lambda x: key(x[-1]))
    print(data)
    treds, accuracies = zip(*data)
    ntreds = len(treds)
    ticks = range(ntreds)
    accuracies = [key(x) for x in accuracies]

    for tick, accuracy in zip(ticks, accuracies):
        ax.plot([tick, tick], [ymin, accuracy], "-", lw=2, color="lightslategray")

    (trios,) = ax.plot(accuracies, "o", mfc="w", mec="b")
    ax.set_title("Mendelian errors based on STR calls in trios in HLI samples")
    ntrios = "Mendelian errors in 802 trios"
    ax.legend([trios], [ntrios], loc="best")

    ax.set_xticks(ticks)
    ax.set_xticklabels(treds, rotation=45, ha="right", size=8)
    ax.set_yticklabels([int(x) for x in ax.get_yticks()], family="Helvetica")
    ax.set_ylabel(r"Mendelian errors (\%)")
    ax.set_ylim(ymin, 100)

    normalize_axes(root)

    image_name = "mendelian_errors2." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def mendelian2(args):
    """
    %prog mendelian2
        XC_kinship_TRIO_annotationed_age_sex_PaternalMaternalAgeWhenChildWasBorn.txt
        hli.20170805.tsv

    Second iteration of Mendelian error calculation. This includes all the read
    counts and gender information to correct error estimate of X-linked loci.
    """
    p = OptionParser(mendelian2.__doc__)
    p.add_option(
        "--treds", default=None, help="Extract specific treds, use comma to separate"
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    triofile, hlitsv = args
    repo = TREDsRepo()
    treds = opts.treds.split(",") if opts.treds else repo.names
    triodata = pd.read_csv(triofile, sep="\t")
    samplekey = lambda x: x.split("_")[1]
    trios = []
    for i, row in triodata.iterrows():
        proband = row["proband"]
        parents = row["parents"]
        proband_sex = row["proband_sex"]
        parents_sex = row["parent1_sex,parent2_sex"]
        proband = samplekey(proband)
        p1, p2 = parents.split(",")
        p1, p2 = samplekey(p1), samplekey(p2)
        p1_sex, p2_sex = parents_sex.split(",")
        if p1_sex == "Male":
            p1, p2 = p2, p1
            p1_sex, p2_sex = p2_sex, p1_sex
        trios.append((proband, proband_sex, p1, p1_sex, p2, p2_sex))

    header = "{0}_ID {0}_Sex {0}_Calls"
    header += " {0}_Full {0}_Partial {0}_Repeat {0}_Paired"
    tredsdata = pd.read_csv(hlitsv, sep="\t", low_memory=False)
    tsvfiles = []
    summary = open("Trios.summary.csv", "w")
    summary_header = (
        "Name,Motif,Inheritance,N_Correct,N_Error,N_missing,"
        "ErrorRate [N_Error / (N_Correct + N_Error))]"
    )
    print(summary_header, file=summary)
    print(summary_header)
    for tred in treds:
        if tred in ("FXS", "AR"):
            continue
        tr = repo[tred]
        tsvfile = "{}.details.tsv".format(tred)
        fw = open(tsvfile, "w")
        td = {}
        for _, row in tredsdata.iterrows():
            s = str(row["SampleKey"])
            inferredGender = row["inferredGender"]
            try:
                calls = row[tred + ".calls"]
                fdp = int(row[tred + ".FDP"])
                pdp = int(row[tred + ".PDP"])
                rdp = int(row[tred + ".RDP"])
                pedp = int(row[tred + ".PEDP"])
                td[s] = [str(x) for x in (inferredGender, calls, fdp, pdp, rdp, pedp)]
            except ValueError:
                logging.error("Invalid row: {}".format(row))
                continue

        h = " ".join((header.format("P1"), header.format("P2"), header.format("Kid")))
        print("\t".join(["MendelianError"] + h.split()), file=fw)
        tredcall = lambda x: td.get(x, ["", "-1|-1", "", "", "", ""])[:]
        counts = defaultdict(int)
        is_xlinked = repo[tred].is_xlinked
        shorten = lambda x: str(int(x[-4:]))  # Simplify SampleKey
        for proband, proband_sex, p1, p1_sex, p2, p2_sex in trios:
            tp1 = tredcall(p1)
            tp2 = tredcall(p2)
            tpp = tredcall(proband)
            m = mendelian_check(tp1, tp2, tpp, is_xlinked=is_xlinked)
            counts[m] += 1
            if is_xlinked:
                for (p, p_sex) in ((tp1, p1_sex), (tp2, p2_sex), (tpp, proband_sex)):
                    if p[1].startswith("-"):
                        p[1] = "n.a."
            cells = [shorten(p1), p1_sex] + tp1[1:]
            cells += [shorten(p2), p2_sex] + tp2[1:]
            cells += [shorten(proband), proband_sex] + tpp[1:]
            print("\t".join([m] + cells), file=fw)
        fw.close()
        tsvfiles.append(tsvfile)

        error_rate = counts["Error"] * 100.0 / (counts["Correct"] + counts["Error"])
        line = ",".join(
            str(x)
            for x in (
                tred,
                tr.motif,
                tr.inheritance,
                counts["Correct"],
                counts["Error"],
                counts["Missing"],
                "{:.1f}%".format(error_rate),
            )
        )
        print(line, file=summary)
        print(line)

    # Combine into a master spreadsheet
    import xlwt

    wb = xlwt.Workbook()
    converter = lambda x: int(x) if is_number(x, cast=int) else x
    header = xlwt.easyxf("font: bold on, name Helvetica; align: horiz center")
    hc = "font: name Helvetica; align: horiz center;"
    horiz_center = xlwt.Style.easyxf(hc)
    correct = xlwt.Style.easyxf(hc + "pattern: pattern solid, fore_colour light_green;")
    error = xlwt.Style.easyxf(hc + "pattern: pattern solid, fore_colour rose;")
    missing = xlwt.Style.easyxf(
        hc + "pattern: pattern solid, fore_colour light_yellow;"
    )
    for tsvfile in tsvfiles:
        sheet = op.basename(tsvfile).split(".", 1)[0]
        ws = wb.add_sheet(sheet)
        fp = open(tsvfile, "rb")
        reader = csv.reader(fp, delimiter="\t")
        for r, row in enumerate(reader):
            style = header if r == 0 else horiz_center
            for c, col in enumerate(row):
                if c == 0 and r:
                    style = {"Correct": correct, "Error": error, "Missing": missing}[
                        col
                    ]
            ws.write(r, c, converter(col), style)
            ws.set_panes_frozen(True)
            ws.set_horz_split_pos(1)

    wb.save("Trios.xls")
    summary.close()


def mendelian_check(tp1, tp2, tpp, is_xlinked=False):
    """
    Compare TRED calls for Parent1, Parent2 and Proband.
    """
    call_to_ints = lambda x: tuple(int(_) for _ in x.split("|") if _ != ".")
    tp1_sex, tp1_call = tp1[:2]
    tp2_sex, tp2_call = tp2[:2]
    tpp_sex, tpp_call = tpp[:2]
    # tp1_evidence = sum(int(x) for x in tp1[2:])
    # tp2_evidence = sum(int(x) for x in tp2[2:])
    # tpp_evidence = sum(int(x) for x in tpp[2:])
    tp1_call = call_to_ints(tp1_call)
    tp2_call = call_to_ints(tp2_call)
    tpp_call = call_to_ints(tpp_call)
    possible_progenies = set(tuple(sorted(x)) for x in product(tp1_call, tp2_call))
    if is_xlinked and tpp_sex == "Male":
        possible_progenies = set(tuple((x,)) for x in tp1_call)
    if -1 in tp1_call or -1 in tp2_call or -1 in tpp_call:
        tag = "Missing"
    else:
        tag = "Correct" if tpp_call in possible_progenies else "Error"
    return tag


def in_region(rname, rstart, target_chr, target_start, target_end):
    """
    Quick check if a point is within the target region.
    """
    return (rname == target_chr) and (target_start <= rstart <= target_end)


def alts(args):
    """
    %prog alts HD

    Build alternative loci based on simulation data.
    """
    import pysam
    from more_itertools import pairwise
    from jcvi.utils.grouper import Grouper

    p = OptionParser(alts.__doc__)
    p.set_outfile(outfile="TREDs.alts.csv")
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    treds = args
    repo = TREDsRepo()
    if "all" in treds:
        treds = repo.names

    pad_left, pad_right = 1000, 10000
    READLEN = 150
    fw = must_open(opts.outfile, "w")
    print("TRED,alts,alts.hg19", file=fw)  # Header
    for tred in treds:
        ref_regions = []

        # Simulate a depth 1000 BAM with 300 repeats
        for ref in ("hg38", "hg19"):

            # This is the region that involves the TRED locus
            repo = TREDsRepo(ref=ref)
            t = repo[tred]
            chr, start, end = t.chr, t.repeat_start, t.repeat_end
            start -= pad_left
            end += pad_right

            tred_ref = "{}_{}".format(tred, ref)
            if not op.isdir(tred_ref):
                simulate(
                    [
                        tred_ref,
                        "300",
                        "300",
                        "--depth=1000",
                        "--ref={}".format(ref),
                        "--tred={}".format(tred),
                    ]
                )
            bamfile = op.join(tred_ref, "300.bam")

            # Parse the BAM file, retrieve all regions
            bamfile = pysam.AlignmentFile(bamfile, "rb")
            nreads = altreads = 0
            alt_points = set()
            for read in bamfile.fetch():
                fname, fstart = (
                    bamfile.getrname(read.reference_id),
                    read.reference_start,
                )
                rname, rstart = (
                    bamfile.getrname(read.next_reference_id),
                    read.next_reference_start,
                )
                f_in_region = in_region(fname, fstart, chr, start, end)
                r_in_region = in_region(rname, rstart, chr, start, end)
                if (not f_in_region) and r_in_region:
                    alt_points.add((fname, fstart))
                    altreads += 1
                if (not r_in_region) and f_in_region:
                    alt_points.add((rname, rstart))
                    altreads += 1
                nreads += 1

            logging.debug(
                "A total of {} reads ({} alts) processed".format(nreads, altreads)
            )
            alt_points = natsorted(alt_points)

            # Chain these points together into regions
            g = Grouper()
            for a in alt_points:
                g.join(a)
            for a, b in pairwise(alt_points):
                achr, apos = a
                bchr, bpos = b
                if achr != bchr:
                    continue
                if (bpos - apos) > READLEN:
                    continue
                g.join(a, b)

            # All regions that contain ALT
            alt_sum = 0
            regions = []
            for c in g:
                chr_min, pos_min = min(c)
                chr_max, pos_max = max(c)
                assert chr_min, chr_max
                pos_min -= READLEN
                pos_max += READLEN
                regions.append((chr_min, pos_min, pos_max))
                alt_sum += pos_max - pos_min

            regions = "|".join(
                [
                    "{}:{}-{}".format(c, start, end)
                    for c, start, end in natsorted(regions)
                ]
            )
            ref_regions.append(regions)

        line = ",".join([tred] + ref_regions)
        print(line, file=sys.stderr)
        print(line, file=fw)
        logging.debug("Alternative region sum: {} bp".format(alt_sum))

    fw.close()


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

    (tsvfile,) = args
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, figsize=(iopts.w, iopts.h)
    )
    plt.tight_layout(pad=6)

    data = pd.read_csv(tsvfile, sep="\t", low_memory=False)

    ids, treds = read_treds()
    for (dp, ax, title) in zip(
        ("FDP", "PDP", "RDP", "PEDP"),
        (ax1, ax2, ax3, ax4),
        ("Spanning reads", "Partial reads", "Repeat-only reads", "Paired-end reads"),
    ):
        logging.debug("Build {}".format(title))
        # Construct related data structure
        xd = []  # (tred, dp)
        mdp = []  # (tred, median_dp)
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
            md = [x for x in data[tred + "." + dp] if x >= 0]
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
        ax.set_yticklabels(yticklabels, family="Helvetica", size=14)

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.04
    panel_labels(
        root,
        (
            (pad, 1 - pad, "A"),
            (1 / 2.0 + pad / 2, 1 - pad, "B"),
            (pad, 0.5 - pad / 2, "C"),
            (1 / 2.0 + pad / 2, 0.5 - pad / 2, "D"),
        ),
    )
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

    (csvfile,) = args
    fig, ax = plt.subplots(ncols=1, nrows=1, figsize=(iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    ymin = -0.2
    df = pd.read_csv(csvfile)
    data = []
    for i, d in df.iterrows():
        if d["TRED"].split()[0] in ignore:
            logging.debug("Ignore {}".format(d["TRED"]))
            continue
        data.append(d)
    treds, duos, trios = zip(*data)
    ntreds = len(treds)
    ticks = range(ntreds)
    treds = [x.split()[0] for x in treds]
    duos = [float(x.rstrip("%")) for x in duos]
    trios = [float(x.rstrip("%")) for x in trios]

    for tick, duo, trio in zip(ticks, duos, trios):
        m = max(duo, trio)
        ax.plot([tick, tick], [ymin, m], "-", lw=2, color="lightslategray")

    (duos,) = ax.plot(duos, "o", mfc="w", mec="g")
    (trios,) = ax.plot(trios, "o", mfc="w", mec="b")
    ax.set_title("Mendelian errors based on trios and duos in HLI samples")
    nduos = "Mendelian errors in 362 duos"
    ntrios = "Mendelian errors in 339 trios"
    ax.legend([trios, duos], [ntrios, nduos], loc="best")

    ax.set_xticks(ticks)
    ax.set_xticklabels(treds, rotation=45, ha="right", size=8)
    yticklabels = [int(x) for x in ax.get_yticks()]
    ax.set_yticklabels(yticklabels, family="Helvetica")
    ax.set_ylabel(r"Mendelian errors (\%)")
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
    parent_keys = [
        x for x in keys if ("parent" in x.lower()) and ("grand" not in x.lower())
    ]
    sib_keys = [
        x for x in keys if ("sibling" in x.lower()) or ("twin" in x.lower())
    ] + self_key
    child_keys = [
        x
        for x in keys
        if ("child" in x.lower())
        and ("grand" not in x.lower())
        and ("self" not in x.lower())
    ]

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
    p.add_option("--tolerance", default=0, type="int", help="Tolernace for differences")
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    triosjson, tredtsv = args
    verbose = opts.verbose
    tolerance = opts.tolerance

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
    # print "\n".join(allterms)
    print("A total of {} families imported".format(len(js)))

    # Read in all data
    df = read_tred_tsv(tredtsv)

    ids, treds = read_treds()
    table = {}
    for tred, inheritance in zip(treds["abbreviation"], treds["inheritance"]):
        x_linked = inheritance[0] == "X"  # X-linked
        name = tred
        if x_linked:
            name += " (X-linked)"
        print("[TRED] {}".format(name))

        n_total = len(duos)
        n_error = 0
        for duo in duos:
            n_error += duo.check_mendelian(
                df, tred, tolerance=tolerance, x_linked=x_linked, verbose=verbose
            )
        tag = "Duos  - Mendelian errors"
        print("{}: {}".format(tag, percentage(n_error, n_total)))
        duo_error = percentage(n_error, n_total, mode=2)
        table[(name, tag)] = "{0:.1f}%".format(duo_error)

        n_total = len(trios)
        n_error = 0
        for trio in trios:
            n_error += trio.check_mendelian(
                df, tred, tolerance=tolerance, x_linked=x_linked, verbose=verbose
            )
        tag = "Trios - Mendelian errors"
        print("{}: {}".format(tag, percentage(n_error, n_total)))
        trio_error = percentage(n_error, n_total, mode=2)
        table[(name, tag)] = "{0:.1f}%".format(trio_error)

    # Summarize
    print(tabulate(table))


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
        print("\t".join(str(x) for x in (c, start - pad, end + pad, td)), file=fw)

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
        print(
            "\t".join(
                str(x)
                for x in (
                    c,
                    start - pad,
                    end + pad,
                    "chrY.unique_ccn.{}".format(nregions),
                )
            ),
            file=fw,
        )
        nregions += 1

    fw.close()
    return filename


def mini(args):
    """
    %prog mini bamfile minibamfile

    Prepare mini-BAMs that contain only the STR loci.
    """
    p = OptionParser(mini.__doc__)
    p.add_option(
        "--pad", default=20000, type="int", help="Add padding to the STR reigons"
    )
    p.add_option(
        "--treds", default=None, help="Extract specific treds, use comma to separate"
    )
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
            i = int(atoms[1].strip("(,"))
            j = int(atoms[2].strip(")"))
            lnL = float(atoms[-1])
            likelihood[(i, j)] = lnL
        if row.startswith("DEBUG:IntegratedCaller:CI(h1)"):
            CI_h1 = [int(x.strip()) for x in row.split("=")[1].split("-")]
        if row.startswith("DEBUG:IntegratedCaller:CI(h2)"):
            CI_h2 = [int(x.strip()) for x in row.split("=")[1].split("-")]
        if row.startswith("DEBUG:IntegratedCaller:ML estimate:"):
            MLE = row.split(":")[3].split("=")[1].split()[:2]
            MLE = [int(x.strip("[],")) for x in MLE]

    return likelihood, CI_h1, CI_h2, MLE


def likelihood(args):
    """
    %prog likelihood

    Plot likelihood surface. Look for two files in the current folder:
    - 100_100.log, haploid model
    - 100_20.log, diploid model
    """
    p = OptionParser(likelihood.__doc__)
    opts, args, iopts = p.set_image_options(
        args, figsize="10x5", style="white", cmap="coolwarm"
    )

    if len(args) != 0:
        sys.exit(not p.print_help())

    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=4)

    # Haploid model
    LL, CI_h1, CI_h2, MLE = parse_log("100_100.log")
    data = []
    for k, v in LL.items():
        data.append((k[0], v))
    data.sort()
    x, y = zip(*data)
    x = np.array(x)
    (curve,) = ax1.plot(x, y, "-", color=lsg, lw=2)
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
    ci = ax1.fill_between(
        x, [ymin] * len(x), y, where=(x >= a) & (x <= b), color=lsg, alpha=0.5
    )
    ax1.legend([curve, ci], ["Likelihood curve", r"95$\%$ CI"], loc="best")

    # Diploid model
    LL, CI_h1, CI_h2, MLE = parse_log("100_20.log")
    _, min_LL = min(data, key=lambda x: x[-1])
    data = np.ones((301, 301)) * min_LL
    for k, v in LL.items():
        a, b = k
        data[a, b] = v
        data[b, a] = v

    data = mask_upper_triangle(data)
    ax_imshow(ax2, data, opts.cmap, LL_label, 20, 104)

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.04
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2.0, 1 - pad, "B")))
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
    ax.plot(x, y, "-", color=lsg, lw=2)
    title = "Marginal distribution for $%s$" % xlabel
    ax.set_title(title)
    if not ticks:
        ax.set_yticks([])

    if a == b:
        ax.plot([h_hat, h_hat], [0, max_P], "-", color=lsg, lw=2)
    else:
        ax.fill_between(
            x, [0] * len(x), y, where=(x >= a) & (x <= b), color=lsg, alpha=0.5
        )
    ax.set_xlim(0, 300)

    ymin, ymax = ax.get_ylim()
    if h_hat < 150:
        ax.text(
            h_hat + 20,
            ymax * 4.0 / 5,
            r"$\hat{%s}=%d$" % (xlabel, h_hat),
            color=lsg,
            va="center",
        )
        ax.text(
            h_hat + 20,
            ymax * 3.0 / 5,
            r"95$\%$ CI" + r"$=%s-%s$" % (a, b),
            color=lsg,
            va="center",
        )
    else:
        ax.text(
            h_hat - 30,
            ymax * 4.0 / 5,
            r"$\hat{%s}=%d$" % (xlabel, h_hat),
            color=lsg,
            ha="right",
            va="center",
        )
        ax.text(
            h_hat - 30,
            ymax * 3.0 / 5,
            r"95$\%$ CI" + r"$=%s-%s$" % (a, b),
            color=lsg,
            ha="right",
            va="center",
        )

    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.05)


def ax_imshow(
    ax,
    P_h1h2,
    cmap,
    label,
    h1_hat,
    h2_hat,
    h1_truth,
    h2_truth,
    r=4,
    draw_circle=True,
    ticks=True,
):
    im = ax.imshow(P_h1h2, cmap=cmap, origin="lower")

    from mpl_toolkits.axes_grid1 import make_axes_locatable

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cb = plt.colorbar(im, cax)
    cb.set_label(label)
    if not ticks:
        cb.set_ticks([])

    if draw_circle:
        circle = plt.Circle((h1_hat, h2_hat), r, ec="w", fill=False)
        ax.add_artist(circle)

    annotation = r"$\hat{h_1}=%d, \hat{h_2}=%d$" % (h1_hat, h2_hat)
    ax.text(200, 100, annotation, color=lsg, ha="center", va="center")

    ax.set_xlabel(r"$h_1$")
    ax.set_ylabel(r"$h_2$")
    title = "Simulated diploid ($h_{1}^{truth}=%d, h_{2}^{truth}=%d$)" % (
        h1_truth,
        h2_truth,
    )
    ax.set_title(title)


def likelihood2(args):
    """
    %prog likelihood2 100_20.json

    Plot the likelihood surface and marginal distributions.
    """
    from matplotlib import gridspec

    p = OptionParser(likelihood2.__doc__)
    opts, args, iopts = p.set_image_options(
        args, figsize="10x5", style="white", cmap="coolwarm"
    )

    if len(args) != 1:
        sys.exit(not p.print_help())

    (jsonfile,) = args
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
    %prog likelihood3 140_20.json 140_70.json

    Plot the likelihood surface and marginal distributions for two settings.
    """
    from matplotlib import gridspec

    p = OptionParser(likelihood3.__doc__)
    opts, args, iopts = p.set_image_options(
        args, figsize="10x10", style="white", cmap="coolwarm"
    )
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
    pad = 0.02
    panel_labels(root, ((pad, 1 - pad, "A"), (pad, 4.0 / 9, "B")))
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
    ax_imshow(
        ax1,
        data,
        cmap,
        label,
        h1_hat,
        h2_hat,
        h1_truth,
        h2_truth,
        draw_circle=False,
        ticks=False,
    )

    CI = calls[tred + ".CI"]
    CI_h1, CI_h2 = CI.split("|")
    CI_h1 = [int(x) for x in CI_h1.split("-")]
    CI_h2 = [int(x) for x in CI_h2.split("-")]
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
    yy = 0.7
    yinterval = 0.1
    height = 0.05
    yp = yy - yinterval - height
    canvas = 0.95
    xstart = 0.025
    convert = lambda x: xstart + x * canvas / 600
    # Symbols
    root.text(
        0.5, 0.9, r"$L$: Read length, $F$: Flank size, $V$: Pair distance", ha="center"
    )
    root.text(0.5, 0.85, r"ex. $L=150bp, F=9bp, V=500bp$", ha="center")
    root.text(
        xstart + canvas,
        yy - height,
        "STR repeat length",
        ha="center",
        color=lsg,
        size=10,
    )

    # Mark the key events
    pad = 0.02
    arrowlen = canvas * 1.05
    arrowprops = dict(
        length_includes_head=True,
        width=0.01,
        fc=lsg,
        lw=0,
        head_length=arrowlen * 0.12,
        head_width=0.04,
    )
    p = FancyArrow(xstart, yy, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    ppad = 30
    keyevents = (
        (0, 0, -1, r"$0$"),
        (150 - 18, 150 - 18 - ppad, 0, r"$L - 2F$"),
        (150 - 9, 150 - 9, 1, r"$L - F$"),
        (150, 150 + ppad, 2, r"$L$"),
        (500 - 9, 500 - 9, 3, r"$V - F$"),
    )
    for event, pos, i, label in keyevents:
        _event = convert(event)
        _pos = convert(pos)
        root.plot((_event, _event), (yy - height / 4, yy + height / 4), "-", color="k")
        root.text(_pos, yy + pad, label, rotation=45, va="bottom", size=8)
        if i < 0:
            continue
        ystart = yp - i * yinterval
        root.plot((_event, _event), (ystart, yy - height / 4), ":", color=lsg)

    # Range on bottom. These are simple 4 rectangles, with the range indicating
    # the predictive range.
    CLOSED, OPEN = range(2)
    ranges = (
        (0, 150 - 18, CLOSED, "Spanning reads"),
        (9, 150 - 9, OPEN, "Partial reads"),
        (150, 500 - 9, CLOSED, "Repeat reads"),
        (0, 500 - 9, CLOSED, "Paired-end reads"),
    )
    for start, end, starttag, label in ranges:
        _start = convert(start)
        _end = convert(end)
        data = (
            [[0.0, 1.0], [0.0, 1.0]] if starttag == OPEN else [[1.0, 0.0], [1.0, 0.0]]
        )
        root.imshow(
            data,
            interpolation="bicubic",
            cmap=plt.cm.Greens,
            extent=[_start, _end, yp, yp + height],
        )
        root.text(_end + pad, yp + height / 2, label, va="center")
        yp -= yinterval

    normalize_axes(root)

    image_name = "diagram." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_allelefreq(ax, df, locus, color="lightslategray"):
    tred = df.ix[locus]
    cnt = af_to_counts(tred["allele_freq"])

    cntx, cnty = zip(*cnt.items())

    motif = tred["motif"]
    cutoff_prerisk = tred["cutoff_prerisk"]
    cutoff_risk = tred["cutoff_risk"]
    npredisease = sum(v for (k, v) in cnt.items() if cutoff_prerisk <= k < cutoff_risk)
    npatients = sum(v for (k, v) in cnt.items() if k >= cutoff_risk)

    ax.bar(cntx, cnty, fc=color)

    ymin, ymax = ax.get_ylim()
    xmax = (cutoff_risk / 10 + 1) * 10 if cutoff_risk > 50 else 50
    pad = xmax * 0.03
    if cutoff_prerisk < cutoff_risk and npredisease:
        ax.axvline(x=cutoff_prerisk, color="k", lw=2)
        ax.text(
            cutoff_prerisk + pad,
            0.5 * ymax,
            r"Pre-disease ($\geq${}$\times${}) - {} alleles".format(
                cutoff_prerisk, motif, npredisease
            ),
            rotation=90,
            color="k",
            ha="center",
            va="center",
        )
    ax.axvline(x=cutoff_risk, color="r", lw=2)

    if locus == "AR":
        npatients = sum(v for (k, v) in cnt.items() if k <= cutoff_risk)
        ax.text(
            cutoff_risk - pad,
            0.5 * ymax,
            r"Disease ($\leq${}$\times${}) - {} alleles".format(
                cutoff_risk, motif, npatients
            ),
            rotation=90,
            color="r",
            ha="center",
            va="center",
        )
    else:
        ax.text(
            cutoff_risk + pad,
            0.5 * ymax,
            r"Disease ($\geq${}$\times${}) - {} alleles".format(
                cutoff_risk, motif, npatients
            ),
            rotation=90,
            color="r",
            ha="center",
            va="center",
        )

    x = []  # All allelels
    for k, v in cnt.items():
        x.extend([k] * v)

    ax.set_xlabel("Number of repeat units")
    ax.set_ylabel("Number of alleles")
    ax.set_xlim(0, xmax)
    ax.set_title(r"{} ({})".format(locus, tred["title"], motif))
    set_helvetica_axis(ax)


def allelefreqall(args):
    """
    %prog allelefreqall HN_Platinum_Gold.20180525.tsv.report.txt

    Plot all 30 STR allele frequencies.
    """
    p = OptionParser(allelefreqall.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (reportfile,) = args
    treds, df = read_treds(reportfile)
    # Prepare 5 pages, each page with 6 distributions
    treds = sorted(treds)
    count = 6
    pdfs = []
    for page in range(len(treds) / count + 1):
        start = page * count
        page_treds = treds[start : start + count]
        if not page_treds:
            break
        allelefreq(
            [
                ",".join(page_treds),
                "--usereport",
                reportfile,
                "--nopanels",
                "--figsize",
                "12x16",
            ]
        )
        outpdf = "allelefreq.{}.pdf".format(page)
        sh("mv allelefreq.pdf {}".format(outpdf))
        pdfs.append(outpdf)

    from jcvi.formats.pdf import cat

    pf = op.basename(reportfile).split(".")[0]
    finalpdf = pf + ".allelefreq.pdf"
    logging.debug("Merging pdfs into `{}`".format(finalpdf))
    cat(pdfs + ["-o", finalpdf, "--cleanup"])


def allelefreq(args):
    """
    %prog allelefreq HD,DM1,SCA1,SCA17,FXTAS,FRAXE

    Plot the allele frequencies of some STRs.
    """
    p = OptionParser(allelefreq.__doc__)
    p.add_option(
        "--nopanels",
        default=False,
        action="store_true",
        help="No panel labels A, B, ...",
    )
    p.add_option("--usereport", help="Use allele frequency in report file")
    opts, args, iopts = p.set_image_options(args, figsize="9x13")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (loci,) = args
    fig, ((ax1, ax2), (ax3, ax4), (ax5, ax6)) = plt.subplots(
        ncols=2, nrows=3, figsize=(iopts.w, iopts.h)
    )
    plt.tight_layout(pad=4)
    if opts.usereport:
        treds, df = read_treds(tredsfile=opts.usereport)
    else:
        treds, df = read_treds()

    df = df.set_index(["abbreviation"])

    axes = (ax1, ax2, ax3, ax4, ax5, ax6)
    loci = loci.split(",")
    for ax, locus in zip(axes, loci):
        plot_allelefreq(ax, df, locus)

    # Delete unused axes
    for ax in axes[len(loci) :]:
        ax.set_axis_off()

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.03
    if not opts.nopanels:
        panel_labels(
            root,
            (
                (pad / 2, 1 - pad, "A"),
                (0.5 + pad, 1 - pad, "B"),
                (pad / 2, 2 / 3.0 - pad / 2, "C"),
                (0.5 + pad, 2 / 3.0 - pad / 2, "D"),
                (pad / 2, 1 / 3.0, "E"),
                (0.5 + pad, 1 / 3.0, "F"),
            ),
        )
    normalize_axes(root)

    image_name = "allelefreq." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def make_fasta(seq, fastafile, id):
    rec = SeqRecord(Seq(seq), description="", id=id)
    fw = open(fastafile, "w")
    SeqIO.write([rec], fw, "fasta")
    fw.close()


def add_simulate_options(p):
    p.add_option("--readlen", default=150, type="int", help="Length of the read")
    p.add_option(
        "--distance",
        default=500,
        type="int",
        help="Outer distance between the two ends",
    )
    p.set_depth(depth=20)


def simulate(args):
    """
    %prog simulate run_dir 1 300

    Simulate BAMs with varying inserts with dwgsim. The above command will
    simulate between 1 to 300 CAGs in the HD region, in a directory called
    `run_dir`.
    """
    p = OptionParser(simulate.__doc__)
    p.add_option(
        "--method", choices=("wgsim", "eagle"), default="eagle", help="Read simulator"
    )
    p.add_option(
        "--ref",
        default="hg38",
        choices=("hg38", "hg19"),
        help="Reference genome version",
    )
    p.add_option("--tred", default="HD", help="TRED locus")
    add_simulate_options(p)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    rundir, startunits, endunits = args
    ref = opts.ref
    ref_fasta = "/mnt/ref/{}.upper.fa".format(ref)
    startunits, endunits = int(startunits), int(endunits)
    basecwd = os.getcwd()
    mkdir(rundir)
    os.chdir(rundir)
    cwd = os.getcwd()

    # TRED region (e.g. Huntington)
    pad_left, pad_right = 1000, 10000
    repo = TREDsRepo(ref=ref)
    tred = repo[opts.tred]
    chr, start, end = tred.chr, tred.repeat_start, tred.repeat_end

    logging.debug("Simulating {}".format(tred))
    fasta = Fasta(ref_fasta)
    seq_left = fasta[chr][start - pad_left : start - 1]
    seq_right = fasta[chr][end : end + pad_right]
    motif = tred.repeat

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
        simulate_method(
            [
                fastafile,
                "--depth={}".format(opts.depth),
                "--readlen={}".format(opts.readlen),
                "--distance={}".format(opts.distance),
                "--outfile={}".format(pf),
            ]
        )

        read1 = pf + ".bwa.read1.fastq"
        read2 = pf + ".bwa.read2.fastq"
        samfile, _ = align([ref_fasta, read1, read2])
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
    p.add_option(
        "--haploid", default="chrY,chrM", help="Use haploid model for these chromosomes"
    )
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (bamlist,) = args
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

    (folder,) = args
    vcf_files = iglob(folder, "*.vcf,*.vcf.gz")
    for vcf_file in vcf_files:
        try:
            p = LobSTRvcf(columnidsfile=None)
            p.parse(vcf_file, filtered=False)
            res = p.items()
            if res:
                k, v = res[0]
                res = v.replace(",", "/")
            else:
                res = "-1/-1"
            num = op.basename(vcf_file).split(".")[0]
            print(num, res)
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
    p.add_option(
        "--csv", default="hli.20170328.tred.tsv", help="TRED csv output to plot"
    )
    opts, args, iopts = p.set_image_options(args, format="pdf")

    if len(args) != 0:
        sys.exit(not p.print_help())

    format = iopts.format

    # Extract sample coverage first
    df = pd.read_csv(
        "qc-export-MeanCoverage.csv",
        header=None,
        names=["Samplekey", "MeanCoverage"],
        index_col=0,
    )

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
    ef.columns = [
        tred + ".1",
        tred + ".2",
        tred + ".FDP",
        tred + ".PDP",
        tred + ".PEDP",
    ]
    ef.index.name = "SampleKey"
    mf = df.merge(ef, how="right", left_index=True, right_index=True)

    # Plot a bunch of figures
    outdir = "output"
    mkdir(outdir)
    xlim = ylim = (0, 100)
    draw_jointplot(
        outdir + "/A",
        "MeanCoverage",
        "HD.FDP",
        data=mf,
        xlim=xlim,
        ylim=ylim,
        format=format,
    )
    draw_jointplot(
        outdir + "/B",
        "MeanCoverage",
        "HD.PDP",
        data=mf,
        color="g",
        xlim=xlim,
        ylim=ylim,
        format=format,
    )
    draw_jointplot(
        outdir + "/C",
        "MeanCoverage",
        "HD.PEDP",
        data=mf,
        color="m",
        xlim=xlim,
        ylim=ylim,
        format=format,
    )

    xlim = (0, 50)
    draw_jointplot(
        outdir + "/D", "HD.2", "HD.FDP", data=mf, xlim=xlim, ylim=ylim, format=format
    )
    draw_jointplot(
        outdir + "/E",
        "HD.2",
        "HD.PDP",
        data=mf,
        color="g",
        xlim=xlim,
        ylim=ylim,
        format=format,
    )
    draw_jointplot(
        outdir + "/F",
        "HD.2",
        "HD.PEDP",
        data=mf,
        color="m",
        xlim=xlim,
        ylim=ylim,
        format=format,
    )


def draw_jointplot(
    figname, x, y, data=None, kind="reg", color=None, xlim=None, ylim=None, format="pdf"
):
    """
    Wraps around sns.jointplot
    """
    import seaborn as sns

    sns.set_context("talk")
    plt.clf()

    register = {
        "MeanCoverage": "Sample Mean Coverage",
        "HD.FDP": "Depth of full spanning reads",
        "HD.PDP": "Depth of partial spanning reads",
        "HD.PEDP": "Depth of paired-end reads",
        "HD.2": "Repeat size of the longer allele",
    }

    g = sns.jointplot(x, y, data=data, kind=kind, color=color, xlim=xlim, ylim=ylim)
    g.ax_joint.set_xlabel(register.get(x, x))
    g.ax_joint.set_ylabel(register.get(y, y))
    savefig(figname + "." + format, cleanup=False)


def long_allele(s, default=19, exclude=None):
    if "_" in s:
        a, b = s.split("_")
    elif "/" in s:
        a, b = s.split("/")
    else:
        raise Exception("Don't know how to split string {}".format(s))

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


def compute_rmsd(truth, a, limit=150):
    truth = truth[:limit]
    a = a[:limit]
    if len(a) > len(truth):
        a = a[: len(truth)]
    return (sum((i - j) ** 2 for (i, j) in zip(truth, a)) / len(truth)) ** 0.5


def compare(args):
    """
    %prog compare Evaluation.csv

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (datafile,) = args
    pf = datafile.rsplit(".", 1)[0]
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, figsize=(iopts.w, iopts.h)
    )
    plt.tight_layout(pad=3)

    bbox = {"facecolor": "tomato", "alpha": 0.2, "ec": "w"}
    pad = 2

    # Read benchmark data
    df = pd.read_csv("Evaluation.csv")
    truth = df["Truth"]
    axes = (ax1, ax2, ax3, ax4)
    progs = ("Manta", "Isaac", "GATK", "lobSTR")
    markers = ("bx-", "yo-", "md-", "c+-")

    for ax, prog, marker in zip(axes, progs, markers):
        ax.plot(truth, df[prog], marker)
        ax.plot(truth, truth, "k--")  # to show diagonal
        ax.axhline(infected_thr, color="tomato")
        ax.text(
            max(truth) - pad,
            infected_thr + pad,
            "Risk threshold",
            bbox=bbox,
            ha="right",
        )
        ax.axhline(ref_thr, color="tomato")
        ax.text(
            max(truth) - pad,
            ref_thr - pad,
            "Reference repeat count",
            bbox=bbox,
            ha="right",
            va="top",
        )
        ax.set_title(SIMULATED_HAPLOID)
        ax.set_xlabel(r"Num of CAG repeats inserted ($\mathit{h}$)")
        ax.set_ylabel("Num of CAG repeats called")
        ax.legend([prog, "Truth"], loc="best")

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.03
    panel_labels(
        root,
        (
            (pad / 2, 1 - pad, "A"),
            (1 / 2.0, 1 - pad, "B"),
            (pad / 2, 1 / 2.0, "C"),
            (1 / 2.0, 1 / 2.0, "D"),
        ),
    )
    normalize_axes(root)

    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def plot_compare(
    ax,
    title,
    tredparse_results,
    lobstr_results,
    pad=8,
    ms=3,
    max_insert=300,
    color="g",
    risk=True,
):
    truth = range(1, max_insert + 1)
    tx, ty, tl, th = zip(*tredparse_results)
    trmsd = compute_rmsd(truth, ty)
    if lobstr_results:
        lx, ly = zip(*lobstr_results)
        lrmsd = compute_rmsd(truth, ly)

    rmsd_tag = "$RMSD_{1:150}$"
    if lobstr_results:
        ax.plot(
            lx, ly, "c+-", ms=ms, label="lobSTR ({}={:.2f})".format(rmsd_tag, lrmsd)
        )
    ax.plot(
        tx,
        ty,
        ".-",
        color=color,
        ms=ms,
        label="TREDPARSE ({}={:.2f})".format(rmsd_tag, trmsd),
    )
    ax.plot(truth, truth, "k--", label="Truth")
    ax.fill_between(
        tx, tl, th, facecolor=color, alpha=0.25, label=r"TREDPARSE 95$\%$ CI"
    )

    ax.set_xlabel(r"Num of CAG repeats inserted ($\mathit{h}$)")
    ax.set_ylabel("Num of CAG repeats called")
    ax.set_title(title)
    ax.legend(loc="best")

    bbox = {"facecolor": "tomato", "alpha": 0.2, "ec": "w"}
    if risk:
        ax.axhline(infected_thr, color="tomato")
        ax.text(
            max(truth) - pad,
            infected_thr + pad,
            "Risk cutoff={}".format(infected_thr) + r"$\times$CAGs",
            bbox=bbox,
            ha="right",
        )
    else:
        readlength, pairdistance = 150 / 3, 500 / 3
        ax.axhline(readlength, color="tomato")
        ax.text(
            max(truth) - pad,
            readlength + pad,
            "Read Length ($L$)",
            bbox=bbox,
            ha="right",
        )
        ax.axhline(pairdistance, color="tomato")
        ax.text(
            max(truth) - pad,
            pairdistance + pad,
            "Paired-end distance($V$)",
            bbox=bbox,
            ha="right",
        )


def compare2(args):
    """
    %prog compare2

    Compare performances of various variant callers on simulated STR datasets.
    """
    p = OptionParser(compare2.__doc__)
    p.add_option(
        "--maxinsert", default=300, type="int", help="Maximum number of repeats"
    )
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x5")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    readlen = opts.readlen
    distance = opts.distance
    max_insert = opts.maxinsert
    fig, (ax1, ax2) = plt.subplots(ncols=2, nrows=1, figsize=(iopts.w, iopts.h))
    plt.tight_layout(pad=2)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo.txt")
    tredparse_results = parse_results("tredparse_results_homo.txt")
    title = SIMULATED_HAPLOID + r" ($D=%s\times, L=%dbp, V=%dbp$)" % (
        depth,
        readlen,
        distance,
    )
    plot_compare(ax1, title, tredparse_results, lobstr_results, max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model
    lobstr_results = parse_results("lobstr_results_het.txt", exclude=20)
    tredparse_results = parse_results("tredparse_results_het.txt", exclude=20)
    title = SIMULATED_DIPLOID + r" ($D=%s\times, L=%dbp, V=%dbp$)" % (
        depth,
        readlen,
        distance,
    )
    plot_compare(ax2, title, tredparse_results, lobstr_results, max_insert=max_insert)

    for ax in (ax1, ax2):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.03
    panel_labels(root, ((pad / 2, 1 - pad, "A"), (1 / 2.0, 1 - pad, "B")))
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def power(args):
    """
    %prog power

    Compare performances of various variant callers on simulated STR datasets.
    This compares the power of various evidence types.
    """
    p = OptionParser(power.__doc__)
    p.add_option(
        "--maxinsert", default=300, type="int", help="Maximum number of repeats"
    )
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10", format="png")

    if len(args) != 0:
        sys.exit(not p.print_help())

    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, figsize=(iopts.w, iopts.h)
    )
    plt.tight_layout(pad=3)

    color = "lightslategray"
    # ax1: Spanning
    tredparse_results = parse_results("tredparse_results_het-spanning.txt")
    title = SIMULATED_DIPLOID + " (Sub-model 1: Spanning reads)"
    plot_compare(
        ax1,
        title,
        tredparse_results,
        None,
        color=color,
        max_insert=max_insert,
        risk=False,
    )

    # ax2: Partial
    tredparse_results = parse_results("tredparse_results_het-partial.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 2: Partial reads)"
    plot_compare(
        ax2,
        title,
        tredparse_results,
        None,
        color=color,
        max_insert=max_insert,
        risk=False,
    )

    # ax3: Repeat
    tredparse_results = parse_results("tredparse_results_het-repeat.txt", exclude=20)
    # HACK (repeat reads won't work under 50)
    tredparse_results = [x for x in tredparse_results if x[0] > 50]
    title = SIMULATED_DIPLOID + " (Sub-model 3: Repeat-only reads)"
    plot_compare(
        ax3,
        title,
        tredparse_results,
        None,
        color=color,
        max_insert=max_insert,
        risk=False,
    )

    # ax4: Pair
    tredparse_results = parse_results("tredparse_results_het-pair.txt", exclude=20)
    title = SIMULATED_DIPLOID + " (Sub-model 4: Paired-end reads)"
    plot_compare(
        ax4,
        title,
        tredparse_results,
        None,
        color=color,
        max_insert=max_insert,
        risk=False,
    )

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.03
    panel_labels(
        root,
        (
            (pad / 2, 1 - pad, "A"),
            (1 / 2.0, 1 - pad, "B"),
            (pad / 2, 1 / 2.0, "C"),
            (1 / 2.0, 1 / 2.0, "D"),
        ),
    )
    normalize_axes(root)

    image_name = "power." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def tredparse(args):
    """
    %prog tredparse

    Compare performances of various variant callers on simulated STR datasets.
    Adds coverage comparisons as panel C and D.
    """
    p = OptionParser(tredparse.__doc__)
    p.add_option(
        "--maxinsert", default=300, type="int", help="Maximum number of repeats"
    )
    add_simulate_options(p)
    opts, args, iopts = p.set_image_options(args, figsize="10x10")

    if len(args) != 0:
        sys.exit(not p.print_help())

    depth = opts.depth
    max_insert = opts.maxinsert
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
        ncols=2, nrows=2, figsize=(iopts.w, iopts.h)
    )
    plt.tight_layout(pad=3)

    # ax1: lobSTR vs TREDPARSE with haploid model
    lobstr_results = parse_results("lobstr_results_homo-20x-150bp-500bp.txt")
    tredparse_results = parse_results("tredparse_results_homo-20x-150bp-500bp.txt")
    title = SIMULATED_HAPLOID + r" (Depth=$%s\times$)" % depth
    plot_compare(ax1, title, tredparse_results, lobstr_results, max_insert=max_insert)

    # ax2: lobSTR vs TREDPARSE with diploid model (depth=20x)
    lobstr_results = parse_results("lobstr_results_het-20x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results(
        "tredparse_results_het-20x-150bp-500bp.txt", exclude=20
    )
    title = SIMULATED_DIPLOID + r" (Depth=$%s\times$)" % depth
    plot_compare(ax2, title, tredparse_results, lobstr_results, max_insert=max_insert)

    # ax3: lobSTR vs TREDPARSE with diploid model (depth=5x)
    lobstr_results = parse_results("lobstr_results_het-5x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results(
        "tredparse_results_het-5x-150bp-500bp.txt", exclude=20
    )
    title = SIMULATED_DIPLOID + r" (Depth=$%s\times$)" % 5
    plot_compare(ax3, title, tredparse_results, lobstr_results, max_insert=max_insert)

    # ax4: lobSTR vs TREDPARSE with diploid model (depth=80x)
    lobstr_results = parse_results("lobstr_results_het-80x-150bp-500bp.txt", exclude=20)
    tredparse_results = parse_results(
        "tredparse_results_het-80x-150bp-500bp.txt", exclude=20
    )
    title = SIMULATED_DIPLOID + r" (Depth=$%s\times$)" % 80
    plot_compare(ax4, title, tredparse_results, lobstr_results, max_insert=max_insert)

    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(0, max_insert)
        ax.set_ylim(0, max_insert)

    root = fig.add_axes([0, 0, 1, 1])
    pad = 0.03
    panel_labels(
        root,
        (
            (pad / 2, 1 - pad, "A"),
            (1 / 2.0, 1 - pad, "B"),
            (pad / 2, 1 / 2.0, "C"),
            (1 / 2.0, 1 / 2.0, "D"),
        ),
    )
    normalize_axes(root)

    image_name = "tredparse." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
