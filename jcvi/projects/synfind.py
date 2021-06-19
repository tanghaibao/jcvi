#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SynFind analyses and visualization.
"""
import os.path as op
import sys
import logging

from copy import deepcopy
from collections import defaultdict
from itertools import groupby

from jcvi.formats.base import get_number, must_open
from jcvi.utils.cbook import SummaryStats, gene_name, percentage
from jcvi.utils.grouper import Grouper
from jcvi.formats.blast import BlastLine
from jcvi.formats.bed import Bed
from jcvi.formats.gff import Gff, load
from jcvi.apps.grid import MakeManager
from jcvi.graphics.base import (
    FancyArrow,
    plt,
    savefig,
    panel_labels,
    markup,
    normalize_axes,
    latex,
)
from jcvi.graphics.glyph import CartoonRegion, RoundRect
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, symlink


def main():

    actions = (
        ("cartoon", "generate cartoon illustration of SynFind"),
        ("ecoli", "gene presence absence analysis in ecoli"),
        ("grass", "validate SynFind pan-grass set against James"),
        ("coge", "prepare coge datasets"),
        # For benchmarking
        ("synfind", "prepare input for SynFind"),
        ("iadhore", "prepare input for iADHoRe"),
        ("mcscanx", "prepare input for MCScanX"),
        ("cyntenator", "prepare input for Cyntenator"),
        ("athalianatruth", "prepare truth pairs for At alpha/beta/gamma"),
        ("yeasttruth", "prepare truth pairs for 14 yeasts"),
        ("grasstruth", "prepare truth pairs for 4 grasses"),
        ("benchmark", "compare SynFind, MCScanX, iADHoRe and OrthoFinder"),
        ("venn", "display benchmark results as Venn diagram"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def grasstruth(args):
    """
    %prog grasstruth james-pan-grass.txt

    Prepare truth pairs for 4 grasses.
    """
    p = OptionParser(grasstruth.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (james,) = args
    fp = open(james)
    pairs = set()
    for row in fp:
        atoms = row.split()
        genes = []
        idx = {}
        for i, a in enumerate(atoms):
            aa = a.split("||")
            for ma in aa:
                idx[ma] = i
            genes.extend(aa)
        genes = [x for x in genes if ":" not in x]
        Os = [x for x in genes if x.startswith("Os")]
        for o in Os:
            for g in genes:
                if idx[o] == idx[g]:
                    continue
                pairs.add(tuple(sorted((o, g))))

    for a, b in sorted(pairs):
        print("\t".join((a, b)))


def synfind(args):
    """
    %prog synfind all.last *.bed

    Prepare input for SynFind.
    """
    p = OptionParser(synfind.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    lastfile = args[0]
    bedfiles = args[1:]
    fp = open(lastfile)
    filteredlast = lastfile + ".filtered"
    fw = open(filteredlast, "w")
    for row in fp:
        b = BlastLine(row)
        if b.query == b.subject:
            continue
        print(b, file=fw)
    fw.close()
    logging.debug("Filtered LAST file written to `{0}`".format(filteredlast))

    allbed = "all.bed"
    fw = open(allbed, "w")
    for i, bedfile in enumerate(bedfiles):
        prefix = chr(ord("A") + i)
        bed = Bed(bedfile)
        for b in bed:
            b.seqid = prefix + b.seqid
            print(b, file=fw)
    fw.close()
    logging.debug("Bed file written to `{0}`".format(allbed))


def yeasttruth(args):
    """
    %prog yeasttruth Pillars.tab *.gff

    Prepare pairs data for 14 yeasts.
    """
    p = OptionParser(yeasttruth.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    pillars = args[0]
    gffiles = args[1:]
    aliases = {}
    pivot = {}
    for gffile in gffiles:
        is_pivot = op.basename(gffile).startswith("Saccharomyces_cerevisiae")
        gff = Gff(gffile)
        for g in gff:
            if g.type != "gene":
                continue
            for a in g.attributes["Alias"]:
                aliases[a] = g.accn
                if is_pivot:
                    pivot[a] = g.accn
    logging.debug("Aliases imported: {0}".format(len(aliases)))
    logging.debug("Pivot imported: {0}".format(len(pivot)))
    fw = open("yeast.aliases", "w")
    for k, v in sorted(aliases.items()):
        print("\t".join((k, v)), file=fw)
    fw.close()

    fp = open(pillars)
    pairs = set()
    fw = must_open(opts.outfile, "w")
    for row in fp:
        atoms = [x for x in row.split() if x != "---"]
        pps = [pivot[x] for x in atoms if x in pivot]
        atoms = [aliases[x] for x in atoms if x in aliases]
        for p in pps:
            for a in atoms:
                if p == a:
                    continue
                pairs.add(tuple(sorted((p, a))))

    for a, b in sorted(pairs):
        print("\t".join((a, b)), file=fw)
    fw.close()


def venn(args):
    """
    %prog venn *.benchmark

    Display benchmark results as Venn diagram.
    """
    from matplotlib_venn import venn2

    p = OptionParser(venn.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="9x9")

    if len(args) < 1:
        sys.exit(not p.print_help())

    bcs = args
    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    pad = 0.02
    ystart = 1
    ywidth = 1.0 / len(bcs)
    tags = ("Bowers", "YGOB", "Schnable")
    for bc, tag in zip(bcs, tags):
        fp = open(bc)
        data = []
        for row in fp:
            prog, pcounts, tcounts, shared = row.split()
            pcounts = int(pcounts)
            tcounts = int(tcounts)
            shared = int(shared)
            data.append((prog, pcounts, tcounts, shared))
        xstart = 0
        xwidth = 1.0 / len(data)
        for prog, pcounts, tcounts, shared in data:
            a, b, c = pcounts - shared, tcounts - shared, shared
            ax = fig.add_axes(
                [
                    xstart + pad,
                    ystart - ywidth + pad,
                    xwidth - 2 * pad,
                    ywidth - 2 * pad,
                ]
            )
            venn2(subsets=(a, b, c), set_labels=(prog, tag), ax=ax)
            message = "Sn={0} Pu={1}".format(
                percentage(shared, tcounts, precision=0, mode=-1),
                percentage(shared, pcounts, precision=0, mode=-1),
            )
            print(message, file=sys.stderr)
            ax.text(
                0.5,
                0.92,
                latex(message),
                ha="center",
                va="center",
                transform=ax.transAxes,
                color="b",
            )
            ax.set_axis_off()
            xstart += xwidth
        ystart -= ywidth

    panel_labels(
        root,
        ((0.04, 0.96, "A"), (0.04, 0.96 - ywidth, "B"), (0.04, 0.96 - 2 * ywidth, "C")),
    )
    panel_labels(
        root,
        (
            (0.5, 0.98, "A. thaliana duplicates"),
            (0.5, 0.98 - ywidth, "14 Yeast genomes"),
            (0.5, 0.98 - 2 * ywidth, "4 Grass genomes"),
        ),
    )
    normalize_axes(root)
    savefig("venn.pdf", dpi=opts.dpi)


def coge(args):
    """
    %prog coge *.gff

    Prepare coge datasets.
    """
    p = OptionParser(coge.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    gffs = args
    for gff in gffs:
        atoms = op.basename(gff).split(".")
        gid = atoms[-2]
        assert gid.startswith("gid")
        gid = get_number(gid)
        genomefasta = "genome_{0}.faa.fasta".format(gid)
        species = "_".join(atoms[0].split("_")[:2])
        cdsfasta = species + ".cds.fasta"
        load(
            [
                gff,
                genomefasta,
                "--id_attribute=Parent",
                "--outfile={0}".format(cdsfasta),
            ]
        )


def calc_sensitivity_specificity(a, truth, tag, fw):
    common = a & truth
    sensitivity = len(common) * 100.0 / len(truth)
    specificity = len(common) * 100.0 / len(a)
    logging.debug("{0}: {1} pairs".format(tag, len(a)))
    logging.debug(
        "{0}: Sensitivity={1:.1f}% Purity={2:.1f}%".format(
            tag, sensitivity, specificity
        )
    )
    print(tag, len(a), len(truth), len(common), file=fw)


def write_pairs(pairs, pairsfile):
    fz = open(pairsfile, "w")
    for a, b in pairs:
        print("\t".join((a, b)), file=fz)
    fz.close()


def benchmark(args):
    """
    %prog benchmark at bedfile

    Compare SynFind, MCScanx, iADHoRe and OrthoFinder against the truth.
    """
    p = OptionParser(benchmark.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    pf, bedfile = args
    truth = pf + ".truth"
    synfind = pf + ".synfind"
    mcscanx = pf + ".mcscanx"
    iadhore = pf + ".iadhore"
    orthofinder = pf + ".orthofinder"
    pivots = set([x.accn for x in Bed(bedfile)])

    fp = open(truth)
    truth = set()
    for row in fp:
        a, b = row.strip().split("\t")[:2]
        pivots.add(a)
        truth.add(tuple(sorted((a, b))))
    logging.debug("Truth: {0} pairs".format(len(truth)))

    fp = open(synfind)
    benchmarkfile = pf + ".benchmark"
    fw = must_open(benchmarkfile, "w")
    synfind = set()
    for row in fp:
        atoms = row.strip().split("\t")
        query, hit, tag = atoms[:3]
        if tag != "S":
            continue
        synfind.add(tuple(sorted((query, hit))))
    calc_sensitivity_specificity(synfind, truth, "SynFind", fw)

    fp = open(mcscanx)
    mcscanx = set()
    for row in fp:
        if row[0] == "#":
            continue
        atoms = row.strip().split(":")[1].split()
        query, hit = atoms[:2]
        mcscanx.add(tuple(sorted((query, hit))))
    calc_sensitivity_specificity(mcscanx, truth, "MCScanX", fw)

    fp = open(iadhore)
    iadhore = set()
    next(fp)
    for row in fp:
        atoms = row.strip().split("\t")
        query, hit = atoms[3:5]
        iadhore.add(tuple(sorted((query, hit))))
    calc_sensitivity_specificity(iadhore, truth, "iADHoRe", fw)

    fp = open(orthofinder)
    orthofinder = set()
    next(fp)
    for row in fp:
        row = row.replace('"', "")
        atoms = row.replace(",", " ").split()
        genes = [x.strip() for x in atoms if not x.startswith("OG")]
        genes = [gene_name(x) for x in genes]
        pps = [x for x in genes if x in pivots]
        for p in pps:
            for g in genes:
                if p == g:
                    continue
                orthofinder.add(tuple(sorted((p, g))))
    # write_pairs(orthofinder, "orthofinder.pairs")
    calc_sensitivity_specificity(orthofinder, truth, "OrthoFinder", fw)
    fw.close()


def write_lst(bedfile):
    pf = op.basename(bedfile).split(".")[0]
    mkdir(pf)
    bed = Bed(bedfile)
    stanza = []
    for seqid, bs in bed.sub_beds():
        fname = op.join(pf, "{0}.lst".format(seqid))
        fw = open(fname, "w")
        for b in bs:
            print("{0}{1}".format(b.accn.replace(" ", ""), b.strand), file=fw)
        stanza.append((seqid, fname))
        fw.close()
    return pf, stanza


def write_txt(bedfile):
    pf = op.basename(bedfile).split(".")[0][:20]
    txtfile = pf + ".txt"
    fw = open(txtfile, "w")
    print("#genome", file=fw)
    bed = Bed(bedfile)
    for b in bed:
        print(
            " ".join(str(x) for x in (b.accn, b.seqid, b.start, b.end, b.strand)),
            file=fw,
        )
    fw.close()
    return txtfile


def cyntenator(args):
    """
    %prog cyntenator athaliana.athaliana.last athaliana.bed

    Prepare input for Cyntenator.
    """
    p = OptionParser(cyntenator.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    lastfile = args[0]
    fp = open(lastfile)
    filteredlastfile = lastfile + ".blast"
    fw = open(filteredlastfile, "w")
    for row in fp:
        b = BlastLine(row)
        if b.query == b.subject:
            continue
        print("\t".join((b.query, b.subject, str(b.score))), file=fw)
    fw.close()

    bedfiles = args[1:]
    fp = open(lastfile)
    b = BlastLine(next(fp))
    subject = b.subject
    txtfiles = []
    for bedfile in bedfiles:
        order = Bed(bedfile).order
        if subject in order:
            db = op.basename(bedfile).split(".")[0][:20]
            logging.debug("Found db: {0}".format(db))
        txtfile = write_txt(bedfile)
        txtfiles.append(txtfile)

    db += ".txt"
    mm = MakeManager()
    for txtfile in txtfiles:
        outfile = txtfile + ".alignment"
        cmd = 'cyntenator -t "({0} {1})" -h blast {2} > {3}'.format(
            txtfile, db, filteredlastfile, outfile
        )
        mm.add((txtfile, db, filteredlastfile), outfile, cmd)
    mm.write()


def iadhore(args):
    """
    %prog iadhore athaliana.athaliana.last athaliana.bed

    Wrap around iADHoRe.
    """
    p = OptionParser(iadhore.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    lastfile = args[0]
    bedfiles = args[1:]
    blast_table = "blast_table.txt"
    fp = open(lastfile)
    seen = set()
    for row in fp:
        c = BlastLine(row)
        a, b = c.query, c.subject
        a, b = gene_name(a), gene_name(b)
        if a > b:
            a, b = b, a
        seen.add((a, b))

    fw = open(blast_table, "w")
    for a, b in seen:
        print("\t".join((a, b)), file=fw)
    fw.close()
    logging.debug(
        "A total of {0} pairs written to `{1}`".format(len(seen), blast_table)
    )

    fw = open("config.txt", "w")
    for bedfile in bedfiles:
        pf, stanza = write_lst(bedfile)
        print("genome={0}".format(pf), file=fw)
        for seqid, fname in stanza:
            print(" ".join((seqid, fname)), file=fw)
        print(file=fw)

    print("blast_table={0}".format(blast_table), file=fw)
    print("cluster_type=colinear", file=fw)
    print("tandem_gap=10", file=fw)
    print("prob_cutoff=0.001", file=fw)
    print("gap_size=20", file=fw)
    print("cluster_gap=20", file=fw)
    print("q_value=0.9", file=fw)
    print("anchor_points=4", file=fw)
    print("alignment_method=gg2", file=fw)
    print("max_gaps_in_alignment=20", file=fw)
    print("output_path=i-adhore_out", file=fw)
    print("number_of_threads=4", file=fw)
    fw.close()


def extract_groups(g, pairs, txtfile):
    register = defaultdict(list)
    fp = open(txtfile)
    next(fp)
    for row in fp:
        if row[0] != ">":
            continue
        track, atg, myname, pairname = row.split()
        pairname = pairname.rstrip("ab").upper()
        register[pairname].append(atg.upper())

    for pairname, genes in register.items():
        tag = pairname[0]
        tag = {"A": "alpha", "B": "beta", "C": "gamma", "S": "others"}[tag]
        pairs.add(tuple(sorted(genes) + [tag]))
        g.join(*genes)


def athalianatruth(args):
    """
    %prog athalianatruth J_a.txt J_bc.txt

    Prepare pairs data for At alpha/beta/gamma.
    """
    p = OptionParser(athalianatruth.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    atxt, bctxt = args
    g = Grouper()
    pairs = set()
    for txt in (atxt, bctxt):
        extract_groups(g, pairs, txt)

    fw = open("pairs", "w")
    for pair in sorted(pairs):
        print("\t".join(pair), file=fw)
    fw.close()

    fw = open("groups", "w")
    for group in list(g):
        print(",".join(group), file=fw)
    fw.close()


def make_gff(bed, prefix, fw):
    bed = Bed(bed)
    nfeats = 0
    for b in bed:
        seqid = prefix + b.seqid
        print("\t".join(str(x) for x in (seqid, b.accn, b.start, b.end)), file=fw)
        nfeats += 1
    logging.debug("A total of {0} features converted to `{1}`".format(nfeats, fw.name))


def mcscanx(args):
    """
    %prog mcscanx athaliana.athaliana.last athaliana.bed

    Wrap around MCScanX.
    """
    p = OptionParser(mcscanx.__doc__)
    opts, args = p.parse_args(args)

    if len(args) < 2:
        sys.exit(not p.print_help())

    blastfile = args[0]
    bedfiles = args[1:]
    prefix = "_".join(op.basename(x)[:2] for x in bedfiles)
    symlink(blastfile, prefix + ".blast")
    allbedfile = prefix + ".gff"
    fw = open(allbedfile, "w")
    for i, bedfile in enumerate(bedfiles):
        prefix = chr(ord("A") + i)
        make_gff(bedfile, prefix, fw)
    fw.close()


def grass(args):
    """
    %prog grass coge_master_table.txt james.txt

    Validate SynFind pan-grass set against James. This set can be generated:

    https://genomevolution.org/r/fhak
    """
    p = OptionParser(grass.__doc__)
    p.set_verbose()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    master, james = args

    fp = open(master)
    next(fp)
    master_store = defaultdict(set)
    for row in fp:
        atoms = row.split()
        s = set()
        for x in atoms[1:6]:
            m = x.split(",")
            s |= set(m)
        if "-" in s:
            s.remove("-")

        a = atoms[1]
        master_store[a] |= set(s)

    fp = open(james)
    next(fp)
    james_store = {}
    tandems = set()
    for row in fp:
        atoms = row.split()
        s = set()
        Os = set()
        for x in atoms[:-1]:
            m = x.split("||")
            if m[0].startswith("Os"):
                Os |= set(m)
            if m[0].startswith("http"):
                continue
            if m[0].startswith("chr"):
                m = ["proxy"]
            if "||" in x:
                tandems |= set(m)
            s |= set(m)

        for x in Os:
            james_store[x] = s

    jaccards = []
    corr_jaccards = []
    perfect_matches = 0
    corr_perfect_matches = 0
    for k, v in james_store.items():
        if k not in master_store:
            continue
        m = master_store[k]
        jaccard = len(v & m) * 100 / len(v | m)
        jaccards.append(jaccard)
        diff = (v ^ m) - tandems
        corr_jaccard = 100 - len(diff) * 100 / len(v | m)
        corr_jaccards.append(corr_jaccard)
        if opts.verbose:
            print(k)
            print(v)
            print(m)
            print(diff)
            print(jaccard)
        if jaccard > 99:
            perfect_matches += 1
        if corr_jaccard > 99:
            corr_perfect_matches += 1

    logging.debug("Perfect matches: {0}".format(perfect_matches))
    logging.debug("Perfect matches (corrected): {0}".format(corr_perfect_matches))
    print("Jaccards:", SummaryStats(jaccards))
    print("Corrected Jaccards:", SummaryStats(corr_jaccards))


def ecoli(args):
    """
    %prog ecoli coge_master_table.txt query.bed

    Perform gene presence / absence analysis in Ecoli master spreadsheet. Ecoli
    spresheets can be downloaded below:

    Ecoli K12 MG1655 (K) as query
    Regenerate this analysis: https://genomevolution.org/r/fggo

    Ecoli O157:H7 EDL933 (O) as query
    Regenerate this analysis: https://genomevolution.org/r/fgt7

    Shigella flexneri 2a 301 (S) as query
    Regenerate this analysis: https://genomevolution.org/r/fgte

    Perform a similar analysis as in:
    Jin et al. (2002) Genome sequence of Shigella flexneri 2a: insights
    into pathogenicity through comparison with genomes of Escherichia
    coli K12 and O157. Nucleic Acid Research.
    """
    p = OptionParser(ecoli.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    master, querybed = args
    fp = open(master)
    header = next(fp)
    assert header[0] == "#"
    qorg = header.strip().split("\t")[1]
    qorg = qorg.split(":")[-1].strip()

    store = {}
    MISSING = ("proxy", "-")
    for row in fp:
        a, b, c = row.strip().split("\t")[1:4]
        store[a] = b in MISSING and c in MISSING

    bed = Bed(querybed)
    tags = []
    for i, b in enumerate(bed):
        accn = b.accn
        if accn not in store:
            logging.warning("missing {0}".format(accn))
            continue
        tags.append((store[accn], accn))

    large = 4  # large segments
    II = []
    II_large = []
    for missing, aa in groupby(tags, key=lambda x: x[0]):
        aa = list(aa)
        if not missing:
            continue
        glist = list(a for missing, a in aa)
        II.append(glist)
        size = len(glist)
        if size >= large:
            II_large.append(glist)

    fw = must_open(opts.outfile, "w")
    for a, t in zip((II, II_large), ("", ">=4 ")):
        nmissing = sum(len(x) for x in a)
        logging.debug(
            "A total of {0} {1}-specific {2}islands found with {3} genes.".format(
                len(a), qorg, t, nmissing
            )
        )

    for x in II:
        print(len(x), ",".join(x), file=fw)


def plot_diagram(ax, x, y, A, B, tag, label):
    ax.text(x, y + 0.14, "{0}: {1}".format(tag, label), ha="center")
    strip = tag != "G"
    A.draw(ax, x, y + 0.06, gene_len=0.02, strip=strip)
    B.draw(ax, x, y, gene_len=0.02, strip=strip)


def cartoon(args):
    """
    %prog synteny.py

    Generate cartoon illustration of SynFind.
    """
    p = OptionParser(cartoon.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="10x7")

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    # Panel A
    A = CartoonRegion(41)
    A.draw(root, 0.35, 0.85, strip=False, color=False)
    x1, x2 = A.x1, A.x2
    lsg = "lightslategray"
    pad = 0.01
    xc, yc = 0.35, 0.88
    arrowlen = x2 - xc - pad
    arrowprops = dict(
        length_includes_head=True,
        width=0.01,
        fc=lsg,
        lw=0,
        head_length=arrowlen * 0.15,
        head_width=0.03,
    )
    p = FancyArrow(xc - pad, yc, -arrowlen, 0, shape="left", **arrowprops)
    root.add_patch(p)
    p = FancyArrow(xc + pad, yc, arrowlen, 0, shape="right", **arrowprops)
    root.add_patch(p)

    yt = yc + 4 * pad
    root.text((x1 + xc) / 2, yt, "20 genes upstream", ha="center")
    root.text((x2 + xc) / 2, yt, "20 genes downstream", ha="center")
    root.plot((xc,), (yc,), "o", mfc="w", mec=lsg, mew=2, lw=2, color=lsg)
    root.text(xc, yt, "Query gene", ha="center")

    # Panel B
    A.draw(root, 0.35, 0.7, strip=False)

    RoundRect(root, (0.07, 0.49), 0.56, 0.14, fc="y", alpha=0.2)
    a = deepcopy(A)
    a.evolve(mode="S", target=10)
    a.draw(root, 0.35, 0.6)
    b = deepcopy(A)
    b.evolve(mode="F", target=8)
    b.draw(root, 0.35, 0.56)
    c = deepcopy(A)
    c.evolve(mode="G", target=6)
    c.draw(root, 0.35, 0.52)

    for x in (a, b, c):
        root.text(0.64, x.y, "Score={0}".format(x.nonwhites), va="center")

    # Panel C
    A.truncate_between_flankers()
    a.truncate_between_flankers()
    b.truncate_between_flankers()
    c.truncate_between_flankers(target=6)

    plot_diagram(root, 0.14, 0.2, A, a, "S", "syntenic")
    plot_diagram(root, 0.37, 0.2, A, b, "F", "missing, with both flankers")
    plot_diagram(root, 0.6, 0.2, A, c, "G", "missing, with one flanker")

    labels = ((0.04, 0.95, "A"), (0.04, 0.75, "B"), (0.04, 0.4, "C"))
    panel_labels(root, labels)

    # Descriptions
    xt = 0.85
    desc = (
        "Extract neighborhood",
        "of *window* size",
        "Count gene pairs within *window*",
        "Find regions above *score* cutoff",
        "Identify flankers",
        "Annotate syntelog class",
    )
    for yt, t in zip((0.88, 0.84, 0.64, 0.6, 0.3, 0.26), desc):
        root.text(xt, yt, markup(t), ha="center", va="center")

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "cartoon"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


if __name__ == "__main__":
    main()
