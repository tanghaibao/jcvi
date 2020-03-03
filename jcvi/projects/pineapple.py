#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the pineapple genome paper.
"""
from __future__ import print_function

import sys
import logging

from jcvi.formats.base import DictFile, LineFile, SetFile, must_open, get_number
from jcvi.formats.bed import Bed
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import Rectangle, panel_labels, plt, savefig
from jcvi.graphics.chromosome import Chromosome
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
from jcvi.graphics.glyph import TextCircle
from jcvi.annotation.ahrd import read_interpro
from jcvi.apps.base import OptionParser, ActionDispatcher


class RegionsLine(object):
    def __init__(self, line):
        args = line.split()
        self.karyotype = args[0][0]
        self.group = args[0][1]
        self.chromosome = int(args[1])
        self.start = int(args[5])
        self.end = int(args[8])
        self.span = abs(self.start - self.end)


class RegionsFile(LineFile):
    def __init__(self, filename):
        super(RegionsFile, self).__init__(filename)
        fp = open(filename)
        next(fp)
        for row in fp:
            self.append(RegionsLine(row))

    @property
    def karyotypes(self):
        return sorted(set(x.karyotype for x in self))

    def get_karyotype(self, k):
        return [x for x in self if x.karyotype == k]


def main():

    actions = (
        # main figures in text
        ("ancestral", "karoytype evolution of pineapple (requires data)"),
        ("ploidy", "plot pineapple macro-synteny (requires data)"),
        # build pseudomolecule
        ("agp", "make agp file"),
        ("breakpoints", "make breakpoints"),
        ("check", "check agreement"),
        # build gene info table
        ("geneinfo", "build gene info table"),
        ("flanking", "extract flanking genes for given SI loci"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def flanking(args):
    """
    %prog flanking SI.ids liftover.bed master.txt master-removed.txt

    Extract flanking genes for given SI loci.
    """
    p = OptionParser(flanking.__doc__)
    p.add_option("-N", default=50, type="int", help="How many genes on both directions")
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    SI, liftover, master, te = args
    N = opts.N
    SI = SetFile(SI, column=0, delimiter=".")
    liftover = Bed(liftover)
    order = liftover.order
    neighbors = set()
    for s in SI:
        si, s = order[s]
        LB = max(si - N, 0)
        RB = min(si + N, len(liftover))
        for j in range(LB, RB + 1):
            a = liftover[j]
            if a.seqid != s.seqid:
                continue
            neighbors.add(a.accn)

    dmain = DictFile(master, keypos=0, valuepos=None, delimiter="\t")
    dte = DictFile(te, keypos=0, valuepos=None, delimiter="\t")
    header = next(open(master))
    print("\t".join(("SI/Neighbor", "Gene/TE", header.strip())))
    for a in liftover:
        s = a.accn
        if s not in neighbors:
            continue

        tag = "SI" if s in SI else "neighbor"
        if s in dmain:
            d = dmain[s]
            print("\t".join([tag, "gene"] + d))
        elif s in dte:
            d = dte[s]
            print("\t".join([tag, "TE"] + d))


def join_nodes_vertical(root, coords, a, b, y, lw=2):
    # Join node a and b to make an internal node
    ax, ay = coords[a]
    bx, by = coords[b]
    nx, ny = (ax + bx) / 2, y
    root.plot((ax, ax), (ay, ny), "k-", lw=lw)
    root.plot((bx, bx), (ay, ny), "k-", lw=lw)
    root.plot((ax, bx), (ny, ny), "k-", lw=lw)
    return nx, ny


def ancestral(args):
    """
    %prog ancestral ancestral.txt assembly.fasta

    Karyotype evolution of pineapple. The figure is inspired by Amphioxus paper
    Figure 3 and Tetradon paper Figure 9.
    """
    p = OptionParser(ancestral.__doc__)
    opts, args, iopts = p.set_image_options(args, figsize="8x7")

    if len(args) != 2:
        sys.exit(not p.print_help())

    regionsfile, sizesfile = args
    regions = RegionsFile(regionsfile)
    sizes = Sizes(sizesfile).mapping
    sizes = dict((k, v) for (k, v) in sizes.iteritems() if k[:2] == "LG")
    maxsize = max(sizes.values())
    ratio = 0.5 / maxsize

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes((0, 0, 1, 1))

    from jcvi.graphics.base import set2

    a, b, c, d, e, f, g = set2[:7]
    set2 = (c, g, b, e, d, a, f)

    # Upper panel is the evolution of segments
    # All segments belong to one of seven karyotypes 1 to 7
    karyotypes = regions.karyotypes
    xgap = 1.0 / (1 + len(karyotypes))
    ygap = 0.05
    mgap = xgap / 4.5
    gwidth = mgap * 0.75
    tip = 0.02
    coords = {}
    for i, k in enumerate(regions.karyotypes):
        x = (i + 1) * xgap
        y = 0.9
        root.text(x, y + tip, "Anc" + k, ha="center")
        root.plot((x, x), (y, y - ygap), "k-", lw=2)
        y -= 2 * ygap
        coords["a"] = (x - 1.5 * mgap, y)
        coords["b"] = (x - 0.5 * mgap, y)
        coords["c"] = (x + 0.5 * mgap, y)
        coords["d"] = (x + 1.5 * mgap, y)
        coords["ab"] = join_nodes_vertical(root, coords, "a", "b", y + ygap / 2)
        coords["cd"] = join_nodes_vertical(root, coords, "c", "d", y + ygap / 2)
        coords["abcd"] = join_nodes_vertical(root, coords, "ab", "cd", y + ygap)
        for n in "abcd":
            nx, ny = coords[n]
            root.text(nx, ny - tip, n, ha="center")
            coords[n] = (nx, ny - ygap / 2)

        kdata = regions.get_karyotype(k)
        for kd in kdata:
            g = kd.group
            gx, gy = coords[g]
            gsize = ratio * kd.span
            gy -= gsize
            p = Rectangle((gx - gwidth / 2, gy), gwidth, gsize, lw=0, color=set2[i])
            root.add_patch(p)
            root.text(
                gx, gy + gsize / 2, kd.chromosome, ha="center", va="center", color="w"
            )
            coords[g] = (gx, gy - tip)

    # Bottom panel shows the location of segments on chromosomes
    # TODO: redundant code, similar to graphics.chromosome
    ystart = 0.54
    chr_number = len(sizes)
    xstart, xend = xgap - 2 * mgap, 1 - xgap + 2 * mgap
    xinterval = (xend - xstart - gwidth) / (chr_number - 1)
    chrpos = {}
    for a, (chr, clen) in enumerate(sorted(sizes.items())):
        chr = get_number(chr)
        xx = xstart + a * xinterval + gwidth / 2
        chrpos[chr] = xx
        root.text(xx, ystart + 0.01, chr, ha="center")
        Chromosome(root, xx, ystart, ystart - clen * ratio, width=gwidth)

    # Start painting
    for r in regions:
        xx = chrpos[r.chromosome]
        yystart = ystart - r.start * ratio
        yyend = ystart - r.end * ratio
        p = Rectangle(
            (xx - gwidth / 2, yystart),
            gwidth,
            yyend - yystart,
            color=set2[int(r.karyotype) - 1],
            lw=0,
        )
        root.add_patch(p)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "pineapple-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def geneinfo(args):
    """
    %prog geneinfo pineapple.20141004.bed liftover.bed pineapple.20150413.bed \
                   note.txt interproscan.txt

    Build gene info table from various sources. The three beds contain
    information on the original scaffolds, linkage groups, and final selected
    loci (after removal of TEs and split loci). The final two text files contain
    AHRD and domain data.
    """
    p = OptionParser(geneinfo.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 5:
        sys.exit(not p.print_help())

    scfbed, liftoverbed, lgbed, note, ipr = args
    note = DictFile(note, delimiter="\t")
    scfbed = Bed(scfbed)
    lgorder = Bed(lgbed).order
    liftover = Bed(liftoverbed).order
    header = (
        "Accession Scaffold-position LG-position "
        "Description Interpro-domain Interpro-description "
        "GO-term KEGG".split()
    )
    ipr = read_interpro(ipr)

    fw_clean = must_open("master.txt", "w")
    fw_removed = must_open("master-removed.txt", "w")

    for fw in (fw_clean, fw_removed):
        print("\t".join(header), file=fw)

    for b in scfbed:
        accession = b.accn
        scaffold_position = b.tag
        if accession in liftover:
            lg_position = liftover[accession][-1].tag
        else:
            lg_position = "split"
        fw = fw_clean if accession in lgorder else fw_removed
        description = note[accession]
        interpro = interpro_description = go = kegg = ""
        if accession in ipr:
            interpro, interpro_description, go, kegg = ipr[accession]
        print(
            "\t".join(
                (
                    accession,
                    scaffold_position,
                    lg_position,
                    description,
                    interpro,
                    interpro_description,
                    go,
                    kegg,
                )
            ),
            file=fw,
        )
    fw.close()


def ploidy(args):
    """
    %prog ploidy seqids karyotype.layout mcscan.out all.bed synteny.layout

    Build a figure that calls graphics.karyotype to illustrate the high ploidy
    of WGD history of pineapple genome. The script calls both graphics.karyotype
    and graphic.synteny.
    """
    p = OptionParser(ploidy.__doc__)
    p.add_option("--switch", help="Rename the seqid with two-column file")
    opts, args, iopts = p.set_image_options(args, figsize="9x7")

    if len(args) != 5:
        sys.exit(not p.print_help())

    seqidsfile, klayout, datafile, bedfile, slayout = args

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])

    Karyotype(fig, root, seqidsfile, klayout)
    Synteny(fig, root, datafile, bedfile, slayout, switch=opts.switch)

    # legend showing the orientation of the genes
    draw_gene_legend(root, 0.27, 0.37, 0.52)

    # annotate the WGD events
    fc = "lightslategrey"
    x = 0.09
    radius = 0.012
    TextCircle(root, x, 0.825, r"$\tau$", radius=radius, fc=fc)
    TextCircle(root, x, 0.8, r"$\sigma$", radius=radius, fc=fc)
    TextCircle(root, x, 0.72, r"$\rho$", radius=radius, fc=fc)
    for ypos in (0.825, 0.8, 0.72):
        root.text(0.12, ypos, r"$\times2$", color=fc, ha="center", va="center")
    root.plot([x, x], [0.85, 0.775], ":", color=fc, lw=2)
    root.plot([x, x], [0.75, 0.675], ":", color=fc, lw=2)

    labels = ((0.04, 0.96, "A"), (0.04, 0.54, "B"))
    panel_labels(root, labels)

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    pf = "pineapple-karyotype"
    image_name = pf + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


scaffold = "scaffold_"


def check(args):
    fp = open("assembly-order.txt")
    next(fp)
    d = {}
    for row in fp:
        atoms = row.split()
        scaf, tag, linkage, no = atoms[:4]
        d[scaf] = tag

    fp = open("chimeric-scaffolds.txt")
    next(fp)
    for row in fp:
        old, new, tag, start, end = row.strip().split("\t")
        if new not in d:
            print(new, "not in sheet1")
            continue
        if d[new] != tag:
            print("{0} => {1} in sheet1 but {2} in sheet2".format(new, d[new], tag))


def agp(args):
    fp = open("assembly-order.txt")
    next(fp)
    sizes = Sizes("SCAFFOLD-SPLIT.fasta").mapping
    for row in fp:
        atoms = row.split()
        assert len(atoms) in (4, 5)
        if len(atoms) == 4:
            atoms.append("?")
        scaf, tag, linkage, no, strand = atoms
        strand = strand.lower()
        strand = {"f": "+", "r": "-", "?": "?"}[strand]
        scaf = "scaffold_" + scaf
        scaf_size = sizes[scaf]
        linkage = "LG{0:02d}".format(ord(linkage.lower()) - ord("a") + 1)
        print("\t".join(str(x) for x in (scaf, 0, scaf_size, linkage, 1000, strand)))


def breakpoints(args):
    fp = open("chimeric-scaffolds.txt")
    next(fp)
    scaffolds = set()
    nbreaks = 0
    for row in fp:
        atoms = row.strip().split("\t")
        if len(atoms) == 3:
            continue
        old, new, tag, start, end = atoms
        old = scaffold + old
        start, end = int(start), int(end)
        if start >= end:
            logging.warn("{0} {1} >= {2}".format(old, start, end))
            start, end = end, start
        print("\t".join(str(x) for x in (old, start - 1, end)))
        nbreaks += 1
        scaffolds.add(old)
    print(
        "{0} breakpoints in total, {1} scaffolds broken".format(
            nbreaks, len(scaffolds)
        ),
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
