#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Scripts for the pineapple genome manuscript (unpublished).
"""

import sys
import logging

from jcvi.formats.base import DictFile, must_open
from jcvi.formats.bed import Bed
from jcvi.formats.sizes import Sizes
from jcvi.graphics.base import plt, savefig, panel_labels
from jcvi.graphics.karyotype import Karyotype
from jcvi.graphics.synteny import Synteny, draw_gene_legend
from jcvi.graphics.glyph import TextCircle
from jcvi.apps.base import OptionParser, ActionDispatcher


def main():

    actions = (
        ('ploidy', 'plot pineapple macro-synteny (requires data)'),
        # build pseudomolecule
        ('agp', 'make agp file'),
        ('breakpoints', 'make breakpoints'),
        ('check', 'check agreement'),
        # build gene info table
        ('geneinfo', 'build gene info table'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def read_interpro(ipr):
    store = {}
    fp = open(ipr)
    # Aco000343.1     0d98a55eb3399a408e06252a2e24efcf        2083    Pfam
    # PF00476 DNA polymerase family A 1685    2075    1.70E-55        T
    # 10-10-2014      IPR001098       "DNA-directed DNA polymerase, family A,
    # palm domain"    GO:0003677|GO:0003887|GO:0006260        KEGG:
    # 00230+2.7.7.7|KEGG: 00240+2.7.7.7
    for row in fp:
        accession, md5, seqlen, analysis, signature, signature_description, \
        start, stop, score, status, date, interpro, interpro_description, GO, \
                        pathway = row.split("\t")
        accession = accession.split(".")[0]
        interpro_description = interpro_description.replace('"', "")
        pathway = pathway.strip()
        if accession not in ipr:
            store[accession] = (interpro, interpro_description, GO, pathway)
    return store


def geneinfo(args):
    """
    %prog geneinfo pineapple.20141007.bed pineapple.20150306.bed \
                   note.txt interproscan.txt

    Build gene info table from various sources.
    """
    p = OptionParser(geneinfo.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 4:
        sys.exit(not p.print_help())

    scfbed, lgbed, note, ipr = args
    note = DictFile(note, delimiter="\t")
    scfbed = Bed(scfbed)
    lgorder = Bed(lgbed).order
    header = "Accession Scaffold-position LG-position "\
             "Description Interpro-domain Interpro-description "\
             "GO-term KEGG".split()
    ipr = read_interpro(ipr)

    fw = must_open(opts.outfile, "w")
    print >> fw, "\t".join(header)

    for b in scfbed:
        accession = b.accn
        scaffold_position = b.tag
        if accession in lgorder:
            lg_position = lgorder[accession][-1].tag
        else:
            lg_position = "split"
        description = note[accession]
        interpro = interpro_description = go = kegg = ""
        if accession in ipr:
            interpro, interpro_description, go, kegg = ipr[accession]
        print >> fw, "\t".join((accession, scaffold_position, lg_position,
                        description, interpro, interpro_description, go, kegg))
    fw.close()


def ploidy(args):
    """
    %prog cotton seqids karyotype.layout mcscan.out all.bed synteny.layout

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
    draw_gene_legend(root, .27, .37, .52)

    # annotate the WGD events
    fc = 'lightslategrey'
    x = .09
    radius = .012
    TextCircle(root, x, .825, r'$\tau$', radius=radius, fc=fc)
    TextCircle(root, x, .8, r'$\sigma$', radius=radius, fc=fc)
    TextCircle(root, x, .72, r'$\rho$', radius=radius, fc=fc)
    for ypos in (.825, .8, .72):
        root.text(.12, ypos, r"$\times2$", color=fc, ha="center", va="center")
    root.plot([x, x], [.85, .775], ":", color=fc, lw=2)
    root.plot([x, x], [.75, .675], ":", color=fc, lw=2)

    labels = ((.04, .96, 'A'), (.04, .54, 'B'))
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
    fp.next()
    d = {}
    for row in fp:
        atoms = row.split()
        scaf, tag, linkage, no = atoms[:4]
        d[scaf] = tag

    fp = open("chimeric-scaffolds.txt")
    fp.next()
    for row in fp:
        old, new, tag, start, end = row.strip().split("\t")
        if new not in d:
            print new, "not in sheet1"
            continue
        if d[new] != tag:
            print "{0} => {1} in sheet1 but {2} in sheet2".format(new, d[new], tag)


def agp(args):
    fp = open("assembly-order.txt")
    fp.next()
    sizes = Sizes("SCAFFOLD-SPLIT.fasta").mapping
    for row in fp:
        atoms = row.split()
        assert len(atoms) in (4, 5)
        if len(atoms) == 4:
            atoms.append('?')
        scaf, tag, linkage, no, strand = atoms
        strand = strand.lower()
        strand = {'f': '+', 'r': '-', '?': '?'}[strand]
        scaf = "scaffold_" + scaf
        scaf_size = sizes[scaf]
        linkage = "LG{0:02d}".format(ord(linkage.lower()) - ord('a') + 1)
        print "\t".join(str(x) for x in \
                    (scaf, 0, scaf_size, linkage, 1000, strand))


def breakpoints(args):
    fp = open("chimeric-scaffolds.txt")
    fp.next()
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
        print "\t".join(str(x) for x in (old, start - 1, end))
        nbreaks += 1
        scaffolds.add(old)
    print >> sys.stderr, "{0} breakpoints in total, {1} scaffolds broken"\
                    .format(nbreaks, len(scaffolds))


if __name__ == '__main__':
    main()
