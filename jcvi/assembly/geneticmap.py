#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Use genetic map to break chimeric scaffolds, order and orient scaffolds onto
chromosomes.
"""
from __future__ import print_function

import os.path as op
import sys
import logging
import numpy as np

from collections import Counter
from functools import lru_cache
from itertools import combinations, groupby

from jcvi.formats.base import BaseFile, LineFile, must_open, read_block
from jcvi.formats.bed import Bed, fastaFromBed
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update


MSTheader = """population_type {0}
population_name LG
distance_function kosambi
cut_off_p_value 0.000001
no_map_dist 10.0
no_map_size 0
missing_threshold {1}
estimation_before_clustering no
detect_bad_data yes
objective_function ML
number_of_loci {2}
number_of_individual {3}
"""


class BinMap(BaseFile, dict):
    def __init__(self, filename):
        super(BinMap, self).__init__(filename)

        fp = open(filename)
        for header, seq in read_block(fp, "group "):
            lg = header.split()[-1]
            self[lg] = []
            for s in seq:
                if s.strip() == "" or s[0] == ";":
                    continue
                marker, pos = s.split()
                pos = int(float(pos) * 1000)
                self[lg].append((marker, pos))

    def print_to_bed(self, filename="stdout", switch=False, sep="."):
        """Print the genetic map in the BED format.

        Args:
            filename (str, optional): Output filename. Defaults to "stdout".
            switch (bool, optional): Use linkage group as seqid. Defaults to False.
            sep (str, optional): Separator that delimits scaffold name and position. Defaults to ".".
        """
        fw = must_open(filename, "w")
        for lg, markers in sorted(self.items()):
            for marker, pos in markers:
                if not switch:
                    line = (lg, pos, pos + 1, marker)
                else:
                    seqid_spos = marker.rsplit(sep, 1)
                    if len(seqid_spos) != 2:
                        logging.error(
                            "Error: `{}` must be in the form e.g. `name{}position`".format(
                                marker, sep
                            )
                        )
                        continue
                    seqid, spos = seqid_spos
                    spos = int(spos)
                    marker = "{0}:{1}".format(lg, pos / 1000.0)
                    line = (seqid, spos - 1, spos, marker)
                print("\t".join(str(x) for x in line), file=fw)
        fw.close()


class MSTMapLine(object):
    def __init__(self, row, startidx=3):
        args = row.split()
        self.id = args[0]
        self.seqid, pos = self.id.split(".")
        self.pos = int(pos)
        self.genotype = "".join(args[startidx:])

    def __len__(self):
        return len(self.genotype)

    def __str__(self):
        return "{0}: {1}".format(self.id, self.genotype)

    @property
    def bedline(self):
        return "\t".join(str(x) for x in (self.seqid, self.pos - 1, self.pos, self.id))


class MSTMap(LineFile):
    def __init__(self, filename):
        super(MSTMap, self).__init__(filename)
        fp = open(filename)
        startidx = 1
        for row in fp:
            if row.startswith("locus_name"):
                if row.split()[1] == "seqid":
                    startidx = 3
                self.header = row.split()
                break

        for row in fp:
            self.append(MSTMapLine(row, startidx=startidx))

        self.nmarkers = len(self)
        self.nind = len(self[0].genotype)
        logging.debug(
            "Map contains {0} markers in {1} individuals".format(
                self.nmarkers, self.nind
            )
        )


class MSTMatrix(object):
    def __init__(self, matrix, markerheader, population_type, missing_threshold):
        self.matrix = matrix
        self.markerheader = markerheader
        self.population_type = population_type
        self.missing_threshold = missing_threshold
        self.ngenotypes = len(matrix)
        self.nind = len(markerheader) - 1
        assert self.nind == len(matrix[0]) - 1
        logging.debug(
            "Imported {0} markers and {1} individuals.".format(
                self.ngenotypes, self.nind
            )
        )

    def write(self, filename="stdout", header=True):
        fw = must_open(filename, "w")
        if header:
            print(
                MSTheader.format(
                    self.population_type,
                    self.missing_threshold,
                    self.ngenotypes,
                    self.nind,
                ),
                file=fw,
            )
        print("\t".join(self.markerheader), file=fw)
        for m in self.matrix:
            print("\t".join(m), file=fw)


def main():

    actions = (
        ("breakpoint", "find scaffold breakpoints using genetic map"),
        ("ld", "calculate pairwise linkage disequilibrium"),
        ("bed", "convert MSTmap output to bed format"),
        ("fasta", "extract markers based on map"),
        ("anchor", "anchor scaffolds based on map"),
        ("rename", "rename markers according to the new mapping locations"),
        ("header", "rename lines in the map header"),
        # Plot genetic map
        ("blat", "make ALLMAPS input csv based on sequences"),
        ("dotplot", "make dotplot between chromosomes and linkage maps"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def blat(args):
    """
    %prog blat map1.txt ref.fasta

    Make ALLMAPS input csv based on sequences. The tab-delimited txt file
    include: name, LG, position, sequence.
    """
    from jcvi.formats.base import is_number
    from jcvi.formats.blast import best as blast_best, bed as blast_bed
    from jcvi.apps.align import blat as blat_align

    p = OptionParser(blat.__doc__)
    _, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    maptxt, ref = args
    pf = maptxt.rsplit(".", 1)[0]
    register = {}
    fastafile = pf + ".fasta"
    fp = open(maptxt)
    fw = open(fastafile, "w")
    for row in fp:
        name, lg, pos, seq = row.split()
        if not is_number(pos):
            continue
        register[name] = (pf + "-" + lg, pos)
        print(">{0}\n{1}\n".format(name, seq), file=fw)
    fw.close()

    blatfile = blat_align([ref, fastafile])
    bestfile = blast_best([blatfile])
    bedfile = blast_bed([bestfile])
    b = Bed(bedfile).order

    pf = ".".join((op.basename(maptxt).split(".")[0], op.basename(ref).split(".")[0]))
    csvfile = pf + ".csv"
    fp = open(maptxt)
    fw = open(csvfile, "w")
    for row in fp:
        name, lg, pos, seq = row.split()
        if name not in b:
            continue
        bbi, bb = b[name]
        scaffold, scaffold_pos = bb.seqid, bb.start
        print(",".join(str(x) for x in (scaffold, scaffold_pos, lg, pos)), file=fw)
    fw.close()


def dotplot(args):
    """
    %prog dotplot map.csv ref.fasta

    Make dotplot between chromosomes and linkage maps.
    The input map is csv formatted, for example:

    ScaffoldID,ScaffoldPosition,LinkageGroup,GeneticPosition
    scaffold_2707,11508,1,0
    scaffold_2707,11525,1,1.2
    """
    from natsort import natsorted
    from jcvi.assembly.allmaps import CSVMapLine
    from jcvi.formats.sizes import Sizes
    from jcvi.graphics.base import shorten
    from jcvi.graphics.dotplot import (
        plt,
        savefig,
        markup,
        normalize_axes,
        downsample,
        plot_breaks_and_labels,
        thousands,
    )

    p = OptionParser(dotplot.__doc__)
    p.set_outfile(outfile=None)
    opts, args, iopts = p.set_image_options(
        args, figsize="8x8", style="dark", dpi=90, cmap="copper"
    )

    if len(args) != 2:
        sys.exit(not p.print_help())

    csvfile, fastafile = args
    sizes = natsorted(Sizes(fastafile).mapping.items())
    seen = set()
    raw_data = []

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])  # the whole canvas
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # the dot plot

    fp = must_open(csvfile)
    for row in fp:
        m = CSVMapLine(row)
        seen.add(m.seqid)
        raw_data.append(m)

    # X-axis is the genome assembly
    ctgs, ctg_sizes = zip(*sizes)
    xsize = sum(ctg_sizes)
    qb = list(np.cumsum(ctg_sizes))
    qbreaks = list(zip(ctgs, [0] + qb, qb))
    qstarts = dict(zip(ctgs, [0] + qb))

    # Y-axis is the map
    key = lambda x: x.lg
    raw_data.sort(key=key)
    ssizes = {}
    for lg, d in groupby(raw_data, key=key):
        ssizes[lg] = max([x.cm for x in d])
    ssizes = natsorted(ssizes.items())
    lgs, lg_sizes = zip(*ssizes)
    ysize = sum(lg_sizes)
    sb = list(np.cumsum(lg_sizes))
    sbreaks = list(zip([("LG" + x) for x in lgs], [0] + sb, sb))
    sstarts = dict(zip(lgs, [0] + sb))

    # Re-code all the scatter dots
    data = [
        (qstarts[x.seqid] + x.pos, sstarts[x.lg] + x.cm, "g")
        for x in raw_data
        if (x.seqid in qstarts)
    ]
    npairs = downsample(data)

    x, y, c = zip(*data)
    ax.scatter(x, y, c=c, edgecolors="none", s=2, lw=0)

    # Flip X-Y label
    gy, gx = op.basename(csvfile).split(".")[:2]
    gx, gy = shorten(gx, maxchar=30), shorten(gy, maxchar=30)
    xlim, ylim = plot_breaks_and_labels(
        fig, root, ax, gx, gy, xsize, ysize, qbreaks, sbreaks
    )
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)

    title = "Alignment: {} vs {}".format(gx, gy)
    title += " ({} markers)".format(thousands(npairs))
    root.set_title(markup(title), x=0.5, y=0.96, color="k")
    logging.debug(title)
    normalize_axes(root)

    image_name = opts.outfile or (csvfile.rsplit(".", 1)[0] + "." + iopts.format)
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)
    fig.clear()


@lru_cache(maxsize=None)
def calc_ldscore(a, b):
    assert len(a) == len(b), "{0}\n{1}".format(a, b)
    # Assumes markers as A/B
    c = Counter(zip(a, b))
    c_aa = c[("A", "A")]
    c_ab = c[("A", "B")]
    c_ba = c[("B", "A")]
    c_bb = c[("B", "B")]
    n = c_aa + c_ab + c_ba + c_bb
    if n == 0:
        return 0

    f = 1.0 / n
    x_aa = c_aa * f
    x_ab = c_ab * f
    x_ba = c_ba * f
    x_bb = c_bb * f
    p_a = x_aa + x_ab
    p_b = x_ba + x_bb
    q_a = x_aa + x_ba
    q_b = x_ab + x_bb
    D = x_aa - p_a * q_a
    denominator = p_a * p_b * q_a * q_b
    if denominator == 0:
        return 0

    r2 = D * D / denominator
    return r2


def ld(args):
    """
    %prog ld map

    Calculate pairwise linkage disequilibrium given MSTmap.
    """
    from random import sample
    from jcvi.algorithms.matrix import symmetrize

    p = OptionParser(ld.__doc__)
    p.add_option(
        "--subsample",
        default=1000,
        type="int",
        help="Subsample markers to speed up",
    )
    opts, args, iopts = p.set_image_options(args, figsize="8x8")

    if len(args) != 1:
        sys.exit(not p.print_help())

    (mstmap,) = args
    subsample = opts.subsample
    data = MSTMap(mstmap)

    markerbedfile = mstmap + ".subsample.bed"
    ldmatrix = mstmap + ".subsample.matrix"
    # Take random subsample while keeping marker order
    if subsample < data.nmarkers:
        data = [data[x] for x in sorted(sample(range(len(data)), subsample))]
    else:
        logging.debug("Use all markers, --subsample ignored")

    nmarkers = len(data)
    if need_update(mstmap, (ldmatrix, markerbedfile)):
        fw = open(markerbedfile, "w")
        print("\n".join(x.bedline for x in data), file=fw)
        logging.debug(
            "Write marker set of size {0} to file `{1}`.".format(
                nmarkers, markerbedfile
            )
        )
        fw.close()

        M = np.zeros((nmarkers, nmarkers), dtype=float)
        for i, j in combinations(range(nmarkers), 2):
            a = data[i]
            b = data[j]
            M[i, j] = calc_ldscore(a.genotype, b.genotype)

        M = symmetrize(M)

        logging.debug("Write LD matrix to file `{0}`.".format(ldmatrix))
        M.tofile(ldmatrix)
    else:
        nmarkers = len(Bed(markerbedfile))
        M = np.fromfile(ldmatrix, dtype="float").reshape(nmarkers, nmarkers)
        logging.debug("LD matrix `{0}` exists ({1}x{1}).".format(ldmatrix, nmarkers))

    from jcvi.graphics.base import plt, savefig, Rectangle, draw_cmap

    plt.rcParams["axes.linewidth"] = 0

    fig = plt.figure(1, (iopts.w, iopts.h))
    root = fig.add_axes([0, 0, 1, 1])
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])  # the heatmap

    ax.matshow(M, cmap=iopts.cmap)

    # Plot chromosomes breaks
    bed = Bed(markerbedfile)
    xsize = len(bed)
    extent = (0, nmarkers)
    chr_labels = []
    ignore_size = 20

    for (seqid, beg, end) in bed.get_breaks():
        ignore = abs(end - beg) < ignore_size
        pos = (beg + end) / 2
        chr_labels.append((seqid, pos, ignore))
        if ignore:
            continue
        ax.plot((end, end), extent, "w-", lw=1)
        ax.plot(extent, (end, end), "w-", lw=1)

    # Plot chromosome labels
    for label, pos, ignore in chr_labels:
        pos = 0.1 + pos * 0.8 / xsize
        if not ignore:
            root.text(
                pos, 0.91, label, ha="center", va="bottom", rotation=45, color="grey"
            )
            root.text(0.09, pos, label, ha="right", va="center", color="grey")

    ax.set_xlim(extent)
    ax.set_ylim(extent)
    ax.set_axis_off()

    draw_cmap(root, "Pairwise LD (r2)", 0, 1, cmap=iopts.cmap)

    root.add_patch(Rectangle((0.1, 0.1), 0.8, 0.8, fill=False, ec="k", lw=2))
    m = mstmap.split(".")[0]
    root.text(
        0.5, 0.06, "Linkage Disequilibrium between {0} markers".format(m), ha="center"
    )

    root.set_xlim(0, 1)
    root.set_ylim(0, 1)
    root.set_axis_off()

    image_name = m + ".subsample" + "." + iopts.format
    savefig(image_name, dpi=iopts.dpi, iopts=iopts)


def header(args):
    """
    %prog header map conversion_table

    Rename lines in the map header. The mapping of old names to new names are
    stored in two-column `conversion_table`.
    """
    from jcvi.formats.base import DictFile

    p = OptionParser(header.__doc__)
    p.add_option("--prefix", default="", help="Prepend text to line number")
    p.add_option("--ids", help="Write ids to file")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mstmap, conversion_table = args
    data = MSTMap(mstmap)
    hd = data.header
    conversion = DictFile(conversion_table)
    newhd = [opts.prefix + conversion.get(x, x) for x in hd]

    print("\t".join(hd))
    print("--->")
    print("\t".join(newhd))

    ids = opts.ids
    if ids:
        fw = open(ids, "w")
        print("\n".join(newhd), file=fw)
        fw.close()


def rename(args):
    """
    %prog rename map markers.bed > renamed.map

    Rename markers according to the new mapping locations.
    """
    p = OptionParser(rename.__doc__)
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mstmap, bedfile = args
    markersbed = Bed(bedfile)
    markers = markersbed.order

    data = MSTMap(mstmap)
    header = data.header
    header = [header[0]] + ["seqid", "start"] + header[1:]
    renamed = []
    for b in data:
        m, geno = b.id, b.genotype
        om = m
        if m not in markers:
            m = m.rsplit(".", 1)[0]
            if m not in markers:
                continue

        i, mb = markers[m]
        renamed.append([om, mb.seqid, mb.start, "\t".join(list(geno))])

    renamed.sort(key=lambda x: (x[1], x[2]))
    fw = must_open(opts.outfile, "w")
    print("\t".join(header), file=fw)
    for d in renamed:
        print("\t".join(str(x) for x in d), file=fw)


def anchor(args):
    """
    %prog anchor map.bed markers.blast > anchored.bed

    Anchor scaffolds based on map.
    """
    from jcvi.formats.blast import bed

    p = OptionParser(anchor.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mapbed, blastfile = args
    bedfile = bed([blastfile])
    markersbed = Bed(bedfile)
    markers = markersbed.order

    mapbed = Bed(mapbed, sorted=False)
    for b in mapbed:
        m = b.accn
        if m not in markers:
            continue

        i, mb = markers[m]
        new_accn = "{0}:{1}-{2}".format(mb.seqid, mb.start, mb.end)
        b.accn = new_accn
        print(b)


def bed(args):
    """
    %prog fasta map.out

    Convert MSTMAP output into bed format.
    """
    p = OptionParser(bed.__doc__)
    p.add_option(
        "--switch",
        default=False,
        action="store_true",
        help="Switch reference and aligned map elements",
    )
    p.add_option(
        "--sep",
        default=".",
        help="Separator that is used to delimit scaffold and position in the marker name",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (mapout,) = args
    pf = mapout.split(".")[0]
    mapbed = pf + ".bed"
    bm = BinMap(mapout)
    bm.print_to_bed(mapbed, switch=opts.switch, sep=opts.sep)

    return mapbed


def fasta(args):
    """
    %prog fasta map.out scaffolds.fasta

    Extract marker sequences based on map.
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(fasta.__doc__)
    p.add_option(
        "--extend",
        default=1000,
        type="int",
        help="Extend seq flanking the gaps",
    )
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    mapout, sfasta = args
    Flank = opts.extend
    pf = mapout.split(".")[0]
    mapbed = pf + ".bed"
    bm = BinMap(mapout)
    bm.print_to_bed(mapbed)

    bed = Bed(mapbed, sorted=False)
    markersbed = pf + ".markers.bed"
    fw = open(markersbed, "w")
    sizes = Sizes(sfasta).mapping
    for b in bed:
        accn = b.accn
        scf, pos = accn.split(".")
        pos = int(pos)
        start = max(0, pos - Flank)
        end = min(pos + Flank, sizes[scf])
        print("\t".join(str(x) for x in (scf, start, end, accn)), file=fw)

    fw.close()

    fastaFromBed(markersbed, sfasta, name=True)


def hamming_distance(a, b, ignore=None):
    dist = 0
    for x, y in zip(a, b):
        if ignore and ignore in (x, y):
            continue
        if x != y:
            dist += 1
    return dist


OK, BREAK, END = range(3)


def check_markers(a, b, maxdiff):

    if a.seqid != b.seqid:
        return END, None
    diff = hamming_distance(a.genotype, b.genotype, ignore="-")
    max_allowed = len(a) * maxdiff
    if diff <= max_allowed:
        return OK, None

    return BREAK, (a.seqid, a.pos, b.pos)


def breakpoint(args):
    """
    %prog breakpoint mstmap.input > breakpoints.bed

    Find scaffold breakpoints using genetic map. Use variation.vcf.mstmap() to
    generate the input for this routine.
    """
    from more_itertools import pairwise

    p = OptionParser(breakpoint.__doc__)
    p.add_option(
        "--diff",
        default=0.1,
        type="float",
        help="Maximum ratio of differences allowed",
    )
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    (mstmap,) = args
    diff = opts.diff
    data = MSTMap(mstmap)

    # Remove singleton markers (avoid double cross-over)
    good = []
    nsingletons = 0
    for i in range(1, len(data) - 1):
        a = data[i]
        left_label, left_rr = check_markers(data[i - 1], a, diff)
        right_label, right_rr = check_markers(a, data[i + 1], diff)

        if left_label == BREAK and right_label == BREAK:
            nsingletons += 1
            continue

        good.append(a)

    logging.debug("A total of {0} singleton markers removed.".format(nsingletons))

    for a, b in pairwise(good):
        label, rr = check_markers(a, b, diff)
        if label == BREAK:
            print("\t".join(str(x) for x in rr))


if __name__ == "__main__":
    main()
