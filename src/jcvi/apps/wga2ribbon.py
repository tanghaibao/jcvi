#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Convert PAF/LASTZ alignments to bed and block files for jcvi.ribbon
using layout file as guide.
"""

import argparse
from collections import namedtuple
import glob
from operator import attrgetter
import os
import shutil
import sys

from jcvi.graphics.base import AbstractLayout


def mainArgs():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert whole-genome alignments to block & bed format \
                    for ribbon.",
        prog="jvci.apps.wga2ribbon",
    )
    parser.add_argument(
        "-f",
        "--format",
        type=str,
        required=True,
        choices=["paf", "lastz"],
        help="Set alignment format as paf or lastz.",
    )
    parser.add_argument(
        "-a",
        "--alignments",
        action="append",
        required=True,
        help='Set arg for each alignment file to be imported.\
                        Has format "[Source chrm Idx],[Hit chrm Idx],\
                        [filename]" i.e. For file with alignmnets between the \
                        first and second tracks in ribbon.layout: \
                        -a "0,1,A2B_Alignment.tab"',
    )
    parser.add_argument(
        "-l", "--layout", type=str, required=True, help="Ribbon plot layout file."
    )
    parser.add_argument(
        "--minID",
        type=float,
        default=50,
        help="Extract alignments with identity => than minID",
    )
    parser.add_argument(
        "--outBed",
        type=str,
        default="ribbons.bed",
        help="Write alignment block coords to this bedfile",
    )
    parser.add_argument(
        "--outBlocks",
        type=str,
        default="ribbons.block",
        help="Write alignment block names to this file.",
    )
    args = parser.parse_args()
    return args


class LayoutLine(object):
    #  0-11: x,y,rotation,ha,va,color,ratio,label,chrmName,rStart,rEnd,chrmMax
    def __init__(self, row, delimiter=","):
        self.hidden = row[0] == "*"
        if self.hidden:
            row = row[1:]
        args = row.rstrip().split(delimiter)
        args = [x.strip() for x in args]
        self.x = float(args[0])
        self.y = float(args[1])
        self.rotation = int(args[2])
        self.ha = args[3]
        self.va = args[4]
        self.color = args[5]
        self.ratio = 1
        if len(args) > 6:
            self.ratio = float(args[6])
        if len(args) > 7:
            self.label = args[7].strip()
        else:
            self.label = None
        if len(args) > 8:
            self.ChrmName = str(args[8])
        if len(args) > 9:
            self.rStart = int(args[9])
            self.rEnd = int(args[10])
        else:
            self.rStart = None
            self.rEnd = None
        if len(args) > 10:
            self.chrmMax = int(args[10])
        else:
            self.chrmMax = None


class Layout(AbstractLayout):
    def __init__(self, filename, delimiter=","):
        super(Layout, self).__init__(filename)
        fp = open(filename)
        self.edges = []
        for row in fp:
            # Skip blank lines
            if row.strip():
                # Skip header rows
                if row[0] == "#":
                    continue
                # Import edges
                if row[0] == "e":
                    args = row.rstrip().split(delimiter)
                    args = [x.strip() for x in args]
                    # From index 1 and 2
                    a, b = args[1:3]
                    a, b = int(a), int(b)
                    assert args[0] == "e"
                    self.edges.append((a, b))
                else:
                    # Parse all other lines as sequence tracks
                    self.append(LayoutLine(row, delimiter=delimiter))

        # Track number of layout lines added
        self.lines = len(self)

        if 3 <= len(self) <= 8:
            self.assign_colors()


def readLASTZ(infile, minID=50):
    """
    Read in LASTZ result file from LASTZ_genome_align.sh
    Populate nested dictionary of hits keyed by Target
    and then Query scaffold names.
    # Check diff between standard LASTZ and mimeo-LASTZ formats
    """
    with open(infile) as f:
        content = f.readlines()
    content = [x.strip().split() for x in content]
    hitsDict = dict()
    # Set named tuple format
    hitTup = namedtuple(
        "Elem",
        [
            "t_start",
            "t_end",
            "t_strand",
            "q_start",
            "q_end",
            "q_strand",
            "idPct",
            "UID",
        ],
    )
    counter = 0
    # Read in split rows
    for row in content:
        # Ignore lines begining with '#'
        if row[0][0] == "#":
            continue
        # Count if hit identity exceeds threshold
        # CS10_Chromosome_08    +    1993504    1993793    supercont1.32    +    15987    16251    8424    72.9
        elif float(row[9]) >= minID:
            counter += 1
            UID = counter
            t_name = str(row[0])
            t_strand = str(row[1])
            t_start = int(row[2])  # - 1 # Convert from idx '1' to idx '0'
            t_end = int(row[3])  # - 1 # Convert from idx '1' to idx '0'
            q_name = str(row[4])
            q_strand = str(row[5])
            idPct = float(row[9])
            # Check that start position < end position
            if int(row[6]) < int(row[7]):
                q_start = int(row[6])  # - 1 # Convert from idx '1' to idx '0'
                q_end = int(row[7])  # - 1 # Convert from idx '1' to idx '0'
            # Correct for inverted coordinates
            else:
                print(("Inverting query sequence coordinates for record number:", UID))
                q_end = int(row[6])  # - 1 # Convert from idx '1' to idx '0'
                q_start = int(row[7])  # - 1 # Convert from idx '1' to idx '0'
            # Create target scaffold dict if not seen
            if t_name not in list(hitsDict.keys()):
                hitsDict[t_name] = dict()
            # Create query scaffold list if not seen
            if q_name not in list(hitsDict[t_name].keys()):
                hitsDict[t_name][q_name] = list()
            # Write record to target:query list as named tuple
            hitsDict[t_name][q_name].append(
                hitTup(t_start, t_end, t_strand, q_start, q_end, q_strand, idPct, UID)
            )
    return hitsDict


def readPAF(infile, minID=50):
    """
    Format info here: https://github.com/lh3/miniasm/blob/master/PAF.md
    From PAF file populate nested dictionary of hits keyed by Target (Ref) and
    then Query scaffold names.
    # Check target and query seq defs for PAF and LASTZ
    # Check what PAF does with rev orientation coords - flipped?
    # Check if ribbon wants intervals indexed from zero
    # PAF is zero based
    # Check how strand orientation is assigned
    """
    # |Col|Type  |Description                               |
    # |--:|:----:|:-----------------------------------------|
    # |0  |string|Query sequence name                       |
    # |1  |int   |Query sequence length                     |
    # |2  |int   |Query start (0-based; BED-like; closed)   |
    # |3  |int   |Query end (0-based; BED-like; open)       |
    # |4  |char  |Relative strand: "+" or "-"               |
    # |5  |string|Target sequence name                      |
    # |6  |int   |Target sequence length                    |
    # |7  |int   |Target start on original strand (0-based) |
    # |8  |int   |Target end on original strand (0-based)   |
    # |9  |int   |Number of residue matches                 |
    # |10 |int   |Alignment block length                    |
    # |11 |int   |Mapping quality (0-255; 255 for missing)  |
    with open(infile) as f:
        content = f.readlines()
    content = [x.strip().split() for x in content]
    hitsDict = dict()
    # Set named tuple format
    hitTup = namedtuple(
        "Elem",
        [
            "t_start",
            "t_end",
            "t_strand",
            "q_start",
            "q_end",
            "q_strand",
            "idPct",
            "UID",
        ],
    )
    counter = 0
    # Read in split rows
    for row in content:
        # Ignore lines begining with '#'
        if row[0][0] == "#":
            continue
        # Count if hit identity exceeds threshold
        # For PAF calc idenity from matching bases and alignment length.
        rowID = float(row[9]) / float(row[10])
        if rowID >= minID:
            counter += 1
            UID = counter
            t_name = str(row[5])
            t_strand = "+"
            t_start = int(row[7])  # + 1 # Convert from idx '0' to idx '1'
            t_end = int(row[8])  # + 1 # Convert from idx '0' to idx '1'
            q_name = str(row[0])
            q_strand = str(row[4])
            idPct = rowID
            # Check that start position < end position
            if int(row[2]) < int(row[3]):
                q_start = int(row[2])  # + 1 # Convert from idx '0' to idx '1'
                q_end = int(row[3])  # + 1 # Convert from idx '0' to idx '1'
            # Correct for inverted coordinates
            else:
                print(("Inverting query sequence coordinates for record number: ", UID))
                q_end = int(row[2])  # + 1 # Convert from idx '0' to idx '1'
                q_start = int(row[3])  # + 1 # Convert from idx '0' to idx '1'
            # Create target scaffold dict if not seen
            if t_name not in list(hitsDict.keys()):
                hitsDict[t_name] = dict()
            # Create query scaffold list if not seen
            if q_name not in list(hitsDict[t_name].keys()):
                hitsDict[t_name][q_name] = list()
            # Write record to target:query list as named tuple
            hitsDict[t_name][q_name].append(
                hitTup(t_start, t_end, t_strand, q_start, q_end, q_strand, idPct, UID)
            )
    return hitsDict


def main():
    # Get cmd line args
    args = mainArgs()
    plotlines = Layout(args.layout)

    # Check number of alignment files = layout lines -1 (i.e. AvB.tab, BvC.tab)
    # if plotlines.lines != len(args.alignments) + 1:
    #    print('Error: Alignment file list must be == number sequence lines in layout - 1')
    #    sys.exit(0)

    # Init dict with layout row indexs as keys and empty list as value
    bedlines = dict()
    for i in range(plotlines.lines):
        bedlines[i] = list()
    # Init blocklines
    blocklines = list()
    for coords in args.alignments:
        assert len(coords.split(",")) == 3
        selfIdx, hitIdx, alnFile = coords.split(",")
        selfIdx = int(selfIdx)
        hitIdx = int(hitIdx)
        hits = readLASTZ(alnFile, minID=args.minID)
        # Extract hits between sequence selfIdx and hitIdx as per layout file.
        for x in sorted(
            hits[plotlines[selfIdx].ChrmName][plotlines[hitIdx].ChrmName],
            key=lambda x: (int(x.t_start), int(x.t_end)),
        ):
            tName = str(plotlines[selfIdx].ChrmName)
            qName = str(plotlines[hitIdx].ChrmName)
            tUID = "_".join(
                [
                    str(selfIdx),
                    str(hitIdx),
                    str(selfIdx),
                    "ID" + str(x.UID),
                    tName,
                    str(x.t_start),
                    str(x.t_end),
                ]
            )
            qUID = "_".join(
                [
                    str(selfIdx),
                    str(hitIdx),
                    str(hitIdx),
                    "ID" + str(x.UID),
                    qName,
                    str(x.q_start),
                    str(x.q_end),
                ]
            )
            bedlines[selfIdx].append(
                tuple(
                    [
                        tName,
                        str(x.t_start),
                        str(x.t_end),
                        tUID,
                        str(x.idPct),
                        str(x.t_strand),
                    ]
                )
            )
            bedlines[hitIdx].append(
                tuple(
                    [
                        qName,
                        str(x.q_start),
                        str(x.q_end),
                        qUID,
                        str(x.idPct),
                        str(x.q_strand),
                    ]
                )
            )
            blockline = ["."] * plotlines.lines
            blockline[selfIdx] = tUID
            blockline[hitIdx] = qUID
            blocklines.append("\t".join(blockline))

    # Sort bedlines for each seq in layout file.
    # Sort by Chrm name, start, end, UID
    with open(args.outBed, "w") as f:
        for layoutIdx in list(bedlines.keys()):
            for line in sorted(
                bedlines[layoutIdx], key=lambda x: (x[0], int(x[1]), int(x[2]), x[3])
            ):
                f.write("\t".join(line) + "\n")

    # Write blocklines
    with open(args.outBlocks, "w") as f:
        for line in blocklines:
            f.write(line + "\n")


if __name__ == "__main__":
    main()
