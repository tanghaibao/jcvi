#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
SSPACE scaffolding-related operations.
"""
from __future__ import print_function

import os.path as op
import sys
import logging

from copy import deepcopy
from more_itertools import pairwise

from jcvi.formats.fasta import gaps
from jcvi.formats.sizes import Sizes
from jcvi.formats.base import BaseFile, read_block, write_file
from jcvi.formats.agp import AGP, AGPLine, reindex, tidy
from jcvi.algorithms.graph import BiGraph
from jcvi.apps.base import OptionParser, ActionDispatcher, download


NO_UPDATE, INSERT_BEFORE, INSERT_AFTER, INSERT_BETWEEN = (
    "NO_UPDATE",
    "INSERT_BEFORE",
    "INSERT_AFTER",
    "INSERT_BETWEEN",
)


class EvidenceLine(object):
    def __init__(self, row, sizes):
        # f_tig3222|size7922|links348|gaps-109|merged16
        args = row.strip().split("|")
        nargs = len(args)

        tig = args[0]
        o, mtig = tig.split("_")
        tig = int(mtig.replace("tig", ""))
        assert o in ("f", "r")
        self.o = ">" if o == "f" else "<"
        self.oo = "+" if o == "f" else "-"

        name, size = sizes[tig]
        self.tig = name
        self.size = int(args[1].replace("size", ""))
        assert self.size == size, "{0} and {1} size mismatch".format(mtig, name)

        if nargs > 2:
            self.links = int(args[2].replace("links", ""))
            self.gaps = int(args[3].replace("gaps", ""))
        if nargs > 4:
            self.merged = int(args[4].replace("merged", ""))


class EvidenceFile(BaseFile):
    def __init__(self, filename, fastafile):
        super(EvidenceFile, self).__init__(filename)
        sz = Sizes(fastafile)
        sizes = [None]  # tig-list starts at 1
        for name, size in sz.iter_sizes():
            sizes.append((name, size))
        self.sizes = sizes
        self.sz = sz.mapping
        self.scf = {}

    def iter_scaffold(self):
        filename = self.filename
        sizes = self.sizes
        fp = open(filename)
        for header, lines in read_block(fp, ">"):
            scaffold, size, tigs = header[1:].split("|")
            lines = [EvidenceLine(x, sizes) for x in lines if x.strip()]
            yield scaffold, lines

    @property
    def graph(self):
        g = BiGraph()
        for scaffold, lines in self.iter_scaffold():
            self.scf[scaffold] = [x.tig for x in lines]

            for a, b in pairwise(lines):
                g.add_edge(a.tig, b.tig, a.o, b.o, length=a.gaps)

            if len(lines) == 1:  # Singleton scaffold
                a = lines[0]
                g.add_node(a.tig)

        return g

    def write_agp(self, filename):
        sizes = self.sz
        agp = []
        for scaffold, lines in self.iter_scaffold():
            for a, b in pairwise(lines):
                cline = AGPLine.cline(scaffold, a.tig, sizes, a.oo)
                gline = AGPLine.gline(scaffold, a.gaps)
                agp.append(cline)
                agp.append(gline)
            a = lines[-1]
            cline = AGPLine.cline(scaffold, a.tig, sizes, a.oo)
            agp.append(cline)

        fw = open(filename, "w")
        for a in agp:
            print(a, file=fw)
        fw.close()

        reindex([filename, "--inplace"])
        return filename


def main():

    actions = (
        ("scaffold", "run SSPACE scaffolding"),
        ("close", "run GapFiller to fill gaps"),
        ("agp", "convert SSPACE scaffold structure to AGP format"),
        ("embed", "embed contigs to upgrade existing structure"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def write_libraries(fastqs, aligner=None):
    from jcvi.assembly.base import get_libs

    libs = get_libs(fastqs)
    assert libs

    libtxt = "libraries.txt"
    contents = []
    for i, (lib, fns) in enumerate(libs):
        fns = " ".join(fns)
        pe = "RF" if lib.read_orientation == "outward" else "FR"
        cc = ["lib{0}".format(i + 1), fns, lib.size, 0.75, pe]
        if aligner:
            cc.insert(1, aligner)
        libline = " ".join(str(x) for x in cc)
        contents.append(libline)

    write_file(libtxt, "\n".join(contents), tee=True)
    return libtxt


def close(args):
    """
    %prog close scaffolds.fasta PE*.fastq

    Run GapFiller to fill gaps.
    """
    p = OptionParser(close.__doc__)
    p.set_home("gapfiller")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    scaffolds = args[0]
    libtxt = write_libraries(args[1:], aligner="bwa")

    cmd = "perl " + op.join(opts.gapfiller_home, "GapFiller.pl")
    cmd += " -l {0} -s {1} -T {2}".format(libtxt, scaffolds, opts.cpus)
    runsh = "run.sh"
    write_file(runsh, cmd)


def scaffold(args):
    """
    %prog scaffold contigs.fasta MP*.fastq

    Run SSPACE scaffolding.
    """
    p = OptionParser(scaffold.__doc__)
    p.set_aligner(aligner="bwa")
    p.set_home("sspace")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    contigs = args[0]
    libtxt = write_libraries(args[1:], aligner=opts.aligner)
    # Requires getopts.pl which may be missing
    download("http://web.vims.edu/bridge/bridge2/aw/lib/getopts.pl")

    cmd = "perl " + op.join(opts.sspace_home, "SSPACE_Standard_v3.0.pl")
    cmd += " -l {0} -s {1} -T {2}".format(libtxt, contigs, opts.cpus)
    runsh = "run.sh"
    write_file(runsh, cmd)


def agp(args):
    """
    %prog agp evidencefile contigs.fasta

    Convert SSPACE scaffold structure to AGP format.
    """
    p = OptionParser(agp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    evidencefile, contigs = args
    ef = EvidenceFile(evidencefile, contigs)

    agpfile = evidencefile.replace(".evidence", ".agp")
    ef.write_agp(agpfile)


def get_target(p, name):
    before, before_tag = p.get_next(name, ">")
    if not before:  # Start of a scaffold
        return (None, ">")
    next, next_tag = p.get_next(name)
    if not next:  # End of a scaffold
        return (None, "<")
    # Internal to a scaffold
    return (next.v, "<")


def get_orientation(o, status):
    o = "+" if o == "<" else "-"
    if status == INSERT_BEFORE:  # Flip orientation for backward traversal
        o = "+" if o == "-" else "-"
    return o


def path_to_agp(g, path, object, sizes, status):
    lines = []
    for (a, ao), (b, bo) in pairwise(path):
        ao = get_orientation(ao, status)
        e = g.get_edge(a.v, b.v)
        cline = AGPLine.cline(object, a.v, sizes, ao)
        gline = AGPLine.gline(object, e.length)
        lines.append(cline)
        lines.append(gline)
    # Do not forget the last one
    z, zo = path[-1]
    zo = get_orientation(zo, status)
    cline = AGPLine.cline(object, z.v, sizes, zo)
    lines.append(cline)

    return lines


def embed(args):
    """
    %prog embed evidencefile scaffolds.fasta contigs.fasta

    Use SSPACE evidencefile to scaffold contigs into existing scaffold
    structure, as in `scaffolds.fasta`. Contigs.fasta were used by SSPACE
    directly to scaffold.

    Rules:
    1. Only update existing structure by embedding contigs small enough to fit.
    2. Promote singleton contigs only if they are big (>= min_length).
    """
    p = OptionParser(embed.__doc__)
    p.set_mingap(default=10)
    p.add_option(
        "--min_length",
        default=200,
        type="int",
        help="Minimum length to consider",
    )
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    evidencefile, scaffolds, contigs = args
    min_length = opts.min_length
    splitfasta, oagp, cagp = gaps(
        [scaffolds, "--split", "--mingap={0}".format(opts.mingap)]
    )

    agp = AGP(cagp)
    p = agp.graph

    ef = EvidenceFile(evidencefile, contigs)
    sizes = ef.sz
    q = ef.graph

    logging.debug("Reference graph: {0}".format(p))
    logging.debug("Patch graph: {0}".format(q))

    newagp = deepcopy(agp)

    seen = set()
    deleted = set()
    for a in agp:
        if a.is_gap:
            continue

        name = a.component_id
        object = a.object
        if name in deleted:
            print("* Skip {0}, already embedded".format(name), file=sys.stderr)
            continue

        seen.add(name)

        target_name, tag = get_target(p, name)
        path = q.get_path(name, target_name, tag=tag)
        path_size = sum([sizes[x.v] for x, t in path]) if path else None
        status = NO_UPDATE

        # Heuristic, the patch must not be too long
        if path and path_size > min_length and len(path) > 3:
            path = None

        if not path:
            print(name, target_name, path, path_size, status, file=sys.stderr)
            continue

        backward = False
        for x, t in path:
            if x.v in seen:
                print(
                    "* Does not allow backward" " patch on {0}".format(x.v),
                    file=sys.stderr,
                )
                backward = True
                break

        if backward:
            continue

        # Build the path plus the ends
        vv = q.get_node(name)
        path.appendleft((vv, tag))
        if tag == ">":
            path.reverse()
            status = INSERT_BEFORE
        elif target_name is None:
            status = INSERT_AFTER
        else:
            target = q.get_node(target_name)
            path.append((target, tag))
            status = INSERT_BETWEEN

        print(name, target_name, path, path_size, status, file=sys.stderr)

        # Trim the ends off from the constructed AGPLines
        lines = path_to_agp(q, path, object, sizes, status)
        if status == INSERT_BEFORE:
            lines = lines[:-1]
            td = newagp.insert_lines(name, lines, delete=True, verbose=True)
        elif status == INSERT_AFTER:
            lines = lines[1:]
            td = newagp.insert_lines(name, lines, after=True, delete=True, verbose=True)
        else:
            lines = lines[1:-1]
            td = newagp.update_between(
                name, target_name, lines, delete=True, verbose=True
            )
        deleted |= td
        seen |= td

    # Recruite big singleton contigs
    CUTOFF = opts.min_length
    for ctg, size in sizes.items():
        if ctg in seen:
            continue
        if size < CUTOFF:
            continue
        newagp.append(AGPLine.cline(ctg, ctg, sizes, "?"))

    # Write a new AGP file
    newagpfile = "embedded.agp"
    newagp.print_to_file(newagpfile, index=True)
    tidy([newagpfile, contigs])


if __name__ == "__main__":
    main()
