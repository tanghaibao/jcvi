#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Create the UCSC chain file which is needed to lift over from one coordinate
system to another.

File format:
<http://genome.ucsc.edu/goldenPath/help/chain.html>

chain 4900 chrY 58368225 + 25985403 25985638 chr5 151006098 - 43257292 43257528 1
  9       1       0
  10      0       5
  48

Header Line:
 chain score tName tSize tStrand tStart tEnd qName qSize qStrand qStart qEnd id
Alignment Data Lines
 size dt dq

NOTE: The last line of the alignment section contains only one number: the ungapped
alignment size of the last block.
"""

import os.path as op
import sys
import logging

from jcvi.formats.base import BaseFile, read_block
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, need_update, \
            which


class ChainLine (object):

    def __init__(self, chain, lines):
        self.chain = chain
        self.blocks = []
        for line in lines:
            atoms = line.split()
            if len(atoms) == 1:
                atoms += [0, 0]
            if len(atoms) == 0:
                continue

            self.blocks.append([int(x) for x in atoms])

        self.ungapped, self.dt, self.dq = zip(*self.blocks)
        self.ungapped = sum(self.ungapped)
        self.dt = sum(self.dt)
        self.dq = sum(self.dq)


class Chain (BaseFile):

    def __init__(self, filename):
        super(Chain, self).__init__(filename)
        self.chains = list(self.iter_chain())

        self.ungapped = sum(x.ungapped for x in self.chains)
        self.dt = sum(x.dt for x in self.chains)
        self.dq = sum(x.dq for x in self.chains)

    def __len__(self):
        return len(self.chains)

    def iter_chain(self):
        fp = open(self.filename)
        for row in fp:
            if row[0] != '#':
                break

        for chain, lines in read_block(fp, "chain"):
            lines = list(lines)
            yield ChainLine(chain, lines)


def main():

    actions = (
        ('blat', 'generate PSL file using BLAT'),
        ('last', 'generate PSL file using LAST'),
        ('frompsl', 'generate chain file from PSL format'),
        ('fromagp', 'generate chain file from AGP format'),
        ('summary', 'provide stats of the chain file'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def summary(args):
    """
    %prog summary old.new.chain old.fasta new.fasta

    Provide stats of the chain file.
    """
    from jcvi.formats.fasta import summary as fsummary
    from jcvi.utils.cbook import percentage, human_size

    p = OptionParser(summary.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    chainfile, oldfasta, newfasta = args
    chain = Chain(chainfile)
    ungapped, dt, dq = chain.ungapped, chain.dt, chain.dq
    print >> sys.stderr, "File `{0}` contains {1} chains.".\
                format(chainfile, len(chain))
    print >> sys.stderr, "ungapped={0} dt={1} dq={2}".\
                format(human_size(ungapped), human_size(dt), human_size(dq))

    oldreal, oldnn, oldlen = fsummary([oldfasta, "--outfile=/dev/null"])
    print >> sys.stderr, "Old fasta (`{0}`) mapped: {1}".\
                format(oldfasta, percentage(ungapped, oldreal))

    newreal, newnn, newlen = fsummary([newfasta, "--outfile=/dev/null"])
    print >> sys.stderr, "New fasta (`{0}`) mapped: {1}".\
                format(newfasta, percentage(ungapped, newreal))


def fromagp(args):
    """
    %prog fromagp agpfile componentfasta objectfasta

    Generate chain file from AGP format. The components represent the old
    genome (target) and the objects represent new genome (query).
    """
    from jcvi.formats.agp import AGP
    from jcvi.formats.sizes import Sizes

    p = OptionParser(fromagp.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    agpfile, componentfasta, objectfasta = args
    chainfile = agpfile.rsplit(".", 1)[0] + ".chain"
    fw = open(chainfile, "w")
    agp = AGP(agpfile)
    componentsizes = Sizes(componentfasta).mapping
    objectsizes = Sizes(objectfasta).mapping
    chain = "chain"
    score = 1000
    tStrand = "+"
    id = 0
    for a in agp:
        if a.is_gap:
            continue

        tName = a.component_id
        tSize = componentsizes[tName]
        tStart = a.component_beg
        tEnd = a.component_end
        tStart -= 1

        qName = a.object
        qSize = objectsizes[qName]
        qStrand = "-" if a.orientation == "-" else "+"
        qStart = a.object_beg
        qEnd = a.object_end
        if qStrand == '-':
            _qStart = qSize - qEnd + 1
            _qEnd = qSize - qStart + 1
            qStart, qEnd = _qStart, _qEnd
        qStart -= 1

        id += 1
        size = a.object_span
        headerline = "\t".join(str(x) for x in (
             chain, score, tName, tSize, tStrand, tStart,
             tEnd, qName, qSize, qStrand, qStart, qEnd, id
        ))
        alignmentline = size
        print >> fw, headerline
        print >> fw, alignmentline
        print >> fw

    fw.close()
    logging.debug("File written to `{0}`.".format(chainfile))


def faToTwoBit(fastafile):
    twobitfile = fastafile.rsplit(".", 1)[0] + ".2bit"
    cmd = "faToTwoBit {0} {1}".format(fastafile, twobitfile)
    if need_update(fastafile, twobitfile):
        sh(cmd)
    return twobitfile


def blat(args):
    """
    %prog blat old.fasta new.fasta

    Generate psl file using blat.
    """
    p = OptionParser(blat.__doc__)
    p.add_option("--minscore", default=100, type="int",
                 help="Matches minus mismatches gap penalty [default: %default]")
    p.add_option("--minid", default=98, type="int",
                 help="Minimum sequence identity [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    oldfasta, newfasta = args
    twobitfiles = []
    for fastafile in args:
        tbfile = faToTwoBit(fastafile)
        twobitfiles.append(tbfile)

    oldtwobit, newtwobit = twobitfiles
    cmd = "pblat -threads={0}".format(opts.cpus) if which("pblat") else "blat"
    cmd += " {0} {1}".format(oldtwobit, newfasta)
    cmd += " -tileSize=12 -minScore={0} -minIdentity={1} ".\
                format(opts.minscore, opts.minid)
    pslfile = "{0}.{1}.psl".format(*(op.basename(x).split('.')[0] \
                for x in (newfasta, oldfasta)))
    cmd += pslfile
    sh(cmd)


def last(args):
    """
    %prog last old.fasta new.fasta

    Generate psl file using last. Calles apps.last() but with special
    parameters: -r5 -q95 -a0 -b95 -e500, which only reports alignments larger
    than 100 bp and >=95 % identity.
    """
    from jcvi.apps.last import main as lastapp

    p = OptionParser(last.__doc__)
    p.add_option("--distant", default=False, action="store_true",
                 help="Assume distant relations")
    p.add_option("--minscore", default=100, type="int",
                 help="Filter alignments by how many bases match [default: %default]")
    p.add_option("--minid", default=95, type="int",
                 help="Minimum sequence identity [default: %default]")
    p.set_outfile()

    opts, args = p.parse_args(args)
    if len(args) != 2:
        sys.exit(not p.print_help())

    oldfasta, newfasta = args
    args = [oldfasta, newfasta, "--format=maf", \
                "--outfile={0}".format(opts.outfile)]

    minscore = opts.minscore
    minid = opts.minid

    r = 100 - minid
    q = minid
    e = minscore * r

    extra = r'--params=-r{0} -q{1} -a0 -b{1} -e{2}'.format(r, q, e)
    if not opts.distant:
        args.append(extra)

    lastapp(args)


def frompsl(args):
    """
    %prog frompsl old.new.psl old.fasta new.fasta

    Generate chain file from psl file. The pipeline is describe in:
    <http://genomewiki.ucsc.edu/index.php/Minimal_Steps_For_LiftOver>
    """
    from jcvi.formats.sizes import Sizes

    p = OptionParser(frompsl.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    pslfile, oldfasta, newfasta = args
    pf = oldfasta.split(".")[0]

    # Chain together alignments from using axtChain
    chainfile = pf + ".chain"
    twobitfiles = []
    for fastafile in (oldfasta, newfasta):
        tbfile = faToTwoBit(fastafile)
        twobitfiles.append(tbfile)
    oldtwobit, newtwobit = twobitfiles

    if need_update(pslfile, chainfile):
        cmd = "axtChain -linearGap=medium -psl {0}".format(pslfile)
        cmd += " {0} {1} {2}".format(oldtwobit, newtwobit, chainfile)
        sh(cmd)

    # Sort chain files
    sortedchain = chainfile.rsplit(".", 1)[0] + ".sorted.chain"
    if need_update(chainfile, sortedchain):
        cmd = "chainSort {0} {1}".format(chainfile, sortedchain)
        sh(cmd)

    # Make alignment nets from chains
    netfile = pf + ".net"
    oldsizes = Sizes(oldfasta).filename
    newsizes = Sizes(newfasta).filename
    if need_update((sortedchain, oldsizes, newsizes), netfile):
        cmd = "chainNet {0} {1} {2}".format(sortedchain, oldsizes, newsizes)
        cmd += " {0} /dev/null".format(netfile)
        sh(cmd)

    # Create liftOver chain file
    liftoverfile = pf + ".liftover.chain"
    if need_update((netfile, sortedchain), liftoverfile):
        cmd = "netChainSubset {0} {1} {2}".\
                format(netfile, sortedchain, liftoverfile)
        sh(cmd)


if __name__ == '__main__':
    main()
