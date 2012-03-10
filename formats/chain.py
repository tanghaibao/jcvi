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

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug, sh, need_update
debug()


def main():

    actions = (
        ('blat', 'generate PSL file using blat'),
        ('frompsl', 'generate chain file from PSL format'),
        ('fromagp', 'generate chain file from AGP format'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    oldfasta, newfasta = args
    twobitfiles = []
    for fastafile in args:
        tbfile = faToTwoBit(fastafile)
        twobitfiles.append(tbfile)

    oldtwobit, newtwobit = twobitfiles
    cmd = "blat {0} {1}".format(oldtwobit, newfasta)
    cmd += " -tileSize=12 -minScore={0} -minIdentity={1} ".\
                format(opts.minscore, opts.minid)
    pslfile = "{0}.{1}.psl".format(*(op.basename(x).split('.')[0] \
                for x in (newfasta, oldfasta)))
    cmd += pslfile
    sh(cmd)


def lastz(args):
    """
    %prog lastz old.fasta new.fasta

    Generate psl file using lastz. To generate psl between two distant
    organisms, turn on --distant. Otherwise, lastz will run using the
    parameters fit for assembly-assembly comparison.

    Bob Harris recommended:
    --seed=match12 --notransition --step=20
    --exact=50 --mismatch=1,100
    --identity=98 --matchcount=100
    """
    p = OptionParser(lastz.__doc__)
    p.add_option("--distant", default=False, action="store_true",
                 help="Assume distant relations")
    p.add_option("--minscore", default=100, type="int",
                 help="Filter alignments by how many bases match [default: %default]")
    p.add_option("--minid", default=98, type="int",
                 help="Minimum sequence identity [default: %default]")


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

    # Use liftUp to change the coordinate system. Requires .lft files
    # This step is skipped as the output psl is the same as input?

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
