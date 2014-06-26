#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Locate IES sequences within MIC genome of tetrahymena.
"""

import sys
import logging

from itertools import groupby

from jcvi.utils.counter import Counter
from jcvi.utils.range import Range, range_interleave, range_chain
from jcvi.formats.bed import Bed, sort, depth, some
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


def main():

    actions = (
        ('deletion', 'find IES based on mapping MAC reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def deletion(args):
    """
    %prog deletion mac.mic.bed mic.gaps.bed

    Find IES based on mapping MAC reads to MIC genome.
    """
    p = OptionParser(deletion.__doc__)
    p.add_option("--mindepth", default=3, type="int",
                 help="Minimum depth to call a deletion")
    p.add_option("--minspan", default=30, type="int",
                 help="Minimum span to call a deletion")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    bedfile, gapsbedfile = args

    pf = bedfile.rsplit(".", 1)[0]
    sortedbedfile = pf + ".sorted.bed"
    if need_update(bedfile, sortedbedfile):
        sort([bedfile, "-u", "--accn"])

    # Find reads that contain multiple matches
    bed = Bed(sortedbedfile, sorted=False)
    ibedfile = pf + ".d.bed"
    fw = open(ibedfile, "w")
    logging.debug("Write deletions to `{0}`.".format(ibedfile))
    for accn, bb in groupby(bed, key=lambda x: x.accn):
        bb = list(bb)
        branges = [(x.seqid, x.start, x.end) for x in bb]
        iranges = range_interleave(branges)
        for seqid, start, end in iranges:
            if end - start + 1 < opts.minspan:
                continue
            print >> fw, "\t".join(str(x) for x in \
                        (seqid, start - 1, end, accn + '-d'))
    fw.close()

    # Uniqify the insertions and count occurrences
    bed = Bed(ibedfile)
    countbedfile = pf + ".uniq.bed"
    fw = open(countbedfile, "w")
    logging.debug("Write counts to `{0}`.".format(countbedfile))
    registry = Counter((x.seqid, x.start, x.end) for x in bed)
    ies_id = 1
    for (seqid, start, end), count in registry.items():
        ies_name = "{0:05d}-r{1}".format(ies_id, count)
        if count < opts.mindepth:
            continue
        print >> fw, "\t".join(str(x) for x in \
                        (seqid, start - 1, end, ies_name))
        ies_id += 1
    fw.close()
    sort([countbedfile, "-i"])

    # Remove deletions that contain some read depth
    depthbedfile = pf + ".depth.bed"
    depth([sortedbedfile, countbedfile, "--outfile={0}".format(depthbedfile)])
    validbedfile = pf + ".valid.bed"
    fw = open(validbedfile, "w")
    logging.debug("Filter valid deletions to `{0}`.".format(validbedfile))
    bed = Bed(depthbedfile)
    for b in bed:
        if float(b.score) >= 1:
            continue
        print >> fw, b
    fw.close()

    # Remove deletions that contain sequencing gaps on its flanks
    flanksbedfile = pf + ".flanks.bed"
    fw = open(flanksbedfile, "w")
    bed = Bed(validbedfile)
    flank = 100
    logging.debug("Write deletion flanks to `{0}`.".format(flanksbedfile))
    for b in bed:
        start, end = b.start, b.end
        b.start, b.end = start, min(start + flank - 1, end)
        print >> fw, b
        b.start, b.end = max(start, end - flank + 1), end
        print >> fw, b
    fw.close()

    intersectidsfile = pf + ".intersect.ids"
    cmd = "intersectBed -a {0} -b {1}".format(flanksbedfile, gapsbedfile)
    cmd += " | cut -f4 | sort -u"
    selectedbedfile = pf + ".selected.bed"
    sh(cmd, outfile=intersectidsfile)
    some([validbedfile, intersectidsfile, "-v",
            "--outfile={0}".format(selectedbedfile)])

    # Find best-scoring non-overlapping set
    bed = Bed(selectedbedfile)
    iesbedfile = pf + ".ies.bed"
    fw = open(iesbedfile, "w")
    logging.debug("Write IES to `{0}`.".format(iesbedfile))
    branges = [Range(x.seqid, x.start, x.end, int(x.accn.rsplit("r")[-1]), i) \
                    for i, x in enumerate(bed)]
    iranges, iscore = range_chain(branges)
    logging.debug("Best chain score: {0} ({1} IES)".\
                    format(iscore, len(iranges)))
    ies_id = 1
    for seqid, start, end, score, id in iranges:
        ies_name = "IES-{0:05d}-r{1}".format(ies_id, score)
        span = end - start + 1
        print >> fw, "\t".join(str(x) for x in \
                        (seqid, start - 1, end, ies_name, span))
        ies_id += 1
    fw.close()


if __name__ == '__main__':
    main()
