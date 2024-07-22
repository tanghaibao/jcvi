from collections import defaultdict
from typing import Dict, Tuple

from ..apps.base import logger
from ..formats.base import BaseFile, read_block, must_open
from ..utils.range import Range


class AnchorFile(BaseFile):
    def __init__(self, filename, minsize=0):
        super().__init__(filename)
        self.blocks = list(self.iter_blocks(minsize=minsize))

    def iter_blocks(self, minsize=0):
        fp = open(self.filename)
        for _, lines in read_block(fp, "#"):
            lines = [x.split() for x in lines]
            if len(lines) >= minsize:
                yield lines

    def iter_pairs(self, minsize=0):
        block_id = -1
        for rows in self.iter_blocks(minsize=minsize):
            block_id += 1
            for row in rows:
                a, b = row[:2]
                yield a, b, block_id

    def make_ranges(self, order, clip=10):
        """Prepare anchors information into a set of ranges for chaining"""
        ranges = []
        block_pairs = defaultdict(dict)
        blocks = self.blocks
        for i, ib in enumerate(blocks):
            q, s, t = zip(*ib)
            if q[0] not in order:
                q, s = s, q

            r = make_range(q, s, t, i, order, block_pairs, clip=clip)
            ranges.append(r)

            assert q[0] in order
            if s[0] not in order:
                continue

            # is_self comparison
            q, s = s, q
            r = make_range(q, s, t, i, order, block_pairs, clip=clip)
            ranges.append(r)
        return ranges, block_pairs

    def filter_blocks(self, accepted: Dict[Tuple[str, str], str]):
        """
        Filter the blocks based on the accepted pairs. This is used to update
        the anchors so that they match the info in the LAST file.
        """
        new_blocks = []
        nremoved = 0
        ncorrected = 0
        nblocks_removed = 0
        for block in self.blocks:
            new_block = []
            for line in block:
                a, b, score = line
                pair = (a, b)
                if pair not in accepted:
                    nremoved += 1
                    continue
                av = accepted[pair]
                if score != av and score != av + "L":
                    score = av
                    ncorrected += 1
                new_block.append((a, b, score))
            if new_block:
                new_blocks.append(new_block)
            else:
                nblocks_removed += 1

        logger.debug("Removed %d existing anchors", nremoved)
        if nblocks_removed:
            logger.debug("Removed %d empty blocks", nblocks_removed)
        logger.debug("Corrected scores for %d anchors", ncorrected)
        self.blocks = new_blocks

    def print_to_file(self, filename="stdout"):
        """
        Print the anchors to a file, optionally filtering them based on the
        accepted pairs.
        """
        fw = must_open(filename, "w")
        for block in self.blocks:
            print("###", file=fw)
            for line in block:
                a, b, score = line
                print("\t".join((a, b, score)), file=fw)
        fw.close()

        logger.debug("Anchors written to `%s`", filename)

    def blast(self, blastfile=None, outfile=None):
        """
        convert anchor file to 12 col blast file
        """
        from ..formats.blast import BlastSlow, BlastLineByConversion

        if not outfile:
            outfile = self.filename + ".blast"

        if blastfile is not None:
            blasts = BlastSlow(blastfile).to_dict()
        else:
            blasts = None

        fw = must_open(outfile, "w", checkexists=True)
        nlines = 0
        for a, b, _ in self.iter_pairs():
            if (a, b) in blasts:
                bline = blasts[(a, b)]
            elif (b, a) in blasts:
                bline = blasts[(b, a)]
            else:
                line = "\t".join((a, b))
                bline = BlastLineByConversion(line, mode="110000000000")

            print(bline, file=fw)
            nlines += 1
        fw.close()

        logger.debug("A total of %d BLAST lines written to `%s`", nlines, outfile)

        return outfile

    @property
    def is_empty(self):
        blocks = self.blocks
        return not blocks or not blocks[0]


def get_best_pair(qs, ss, ts):
    pairs = {}
    for q, s, t in zip(qs, ss, ts):
        t = int(t[:-1]) if t[-1] == "L" else int(t)
        if q not in pairs or pairs[q][1] < t:
            pairs[q] = (s, t)

    # Discard score
    spairs = dict((q, s) for q, (s, t) in pairs.items())
    return spairs


def make_range(q, s, t, i, order, block_pairs, clip=10):
    pairs = get_best_pair(q, s, t)
    score = len(pairs)
    block_pairs[i].update(pairs)

    q = [order[x][0] for x in q]
    q.sort()
    qmin = q[0]
    qmax = q[-1]
    if qmax - qmin >= 2 * clip:
        qmin += clip / 2
        qmax -= clip / 2

    return Range("0", qmin, qmax, score=score, id=i)
