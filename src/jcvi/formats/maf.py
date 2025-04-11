#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
MAF format specification:
<http://genome.ucsc.edu/FAQ/FAQformat#format5>
"""
from bisect import bisect
from dataclasses import dataclass
import sys

from Bio import AlignIO, SeqIO
from bx import interval_index_file
from bx.align import maf

from ..apps.base import ActionDispatcher, OptionParser, need_update
from ..apps.lastz import blastz_score_to_ncbi_bits, blastz_score_to_ncbi_expectation
from .base import BaseFile, logger

FLANK = 60


class Maf(BaseFile, dict):
    def __init__(self, filename, index=False):
        super().__init__(filename)

        indexfile = filename + ".idx"
        if index:
            if need_update(filename, indexfile):
                self.build_index(filename, indexfile)

            self.index = maf.Index(filename, indexfile)

        fp = open(filename)
        self.reader = maf.Reader(fp)

    def build_index(self, filename, indexfile):
        """
        Recipe from Brad Chapman's blog
        <http://bcbio.wordpress.com/2009/07/26/sorting-genomic-alignments-using-python/>
        """
        indexes = interval_index_file.Indexes()
        in_handle = open(filename)

        reader = maf.Reader(in_handle)
        while True:
            pos = reader.file.tell()
            rec = next(reader)
            if rec is None:
                break
            for c in rec.components:
                indexes.add(
                    c.src,
                    c.forward_strand_start,
                    c.forward_strand_end,
                    pos,
                    max=c.src_size,
                )

        index_handle = open(indexfile, "w")
        indexes.write(index_handle)
        index_handle.close()


@dataclass
class Breakpoint:
    arec: str
    astart: int
    brec: str
    bstart: int

    def __str__(self):
        return f"{self.arec}:{self.astart}-{self.brec}:{self.bstart}"


def main():

    actions = (
        ("bed", "convert MAF to BED format"),
        ("blast", "convert MAF to BLAST tabular format"),
        ("breakpoints", "find breakpoints in MAF and 'simulate' chimeric contigs"),
    )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def breakpoints(args):
    """
    %prog breakpoints A.B.maf A.fa B.fa AB 1000000 2000000

    Find breakpoints in MAF and 'simulate' chimeric contigs in `AB.fa`.
    Breakpoints are 'roughly' around the user defined positions. The idea is
    to simulate chimeric contigs, which are useful for testing algorithms,
    e.g. klassify.
    """
    p = OptionParser(breakpoints.__doc__)
    p.add_argument(
        "--minsize",
        default=10000,
        type=int,
        help="Minimum size of alignment to consider",
    )
    opts, args = p.parse_args(args)

    if len(args) not in (5, 6):
        sys.exit(not p.print_help())

    maf_file, a_fasta, b_fasta, ab = args[:4]
    bps = sorted(int(x) for x in args[4:])
    minsize = opts.minsize

    filtered_msa = []
    for msa in AlignIO.parse(maf_file, "maf"):
        arec, brec = msa
        if brec.annotations["size"] < minsize:
            continue
        filtered_msa.append((brec.annotations["start"], arec, brec))
    logger.info("Total alignments: %d", len(filtered_msa))

    final = []
    # Load the sequences
    ar = next(SeqIO.parse(a_fasta, "fasta"))
    br = next(SeqIO.parse(b_fasta, "fasta"))
    for bp in bps:
        i = bisect(filtered_msa, (bp,))
        _, arec, brec = filtered_msa[i]
        logger.info("%s", arec)
        logger.info("%s", brec)
        assert len(arec) == len(brec)
        # Find the midpoint, safe to crossover there
        midpoint = len(arec) // 2
        aseq = arec.seq[:midpoint]
        astart = arec.annotations["start"] + len(aseq) - aseq.count("-")
        logger.info("%s|%s", aseq[-FLANK:], arec.seq[midpoint:][:FLANK])
        bseq = brec.seq[:midpoint]
        bstart = brec.annotations["start"] + len(bseq) - bseq.count("-")
        logger.info("%s|%s", bseq[-FLANK:], brec.seq[midpoint:][:FLANK])
        bpt = Breakpoint(arec.id, astart, brec.id, bstart)
        logger.info("-" * FLANK * 2 + ">")
        logger.info("%s|%s", ar.seq[:astart][-FLANK:], br.seq[bstart:][:FLANK])
        final.append(bpt)

    logger.info("Breakpoints found: %s", final)
    if len(final) == 2:
        bp1, bp2 = final[:2]
        # ====-------=======
        #   bp1     bp2
        abseq = (
            ar.seq[: bp1.astart]
            + br.seq[bp1.bstart : bp2.bstart]
            + ar.seq[bp2.astart :]
        )
    elif len(final) == 1:
        bp = final[0]
        abseq = ar.seq[: bp.astart] + br.seq[bp.bstart :]
    abrec = SeqIO.SeqRecord(abseq, id=ab, description="")
    ab_fasta = f"{ab}.fa"
    SeqIO.write([abrec], ab_fasta, "fasta")
    logger.info("Writing to %s", ab_fasta)


def bed(args):
    """
    %prog bed maffiles > out.bed

    Convert a folder of maf alignments to the bed features
    then useful to check coverage, etc.
    """
    p = OptionParser(bed.__doc__)
    _, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(p.print_help())

    flist = args
    prefix = flist[0].split(".")[0]

    j = 0
    for f in flist:
        reader = Maf(f).reader
        for rec in reader:
            a, b = rec.components

            for a, tag in zip((a, b), "ab"):
                name = "{0}_{1:07d}{2}".format(prefix, j, tag)
                print(
                    "\t".join(
                        str(x)
                        for x in (
                            a.src,
                            a.forward_strand_start,
                            a.forward_strand_end,
                            name,
                        )
                    )
                )

            j += 1


def alignment_details(a, b):
    nmatch = 0
    nmismatch = 0
    ngaps = 0

    assert len(a) == len(b)
    l = len(a)

    for i in range(l):
        if a[i] == b[i]:
            nmatch += 1
        elif a[i] == "-" or b[i] == "-":
            ngaps += 1
        else:
            nmismatch += 1

    pctid = 100.0 * nmatch / l
    return pctid, nmismatch, ngaps


def maf_to_blast8(f):
    """
    Convert a MAF file to BLAST tabular format.
    """
    reader = Maf(f).reader
    for rec in reader:
        a, b = rec.components
        query = a.src
        subject = b.src
        qstart = a.forward_strand_start
        qstop = a.forward_strand_end
        sstart = b.forward_strand_start
        sstop = b.forward_strand_end
        score = rec.score

        evalue = blastz_score_to_ncbi_expectation(score)
        score = blastz_score_to_ncbi_bits(score)
        evalue, score = "{0:.2g}".format(evalue), "{0:.1f}".format(score)
        hitlen = len(a.text)

        pctid, nmismatch, ngaps = alignment_details(a.text, b.text)
        print(
            "\t".join(
                str(x)
                for x in (
                    query,
                    subject,
                    pctid,
                    hitlen,
                    nmismatch,
                    ngaps,
                    qstart,
                    qstop,
                    sstart,
                    sstop,
                    evalue,
                    score,
                )
            )
        )


def blast(args):
    """
    %prog blast maffiles > out.blast

    From a folder of .maf files, generate .blast file with tabular format.
    """
    p = OptionParser(blast.__doc__)
    _, args = p.parse_args(args)

    if len(args) == 0:
        sys.exit(p.print_help())

    flist = args

    for f in flist:
        maf_to_blast8(f)


if __name__ == "__main__":
    main()
