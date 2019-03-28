#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Script to write and assist SOAPdenovo assembly.
"""
from __future__ import print_function

import os.path as op
import sys

from jcvi.formats.fastq import guessoffset, readlen, is_fastq
from jcvi.assembly.base import FastqNamings, Library, get_libs
from jcvi.apps.base import OptionParser, ActionDispatcher, need_update, sh


class FillLine (object):

    def __init__(self, row):
        args = row.split()
        self.start = int(args[0])
        self.end = int(args[1])
        self.leftextend = int(args[2])
        self.rightextend = int(args[3])
        self.closed = (int(args[4]) == 1)
        self.extendlength = int(args[5])
        self.before = int(args[6])
        self.after = int(args[7])
        # Convert from unsigned to signed
        # <http://stackoverflow.com/questions/1375897/how-to-get-the-signed-integer-value-of-a-long-in-python>
        if self.after > 0 and (self.after & 0x80000000):
            self.after += -0x100000000

    @property
    def delta(self):
        return self.after - self.before


def main():

    actions = (
        ('clean', 'clean and dedup paired FASTQ files'),
        ('correct', 'correct reads using ErrorCorrection'),
        ('prepare', 'prepare SOAP config files and run script'),
        ('fillstats', 'build stats on .fill file from GapCloser'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


SOAPHEADER = """
P={0}
K={1}
S=soap.config
G=soap.gc.config
C={2}
A=asm$K
"""

GCRUN = "GapCloser_v1.12 -a ${A}.scafSeq -b $G -l 155 -o ${A}.closed.scafSeq -p 31 -t $P"
GCRUNG = "GapCloser_v1.12 -a {0} -b $G -l 155 -o {1} -p 31 -t $P"

SOAPRUN = """
$C pregraph -s $S -d 1 -K $K -o $A -R -p $P
$C contig -s $S -g $A -M 1 -R -p $P
$C map -s $S -g $A -p $P
$C scaff -g $A -F -p $P
""" + GCRUN

SCFRUN = """
prepare -K $K -c %s -g $A
$C map -s $S -g $A -p $P
$C scaff -z -g $A -F -p $P
""" + GCRUN


def get_size(filename):

    library_name = lambda x: "-".join(\
                op.basename(x).split(".")[0].split("-")[:2])

    lib = Library(library_name(filename))
    return lib.size


def correct(args):
    """
    %prog correct *.fastq

    Correct reads using ErrorCorrection. Only PE will be used to build the K-mer
    table.
    """
    p = OptionParser(correct.__doc__)
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    lstfile = "reads2cor.lst"
    fw = open(lstfile, "w")
    print("\n".join(x for x in args if x[:2] == "PE"), file=fw)
    fw.close()

    p1 = args[0]
    offset = guessoffset([p1])
    cpus = opts.cpus

    freq = "output.freq.cz"
    freqlen = freq + ".len"
    if need_update(args, (freq, freqlen)):
        cmd = "KmerFreq_AR_v2.0 -k 17 -c -1 -q {0}".format(offset)
        cmd += " -m 1 -t {0}".format(cpus)
        cmd += " -p output {0}".format(lstfile)
        sh(cmd)

    fw = open(lstfile, "w")
    print("\n".join(args), file=fw)
    fw.close()

    cmd = "Corrector_AR_v2.0 -k 17 -l 3 -m 5 -c 5 -a 0 -e 1 -w 0 -r 45"
    cmd += " -Q {0} -q 30 -x 8 -t {1} -o 1 ".format(offset, cpus)
    cmd += " {0} {1} {2}".format(freq, freqlen, lstfile)
    sh(cmd)


def clean(args):
    """
    %prog clean 1.fastq 2.fastq [insertsize]

    Clean and dedup paired FASTQ files.
    """
    p = OptionParser(clean.__doc__)
    p.add_option("-a", default=0, type="int",
                 help="Trim length at 5' end [default: %default]")
    p.add_option("-b", default=50, type="int",
                 help="Trim length at 3' end [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) == 2:
        p1, p2 = args
        size = get_size(p1)
    elif len(args) == 3:
        p1, p2, size = args
        size = int(size)
    else:
        sys.exit(not p.print_help())

    pf = p1.split(".")[0]
    cpus = opts.cpus

    offset = guessoffset([p1])
    a, b = opts.a, opts.b

    p1_clean = p1 + ".clean"
    p1_cleangz = p1_clean + ".gz"
    p2_clean = p2 + ".clean"
    p2_cleangz = p2_clean + ".gz"
    if need_update([p1, p2], [p1_cleangz, p2_cleangz]):
        cmd = "SOAPfilter_v2.0 -t {0} -m 2000000 -p -y -z -g".format(cpus)
        cmd += " -q {0} -w 10 -B 50 -f 0".format(offset)
        cmd += " -l {0} -a {1} -b {2} -c {1} -d {2}".format(size, a, b, a, b)
        cmd += " {0} {1} {2}.clean.stat {3} {4}".\
                    format(p1, p2, pf, p1_clean, p2_clean)
        sh(cmd)


def fillstats(args):
    """
    %prog fillstats genome.fill

    Build stats on .fill file from GapCloser.
    """
    from jcvi.utils.cbook import SummaryStats, percentage, thousands

    p = OptionParser(fillstats.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fillfile, = args
    fp = open(fillfile)
    scaffolds = 0
    gaps = []
    for row in fp:
        if row[0] == ">":
            scaffolds += 1
            continue
        fl = FillLine(row)
        gaps.append(fl)

    print("{0} scaffolds in total".format(scaffolds), file=sys.stderr)

    closed = [x for x in gaps if x.closed]
    closedbp = sum(x.before for x in closed)
    notClosed = [x for x in gaps if not x.closed]
    notClosedbp = sum(x.before for x in notClosed)

    totalgaps = len(closed) + len(notClosed)

    print("Closed gaps: {0} size: {1} bp".\
                        format(percentage(len(closed), totalgaps), thousands(closedbp)), file=sys.stderr)
    ss = SummaryStats([x.after for x in closed])
    print(ss, file=sys.stderr)

    ss = SummaryStats([x.delta for x in closed])
    print("Delta:", ss, file=sys.stderr)

    print("Remaining gaps: {0} size: {1} bp".\
                        format(percentage(len(notClosed), totalgaps), thousands(notClosedbp)), file=sys.stderr)
    ss = SummaryStats([x.after for x in notClosed])
    print(ss, file=sys.stderr)


def prepare(args):
    """
    %prog prepare *.fastq

    Scan input fastq files (see below) and write SOAP config files based
    on inputfiles. Use "--scaffold contigs.fasta" to perform scaffolding.
    """
    from jcvi.formats.base import write_file

    p = OptionParser(prepare.__doc__ + FastqNamings)
    p.add_option("-K", default=45, type="int",
                 help="K-mer size [default: %default]")
    p.add_option("--assemble_1st_rank_only", default=False, action="store_true",
                 help="Assemble the first rank only, other libs asm_flags=2 [default: %default]")
    p.add_option("--scaffold",
                 help="Only perform scaffolding [default: %default]")
    p.add_option("--gapclose",
                 help="Only perform gap closure [default: %default]")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) < 1:
        sys.exit(not p.print_help())

    fnames = args
    K = opts.K
    for x in fnames:
        assert op.exists(x), "File `{0}` not found.".format(x)

    a1st = opts.assemble_1st_rank_only

    cfgfile = "soap.config"
    gc_cfgfile = "soap.gc.config"
    fw = open(cfgfile, "w")
    fw_gc = open(gc_cfgfile, "w")

    libs = get_libs(fnames)
    rank = 0
    singletons = []
    max_rd_len = max(readlen([f]) for f in fnames)

    block = "max_rd_len={0}\n".format(max_rd_len)
    for stream in (sys.stderr, fw, fw_gc):
        print(block, file=stream)

    # Collect singletons first
    singletons = []
    for lib, fs in libs:
        if lib.size == 0:
            singletons += fs
            continue

    for lib, fs in libs:
        size = lib.size
        if size == 0:
            continue

        rank += 1
        block = "[LIB]\n"
        block += "avg_ins={0}\n".format(size)
        f = fs[0]
        block += "reverse_seq={0}\n".format(lib.reverse_seq)
        asm_flags = 2 if (rank > 1 and a1st) else lib.asm_flags
        block += "asm_flags={0}\n".format(asm_flags)
        block += "rank={0}\n".format(rank)
        if lib.reverse_seq:
            pair_num_cutoff = 3
            block += "pair_num_cutoff={0}\n".format(pair_num_cutoff)
        block += "map_len=35\n"

        for f in fs:
            if ".1." in f:
                tag = "q1"
            elif ".2." in f:
                tag = "q2"
            block += "{0}={1}\n".format(tag, f)

        if rank == 1:
            for s in singletons:
                tag = "q" if is_fastq(s) else "f"
                block += tag + "={0}\n".format(s)

        print(block, file=sys.stderr)
        print(block, file=fw)

        if asm_flags > 2:
            print(block, file=fw_gc)

    runfile = "run.sh"
    scaffold = opts.scaffold
    bb = 63 if K <= 63 else 127
    binary = "SOAPdenovo-{0}mer".format(bb)
    header = SOAPHEADER.format(opts.cpus, K, binary)
    if opts.gapclose:
        gapclose = opts.gapclose
        outfile = gapclose.rsplit(".", 1)[0] + ".closed.fasta"
        template = header + GCRUNG.format(gapclose, outfile)
    else:
        template = header + (SCFRUN % scaffold if scaffold else SOAPRUN)

    write_file(runfile, template)
    fw.close()
    fw_gc.close()


if __name__ == '__main__':
    main()
