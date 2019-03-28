#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run PBJelly on reference FASTA and a set of patching reads. Sometimes patching
'reads' can also be contigs, which require setting the blasr criteria to be
higher in `Protocol.xml`.
"""
from __future__ import print_function

import os
import os.path as op
import sys
import logging

from collections import defaultdict

from jcvi.formats.base import must_open
from jcvi.utils.cbook import percentage
from jcvi.apps.base import OptionParser, ActionDispatcher, sh, which


class Protocol (object):

    def __init__(self, outputDir, reference, reads, highqual=False):
        self.outputDir = outputDir
        self.reference = reference
        self.reads = reads
        oblasr = (16, 98) if highqual else (8, 75)
        self.blasr = "-minMatch {0} -sdpTupleSize 8 -minPctSimilarity {1}".format(*oblasr)
        self.blasr += " -bestn 1 -nCandidates 10 -maxScore -500"
        self.blasr += " -nproc 64 -noSplitSubreads"

    def write_xml(self, filename="Protocol.xml"):
        import xml.etree.cElementTree as ET
        import xml.dom.minidom as XD

        jellyProtocol = ET.Element("jellyProtocol")
        reference = ET.SubElement(jellyProtocol, "reference")
        reference.text = self.reference
        outputDir = ET.SubElement(jellyProtocol, "outputDir")
        outputDir.text = self.outputDir
        blasr = ET.SubElement(jellyProtocol, "blasr")
        blasr.text = self.blasr
        inp = ET.SubElement(jellyProtocol, "input")
        baseDir, readsfile = op.split(self.reads)
        inp.set("baseDir", baseDir)
        job = ET.SubElement(inp, "job")
        job.text = readsfile

        s = ET.tostring(jellyProtocol)
        s = XD.parseString(s)

        fw = open(filename, "w")
        print(s.toprettyxml(), file=sys.stderr)
        print(s.toprettyxml(), file=fw)
        logging.debug("XML configuration written to `{0}`".format(filename))
        fw.close()


class M4Line (object):
    """
    See doc:

    https://github.com/mchaisso/blasr
    """
    def __init__(self, sline):
        args = sline.split()
        self.query = args[0]
        self.subject = args[1]
        self.score = int(args[2])
        self.qstrand = args[3]
        self.qstart = int(args[4])
        self.qstop = int(args[5])
        self.qlen = int(args[6])
        self.sstrand = args[7]
        self.sstart = int(args[8])
        self.sstop = int(args[9])
        self.slen = int(args[10])
        self.spaceusage = int(args[11])
        self.quality = int(args[12])


def main():

    actions = (
        ('patch', 'run PBJelly with reference and reads'),
        ('spancount', 'count support for each gap'),
        ('filterm4', 'filter .m4 file after blasr is run'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def filterm4(args):
    """
    %prog filterm4 sample.m4 > filtered.m4

    Filter .m4 file after blasr is run. As blasr takes a long time to run,
    changing -bestn is undesirable. This screens the m4 file to retain top hits.
    """
    p = OptionParser(filterm4.__doc__)
    p.add_option("--best", default=1, type="int", help="Only retain best N hits")
    p.set_outfile()
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    m4file, = args
    best = opts.best
    fp = open(m4file)
    fw = must_open(opts.outfile, "w")
    seen = defaultdict(int)
    retained = total = 0
    for row in fp:
        r = M4Line(row)
        total += 1
        if total % 100000 == 0:
            logging.debug("Retained {0} lines".\
                            format(percentage(retained, total)))
        if seen.get(r.query, 0) < best:
            fw.write(row)
            seen[r.query] += 1
            retained += 1
    fw.close()


def spancount(args):
    """
    %prog spancount list_of_fillingMetrics

    Count span support for each gap. A file with paths of all fillingMetrics can
    be built with Linux `find`.

    $ (find assembly -name "fillingMetrics.json" -print > list_of_fillMetrics 2>
    /dev/null &)
    """
    import json

    p = OptionParser(spancount.__doc__)
    opts, args = p.parse_args(args)

    if len(args) != 1:
        sys.exit(not p.print_help())

    fof, = args
    fp = open(fof)
    flist = [row.strip() for row in fp]
    spanCount = "spanCount"
    avgSpanBases = "avgSpanBases"
    fw = open(spanCount, "w")
    for f in flist:
        fp = open(f)
        j = json.load(fp)
        sc = j.get(spanCount, None)
        asb = j.get(avgSpanBases, None)
        print(f, asb, sc, file=fw)
        fw.flush()
    fw.close()


def fake_quals(fa):
    faq = fa.rsplit(".", 1)[0] + ".qual"
    if op.exists(faq):
        logging.debug("Qual file `{0}` found.".format(faq))
    else:
        sh("fakeQuals.py {0} {1}".format(fa, faq))
    return fa, faq


def patch(args):
    """
    %prog patch reference.fasta reads.fasta

    Run PBJelly with reference and reads.
    """
    from jcvi.formats.base import write_file
    from jcvi.formats.fasta import format

    p = OptionParser(patch.__doc__)
    p.add_option("--cleanfasta", default=False, action="store_true",
                 help="Clean FASTA to remove description [default: %default]")
    p.add_option("--highqual", default=False, action="store_true",
                 help="Reads are of high quality [default: %default]")
    p.set_home("pbjelly")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ref, reads = args
    cpus = opts.cpus
    cmd = op.join(opts.pbjelly_home, "setup.sh")
    setup = "source {0}".format(cmd)
    if not which("fakeQuals.py"):
        sh(setup)

    pf = ref.rsplit(".", 1)[0]
    pr, px = reads.rsplit(".", 1)
    # Remove description line
    if opts.cleanfasta:
        oref = pf + ".f.fasta"
        oreads = pr + ".f.fasta"
        format([ref, oref])
        format([reads, oreads])
        ref, reads = oref, oreads

    # Check if the FASTA has qual
    ref, refq = fake_quals(ref)
    convert_reads = not px in ("fq", "fastq", "txt")
    if convert_reads:
        reads, readsq = fake_quals(reads)
        readsfiles = " ".join((reads, readsq))
    else:
        readsfiles = reads

    # Make directory structure
    dref, dreads = "data/reference", "data/reads"
    cwd = os.getcwd()
    reference = op.join(cwd, "{0}/{1}".format(dref, ref))
    reads = op.join(cwd, "{0}/{1}".format(dreads, reads))
    if not op.exists(reference):
        sh("mkdir -p {0}".format(dref))
        sh("cp {0} {1}/".format(" ".join((ref, refq)), dref))
    if not op.exists(reads):
        sh("mkdir -p {0}".format(dreads))
        sh("cp {0} {1}/".format(readsfiles, dreads))

    outputDir = cwd
    p = Protocol(outputDir, reference, reads, highqual=opts.highqual)
    p.write_xml()

    # Build the pipeline
    runsh = [setup]
    for action in "setup|mapping|support|extraction".split("|"):
        runsh.append("Jelly.py {0} Protocol.xml".format(action))

    runsh.append('Jelly.py assembly Protocol.xml -x "--nproc={0}"'.format(cpus))
    runsh.append("Jelly.py output Protocol.xml")

    runfile = "run.sh"
    contents = "\n".join(runsh)
    write_file(runfile, contents)


if __name__ == '__main__':
    main()
