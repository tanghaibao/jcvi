#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run PBJelly on reference FASTA and a set of patching reads. Sometimes patching
'reads' can also be contigs, which require setting the blasr criteria to be
higher in `Protocol.xml`.
"""

import os
import os.path as op
import sys
import logging

from jcvi.apps.base import OptionParser, ActionDispatcher, sh, which


class Protocol (object):

    def __init__(self, outputDir, reference, reads, highqual=False):
        self.outputDir = outputDir
        self.reference = reference
        self.reads = reads
        oblasr = (24, 98) if highqual else (12, 75)
        self.blasr = "-minMatch {0} -minPctSimilarity {1}".format(*oblasr)
        self.blasr += " -bestn 8 -maxScore -500 -nproc 64 -noSplitSubreads"

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
        print >> sys.stderr, s.toprettyxml()
        print >> fw, s.toprettyxml()
        logging.debug("XML configuration written to `{0}`".format(filename))
        fw.close()


def main():

    actions = (
        ('patch', 'run PBJelly with reference and reads'),
        ('spancount', 'count support for each gap'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
        print >> fw, f, asb, sc
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
