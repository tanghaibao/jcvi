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

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug, sh, which
debug()


class Protocol (object):

    def __init__(self, outputDir, reference, reads, highqual=False):
        self.outputDir = outputDir
        self.reference = reference
        self.reads = reads
        oblasr = (20, 98) if highqual else (8, 70)
        self.blasr = "-minMatch {0} -minPctIdentity {1}".format(*oblasr)
        self.blasr += " -bestn 8 -nCandidates 30 -maxScore -500 -nproc 64 -noSplitSubreads"

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
        print >> fw, s.toprettyxml()
        logging.debug("XML configuration written to `{0}`".format(filename))
        fw.close()


def main():

    actions = (
        ('patch', 'run PBJelly with reference and reads'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


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
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    ref, reads = args
    pf = ref.rsplit(".", 1)[0]
    pr = reads.rsplit(".", 1)[0]
    # Remove description line
    if opts.cleanfasta:
        oref = pf + ".f.fasta"
        oreads = pr + ".f.fasta"
        format([ref, oref])
        format([reads, oreads])
        ref, reads = oref, oreads

    # Check if the FASTA has qual
    ref, refq = fake_quals(ref)
    reads, readsq = fake_quals(reads)

    # Make directory structure
    dref, dreads = "data/reference", "data/reads"
    sh("mkdir -p {0}".format(dref))
    sh("mkdir -p {0}".format(dreads))
    sh("mv {0} {1}/".format(" ".join((ref, refq)), dref))
    sh("mv {0} {1}/".format(" ".join((reads, readsq)), dreads))
    cwd = os.getcwd()

    outputDir = cwd
    reference = op.join(cwd, "{0}/{1}".format(dref, ref))
    reads = op.join(cwd, "{0}/{1}".format(dreads, reads))
    p = Protocol(outputDir, reference, reads, highqual=opts.highqual)
    p.write_xml()

    # Make sure we have the patched version of Extraction.py
    # See discussion <http://seqanswers.com/forums/showthread.php?t=27599>
    extpy = which("Extraction.py")
    lines = open(extpy).readlines()
    patchline = lines[191].strip()
    assert patchline.split()[0] == "return"

    # Build the pipeline
    runsh = []
    for action in "setup|mapping|support|extraction".split("|"):
        runsh.append("Jelly.py {0} Protocol.xml".format(action))

    pcmds = """find assembly -name "ref*" -exec echo \\
        "WrapAssembly.py {{}} {0}/{1}/{2}.gapInfo.bed \\
        > {{}}/assembly.out 2> {{}}/assembly.err" \; > commands.list""".\
        format(outputDir, dref, pf)
    runsh.append(pcmds)

    runsh.append("parallel < commands.list")
    runsh.append("Jelly.py output Protocol.xml")

    runfile = "run.sh"
    contents = "\n".join(runsh)
    write_file(runfile, contents, meta="run script")


if __name__ == '__main__':
    main()
