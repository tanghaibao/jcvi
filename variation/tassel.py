#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Run TASSEL GBS module, following manual here:
http://www.maizegenetics.net/tassel/docs/TasselPipelineGBS.pdf
"""

import os.path as op
import sys

from jcvi.formats.base import write_file
from jcvi.apps.base import OptionParser, ActionDispatcher, mkdir, get_abs_path


def main():

    actions = (
        ('prepare', 'prepare TASSEL pipeline'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def run_pipeline(thome, plugin, options, log=False):
    cmd = op.join(thome, "run_pipeline.pl")
    cmd += " -Xmx16g -fork1 -{0} {1} -endPlugin -runfork1".format(plugin, options)
    if log:
        cmd += " | tee {0}.log".format(plugin)
    return cmd


def prepare(args):
    """
    %prog prepare barcode_key.csv reference.fasta

    Prepare TASSEL pipeline.
    """
    valid_enzymes = "ApeKI|ApoI|BamHI|EcoT22I|HinP1I|HpaII|MseI|MspI|" \
                    "NdeI|PasI|PstI|Sau3AI|SbfI|AsiSI-MspI|BssHII-MspI|" \
                    "FseI-MspI|PaeR7I-HhaI|PstI-ApeKI|PstI-EcoT22I|PstI-MspI" \
                    "PstI-TaqI|SalI-MspI|SbfI-MspI".split("|")
    p = OptionParser(prepare.__doc__)
    p.add_option("--enzyme", default="ApeKI", choices=valid_enzymes,
                 help="Restriction enzyme used [default: %default]")
    p.set_home("tassel")
    p.set_aligner(aligner="bwa")
    p.set_cpus()
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    barcode, reference = args
    thome = opts.tassel_home
    reference = get_abs_path(reference)
    folders = ("fastq", "tagCounts", "mergedTagCounts", "topm",
               "tbt", "mergedTBT", "hapmap", "hapmap/raw",
               "hapmap/mergedSNPs", "hapmap/filt", "hapmap/bpec")
    for f in folders:
        mkdir(f)

    # Build the pipeline
    runsh = []
    o = "-i fastq -k {0} -e {1} -o tagCounts".format(barcode, opts.enzyme)
    cmd = run_pipeline(thome, "FastqToTagCountPlugin", o)
    runsh.append(cmd)

    o = "-i tagCounts -o mergedTagCounts/myMasterTags.cnt"
    o += " -c 5 -t mergedTagCounts/myMasterTags.cnt.fq"
    cmd = run_pipeline(thome, "MergeMultipleTagCountPlugin", o)
    runsh.append(cmd)
    runsh.append("cd mergedTagCounts")

    cmd = "python -m jcvi.apps.{0} align --cpus {1}".\
                format(opts.aligner, opts.cpus)
    cmd += " {0} myMasterTags.cnt.fq".format(reference)
    runsh.append(cmd)
    runsh.append("cd ..")

    o = "-i mergedTagCounts/*.sam -o topm/myMasterTags.topm"
    cmd = run_pipeline(thome, "SAMConverterPlugin", o)
    runsh.append(cmd)

    o = "-i mergedTBT/myStudy.tbt.byte -y -m topm/myMasterTags.topm"
    o += " -mUpd topm/myMasterTagsWithVariants.topm"
    o += " -o hapmap/raw/myGBSGenos_chr+.hmp.txt"
    o += " -mnF 0.8 -p myPedigreeFile.ped -mnMAF 0.02 -mnMAC 100000"
    o += " -ref {0} -sC 1 -eC 10".format(reference)
    cmd = run_pipeline(thome, "TagsToSNPByAlignmentPlugin", o)
    runsh.append(cmd)

    o = "-hmp hapmap/raw/myGBSGenos_chr+.hmp.txt"
    o += " -o hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.hmp.txt"
    o += " -misMat 0.1 -p myPedigreeFile.ped -callHets -sC 1 -eC 10"
    cmd = run_pipeline(thome, "MergeDuplicateSNPsPlugin", o)
    runsh.append(cmd)

    o = "-hmp hapmap/mergedSNPs/myGBSGenos_mergedSNPs_chr+.hmp.txt"
    o += " -o hapmap/filt/myGBSGenos_mergedSNPsFilt_chr+.hmp.txt"
    o += " -mnTCov 0.01 -mnSCov 0.2 -mnMAF 0.01 -sC 1 -eC 10"
    #o += "-hLD -mnR2 0.2 -mnBonP 0.005"
    cmd = run_pipeline(thome, "GBSHapMapFiltersPlugin", o)
    runsh.append(cmd)

    runfile = "run.sh"
    write_file(runfile, "\n".join(runsh))


if __name__ == '__main__':
    main()
