#!/usr/bin/env python
# -*- coding: UTF-8 -*-

"""
Train ab initio gene predictors.
"""

import os
import sys
import logging

from optparse import OptionParser

from jcvi.apps.base import ActionDispatcher, debug, mkdir, sh, need_update
debug()


def main():

    actions = (
        ('pasa', 'extract pasa training models'),
        ('snap', 'train snap model'),
        ('augustus', 'train augustus model'),
            )
    p = ActionDispatcher(actions)
    p.dispatch(globals())


def pasa(args):
    """
    %prog ${pasadb}.assemblies.fasta ${pasadb}.pasa_assemblies.gff3

    Wraps `pasa_asmbls_to_training_set.dbi`.
    """
    from jcvi.formats.base import SetFile
    from jcvi.formats.gff import Gff

    p = OptionParser(pasa.__doc__)
    p.add_option("--pasa_home", default="~/htang/export/PASA2-r20130605p1",
                 help="Home directory for PASA [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 2:
        sys.exit(not p.print_help())

    fastafile, gffile = args
    transcodergff = fastafile + ".transdecoder.gff3"
    transcodergenomegff = fastafile + ".transdecoder.genome.gff3"
    if need_update((fastafile, gffile), (transcodergff, transcodergenomegff)):
        cmd = "{0}/scripts/pasa_asmbls_to_training_set.dbi".format(opts.pasa_home)
        cmd += " --pasa_transcripts_fasta {0} --pasa_transcripts_gff3 {1}".\
                format(fastafile, gffile)
        sh(cmd)

    completeids = fastafile.rsplit(".", 1)[0] + ".complete.ids"
    if need_update(transcodergff, completeids):
        cmd = "grep complete {0} | cut -f1 | sort -u".format(transcodergff)
        sh(cmd, outfile=completeids)

    complete = SetFile(completeids)
    seen = set()
    completegff = transcodergenomegff.rsplit(".", 1)[0] + ".complete.gff3"
    fw = open(completegff, "w")
    gff = Gff(transcodergenomegff)
    for g in gff:
        a = g.attributes
        if "Parent" in a:
            id = a["Parent"][0]
        else:
            id = a["ID"][0]
        asmbl_id = id.split("|")[0]
        if asmbl_id not in complete:
            continue
        print >> fw, g
        if g.type == "gene":
            seen.add(id)

    fw.close()
    logging.debug("A total of {0} complete models extracted to `{1}`.".\
                    format(len(seen), completegff))


def snap(args):
    """
    %prog species gffile fastafile

    Train SNAP model given gffile and fastafile. Whole procedure taken from:
    <http://gmod.org/wiki/MAKER_Tutorial_2012>
    """
    p = OptionParser(snap.__doc__)
    p.add_option("--maker_home", default="~/htang/export/maker",
                 help="Home directory for MAKER [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    species, gffile, fastafile = args
    mhome = opts.maker_home
    snapdir = "snap"
    mkdir(snapdir)

    cwd = os.getcwd()
    os.chdir(snapdir)

    newgffile = "training.gff3"
    logging.debug("Construct GFF file combined with sequence ...")
    sh("cat ../{0} > {1}".format(gffile, newgffile))
    sh('echo "##FASTA" >> {0}'.format(newgffile))
    sh("cat ../{0} >> {1}".format(fastafile, newgffile))

    logging.debug("Make models ...")
    sh("{0}/bin/maker2zff training.gff3".format(mhome))
    sh("{0}/exe/snap/fathom -categorize 1000 genome.ann genome.dna".format(mhome))
    sh("{0}/exe/snap/fathom -export 1000 -plus uni.ann uni.dna".format(mhome))
    sh("{0}/exe/snap/forge export.ann export.dna".format(mhome))
    sh("{0}/exe/snap/hmm-assembler.pl {1} . > {1}.hmm".format(mhome, species))

    os.chdir(cwd)
    logging.debug("SNAP matrix written to `{0}/{1}.hmm`".format(snapdir, species))


def augustus(args):
    """
    %prog species gffile fastafile

    Train AUGUSTUS model given gffile and fastafile. Whole procedure taken from:
    <http://www.molecularevolution.org/molevolfiles/exercises/augustus/training.html>
    """
    p = OptionParser(snap.__doc__)
    p.add_option("--augustus_home", default="~/htang/export/augustus.2.5.5",
                 help="Home directory for AUGUSTUS [default: %default]")
    opts, args = p.parse_args(args)

    if len(args) != 3:
        sys.exit(not p.print_help())

    species, gffile, fastafile = args
    mhome = opts.augustus_home
    augdir = "augustus"
    mkdir(augdir)
    os.chdir(augdir)

    sh("{0}/scripts/new_species.pl --species={1}".format(mhome, species))
    sh("{0}/scripts/gff2gbSmallDNA.pl ../{1} ../{2} 1000 raw.gb".\
            format(mhome, gffile, fastafile))
    sh("{0}/bin/etraining --species={1} raw.gb 2> train.err".\
            format(mhome, species))
    sh("cat train.err | perl -pe 's/.*in sequence (\S+): .*/$1/' > badgenes.lst")
    sh("{0}/scripts/filterGenes.pl badgenes.lst raw.gb > training.gb".\
            format(mhome))
    sh("grep -c LOCUS raw.gb training.gb")
    sh("{0}/scripts/randomSplit.pl training.gb 100".format(mhome))
    sh("{0}/bin/etraining --species={1} training.gb.train".format(mhome, species))

    msg = """Now ready to run:
    $ {0}/scripts/optimize_augustus.pl --species={1} training.gb.train
    """.format(mhome, species)
    logging.debug(msg)


if __name__ == '__main__':
    main()
