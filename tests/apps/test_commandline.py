"""
Command lines built by jcvi.apps.commandline.

Expected strings were produced by the equivalent Biopython 1.85 wrappers
(Bio.Emboss/Phylo/Align.Applications) before Biopython removed them.
"""

import pytest

from jcvi.apps.commandline import (
    ClustalwCommandline,
    FConsenseCommandline,
    FDNADistCommandline,
    FNeighborCommandline,
    FSeqBootCommandline,
    MuscleCommandline,
    NeedleCommandline,
    PhymlCommandline,
    RaxmlCommandline,
)


@pytest.mark.parametrize(
    "cline,expected",
    [
        (
            FSeqBootCommandline(
                "fseqboot",
                sequence="w/aln.phy",
                outfile="w/aln.fseqboot",
                seqtype="d",
                reps=100,
                seed=12345,
            ),
            "fseqboot -outfile=w/aln.fseqboot -sequence=w/aln.phy -seqtype=d -reps=100 -seed=12345",
        ),
        (
            FDNADistCommandline(
                "fdnadist",
                sequence="w/aln.fseqboot",
                outfile="w/aln.fdnadist",
                method="f",
            ),
            "fdnadist -outfile=w/aln.fdnadist -sequence=w/aln.fseqboot -method=f",
        ),
        (
            FNeighborCommandline(
                "fneighbor",
                datafile="w/aln.fdnadist",
                outfile="w/aln.fneighbor",
                outtreefile="w/aln.njtree",
            ),
            "fneighbor -outfile=w/aln.fneighbor -datafile=w/aln.fdnadist -outtreefile=w/aln.njtree",
        ),
        (
            FConsenseCommandline(
                "fconsense",
                intreefile="w/aln.njtree",
                outfile="w/aln.fconsense",
                outtreefile="w/aln.ct",
            ),
            "fconsense -outfile=w/aln.fconsense -intreefile=w/aln.njtree -outtreefile=w/aln.ct",
        ),
        (
            NeedleCommandline(
                asequence="a.fa",
                bsequence="b.fa",
                gapopen=10,
                gapextend=0.5,
                outfile="ab.needle",
            ),
            "needle -outfile=ab.needle -asequence=a.fa -bsequence=b.fa -gapopen=10 -gapextend=0.5",
        ),
        (PhymlCommandline(cmd="phyml", input="w/aln.phy"), "phyml -i w/aln.phy"),
        (
            PhymlCommandline(
                cmd="phyml",
                input="w/aln.phy",
                datatype="nt",
                bootstrap=100,
                model="GTR",
            ),
            "phyml -i w/aln.phy -d nt -b 100 -m GTR",
        ),
        (
            RaxmlCommandline(
                cmd="raxmlHPC",
                sequences="w/aln.phy",
                algorithm="a",
                model="GTRGAMMA",
                parsimony_seed=12345,
                rapid_bootstrap_seed=12345,
                num_replicates=100,
                name="aln",
                working_dir="/abs/raxml_work",
            ),
            "raxmlHPC -f a -m GTRGAMMA -n aln -p 12345 -s w/aln.phy -w /abs/raxml_work -x 12345 -N 100",
        ),
        (
            RaxmlCommandline(
                cmd="raxmlHPC",
                sequences="w/aln.phy",
                algorithm="h",
                model="GTRGAMMA",
                name="SH",
                starting_tree="ref.dnd",
                bipartition_filename="q.dnd",
                working_dir="/abs/raxml_work",
            ),
            "raxmlHPC -f h -m GTRGAMMA -n SH -p 10000 -s w/aln.phy -t ref.dnd -w /abs/raxml_work -z q.dnd",
        ),
        (
            ClustalwCommandline(
                cmd="clustalw2",
                infile="w/p.fasta",
                outfile="w/p.aln",
                outorder="INPUT",
                type="PROTEIN",
            ),
            "clustalw2 -infile=w/p.fasta -type=PROTEIN -outfile=w/p.aln -outorder=INPUT",
        ),
        (
            MuscleCommandline(
                cmd="muscle",
                input="w/p.fasta",
                out="w/p.aln",
                seqtype="protein",
                clwstrict=True,
            ),
            "muscle -in w/p.fasta -out w/p.aln -seqtype protein -clwstrict",
        ),
    ],
)
def test_command_line_matches_biopython(cline, expected):
    assert str(cline) == expected


def test_attributes_and_aliases():
    cl = MuscleCommandline(input="in.fa", out="out.aln")
    assert cl.input == cl.__getattr__("in") == "in.fa"
    assert cl.out == "out.aln"
    assert cl.seqtype is None
    cl.seqtype = "protein"
    assert str(cl) == "muscle -in in.fa -out out.aln -seqtype protein"
    assert ClustalwCommandline(INFILE="x.fa").infile == "x.fa"


def test_unknown_parameter_raises():
    with pytest.raises(ValueError, match="no parameter 'bogus'"):
        NeedleCommandline(bogus=1)


def test_values_are_shell_quoted():
    cl = FDNADistCommandline(sequence="my aln.phy")
    assert str(cl) == "fdnadist -sequence='my aln.phy'"


def test_call_returns_output():
    stdout, stderr = PhymlCommandline(cmd="echo", input="hello")()
    assert stdout.strip() == "-i hello"
    assert stderr == ""
