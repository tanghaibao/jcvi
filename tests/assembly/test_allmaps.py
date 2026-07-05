import os
import os.path as op
from jcvi.apps.base import mkdir
from jcvi.assembly.allmaps import path


def datafile(path: str) -> str:
    """Returns the path to a file in the designated test directory."""
    return op.join(op.dirname(__file__), "allmaps.py", path)


def test_weights():
    from jcvi.assembly.allmaps import Weights

    weightsfile = datafile("inputs/weights.txt")
    weights = Weights(weightsfile, ["JMMale", "JMFemale"])
    assert weights.maps == ["JMMale", "JMFemale"]

    weights = Weights(weightsfile, ["JMMale", "b"])
    assert weights.maps == ["JMMale", "JMFemale"]


def test_liftover(tmp_path, monkeypatch):
    from jcvi.assembly.allmaps import liftover
    from ..config import compare_line_by_line

    monkeypatch.chdir(tmp_path)
    chainfile = datafile("inputs/JM-2.chain")
    bedfile = datafile("inputs/JM-2.bed")
    liftedbedfile = "JM-2.lifted.bed"
    expected = datafile("references/JM-2.lifted.bed")
    liftover(chainfile, bedfile, liftedbedfile, unmapfile="unmapped", cstyle="l")
    compare_line_by_line(liftedbedfile, expected)


def test_path(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    bedfile = datafile("inputs/JM-2.bed")
    fastafile = datafile("inputs/scaffolds.fasta.gz")
    weightsfile = datafile("inputs/weights.txt")
    output_image = "chr23.pdf"
    testdir = str(tmp_path / "chr23")
    mkdir(testdir)
    os.chdir(testdir)
    path(
        [
            bedfile,
            fastafile,
            "-w",
            weightsfile,
        ]
    )
    assert op.exists(output_image)
