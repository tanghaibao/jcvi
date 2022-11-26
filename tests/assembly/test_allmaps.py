import os.path as op


def test_weights():
    from jcvi.assembly.allmaps import Weights

    weightsfile = op.join(op.dirname(__file__), "weights.txt")
    weights = Weights(weightsfile, ["JMMale", "JMFemale"])
    assert weights.maps == ["JMMale", "JMFemale"]

    weights = Weights(weightsfile, ["JMMale", "b"])
    assert weights.maps == ["JMMale", "JMFemale"]


def test_liftover():
    from jcvi.apps.base import cleanup
    from jcvi.assembly.allmaps import liftover
    from ..config import compare_line_by_line

    chainfile = op.join(op.dirname(__file__), "allmaps.py", "inputs", "JM-2.chain")
    bedfile = op.join(op.dirname(__file__), "allmaps.py", "inputs", "JM-2.bed")
    liftedbedfile = "JM-2.lifted.bed"
    expected = op.join(
        op.dirname(__file__), "allmaps.py", "references", "JM-2.lifted.bed"
    )
    liftover(chainfile, bedfile, liftedbedfile, unmapfile="unmapped", cstyle="l")
    compare_line_by_line(liftedbedfile, expected)
    cleanup(liftedbedfile, "unmapped")
