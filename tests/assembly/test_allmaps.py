import pytest
import os.path as op


def test_weights():
    from jcvi.assembly.allmaps import Weights

    weightsfile = op.join(op.dirname(__file__), "weights.txt")
    weights = Weights(weightsfile, ["JMMale", "JMFemale"])
    assert weights.maps == ["JMMale", "JMFemale"]

    weights = Weights(weightsfile, ["JMMale", "b"])
    assert weights.maps == ["JMMale", "JMFemale"]
