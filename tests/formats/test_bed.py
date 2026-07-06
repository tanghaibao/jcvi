import os.path as op

from jcvi.formats.bed import summary


def test_summary():
    data_dir = op.join(op.dirname(__file__), "data")
    summary([op.join(data_dir, "custom.bed")])
