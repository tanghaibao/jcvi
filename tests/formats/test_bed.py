import os
import os.path as op

from jcvi.formats.bed import summary


def test_summary():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    summary(["custom.bed"])
    os.chdir(cwd)
