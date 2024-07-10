import os
import os.path as op

from jcvi.apps.base import cleanup
from jcvi.graphics.landscape import stack


def test_stack():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    cleanup("TAIR10_organelles.fas.pdf")
    image_name = stack(
        ["TAIR10_organelles.fas.gz", "--stacks=exons", "--window=25000", "--shift=5000"]
    )
    assert op.exists(image_name)
    os.chdir(cwd)
