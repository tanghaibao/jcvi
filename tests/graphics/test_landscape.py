import os
import os.path as op

from jcvi.apps.base import cleanup
from jcvi.graphics.landscape import depth, stack


def test_depth():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    cleanup("depth.pdf")
    image_name = depth(
        [
            "VAR0_srtd.wgs.regions.bed.gz",
            "VAR1_srtd.wgs.regions.bed.gz",
            "VAR2_srtd.wgs.regions.bed.gz",
            "--chrinfo=chrinfo.txt",
            "--titleinfo=titleinfo.txt",
            "--figsize=11x7",
        ]
    )
    assert op.exists(image_name)
    os.chdir(cwd)


def test_stack():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    cleanup("TAIR10_organelles.fas.pdf")
    image_name = stack(
        ["TAIR10_organelles.fas.gz", "--stacks=exons", "--window=25000", "--shift=5000"]
    )
    assert op.exists(image_name)
    os.chdir(cwd)
