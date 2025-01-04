import os
import os.path as op

from jcvi.apps.base import cleanup
from jcvi.graphics.grabseeds import calibrate, seeds


def test_main():
    cwd = os.getcwd()
    os.chdir(op.join(op.dirname(__file__), "data"))
    output_image = "test.pdf"
    json_file = "calibrate.json"
    cleanup(output_image)
    seeds(["test.JPG"])
    assert op.exists(output_image)

    # Test calibrate
    json_file = "calibrate.json"
    cleanup(json_file, output_image)
    json_file = calibrate(["calibrate.JPG", "1"])  # `1` for the boxsize arg
    assert op.exists(json_file)
    seeds(["test.JPG", "--calibrate", json_file])

    os.chdir(cwd)
