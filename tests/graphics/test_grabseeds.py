import os.path as op

from jcvi.graphics.grabseeds import calibrate, seeds


def test_main(copy_data):
    copy_data(op.join(op.dirname(__file__), "data"))

    output_image = "test.pdf"
    seeds(["test.JPG"])
    assert op.exists(output_image)

    # Test calibrate
    json_file = calibrate(["calibrate.JPG", "1"])  # `1` for the boxsize arg
    assert op.exists(json_file)
    seeds(["test.JPG", "--calibrate", json_file])
