import os
import os.path as op

from jcvi.graphics.grabseeds import calibrate, seeds


def test_main(tmp_path, monkeypatch):
    data_dir = op.join(op.dirname(__file__), "data")
    for fname in os.listdir(data_dir):
        os.symlink(op.join(data_dir, fname), op.join(str(tmp_path), fname))
    monkeypatch.chdir(tmp_path)

    output_image = "test.pdf"
    seeds(["test.JPG"])
    assert op.exists(output_image)

    # Test calibrate
    json_file = calibrate(["calibrate.JPG", "1"])  # `1` for the boxsize arg
    assert op.exists(json_file)
    seeds(["test.JPG", "--calibrate", json_file])
