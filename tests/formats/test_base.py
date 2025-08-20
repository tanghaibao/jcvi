import os
import gzip
import pytest

from pathlib import Path

from jcvi.formats.base import FileMerger  # adjust import


def w(path: Path, data: bytes):
    path.write_bytes(data)


def wz(path: Path, data: bytes):
    with gzip.open(path, "wb") as f:
        f.write(data)


def rz(path: Path) -> bytes:
    with gzip.open(path, "rb") as f:
        return f.read()


def test_plain_to_plain_merge(tmp_path: Path):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    w(a, b"A")
    w(b, b"B")
    out = tmp_path / "out.txt"
    FileMerger([a, b], out).merge()
    assert out.read_bytes() == b"AB"


def test_gz_to_gz_fast_concat(tmp_path: Path):
    a, b = tmp_path / "a.gz", tmp_path / "b.gz"
    wz(a, b"A")
    wz(b, b"B")
    out = tmp_path / "out.gz"
    FileMerger([a, b], out).merge()
    # The concatenated .gz should decompress to AB
    assert rz(out) == b"AB"


def test_gz_to_plain_decompress(tmp_path: Path):
    a, b = tmp_path / "a.gz", tmp_path / "b.gz"
    wz(a, b"foo")
    wz(b, b"bar")
    out = tmp_path / "out.txt"
    FileMerger([a, b], out).merge()
    assert out.read_bytes() == b"foobar"


def test_plain_to_gz_compress(tmp_path: Path):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    w(a, b"foo")
    w(b, b"bar")
    out = tmp_path / "out.gz"
    FileMerger([a, b], out).merge()
    assert rz(out) == b"foobar"


def test_mixed_inputs_to_plain(tmp_path: Path):
    a, b = tmp_path / "a.txt", tmp_path / "b.gz"
    w(a, b"X")
    wz(b, b"Y")
    out = tmp_path / "out.txt"
    FileMerger([a, b], out).merge()
    assert out.read_bytes() == b"XY"


def test_mixed_inputs_to_gz(tmp_path: Path):
    a, b = tmp_path / "a.txt", tmp_path / "b.gz"
    w(a, b"X")
    wz(b, b"Y")
    out = tmp_path / "out.gz"
    FileMerger([a, b], out).merge()
    assert rz(out) == b"XY"


def test_single_file_copy_same_compression(tmp_path: Path):
    a = tmp_path / "a.txt"
    w(a, b"hello")
    out = tmp_path / "out.txt"
    FileMerger([a], out).merge()
    assert out.read_bytes() == b"hello"


def test_single_file_transform_to_gz(tmp_path: Path):
    a = tmp_path / "a.txt"
    w(a, b"hello")
    out = tmp_path / "out.gz"
    FileMerger([a], out).merge()
    assert rz(out) == b"hello"


def test_skip_if_uptodate(tmp_path: Path, monkeypatch):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    w(a, b"A")
    w(b, b"B")
    out = tmp_path / "out.txt"
    FileMerger([a, b], out).merge()
    before = out.read_bytes()
    # Make outfile newer than inputs
    os.utime(out, None)
    FileMerger([a, b], out).merge(skip_if_uptodate=True)
    assert out.read_bytes() == before  # unchanged


def test_outfile_same_as_input_noop(tmp_path: Path):
    a = tmp_path / "a.txt"
    w(a, b"one")
    # outfile == input -> no-op
    FileMerger([a], a).merge()
    assert a.read_bytes() == b"one"


def test_error_on_missing_input(tmp_path: Path):
    out = tmp_path / "out.txt"
    with pytest.raises(FileNotFoundError):
        FileMerger([tmp_path / "missing.txt"], out).merge()


def test_overwrite_false_raises(tmp_path: Path):
    a, b = tmp_path / "a.txt", tmp_path / "b.txt"
    w(a, b"A")
    w(b, b"B")
    out = tmp_path / "out.txt"
    FileMerger([a], out).merge()
    with pytest.raises(FileExistsError):
        FileMerger([b], out).merge(overwrite=False)
