import os
import io
import shutil
import pytest
from pathlib import Path

from pypdf import PdfReader, PdfWriter

# Import the function under test and the module for patching
from jcvi.formats.pdf import cat as pdf_cat
import jcvi.formats.pdf as pdfmod


def _make_pdf(path: Path, sizes):
    """
    Create a PDF at `path` with one blank page per (w,h) in `sizes`.
    """
    writer = PdfWriter()
    for w, h in sizes:
        writer.add_blank_page(width=w, height=h)
    with open(path, "wb") as fh:
        writer.write(fh)


def _page_sizes(pdf_path: Path):
    r = PdfReader(str(pdf_path))
    out = []
    for p in r.pages:
        mb = p.mediabox
        out.append((float(mb.width), float(mb.height)))
    return out


def test_cat_basic_merge(tmp_path: Path):
    a = tmp_path / "a.pdf"
    b = tmp_path / "b.pdf"
    out = tmp_path / "out.pdf"

    _make_pdf(a, [(200, 200), (200, 200)])  # 2 pages
    _make_pdf(b, [(300, 300), (300, 300), (300, 300)])  # 3 pages

    # Merge all pages from both
    pdf_cat(["-o", str(out), str(a), ":", str(b), ":"])
    assert out.exists()
    r = PdfReader(str(out))
    assert len(r.pages) == 5


def test_cat_page_ranges(tmp_path: Path):
    x = tmp_path / "x.pdf"
    out = tmp_path / "out.pdf"

    # 5 pages of varying sizes so we can validate slicing easily
    sizes = [(200, 200), (210, 210), (220, 220), (230, 230), (240, 240)]
    _make_pdf(x, sizes)

    # Take pages [1:4] -> indices 1,2,3
    pdf_cat(["-o", str(out), str(x), "1:4"])
    assert _page_sizes(out) == sizes[1:4]


def test_cat_nosort_and_default_sort(tmp_path: Path):
    # Create two one-page PDFs with distinct sizes
    chap10 = tmp_path / "chapter10.pdf"
    chap2 = tmp_path / "chapter2.pdf"
    out_sorted = tmp_path / "sorted.pdf"
    out_nosort = tmp_path / "nosort.pdf"

    _make_pdf(chap10, [(300, 300)])  # intended to come after chap2 when sorted
    _make_pdf(chap2, [(200, 200)])

    # Default behavior: natural sort -> chapter2 then chapter10
    pdf_cat(["-o", str(out_sorted), str(chap10), str(chap2)])
    assert _page_sizes(out_sorted)[0] == (200.0, 200.0)

    # --nosort keeps the provided order -> chapter10 first
    pdf_cat(["--nosort", "-o", str(out_nosort), str(chap10), str(chap2)])
    assert _page_sizes(out_nosort)[0] == (300.0, 300.0)


def test_cat_single_file_copy_fastpath(tmp_path: Path, monkeypatch):
    src = tmp_path / "src.pdf"
    out = tmp_path / "out.pdf"
    _make_pdf(src, [(222, 222), (333, 333)])

    calls = {"count": 0}
    real_copy = shutil.copyfile

    def _spy_copy(src_path, dst_path):
        calls["count"] += 1
        return real_copy(src_path, dst_path)

    # Patch the module-level shutil.copyfile used inside cat()
    monkeypatch.setattr(pdfmod.shutil, "copyfile", _spy_copy)

    pdf_cat(["-o", str(out), str(src)])
    assert out.exists()
    # Fast-path must have used shutil.copyfile exactly once
    assert calls["count"] == 1
    # And the bytes should represent the same two pages
    assert _page_sizes(out) == _page_sizes(src)


def test_cat_ignores_outfile_if_present_in_args(tmp_path: Path):
    a = tmp_path / "a.pdf"
    b = tmp_path / "b.pdf"
    out = tmp_path / "out.pdf"

    _make_pdf(a, [(200, 200)])
    _make_pdf(b, [(300, 300)])
    # Pre-create an 'out.pdf' that will also (wrongly) appear in input args
    _make_pdf(out, [(999, 999)])

    # out.pdf appears in args but should be removed by cat() before merging
    pdf_cat(["-o", str(out), str(a), str(out), str(b)])
    assert _page_sizes(out) == [(200.0, 200.0), (300.0, 300.0)]


def test_cat_cleanup_removes_inputs(tmp_path: Path):
    x = tmp_path / "x.pdf"
    out = tmp_path / "out.pdf"
    _make_pdf(x, [(200, 200)])

    pdf_cat(["--cleanup", "-o", str(out), str(x)])
    assert out.exists()
    assert not x.exists()  # input removed


def test_cat_missing_input_raises_systemexit(tmp_path: Path):
    out = tmp_path / "out.pdf"
    bad = tmp_path / "nope.pdf"  # does not exist

    with pytest.raises(SystemExit):
        pdf_cat(["-o", str(out), str(bad)])
