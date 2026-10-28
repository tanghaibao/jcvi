#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import io


def _file(name, status="RESTORED", added="2020-01-01", **meta):
    return {
        "_id": name,
        "file_name": name,
        "file_status": status,
        "added_date": added,
        "metadata": meta,
    }


def _organism(genome_id, version, files, policy="Unrestricted"):
    return {
        "id": f"Phytozome-{genome_id}",
        "name": f"Species {version}",
        "data_utilization_status": policy,
        "top_hit": {"metadata": {"phytozome": {"version": version}}},
        "files": files,
    }


def test_classify_phytozome_file():
    from jcvi.apps.fetch import classify_phytozome_file

    assert classify_phytozome_file(_file("Vv_457_v2.1.gene.gff3.gz")) == "gff"
    assert classify_phytozome_file(_file("Vv_145_gene.gff3.gz")) == "gff"
    assert classify_phytozome_file(_file("Vv_457_v2.1.gene_exons.gff3.gz")) is None
    assert classify_phytozome_file(_file("Vv_457_v2.1.cds.fa.gz")) == "cds"
    assert classify_phytozome_file(_file("Vv_457_v2.1.lowConfidence.cds.fa.gz")) is None
    assert (
        classify_phytozome_file(_file("Vv_457_v2.1.cds_primaryTranscriptOnly.fa.gz"))
        is None
    )
    assert classify_phytozome_file(_file("Vv_457_12X.fa.gz", type="assembly")) == "asm"
    masked = _file("Vv_457_12X.softmasked.fa.gz", type="assembly")
    assert classify_phytozome_file(masked) is None
    assert classify_phytozome_file(_file("Vv_457_v2.1.protein.fa.gz")) is None


def test_select_phytozome_files():
    from jcvi.apps.fetch import select_phytozome_files

    old = _organism(
        145, "12X", [_file("Vv_145_gene.gff3.gz"), _file("Vv_145_cds.fa.gz")]
    )
    new = _organism(
        457,
        "v2.1",
        [
            _file("Vv_457_v2.0.gene.gff3.gz", added="2024-01-01"),
            _file("Vv_457_v2.1.gene.gff3.gz"),
            _file("Vv_457_v2.1.cds.fa.gz", status="PURGED", added="2024-01-01"),
            _file("Vv_457_v2.1.cds.fa.gz"),
        ],
    )
    # Newest genome wins, on-disk and version-matched files are preferred
    res = select_phytozome_files([old, new])
    assert res["gff"]["file_name"] == "Vv_457_v2.1.gene.gff3.gz"
    assert res["cds"]["file_status"] == "RESTORED"
    # Fall back to older genome if the newest lacks requested files
    new["files"] = new["files"][:2]
    res = select_phytozome_files([old, new])
    assert res["cds"]["file_name"] == "Vv_145_cds.fa.gz"
    assert select_phytozome_files([old], assembly=True) == {}
    assert select_phytozome_files([]) == {}


def test_get_jgi_token(monkeypatch):
    from jcvi.apps.fetch import JGI_TOKEN_ENV, get_jgi_token

    monkeypatch.setenv(JGI_TOKEN_ENV, " /api/sessions/abc\n")
    assert get_jgi_token() == "Bearer /api/sessions/abc"
    monkeypatch.setenv(JGI_TOKEN_ENV, "Bearer /api/sessions/abc")
    assert get_jgi_token() == "Bearer /api/sessions/abc"


def test_usage_with_percent_prog():
    """Test that usage strings with %prog don't cause ValueError in argparse.

    This is a regression test for the issue where docstrings containing %prog
    were set to p.usage after OptionParser initialization, bypassing the
    %prog -> %(prog)s conversion that OptionParser.__init__ performs.
    """
    from jcvi.apps.base import OptionParser

    # Create an OptionParser with a docstring containing %prog
    doc = """
    %prog test species

    Test description here. Example:

    $ %prog test species_name
    """
    p = OptionParser(doc)

    # Simulate what fetch.py does: update usage with additional content
    additional_info = "species1  species2  species3"
    updated_doc = "\n".join((doc, additional_info))
    # The fix: replace %prog with %(prog)s
    p.usage = updated_doc.replace("%prog", "%(prog)s")

    # This should not raise ValueError
    output = io.StringIO()
    p.print_help(file=output)
    help_text = output.getvalue()

    # Verify the help text contains the expected content
    assert "test species" in help_text
    assert "species1" in help_text
