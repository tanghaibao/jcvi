#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import io
import os
from unittest.mock import patch


@patch("builtins.input", return_value="username")
@patch("getpass.getpass", return_value="password")
def test_get_cookies(mock_username, mock_password):
    from jcvi.apps.fetch import get_cookies, PHYTOZOME_COOKIES
    from jcvi.apps.base import cleanup, which

    cleanup(PHYTOZOME_COOKIES)
    if which("curl"):
        assert get_cookies() == PHYTOZOME_COOKIES
    else:
        assert get_cookies() is None  # errored out with "curl not found"
    cleanup(PHYTOZOME_COOKIES)


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
