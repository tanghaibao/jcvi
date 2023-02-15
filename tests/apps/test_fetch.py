#!/usr/bin/env python
# -*- coding: UTF-8 -*-

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
