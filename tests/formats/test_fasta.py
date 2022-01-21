#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# test_gff.py
# formats
#
# Created by Haibao Tang on 01/21/22
# Copyright © 2022 Haibao Tang. All rights reserved.
#

from unittest.mock import patch, mock_open


def test_iter_clean_fasta():
    from jcvi.formats.fasta import iter_clean_fasta

    with patch(
        "builtins.open", mock_open(read_data=">seq1\nACGTxyz*\n>seq2\nG^CC@")
    ) as m:
        records = []
        for header, seq in iter_clean_fasta("test.fasta"):
            records.append((header, seq))
        assert records == [("seq1", "ACGTxyz*"), ("seq2", "GCC")]

    m.assert_called_once_with("test.fasta", "r")
