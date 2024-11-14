#!/usr/bin/env python3
# -*- coding:utf-8 -*-
#
# webcolors.py
# utils
#
# Created by Haibao Tang on 01/28/20
# Copyright © 2021 Haibao Tang. All rights reserved.
#
import logging

import numpy as np

from skimage.color import rgb2lab, deltaE_cmc
from webcolors._definitions import _CSS3_NAMES_TO_HEX
from webcolors import hex_to_rgb


def color_diff(rgb1, rgb2):
    """
    Calculate distance between two RGB colors. See discussion:

    http://stackoverflow.com/questions/8863810/python-find-similar-colors-best-way

    - for basic / fast calculations, you can use dE76 but beware of its problems
    - for graphics arts use we recommend dE94 and perhaps dE-CMC 2:1
    - for textiles use dE-CMC
    """
    rgb1 = np.array(rgb1, dtype="float64").reshape(1, 1, 3) / 255.0
    rgb2 = np.array(rgb2, dtype="float64").reshape(1, 1, 3) / 255.0
    lab1 = rgb2lab(rgb1)
    lab2 = rgb2lab(rgb2)
    return deltaE_cmc(lab1, lab2, kL=2, kC=1)[0, 0]


def closest_color(requested_color):
    """
    Find closest color name for the request RGB tuple.
    """
    logging.disable(logging.DEBUG)
    colors = []
    for name, hex in _CSS3_NAMES_TO_HEX.items():
        diff = color_diff(hex_to_rgb(hex), requested_color)
        colors.append((diff, name))
    logging.disable(logging.NOTSET)
    _, min_color = min(colors)

    return min_color


if __name__ == "__main__":
    import doctest

    doctest.testmod()
