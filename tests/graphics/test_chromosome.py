#!/usr/bin/env python
# -*- coding: UTF-8 -*-

import matplotlib

matplotlib.use("Agg")
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
from matplotlib.transforms import Affine2D
import numpy as np

from jcvi.graphics.chromosome import (
    Chromosome,
    ChromosomeWithCentromere,
    HorizontalChromosome,
)


def _canvas():
    fig = plt.figure(figsize=(2, 2), dpi=100)
    ax = fig.add_axes((0, 0, 1, 1))
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_axis_off()
    return fig, ax


def _render(fig):
    fig.canvas.draw()
    rgb = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    plt.close(fig)
    return rgb


def _is_red(rgb, x, y):
    # Image rows run top-to-bottom; convert from axes fraction
    h, w = rgb.shape[:2]
    r, g, b = rgb[int((1 - y) * (h - 1)), int(x * (w - 1))]
    return r > 200 and g < 80 and b < 80


def test_chromosome_clip_keeps_painted_region_inside_caps():
    fig, ax = _canvas()
    chrom = Chromosome(ax, 0.5, 0.1, 0.9, width=0.2)
    painted = Rectangle((0.4, 0.8), 0.2, 0.1, fc="red", lw=0)
    chrom.clip(painted)
    ax.add_patch(painted)
    rgb = _render(fig)
    assert _is_red(rgb, 0.5, 0.85)  # inside the cap
    assert not _is_red(rgb, 0.41, 0.895)  # corner outside the rounded cap


def test_chromosome_patch_option_is_clipped():
    fig, ax = _canvas()
    Chromosome(
        ax, 0.5, 0.1, 0.9, width=0.2, patch=[0.1, 0.2, 0.8, 0.9], patchcolor="red"
    )
    rgb = _render(fig)
    assert _is_red(rgb, 0.5, 0.85)
    assert not _is_red(rgb, 0.42, 0.895)


def test_horizontal_chromosome_patch_option_is_clipped():
    fig, ax = _canvas()
    HorizontalChromosome(
        ax, 0.1, 0.9, 0.5, height=0.2, patch=[0.1, 0.2, 0.8, 0.9], patchcolor="red"
    )
    rgb = _render(fig)
    assert _is_red(rgb, 0.85, 0.5)
    assert not _is_red(rgb, 0.895, 0.42)


def test_chromosome_with_centromere_clips_to_both_arms():
    fig, ax = _canvas()
    chrom = ChromosomeWithCentromere(ax, 0.5, 0.9, 0.5, 0.1, width=0.2)
    painted = Rectangle((0.4, 0.1), 0.2, 0.8, fc="red", lw=0)
    chrom.clip(painted)
    ax.add_patch(painted)
    rgb = _render(fig)
    assert _is_red(rgb, 0.5, 0.7)  # upper arm
    assert _is_red(rgb, 0.5, 0.3)  # lower arm
    assert not _is_red(rgb, 0.41, 0.895)  # outside the top cap
    assert not _is_red(rgb, 0.41, 0.105)  # outside the bottom cap


def test_horizontal_chromosome_patch_follows_set_transform():
    fig, ax = _canvas()
    hc = HorizontalChromosome(
        ax, 0.1, 0.9, 0.3, height=0.2, patch=[0.1, 0.2, 0.8, 0.9], patchcolor="red"
    )
    hc.set_transform(Affine2D().translate(0, 0.4) + ax.transData)
    rgb = _render(fig)
    assert _is_red(rgb, 0.85, 0.7)  # inside the shifted cap
    assert not _is_red(rgb, 0.895, 0.62)  # corner outside the shifted cap
    assert not _is_red(rgb, 0.85, 0.3)  # nothing left at the old position
