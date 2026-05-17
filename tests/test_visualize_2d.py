"""Smoke tests for the 2D Manna visualization module.

These verify that the animation and snapshot functions produce a file
at the requested path with non-zero size. They do not inspect the
contents of the output (a real image-content test belongs at a higher
level). matplotlib is run in the non-interactive 'Agg' backend so the
test does not require a display.
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from void_cascade.visualize_2d import record_animation, snapshot_large_avalanche


def test_record_animation_writes_gif(tmp_path):
    out = tmp_path / "anim.gif"
    path = record_animation(L=12, max_frames=40, seed=3, out_path=out, fps=10)
    assert path == out
    assert out.exists()
    assert out.stat().st_size > 0


def test_snapshot_large_avalanche_writes_png(tmp_path):
    out = tmp_path / "snap.png"
    # L=16 + short burn-in + tiny min_area so the search terminates quickly.
    path = snapshot_large_avalanche(
        L=16,
        burn_in=200,
        min_area=8,
        max_search_drops=2_000,
        seed=5,
        out_path=out,
    )
    assert path == out
    assert out.exists()
    assert out.stat().st_size > 0
