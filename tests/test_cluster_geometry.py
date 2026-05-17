"""Tests for fractal-dimension estimators.

We hit each estimator with a synthetic mask of known dimension so the
tests fail loudly if a refactor breaks the relation a ~ R^{D_f} or the
box-count log-log fit.
"""

from __future__ import annotations

import numpy as np

from void_cascade.cluster_geometry import box_count_dimension, fit_fractal_dimension


def _disk_mask(L: int, radius: float) -> np.ndarray:
    cy, cx = (L - 1) / 2.0, (L - 1) / 2.0
    yy, xx = np.mgrid[0:L, 0:L]
    return (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2


def test_box_count_disk_dimension_near_two():
    # A solid disk is a 2D set; box-counting should give D ~ 2.
    mask = _disk_mask(L=256, radius=100)
    D, D_err, sizes, counts = box_count_dimension(mask)
    assert sizes.size >= 4
    assert 1.85 < D < 2.05, f"disk D_box = {D:.3f}, expected ~2"


def test_box_count_line_dimension_near_one():
    # A horizontal line in a 2D lattice has box-counting dimension 1.
    L = 256
    mask = np.zeros((L, L), dtype=bool)
    mask[L // 2, :] = True
    D, _, _, _ = box_count_dimension(mask)
    assert 0.85 < D < 1.15, f"line D_box = {D:.3f}, expected ~1"


def test_box_count_rejects_non_square():
    mask = np.ones((8, 16), dtype=bool)
    try:
        box_count_dimension(mask)
    except ValueError:
        return
    raise AssertionError("Expected ValueError on non-square mask")


def test_box_count_3d_solid_ball_near_three():
    # Solid 3D ball is a compact 3-dimensional set; D_box -> 3 in the
    # limit of infinite resolution. On a finite lattice, the standard
    # boundary-box-overlap bias drags the apparent slope down a few
    # percent; 2.7 is the realistic floor for R=40 on a 128^3 lattice.
    L = 128
    radius = 40
    cz = cy = cx = (L - 1) / 2.0
    zz, yy, xx = np.mgrid[0:L, 0:L, 0:L]
    mask = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
    D, _, _, _ = box_count_dimension(mask)
    assert 2.7 < D < 3.1, f"3D ball D_box = {D:.3f}, expected ~3"


def test_box_count_3d_plane_near_two():
    # An axis-aligned plane in 3D has box-counting dimension 2.
    L = 64
    mask = np.zeros((L, L, L), dtype=bool)
    mask[L // 2, :, :] = True
    D, _, _, _ = box_count_dimension(mask)
    assert 1.85 < D < 2.15, f"3D plane D_box = {D:.3f}, expected ~2"


def test_box_count_3d_line_near_one():
    # A 1D line embedded in 3D has D ~ 1.
    L = 64
    mask = np.zeros((L, L, L), dtype=bool)
    mask[L // 2, L // 2, :] = True
    D, _, _, _ = box_count_dimension(mask)
    assert 0.85 < D < 1.15, f"3D line D_box = {D:.3f}, expected ~1"


def test_box_count_rejects_4d():
    mask = np.ones((4, 4, 4, 4), dtype=bool)
    try:
        box_count_dimension(mask)
    except ValueError:
        return
    raise AssertionError("Expected ValueError on 4D mask")


def test_box_count_rejects_empty():
    mask = np.zeros((16, 16), dtype=bool)
    try:
        box_count_dimension(mask)
    except ValueError:
        return
    raise AssertionError("Expected ValueError on empty mask")


def test_fit_fractal_dimension_compact_clusters():
    # Synthesize avalanches as filled disks of varying radius: area = pi R^2,
    # R_g of a uniform disk = R / sqrt(2). The slope of log(a) vs log(R_g) is 2.
    rng = np.random.default_rng(0)
    radii = rng.uniform(2, 50, size=2000)
    areas = np.pi * radii ** 2
    rgyrs = radii / np.sqrt(2.0)
    D_f, D_f_err, _, _ = fit_fractal_dimension(
        areas, rgyrs, r_min=2.0, n_bins=15
    )
    assert abs(D_f - 2.0) < 0.05, f"compact-disk D_f = {D_f:.3f}, expected 2"


def test_fit_fractal_dimension_line_clusters():
    # Linear clusters: area = L, R_g of uniform segment of length L = L/sqrt(12).
    # Slope of log(area) vs log(R_g) is 1.
    rng = np.random.default_rng(0)
    lengths = rng.uniform(4, 200, size=2000)
    areas = lengths
    rgyrs = lengths / np.sqrt(12.0)
    D_f, _, _, _ = fit_fractal_dimension(areas, rgyrs, r_min=2.0, n_bins=15)
    assert abs(D_f - 1.0) < 0.05, f"line D_f = {D_f:.3f}, expected 1"
