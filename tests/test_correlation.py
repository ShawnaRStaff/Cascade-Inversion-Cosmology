"""Tests for the two-point correlation function on 3D bool masks.

Each test hits the estimator with a synthetic field whose xi(r) is
predictable, so the test fails loudly if a refactor breaks the
normalization or the radial binning.
"""

from __future__ import annotations

import numpy as np

from void_cascade.correlation import (
    power_law_galaxy_xi,
    two_point_correlation_3d,
)


def test_uncorrelated_random_xi_near_zero():
    # IID Bernoulli(p) field: xi(r) ~ 0 at every r > 0, with statistical
    # scatter that shrinks with lattice size and number of pair samples.
    rng = np.random.default_rng(0)
    L = 64
    p = 0.2
    mask = rng.random((L, L, L)) < p
    r, xi, n_pairs = two_point_correlation_3d(mask, n_bins=10)
    # Drop empty bins and the tiniest-r bin where shot noise is worst.
    valid = ~np.isnan(xi) & (n_pairs > 1000)
    assert valid.sum() >= 5
    # All values within +/- 0.02 of zero.
    assert np.all(np.abs(xi[valid]) < 0.02), (
        f"uncorrelated xi values: {xi[valid]}"
    )


def test_full_mask_xi_is_zero():
    # n(x) = 1 everywhere: <nn> = 1 = <n>^2, xi = 0 exactly.
    mask = np.ones((32, 32, 32), dtype=bool)
    _r, xi, _np = two_point_correlation_3d(mask, n_bins=5)
    valid = ~np.isnan(xi)
    assert np.allclose(xi[valid], 0.0, atol=1e-9)


def test_single_blob_clustered_at_small_r():
    # A solid sphere has positive xi at separations smaller than the
    # blob diameter (pairs inside the blob are overrepresented relative
    # to the uniform expectation), and falls toward -<n> at large r
    # where the only possible pair-occupancy structure is the blob's
    # geometry.
    L = 64
    radius = 8
    cz = cy = cx = L // 2
    zz, yy, xx = np.mgrid[0:L, 0:L, 0:L]
    mask = (zz - cz) ** 2 + (yy - cy) ** 2 + (xx - cx) ** 2 <= radius ** 2
    r, xi, n_pairs = two_point_correlation_3d(mask, n_bins=12, log_bins=True)
    valid = ~np.isnan(xi) & (n_pairs > 100)
    # Small-r values strongly positive (inside-blob pair probability is
    # huge relative to the overall sparse density of the field).
    small_r = (r < radius) & valid
    assert small_r.any()
    assert np.all(xi[small_r] > 1.0)


def test_rejects_non_cubic():
    try:
        two_point_correlation_3d(np.ones((8, 8, 16), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on non-cubic mask")


def test_rejects_empty():
    try:
        two_point_correlation_3d(np.zeros((8, 8, 8), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on empty mask")


def test_rejects_non_3d():
    try:
        two_point_correlation_3d(np.ones((8, 8), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on 2D mask")


def test_power_law_reference_form():
    r = np.array([1.0, 2.0, 5.0, 10.0])
    xi = power_law_galaxy_xi(r, r0=5.0, gamma=1.8)
    # By construction, xi(r0) = 1.
    assert abs(xi[2] - 1.0) < 1e-12
    # And xi(r > r0) < 1.
    assert xi[3] < 1.0
    # And xi(r < r0) > 1.
    assert xi[0] > 1.0
    assert xi[1] > 1.0
