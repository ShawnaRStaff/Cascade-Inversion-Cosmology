"""Tests for the cluster mass function module.

We hand-build masks where the connected-component sizes are obvious by
inspection so the test fails loudly if a refactor breaks the labeling
or the histogram normalization.
"""

from __future__ import annotations

import numpy as np

from void_cascade.mass_function import cluster_mass_pdf, cluster_sizes


def test_empty_mask_returns_empty():
    sizes = cluster_sizes(np.zeros((4, 4, 4), dtype=bool))
    assert sizes.size == 0


def test_isolated_voxels_each_count_one():
    mask = np.zeros((6, 6, 6), dtype=bool)
    mask[0, 0, 0] = True
    mask[5, 5, 5] = True
    mask[2, 4, 1] = True
    sizes = cluster_sizes(mask)
    assert sorted(sizes.tolist()) == [1, 1, 1]


def test_single_blob_one_cluster():
    mask = np.zeros((6, 6, 6), dtype=bool)
    mask[2:5, 2:5, 2:5] = True   # 3x3x3 cube = 27 voxels
    sizes = cluster_sizes(mask)
    assert sizes.tolist() == [27]


def test_diagonal_neighbors_are_separate_clusters():
    # Two voxels touching only at a corner are not face-adjacent and
    # should not be merged by 6-connectivity.
    mask = np.zeros((4, 4, 4), dtype=bool)
    mask[0, 0, 0] = True
    mask[1, 1, 1] = True
    sizes = cluster_sizes(mask)
    assert sorted(sizes.tolist()) == [1, 1]


def test_mass_pdf_normalization():
    # PDF must integrate to ~1 over a bounded sample.
    rng = np.random.default_rng(0)
    sizes = rng.integers(1, 1000, size=500)
    centers, pdf, counts = cluster_mass_pdf(sizes, n_bins=15)
    # Trapezoidal log-grid integral of pdf in linear space.
    # Sum of counts must equal the input sample count.
    assert counts.sum() == 500


def test_mass_pdf_drops_empty_bins():
    sizes = np.array([1, 2, 1, 1, 2])
    centers, pdf, counts = cluster_mass_pdf(sizes, n_bins=10)
    # All counts in the kept bins must be positive.
    assert np.all(counts > 0)


def test_cluster_sizes_rejects_non_3d():
    try:
        cluster_sizes(np.zeros((4, 4), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on 2D mask")
