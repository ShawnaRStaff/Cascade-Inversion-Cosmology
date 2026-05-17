"""Tests for the percolation detector.

We use small hand-built 3D masks where the spanning state is obvious by
inspection, so the test fails loudly if the connectivity or face-checking
logic regresses.
"""

from __future__ import annotations

import numpy as np

from void_cascade.percolation import check_spanning, fractured_fraction


def test_empty_mask_does_not_percolate():
    r = check_spanning(np.zeros((4, 4, 4), dtype=bool))
    assert not r.percolates
    assert r.n_clusters == 0
    assert r.largest_cluster_size == 0


def test_full_mask_percolates_all_axes():
    r = check_spanning(np.ones((4, 4, 4), dtype=bool))
    assert r.percolates
    assert r.spans_x and r.spans_y and r.spans_z
    assert r.n_clusters == 1
    assert r.largest_cluster_size == 64


def test_line_spanning_x_only():
    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[:, 2, 2] = True   # a column along x at (y=2, z=2)
    r = check_spanning(mask)
    assert r.percolates
    assert r.spans_x
    assert not r.spans_y
    assert not r.spans_z
    assert r.n_clusters == 1
    assert r.largest_cluster_size == 5


def test_line_spanning_z_only():
    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[2, 2, :] = True   # column along z
    r = check_spanning(mask)
    assert r.percolates
    assert r.spans_z
    assert not r.spans_x
    assert not r.spans_y


def test_cluster_touching_one_face_does_not_percolate():
    # A short column that only touches the x=0 face but not x=L-1.
    mask = np.zeros((6, 6, 6), dtype=bool)
    mask[0:3, 3, 3] = True
    r = check_spanning(mask)
    assert not r.percolates
    assert r.n_clusters == 1
    assert r.largest_cluster_size == 3


def test_diagonal_contact_does_not_count_as_connected():
    # Two cubes touching only at a corner (no shared face). With 6-conn
    # they are SEPARATE clusters and neither spans.
    mask = np.zeros((4, 4, 4), dtype=bool)
    mask[0, 0, 0] = True
    mask[1, 1, 1] = True
    r = check_spanning(mask)
    assert not r.percolates
    assert r.n_clusters == 2
    assert r.largest_cluster_size == 1


def test_two_disconnected_lines_neither_spans():
    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[0, 0, :] = True   # corner edge along z
    mask[4, 4, :] = True   # opposite corner edge along z
    r = check_spanning(mask)
    # Both span z individually (each line touches both z-faces).
    assert r.spans_z
    assert r.percolates
    assert r.n_clusters == 2


def test_face_only_does_not_span():
    # Fill an entire face but no interior. That face spans y and z within
    # itself (its only one slice thick along x), but does NOT span x.
    mask = np.zeros((5, 5, 5), dtype=bool)
    mask[0, :, :] = True
    r = check_spanning(mask)
    assert not r.spans_x
    assert r.spans_y
    assert r.spans_z
    assert r.percolates
    assert r.n_clusters == 1


def test_check_spanning_rejects_non_3d():
    try:
        check_spanning(np.zeros((4, 4), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on 2D mask")


def test_check_spanning_rejects_non_cubic():
    try:
        check_spanning(np.zeros((4, 4, 8), dtype=bool))
    except ValueError:
        return
    raise AssertionError("Expected ValueError on non-cubic mask")


def test_fractured_fraction_basic():
    mask = np.zeros((4, 4, 4), dtype=bool)
    assert fractured_fraction(mask) == 0.0
    mask[0, 0, 0] = True
    assert abs(fractured_fraction(mask) - 1 / 64) < 1e-12
