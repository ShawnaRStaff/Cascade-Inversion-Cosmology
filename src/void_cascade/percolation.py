"""Percolation detection on a 3D boolean lattice.

The cosmologically meaningful event in the void-cascade model is the
moment when the cumulative fractured region first becomes globally
connected - that is, when a single connected cluster of toppled sites
spans the box from one face to the opposite face.

This module wraps `scipy.ndimage.label` with a 6-connectivity structuring
element (faces only, not edges or corners) and reports:

- the spanning state (does any cluster touch two opposite faces?),
- which axis spans (x, y, z, or any subset),
- the size of the largest cluster.

We use 6-connectivity rather than 26-connectivity because the underlying
lattice topology is cubic with face-sharing neighbors - that matches the
toppling rule in sandpile_3d.py. Choosing 26-connectivity would let
clusters connect through diagonal contacts that the sandpile dynamics
cannot produce; we'd then count percolation transitions that are
artifacts of the connectivity convention, not the dynamics.

Spanning convention: a cluster "spans the x-axis" if and only if at least
one of its sites lies on the i=0 face AND at least one lies on the i=L-1
face. Same for y, z. Total percolation = "spans at least one axis."
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import ndimage


# 6-connectivity in 3D: face neighbors only.
_STRUCT_6 = ndimage.generate_binary_structure(rank=3, connectivity=1)


@dataclass
class PercolationResult:
    """Outcome of one spanning check.

    Attributes
    ----------
    percolates : bool
        True iff at least one cluster spans some axis.
    spans_x, spans_y, spans_z : bool
        Per-axis spanning flags.
    largest_cluster_size : int
        Number of sites in the largest connected component. 0 if the
        mask is empty.
    n_clusters : int
        Number of distinct connected components.
    """

    percolates: bool
    spans_x: bool
    spans_y: bool
    spans_z: bool
    largest_cluster_size: int
    n_clusters: int


def check_spanning(mask: np.ndarray) -> PercolationResult:
    """Run connected-component analysis on a 3D bool mask and report spanning.

    Parameters
    ----------
    mask : 3D bool ndarray, shape (L, L, L)
        The cumulative ever-toppled set (or any binary 3D field you want
        to test for spanning).

    Returns
    -------
    PercolationResult
    """
    if mask.ndim != 3:
        raise ValueError("mask must be 3D")
    L = mask.shape[0]
    if mask.shape[1] != L or mask.shape[2] != L:
        raise ValueError("mask must be cubic (L, L, L)")
    if not mask.any():
        return PercolationResult(
            percolates=False,
            spans_x=False,
            spans_y=False,
            spans_z=False,
            largest_cluster_size=0,
            n_clusters=0,
        )

    labels, n_clusters = ndimage.label(mask, structure=_STRUCT_6)

    # Cluster IDs present on each face. Label 0 is the background; exclude.
    def _ids_on_face(face: np.ndarray) -> set[int]:
        ids = np.unique(face)
        return {int(v) for v in ids if v != 0}

    x0 = _ids_on_face(labels[0, :, :])
    x1 = _ids_on_face(labels[-1, :, :])
    y0 = _ids_on_face(labels[:, 0, :])
    y1 = _ids_on_face(labels[:, -1, :])
    z0 = _ids_on_face(labels[:, :, 0])
    z1 = _ids_on_face(labels[:, :, -1])

    spans_x = len(x0 & x1) > 0
    spans_y = len(y0 & y1) > 0
    spans_z = len(z0 & z1) > 0

    # Largest cluster size via bincount over labels (label 0 = background).
    sizes = np.bincount(labels.ravel())
    largest = int(sizes[1:].max()) if sizes.size > 1 else 0

    return PercolationResult(
        percolates=spans_x or spans_y or spans_z,
        spans_x=spans_x,
        spans_y=spans_y,
        spans_z=spans_z,
        largest_cluster_size=largest,
        n_clusters=int(n_clusters),
    )


def fractured_fraction(mask: np.ndarray) -> float:
    """Return the fraction of sites that are True (i.e. have toppled).

    Useful as the running "occupation probability" for plotting the
    largest-cluster size or the spanning state against the standard
    percolation control parameter p.
    """
    return float(mask.mean())
