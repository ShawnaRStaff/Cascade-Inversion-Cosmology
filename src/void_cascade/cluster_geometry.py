"""Fractal-dimension estimators for avalanche clusters and binary masks.

Two complementary measurements:

1. **Population fit** `fit_fractal_dimension(areas, rgyrs)`
   Treats the avalanche ensemble as samples of a self-similar object.
   The relation a ~ R_g^{D_f} holds when avalanches at different sizes
   are statistically similar, which is the case in the SOC steady state
   away from finite-size cutoffs. Slope of log(a) vs log(R_g) over the
   scaling range is D_f.

2. **Single-cluster box count** `box_count_dimension(mask)`
   Classical box-counting on one binary lattice snapshot. Counts boxes
   of side b that contain at least one lit cell, fits log(N) vs log(1/b).
   For a compact 2D blob D_f -> 2; for a line D_f -> 1; for a Sierpinski
   triangle D_f ~ 1.585.

Both are deliberately simple. Mandelbrot-style sandbox counting and
multifractal spectra are out of scope for Milestone 2 — the goal is to
verify that avalanche clusters in 2D Manna come out with a fractal
dimension we can compare to literature (~2.0 for the support, near-compact).
"""

from __future__ import annotations

import numpy as np


def fit_fractal_dimension(
    areas: np.ndarray,
    rgyrs: np.ndarray,
    r_min: float = 2.0,
    r_max: float | None = None,
    n_bins: int = 20,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Fit a ~ R_g^{D_f} from a population of avalanche clusters.

    To reduce noise we bin avalanches by R_g (log-spaced) and fit the
    mean area per bin. Equal weighting in log-R_g space avoids the
    dominance of the abundant small-cluster bins.

    Parameters
    ----------
    areas, rgyrs : 1D arrays of equal length
        Per-avalanche area and radius of gyration. NaN rgyrs and
        zero-area avalanches are dropped.
    r_min : float
        Lower fit cutoff. R_g below this is too small to resolve a
        meaningful shape (one or two sites).
    r_max : float or None
        Upper fit cutoff. Defaults to no cutoff; callers should set
        this below the finite-size knee R_g ~ L/2 if visible.
    n_bins : int
        Number of log-spaced R_g bins.

    Returns
    -------
    D_f : float
    D_f_err : float
        Standard error of the slope.
    centers : ndarray
        Bin centers (geometric mean of bin edges) used in the fit.
    mean_area : ndarray
        Mean area within each bin.
    """
    areas = np.asarray(areas, dtype=np.float64)
    rgyrs = np.asarray(rgyrs, dtype=np.float64)
    mask = np.isfinite(rgyrs) & (areas > 0) & (rgyrs > 0)
    a = areas[mask]
    r = rgyrs[mask]

    fit_mask = r >= r_min
    if r_max is not None:
        fit_mask &= r <= r_max
    a = a[fit_mask]
    r = r[fit_mask]
    if a.size < n_bins:
        raise ValueError(
            f"Not enough samples ({a.size}) for {n_bins} bins; "
            "lower n_bins or run a longer simulation."
        )

    edges = np.logspace(np.log10(r.min()), np.log10(r.max()) + 1e-9, n_bins + 1)
    bin_idx = np.digitize(r, edges) - 1
    bin_idx = np.clip(bin_idx, 0, n_bins - 1)

    centers = np.sqrt(edges[:-1] * edges[1:])
    mean_area = np.zeros(n_bins)
    counts = np.zeros(n_bins, dtype=np.int64)
    for b in range(n_bins):
        sel = bin_idx == b
        if sel.any():
            mean_area[b] = float(a[sel].mean())
            counts[b] = int(sel.sum())

    ok = counts >= 3
    if ok.sum() < 2:
        raise ValueError("Too few populated bins to fit.")

    log_r = np.log10(centers[ok])
    log_a = np.log10(mean_area[ok])
    coeffs, cov = np.polyfit(log_r, log_a, 1, cov=True)
    D_f = float(coeffs[0])
    D_f_err = float(np.sqrt(cov[0, 0]))
    return D_f, D_f_err, centers[ok], mean_area[ok]


def box_count_dimension(
    mask: np.ndarray,
    box_sizes: list[int] | None = None,
) -> tuple[float, float, np.ndarray, np.ndarray]:
    """Box-counting fractal dimension of a hypercubic binary mask (2D or 3D).

    For each box side b, count N(b) = number of b^d tiles (d = mask.ndim)
    that contain at least one True cell. Fit log(N) vs log(1/b) over the
    box sizes that are well inside (cell size, lattice size).

    For a compact d-dimensional object D_box -> d. For a (d-1)-dimensional
    sheet D_box -> d-1. For fractal sets, D_box is the non-integer
    Hausdorff dimension (in the usable scaling regime).

    Parameters
    ----------
    mask : bool ndarray, must be hypercubic (all shape entries equal) in
        2 or 3 dimensions.
    box_sizes : list of int, optional
        Box side lengths in lattice units. Defaults to powers of two from
        1 up to bbox_max / 4 of the cluster's bounding box.

    Returns
    -------
    D_box : float
    D_box_err : float
    box_sizes_used : ndarray
    counts : ndarray
        N(b) for each box size used.
    """
    mask = np.asarray(mask, dtype=bool)
    if mask.ndim not in (2, 3):
        raise ValueError(f"mask must be 2D or 3D, got {mask.ndim}D")
    L = mask.shape[0]
    if any(s != L for s in mask.shape):
        raise ValueError(f"mask must be hypercubic (got shape {mask.shape})")
    if not mask.any():
        raise ValueError("mask is empty; cannot box-count an empty set")

    if box_sizes is None:
        # Restrict to the scaling regime where the box is small compared
        # to the cluster's bounding box along its longest axis. Above
        # ~1/4 of the bbox the count saturates near 1 and biases the slope.
        bbox_max = 0
        for axis in range(mask.ndim):
            other = tuple(a for a in range(mask.ndim) if a != axis)
            idx = np.argwhere(np.any(mask, axis=other))
            extent = int(idx[-1, 0] - idx[0, 0] + 1)
            bbox_max = max(bbox_max, extent)
        b_max = max(bbox_max // 4, 2)
        bs = []
        b = 1
        while b <= b_max:
            bs.append(b)
            b *= 2
        box_sizes = bs
    box_sizes_arr = np.array(box_sizes, dtype=int)

    counts = np.zeros(len(box_sizes_arr), dtype=np.int64)
    for k, b in enumerate(box_sizes_arr):
        # Pad each axis with False so its length is a multiple of b, then
        # reshape into (n_tiles_axis_0, b, n_tiles_axis_1, b, ...) and
        # OR-reduce over the within-box axes.
        pad_width = tuple((0, (-L) % b) for _ in range(mask.ndim))
        if any(p[1] > 0 for p in pad_width):
            m = np.pad(mask, pad_width)
        else:
            m = mask
        n_tiles = [m.shape[axis] // b for axis in range(mask.ndim)]
        # Build interleaved shape: (n0, b, n1, b, ...).
        reshape_dims = []
        for n in n_tiles:
            reshape_dims.extend([n, b])
        tiled = m.reshape(reshape_dims)
        # Reduce over every odd-numbered axis (the within-box ones).
        within_box_axes = tuple(2 * i + 1 for i in range(mask.ndim))
        present = tiled.any(axis=within_box_axes)
        counts[k] = int(present.sum())

    log_inv_b = np.log10(1.0 / box_sizes_arr.astype(np.float64))
    log_N = np.log10(counts.astype(np.float64))
    if log_inv_b.size < 2:
        raise ValueError("Need at least 2 box sizes to fit.")
    coeffs, cov = np.polyfit(log_inv_b, log_N, 1, cov=True)
    return float(coeffs[0]), float(np.sqrt(cov[0, 0])), box_sizes_arr, counts
