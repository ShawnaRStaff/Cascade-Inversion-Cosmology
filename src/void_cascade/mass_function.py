"""Connected-cluster mass function n(s) for a 3D bool occupation field.

For a percolation-like system, the distribution of connected-cluster
sizes (or "masses") is a standard observable. At criticality random
percolation gives n(s) ~ s^{-tau_s} with tau_s ~ 2.189 in 3D (Stauffer
& Aharony). Away from p_c, the distribution acquires a finite-size
cutoff: above the typical cluster mass, n(s) falls off (exponentially
or super-exponentially depending on the regime).

For galaxy/halo cosmology the standard reference forms are
Press-Schechter (1974) and its descendants (Sheth-Tormen 1999):

    n(M) ~ M^{-alpha} exp(-(M / M_star)^{beta})

with alpha ~ 2 on the small-mass end and an exponential cutoff at M_star.
The shape, dimensionless on log axes, is what we compare to.

This module:
1. Labels connected components in a 3D bool mask using the same
   6-connectivity as the percolation detector (matching the toppling
   topology of the Manna model).
2. Histograms the resulting cluster sizes into log-spaced bins.
3. Returns the histogram for log-log fitting.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage


_STRUCT_6 = ndimage.generate_binary_structure(rank=3, connectivity=1)


def cluster_sizes(mask: np.ndarray) -> np.ndarray:
    """Return the size of every connected component in a 3D bool mask.

    6-connectivity, matching the Manna toppling topology and the
    percolation detector. Returns an ndarray of ints in arbitrary order;
    cluster size 0 (background) is excluded.
    """
    if mask.ndim != 3:
        raise ValueError(f"mask must be 3D (got {mask.ndim}D)")
    if not mask.any():
        return np.array([], dtype=np.int64)
    labels, _n = ndimage.label(mask, structure=_STRUCT_6)
    sizes = np.bincount(labels.ravel())
    return sizes[1:].astype(np.int64)   # drop background (label 0)


def cluster_mass_pdf(
    sizes: np.ndarray, n_bins: int = 20, min_size: int = 1
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Log-spaced histogram of cluster sizes.

    Returns (centers, pdf, counts) where:
      centers : geometric-mean cluster size per bin
      pdf     : counts / (n_clusters * bin_width)
      counts  : raw cluster count per bin

    Empty bins are dropped from the output. Use the same convention as
    the avalanche-PDF estimator in analysis.py so the existing log-log
    fitter applies unchanged.
    """
    sizes = np.asarray(sizes, dtype=np.int64)
    sizes = sizes[sizes >= min_size]
    if sizes.size == 0:
        return np.array([]), np.array([]), np.array([])
    edges = np.logspace(np.log10(min_size), np.log10(sizes.max()) + 1e-9, n_bins + 1)
    counts, _ = np.histogram(sizes, bins=edges)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / (sizes.size * widths)
    keep = counts > 0
    return centers[keep], pdf[keep], counts[keep]
