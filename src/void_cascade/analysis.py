"""Analysis helpers for sandpile avalanche distributions.

SOC avalanche statistics typically have power-law tails P(s) ~ s^{-tau} with
a finite-size upper cutoff. Two operations recur throughout this project:

1. Log-binning the empirical histogram of s. Linear bins waste resolution at
   large s because each bin then contains very few events. Logarithmic bins
   yield roughly equal statistical noise per bin across decades.
2. Straight-line fitting on log-log axes to extract tau over the scaling
   range, picking the lower and upper cutoffs by eye to exclude the small-s
   rollover and the finite-size knee.

These are intentionally simple. For a publication-grade estimate of tau the
right tool is maximum-likelihood with goodness-of-fit (Clauset, Shalizi,
Newman 2009), but the goal at this milestone is just to confirm that a power
law exists over a couple of decades.
"""

from __future__ import annotations

import numpy as np


def log_binned_pdf(
    s: np.ndarray, n_bins: int = 30
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Probability density of s over logarithmically spaced bins.

    Zero-sized avalanches (no toppling) are dropped because they are not
    power-law-distributed and would dominate the smallest bin.

    Returns
    -------
    centers : ndarray
        Geometric-mean bin centers (the natural choice for log bins).
    pdf : ndarray
        Normalized density: counts_in_bin / (n_samples * bin_width).
    counts : ndarray
        Raw event counts per bin. Useful for assessing noise in the tail.
    """
    s = np.asarray(s)
    s = s[s > 0]
    if s.size == 0:
        return np.array([]), np.array([]), np.array([])

    edges = np.logspace(0.0, np.log10(s.max()) + 1e-9, n_bins + 1)
    counts, _ = np.histogram(s, bins=edges)
    widths = np.diff(edges)
    centers = np.sqrt(edges[:-1] * edges[1:])
    pdf = counts / (s.size * widths)
    mask = counts > 0
    return centers[mask], pdf[mask], counts[mask]


def fit_power_law(
    s: np.ndarray,
    pdf: np.ndarray,
    s_min: float = 1.0,
    s_max: float | None = None,
) -> tuple[float, float]:
    """Least-squares fit of log10(pdf) = -tau * log10(s) + c.

    The fit range [s_min, s_max] should exclude the small-s rollover and the
    finite-size cutoff; pick it by inspection of a log-log plot.

    Returns (tau, c) where tau is the power-law exponent and c is the log10
    of the prefactor.
    """
    mask = s >= s_min
    if s_max is not None:
        mask &= s <= s_max
    if mask.sum() < 2:
        raise ValueError("Need at least two points in the fit range.")
    slope, intercept = np.polyfit(np.log10(s[mask]), np.log10(pdf[mask]), 1)
    return float(-slope), float(intercept)
