"""Two-point correlation function xi(r) for a 3D boolean occupation field.

The two-point correlation function is the standard cosmological-structure
observable. For an occupation field n(x) with mean density <n>:

    xi(r) = <n(x) n(x+r)> / <n>^2 - 1

xi(r) > 0 means sites at separation r are clustered (more likely to be
co-occupied than chance); xi(r) < 0 means they are anti-correlated; xi(r)
= 0 means uncorrelated. For galaxies, published fits are approximately
xi(r) ~ (r / r_0)^(-gamma) with r_0 ~ 5 h^-1 Mpc and gamma ~ 1.8 on small
scales (Peebles 1980; Davis & Peebles 1983; Zehavi et al. 2011 SDSS).

This module computes xi(r) for the union-of-ever-toppled set from the 3D
Manna simulation. The "lattice spacing in Mpc" is a free parameter we fit
to data; absolute scales depend on it, but the *shape* of xi(r) on a
dimensionless r / L_box axis does not.

Method
------
FFT-based autocorrelation:

    F     = fftn(n)
    pairs = ifftn(|F|^2).real           # this gives sum_x n(x) n(x+r)
    <nn>  = pairs / N
    xi(r) = <nn> / <n>^2 - 1

then radially average over lattice sites with the same |r|. We use the
shortest-displacement metric on the torus (the "wrap" minimum) which
matches what the FFT produces.

Caveats
-------
The FFT assumes periodic boundary conditions; the underlying Manna
simulation has open boundaries. This introduces a boundary artifact:
sites near a face are less likely to be in the ever-toppled set because
their grains can leave the system. The artifact is bounded and
concentrated near r ~ L/2 in this estimator. Small-r structure (the
cosmologically interesting range) is robust. For publication-grade work
on open-boundary systems we would mask to a central sub-volume and use
a Landy-Szalay estimator with explicit random catalogs, but for a first
cut this FFT estimator is the right tool.
"""

from __future__ import annotations

import numpy as np


def two_point_correlation_3d(
    mask: np.ndarray,
    r_max: float | None = None,
    n_bins: int = 20,
    log_bins: bool = True,
    min_r: float = 1.0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute xi(r) for a 3D bool occupation mask.

    Parameters
    ----------
    mask : 3D bool ndarray of shape (L, L, L)
    r_max : float, optional
        Upper limit of the radial range. Defaults to L/2, which is the
        largest separation that exists on a periodic lattice.
    n_bins : int
    log_bins : bool
        If True, log-spaced bins (right for power-law tails). If False,
        linear.
    min_r : float
        Lower bin edge in lattice units.

    Returns
    -------
    r_centers : ndarray
        Mean separation within each populated bin.
    xi : ndarray
        Two-point correlation function value in each bin. NaN where the
        bin is empty (no lattice pair has |r| in that bin).
    n_pairs : ndarray
        Number of lattice pairs contributing to each bin. Use this as a
        weight if fitting; small-pair bins are noise-dominated.
    """
    if mask.ndim != 3:
        raise ValueError(f"mask must be 3D (got {mask.ndim}D)")
    L = mask.shape[0]
    if any(s != L for s in mask.shape):
        raise ValueError(f"mask must be cubic (got {mask.shape})")
    if not mask.any():
        raise ValueError("mask is empty")

    mask_f = mask.astype(np.float64)
    N = mask_f.size
    mean_n = float(mask_f.mean())

    F = np.fft.fftn(mask_f)
    pairs = np.fft.ifftn(F * F.conj()).real
    pair_prob = pairs / N  # = <n(x) n(x+r)> with periodic wrap

    # Shortest-displacement metric on the torus: r along each axis is
    # min(i, L - i), so r = 0 at the origin and grows to L//2 at the
    # antipode. Matches the FFT output's natural indexing.
    coords = np.arange(L)
    delta = np.minimum(coords, L - coords)
    dx, dy, dz = np.meshgrid(delta, delta, delta, indexing="ij")
    r_lattice = np.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

    if r_max is None:
        r_max = float(L) / 2.0

    if log_bins:
        edges = np.logspace(np.log10(min_r), np.log10(r_max), n_bins + 1)
    else:
        edges = np.linspace(min_r, r_max, n_bins + 1)

    r_flat = r_lattice.ravel()
    pair_prob_flat = pair_prob.ravel()

    xi = np.full(n_bins, np.nan)
    n_pairs = np.zeros(n_bins, dtype=np.int64)
    r_centers = np.zeros(n_bins)

    for i in range(n_bins):
        sel = (r_flat >= edges[i]) & (r_flat < edges[i + 1])
        n_sel = int(sel.sum())
        if n_sel == 0:
            r_centers[i] = np.sqrt(edges[i] * edges[i + 1])
            continue
        r_centers[i] = float(r_flat[sel].mean())
        mean_pp = float(pair_prob_flat[sel].mean())
        xi[i] = mean_pp / (mean_n ** 2) - 1.0
        n_pairs[i] = n_sel

    return r_centers, xi, n_pairs


def power_law_with_cutoff(
    r: np.ndarray, A: float, alpha: float, xi_corr: float
) -> np.ndarray:
    """Power-law-with-exponential-cutoff functional form.

        xi(r) = A * r^(-alpha) * exp(-r / xi_corr)

    This is the standard finite-correlation-length form: a self-similar
    underlying correlation modulated by an exponential damping above
    the correlation length xi_corr. As xi_corr -> infinity, recovers a
    pure power law. As xi_corr -> 0, the cutoff dominates and the
    apparent power-law slope (from a pure-power-law fit) becomes
    arbitrarily steep.

    The hypothesis (Speculation #1 from the M4 design discussion): a
    single (alpha, A) pair plus a p-dependent xi_corr(p) fits all
    Manna xi(r) measurements. If true, the apparent gamma we measured
    is a fitting artifact and the underlying correlation has a fixed
    universal slope.
    """
    return A * r ** (-alpha) * np.exp(-r / xi_corr)


def fit_power_law_with_cutoff(
    r: np.ndarray,
    xi: np.ndarray,
    sigma: np.ndarray | None = None,
    r_min: float = 1.0,
    r_max: float | None = None,
) -> tuple[float, float, float, np.ndarray] | None:
    """Fit xi(r) = A r^(-alpha) exp(-r / xi_corr).

    Returns (A, alpha, xi_corr, perr) where perr is the 3-element array
    of standard errors on (A, alpha, xi_corr). None if the fit fails or
    too few valid points.
    """
    from scipy.optimize import curve_fit

    mask = np.isfinite(xi) & (xi > 0)
    if r_min is not None:
        mask &= r >= r_min
    if r_max is not None:
        mask &= r <= r_max
    if mask.sum() < 4:
        return None
    r_fit = r[mask]
    xi_fit = xi[mask]
    sig_fit = sigma[mask] if sigma is not None else None

    # Initial guesses: alpha ~ 1, xi_corr ~ half the fit range, A so that
    # xi(r=1) ~ xi[smallest r] (roughly).
    p0 = [
        float(xi_fit[0]) * r_fit[0] ** 1.0,
        1.0,
        max(r_fit.max() / 2.0, 2.0),
    ]
    try:
        popt, pcov = curve_fit(
            power_law_with_cutoff, r_fit, xi_fit, p0=p0,
            sigma=sig_fit, absolute_sigma=sig_fit is not None,
            maxfev=10_000,
            bounds=([1e-9, 0.0, 1e-3], [np.inf, 10.0, np.inf]),
        )
    except Exception:
        return None
    perr = np.sqrt(np.diag(pcov))
    return float(popt[0]), float(popt[1]), float(popt[2]), perr


def power_law_galaxy_xi(r: np.ndarray, r0: float = 5.0, gamma: float = 1.8) -> np.ndarray:
    """Standard galaxy two-point correlation reference form.

        xi(r) = (r / r_0) ** -gamma

    Defaults are the canonical Peebles 1980 / Davis & Peebles 1983 fit
    for SDSS-like galaxy samples on small scales (r_0 ~ 5 h^-1 Mpc,
    gamma ~ 1.8). r is in the same units as r_0; for a dimensionless
    comparison to our lattice xi(r), we fit r_0 as a free scale.
    """
    return (r / r0) ** (-gamma)
