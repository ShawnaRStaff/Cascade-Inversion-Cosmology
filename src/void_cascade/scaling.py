"""Finite-size scaling for SOC avalanche statistics.

Standard scaling ansatz for an SOC system of linear size L:

    P(x; L) = x^{-a} f(x / L^B)

where x is either avalanche size s (with exponents a=tau, B=D) or duration T
(with exponents a=alpha, B=z, the dynamical exponent). f is a universal
scaling function.

We extract (a, B) by the moment method (Tebaldi, De Menech & Stella 1999;
Pruessner, "Self-Organised Criticality", 2012, Chapter 7):

    <x^k> ~ L^{B * (1 + k - a)}    for k > a - 1

Fitting log<x^k> vs log L for several k gives slopes sigma_k = B(1+k-a). The
slopes are linear in k: sigma_k = B*k - B*(a-1). A single straight-line fit
of sigma_k vs k then yields B (slope) and a = 1 - intercept/B.

For 1D Oslo the reference values are tau ~ 1.557, D ~ 2.245 (Pruessner &
Jensen 2003), alpha ~ 1.76, z ~ 1.42. The same code path is used to fit
either the (tau, D) pair from sizes or (alpha, z) from durations.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .analysis import log_binned_pdf


@dataclass
class MomentScalingResult:
    """Result of a moment-scaling fit for a single observable.

    Attributes
    ----------
    k_values : ndarray of float
        Moment orders used in the fit.
    L_values : ndarray of float
        Lattice sizes used, in ascending order.
    moments : ndarray, shape (n_k, n_L)
        Empirical <x^k> for each (k, L). Computed from positive samples only.
    sigma_k : ndarray, shape (n_k,)
        Fitted slopes of log10<x^k> vs log10(L), one per moment.
    sigma_k_err : ndarray, shape (n_k,)
        Standard errors on sigma_k from the per-moment linear fits.
    cutoff_exponent : float
        Either D (for sizes) or z (for durations). The slope of sigma_k vs k.
    cutoff_exponent_err : float
    tail_exponent : float
        Either tau (for sizes) or alpha (for durations). Inferred via
        tail_exponent = 1 - intercept / cutoff_exponent.
    tail_exponent_err : float
    """

    k_values: np.ndarray
    L_values: np.ndarray
    moments: np.ndarray
    sigma_k: np.ndarray
    sigma_k_err: np.ndarray
    cutoff_exponent: float
    cutoff_exponent_err: float
    tail_exponent: float
    tail_exponent_err: float


def moment_scaling(
    samples_by_L: dict[int, np.ndarray],
    k_values: list[int] | None = None,
) -> MomentScalingResult:
    """Extract (tail_exponent, cutoff_exponent) from samples at several L.

    Parameters
    ----------
    samples_by_L : dict[L -> array]
        Steady-state samples (sizes or durations) for each lattice size L.
        Zero-valued samples are discarded; their fraction is not part of the
        scaling ansatz and would only contaminate moments at small k.
    k_values : list of int, optional
        Moment orders. Defaults to [1, 2, 3, 4]. Must satisfy k > a-1 for each
        k, so the moment is dominated by the cutoff region rather than the
        lower bound; for tau ~ 1.55 this is satisfied by k >= 1.

    Raises
    ------
    ValueError
        If fewer than three lattice sizes are provided (the linear fit
        sigma_k vs L is underconstrained otherwise), or if any L has no
        positive samples.
    """
    if k_values is None:
        k_values_arr = np.array([1, 2, 3, 4], dtype=float)
    else:
        k_values_arr = np.array(k_values, dtype=float)

    if len(samples_by_L) < 3:
        raise ValueError("Need samples for at least 3 lattice sizes.")

    L_values = np.array(sorted(samples_by_L.keys()), dtype=float)
    log_L = np.log10(L_values)

    moments = np.zeros((len(k_values_arr), len(L_values)))
    for j, L in enumerate(L_values):
        x = np.asarray(samples_by_L[int(L)], dtype=np.float64)
        x = x[x > 0]
        if x.size == 0:
            raise ValueError(f"No positive samples for L={int(L)}")
        for i, k in enumerate(k_values_arr):
            moments[i, j] = float(np.mean(x ** k))

    log_moments = np.log10(moments)

    # Per-moment slope sigma_k = d log<x^k> / d log L.
    sigma_k = np.zeros(len(k_values_arr))
    sigma_k_err = np.zeros(len(k_values_arr))
    for i in range(len(k_values_arr)):
        coeffs, cov = np.polyfit(log_L, log_moments[i], 1, cov=True)
        sigma_k[i] = float(coeffs[0])
        sigma_k_err[i] = float(np.sqrt(cov[0, 0]))

    # sigma_k is linear in k: sigma_k = B*k - B*(a-1).
    # numpy's polyfit takes weights as 1/sigma (not 1/sigma^2).
    weights = 1.0 / sigma_k_err
    coeffs, cov = np.polyfit(k_values_arr, sigma_k, 1, w=weights, cov=True)
    B = float(coeffs[0])
    intercept = float(coeffs[1])
    B_err = float(np.sqrt(cov[0, 0]))
    intercept_err = float(np.sqrt(cov[1, 1]))
    cov_B_intercept = float(cov[0, 1])

    # a = 1 - intercept / B. Standard propagation, keeping the off-diagonal.
    a = 1.0 - intercept / B
    a_var = (
        (1.0 / B) ** 2 * intercept_err ** 2
        + (intercept / B ** 2) ** 2 * B_err ** 2
        - 2.0 * (intercept / B ** 3) * cov_B_intercept
    )
    a_err = float(np.sqrt(max(a_var, 0.0)))

    return MomentScalingResult(
        k_values=k_values_arr,
        L_values=L_values,
        moments=moments,
        sigma_k=sigma_k,
        sigma_k_err=sigma_k_err,
        cutoff_exponent=B,
        cutoff_exponent_err=B_err,
        tail_exponent=a,
        tail_exponent_err=a_err,
    )


def collapsed_distribution(
    x: np.ndarray,
    L: float,
    tail_exponent: float,
    cutoff_exponent: float,
    n_bins: int = 30,
) -> tuple[np.ndarray, np.ndarray]:
    """Return (x / L^B, x^a P(x)) for a finite-size data-collapse plot.

    If the scaling ansatz P(x; L) = x^{-a} f(x/L^B) holds, plotting the second
    component vs the first should produce the same curve f(.) for every L,
    visually confirming the fitted exponents.
    """
    centers, pdf, _ = log_binned_pdf(x, n_bins=n_bins)
    x_scaled = centers / (L ** cutoff_exponent)
    y_scaled = centers ** tail_exponent * pdf
    return x_scaled, y_scaled
