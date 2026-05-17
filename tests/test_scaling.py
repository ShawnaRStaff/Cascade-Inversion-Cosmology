"""Tests for moment-scaling exponent extraction.

The Oslo simulation is the realistic test, but it's slow and noisy. These
tests use synthetic truncated power-law samples with known (a, B) to verify
that moment_scaling recovers the correct values within statistical error.
"""

from __future__ import annotations

import numpy as np
import pytest

from void_cascade.scaling import collapsed_distribution, moment_scaling


def _sample_truncated_power_law(
    a: float, cutoff: float, n: int, rng: np.random.Generator
) -> np.ndarray:
    """Draw n samples from P(x) ~ x^{-a} on [1, cutoff], a > 1.

    Inverse-CDF sampling: F(x) = (1 - x^{1-a}) / (1 - cutoff^{1-a}) for
    x in [1, cutoff].
    """
    u = rng.uniform(size=n)
    norm = 1.0 - cutoff ** (1.0 - a)
    return (1.0 - u * norm) ** (1.0 / (1.0 - a))


def test_moment_scaling_recovers_synthetic_exponents():
    a_true = 1.55
    B_true = 2.0
    rng = np.random.default_rng(42)
    samples = {}
    for L in (32, 64, 128, 256):
        cutoff = float(L) ** B_true
        samples[L] = _sample_truncated_power_law(a_true, cutoff, 200_000, rng)

    r = moment_scaling(samples, k_values=[1, 2, 3, 4])

    assert abs(r.tail_exponent - a_true) < 0.05, (
        f"recovered tail exponent {r.tail_exponent:.3f} vs true {a_true}"
    )
    assert abs(r.cutoff_exponent - B_true) < 0.1, (
        f"recovered cutoff exponent {r.cutoff_exponent:.3f} vs true {B_true}"
    )


def test_moment_scaling_requires_three_or_more_L():
    rng = np.random.default_rng(0)
    samples = {
        L: _sample_truncated_power_law(1.5, float(L) ** 2, 1000, rng)
        for L in (32, 64)
    }
    with pytest.raises(ValueError):
        moment_scaling(samples)


def test_moment_scaling_rejects_empty_samples():
    samples = {
        32: np.zeros(100, dtype=np.int64),
        64: np.array([1, 2, 3]),
        128: np.array([1, 2, 3]),
    }
    with pytest.raises(ValueError):
        moment_scaling(samples)


def test_collapsed_distribution_shapes():
    rng = np.random.default_rng(0)
    x = _sample_truncated_power_law(1.5, 1000.0, 5000, rng)
    x_scaled, y_scaled = collapsed_distribution(
        x, L=128.0, tail_exponent=1.5, cutoff_exponent=2.0
    )
    assert x_scaled.shape == y_scaled.shape
    assert x_scaled.size > 0
