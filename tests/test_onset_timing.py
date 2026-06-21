"""Tests for onset-timing observables.

Three onset signals clocked on a shared time axis:
  (1) n_over:      cells above load threshold -- proxy for avalanche activity
  (2) corr_length: spatial correlation length of the load field -- spatial structure
  (3) n_hot:       cells above ignition temperature -- ignition / plasma

correlation_length() is a pure function (FFT ACF, no state).
run_onset_measurement() drives the simulation and samples all three.
"""
import numpy as np
import pytest

from void_cascade.cascade_breakdown import (
    BreakdownParams,
    correlation_length,
    run_onset_measurement,
)
# Note: breakdown always clears cells above threshold, so n_over is always 0
# post-step. The honest proxy for avalanche activity is load.std() (stress
# clustering). run_onset_measurement returns load_std, not n_over.

_FAST = BreakdownParams(hpc=0.2, drive_amount=2.0, drive_sites=8)


# ---------------------------------------------------------------------------
# correlation_length
# ---------------------------------------------------------------------------


def test_correlation_length_uniform_field():
    """Uniform field has no spatial variation: correlation length = 0."""
    assert correlation_length(np.ones((16, 16))) == 0.0


def test_correlation_length_broad_greater_than_sharp():
    """A smooth broad blob has a longer correlation length than a sharp peak."""
    L = 32
    cy, cx = np.mgrid[0:L, 0:L] - L // 2
    sharp = np.zeros((L, L)); sharp[L//2, L//2] = 1.0
    broad = np.exp(-(cy**2 + cx**2) / 25.0)
    assert correlation_length(broad) > correlation_length(sharp)


def test_correlation_length_grows_with_blob_width():
    """Wider blobs have longer correlation lengths."""
    L = 32
    cy, cx = np.mgrid[0:L, 0:L] - L // 2
    narrow = np.exp(-(cy**2 + cx**2) / 4.0)
    wide   = np.exp(-(cy**2 + cx**2) / 36.0)
    assert correlation_length(wide) > correlation_length(narrow)


def test_correlation_length_non_negative():
    """Correlation length is always >= 0 for any field."""
    rng = np.random.default_rng(0)
    for _ in range(5):
        f = rng.random((20, 20))
        assert correlation_length(f) >= 0.0


# ---------------------------------------------------------------------------
# run_onset_measurement
# ---------------------------------------------------------------------------


def test_run_onset_required_keys():
    r = run_onset_measurement(L=16, steps=80, params=_FAST, seed=0)
    for key in ("tip_step", "t_axis", "load_std", "corr_lengths", "n_hot"):
        assert key in r


def test_run_onset_series_lengths_match():
    r = run_onset_measurement(L=16, steps=80, params=_FAST, seed=0)
    n = len(r["t_axis"])
    assert len(r["load_std"])     == n
    assert len(r["corr_lengths"]) == n
    assert len(r["n_hot"])        == n


def test_run_onset_cold_control_no_tip():
    """hpc=0: no stress->heat, no tip, n_hot stays zero."""
    p = BreakdownParams(hpc=0.0, drive_amount=2.0, drive_sites=8)
    r = run_onset_measurement(L=16, steps=100, params=p, seed=0)
    assert r["tip_step"] is None
    assert all(v == 0 for v in r["n_hot"])


def test_run_onset_load_std_nonzero():
    """load.std() > 0 at some point during the buildup -- stress clusters before tip.

    On very small grids with aggressive drive, combust may clear all load post-tip
    so load_std at the final sample can be 0.  The honest check is that it was
    nonzero SOMEWHERE during the run (stress clustering did happen).
    On L=40 with default slow drive the system stays in buildup for all 300 steps
    and load_std rises monotonically.
    """
    r = run_onset_measurement(L=40, steps=300, seed=0)  # default (slow) params
    assert any(v > 0 for v in r["load_std"])
