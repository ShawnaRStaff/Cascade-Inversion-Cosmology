"""Tests for the buildup -> breakdown -> tip engine (cascade_breakdown).

These pin the behaviour we earned through the probes, so it can never silently
regress: exact conservation, a cold control that never tips (not rigged), a slow
drive that can fizzle, and a fast drive that tips to plasma AFTER a quiet period.
"""
import numpy as np

from void_cascade.cascade_breakdown import (
    BreakdownParams, breakdown, combust, drive, run_buildup_tip, total_energy,
)


def test_drive_adds_exact_stress_and_is_pure():
    load = np.zeros((4, 4))
    rng = np.random.default_rng(0)
    out = drive(load, rng, n_sites=3, amount=2.0)
    assert np.isclose(out.sum(), 6.0)      # 3 sites * 2.0 stress
    assert load.sum() == 0.0               # original untouched (pure)


def test_breakdown_conserves_and_relaxes():
    L = 8
    load = np.zeros((L, L)); load[4, 4] = 5.0
    rho = np.ones((L, L)); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)
    before = total_energy(E, load)
    load2, E2 = breakdown(load, rho, momx, momy, E, thr=2.0, hpc=0.1, e_ign=2.5)
    assert np.isclose(before, total_energy(E2, load2))   # load+E conserved
    assert load2[4, 4] < 5.0                              # it shed stress
    assert load2.max() < 2.0                              # relaxed below threshold


def test_breakdown_zero_hpc_makes_no_heat():
    """Cold control at the unit level: no stress->heat => fluid energy untouched."""
    L = 8
    load = np.zeros((L, L)); load[4, 4] = 5.0
    rho = np.ones((L, L)); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)
    load2, E2 = breakdown(load, rho, momx, momy, E, thr=2.0, hpc=0.0, e_ign=2.5)
    assert np.allclose(E2, E)                       # no heat created
    assert np.isclose(load2.sum(), load.sum())      # stress only redistributed


def test_combust_releases_when_hot_conserving():
    L = 4
    rho = np.ones((L, L)); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), 5.0)              # specific internal energy 5 > e_ign
    load = np.full((L, L), 1.0)
    E2, load2, ig = combust(rho, momx, momy, E, load, e_ign=2.5)
    assert ig is True
    assert np.isclose(load2.sum(), 0.0)
    assert np.isclose(E2.sum(), E.sum() + load.sum())   # conserved


def test_combust_quiet_when_cold():
    L = 4
    rho = np.ones((L, L)); momx = np.zeros((L, L)); momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)             # cold
    load = np.full((L, L), 1.0)
    E2, load2, ig = combust(rho, momx, momy, E, load, e_ign=2.5)
    assert ig is False
    assert np.allclose(E2, E) and np.allclose(load2, load)


def test_run_conserves_energy_exactly():
    r = run_buildup_tip(L=40, steps=200, seed=0)
    assert abs(r["energy_residual"]) < 1e-6


def test_cold_control_never_tips():
    """No stress->heat (hpc=0) => can never reach plasma, however hard we drive."""
    p = BreakdownParams(hpc=0.0, drive_amount=1.0, drive_sites=4)
    r = run_buildup_tip(L=40, steps=400, params=p, seed=0)
    assert r["tip_step"] is None
    assert r["max_temp"] < 2.5


def test_slow_drive_can_fizzle():
    """Guardrail: a slow enough drive must be ABLE to stay cold (not tuned to win)."""
    p = BreakdownParams(hpc=0.1, drive_amount=0.1, drive_sites=1)
    r = run_buildup_tip(L=40, steps=400, params=p, seed=0)
    assert r["tip_step"] is None


def test_fast_drive_tips_to_plasma_after_a_quiet():
    """Fast enough drive: sits quiet, then tips itself to plasma -- force-free."""
    p = BreakdownParams(hpc=0.2, drive_amount=2.0, drive_sites=8)
    r = run_buildup_tip(L=40, steps=800, params=p, seed=0)
    assert r["tip_step"] is not None        # it tips
    assert r["max_temp"] >= 2.5             # reaches plasma
    assert r["tip_step"] > 20               # only AFTER a quiet period, not instant
