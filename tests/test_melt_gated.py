"""Tests for melt-gated rigidity: cold cells absorb kinetic energy as load.

Mechanism: after the LF fluid step, cold cells convert their kinetic energy
into mechanical stress (load) rather than keeping it as motion. Hot cells flow
freely. Cold cells are stiff walls that get loaded up — which can then trigger
more breakdown in the next step, potentially intensifying the cascade.

Conservation is exact: KE removed from E, added to load, E+load unchanged.
melt_frac=0  -> identical to step().
melt_frac=1  -> cold cells fully absorb all post-LF kinetic energy into load.
"""
import numpy as np
import pytest

from void_cascade.cascade_breakdown import (
    BreakdownParams,
    step,
    step_melt,
    run_melt_arc,
    total_energy,
)
from void_cascade.material_motion_2d import GAMMA, internal_energy


# --- helpers ---

def _hot_cold_fields(L=8, hot_val=5.0, cold_val=1.0, P0=1.0):
    """Hot centre cell, everything else cold and still."""
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), P0 / (GAMMA - 1.0))
    load = np.zeros((L, L))
    E[L // 2, L // 2] = hot_val        # one hot cell in the centre
    return rho, momx, momy, E, load


# --- step_melt: basic contracts ---

def test_step_melt_zero_identical_to_step():
    """melt_frac=0 must be numerically identical to step()."""
    p   = BreakdownParams(drive_sites=2, drive_amount=0.5)
    rng_a = np.random.default_rng(42)
    rng_b = np.random.default_rng(42)
    rho, momx, momy, E, load = _hot_cold_fields()
    r_s  = step(rho, momx, momy, E, load, p, rng_a)
    r_m  = step_melt(rho, momx, momy, E, load, p, rng_b, melt_frac=0.0)
    for a, b in zip(r_s[:5], r_m[:5]):        # rho, momx, momy, E, load
        assert np.allclose(a, b), "melt_frac=0 must equal step()"


def test_step_melt_conserves_E_plus_load_zero():
    p   = BreakdownParams()
    rng = np.random.default_rng(7)
    rho, momx, momy, E, load = _hot_cold_fields()
    before = total_energy(E, load)
    rho2, momx2, momy2, E2, load2, _ = step_melt(
        rho, momx, momy, E, load, p, rng, melt_frac=0.0)
    after = total_energy(E2, load2)
    assert abs(after - before - p.drive_sites * p.drive_amount) < 1e-8


def test_step_melt_conserves_E_plus_load_half():
    p   = BreakdownParams()
    rng = np.random.default_rng(7)
    rho, momx, momy, E, load = _hot_cold_fields()
    before = total_energy(E, load)
    rho2, momx2, momy2, E2, load2, _ = step_melt(
        rho, momx, momy, E, load, p, rng, melt_frac=0.5)
    after = total_energy(E2, load2)
    assert abs(after - before - p.drive_sites * p.drive_amount) < 1e-8


def test_step_melt_conserves_E_plus_load_full():
    p   = BreakdownParams()
    rng = np.random.default_rng(7)
    rho, momx, momy, E, load = _hot_cold_fields()
    before = total_energy(E, load)
    rho2, momx2, momy2, E2, load2, _ = step_melt(
        rho, momx, momy, E, load, p, rng, melt_frac=1.0)
    after = total_energy(E2, load2)
    assert abs(after - before - p.drive_sites * p.drive_amount) < 1e-8


def test_step_melt_cold_cells_momentum_reduced():
    """melt_frac=1.0: cold cells must end up with less kinetic energy than melt_frac=0."""
    p    = BreakdownParams(drive_sites=0, drive_amount=0.0)   # no drive, isolate effect
    L    = 8
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    # Mild hot cell: E=3.0 (just above e_ign=2.5) on cool background E=0.5.
    # After LF averaging, the neighbour gets E_avg=(3.0+0.5+0.5+0.5)/4=1.125 < 2.5
    # so it stays cold and acquires nonzero KE that the melt gate can absorb.
    E    = np.full((L, L), 0.5)
    E[L//2, L//2] = 3.0
    load = np.zeros((L, L))

    rng_a = np.random.default_rng(0)
    rng_b = np.random.default_rng(0)
    _, mx0, my0, E0, ld0, _ = step_melt(rho, momx, momy, E, load,
                                          p, rng_a, melt_frac=0.0)
    _, mx1, my1, E1, ld1, _ = step_melt(rho, momx, momy, E, load,
                                          p, rng_b, melt_frac=1.0)

    e_spec0 = internal_energy(rho, mx0, my0, E0) / rho
    cold0   = e_spec0 < 2.5
    KE0_cold = (0.5 * (mx0**2 + my0**2) / rho * cold0).sum()

    e_spec1 = internal_energy(rho, mx1, my1, E1) / rho
    cold1   = e_spec1 < 2.5
    KE1_cold = (0.5 * (mx1**2 + my1**2) / rho * cold1).sum()

    assert KE1_cold < KE0_cold, (
        "melt_frac=1 cold cells should have less KE than melt_frac=0"
    )


def test_step_melt_cold_KE_transfers_to_load():
    """Cold cells' absorbed KE ends up in load (not lost)."""
    p    = BreakdownParams(drive_sites=0, drive_amount=0.0)
    L    = 8
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), 0.5)     # cool background stays cold after LF averaging
    E[L//2, L//2] = 3.0             # mild hot centre (just above e_ign=2.5)
    load = np.zeros((L, L))

    rng_a = np.random.default_rng(1)
    rng_b = np.random.default_rng(1)
    _, mx0, my0, E0, ld0, _ = step_melt(rho, momx, momy, E, load,
                                          p, rng_a, melt_frac=0.0)
    _, mx1, my1, E1, ld1, _ = step_melt(rho, momx, momy, E, load,
                                          p, rng_b, melt_frac=1.0)

    # melt_frac=1 should have more load than melt_frac=0
    assert ld1.sum() > ld0.sum(), "absorbed cold KE should increase total load"
    # And E should be lower by the same amount
    delta_load = ld1.sum() - ld0.sum()
    delta_E    = E1.sum() - E0.sum()
    assert abs(delta_load + delta_E) < 1e-8, "load gain == E loss"


def test_step_melt_E_never_negative():
    """E must stay >= 0 in every cell after step_melt (no overcooling)."""
    p   = BreakdownParams()
    rng = np.random.default_rng(99)
    rho, momx, momy, E, load = _hot_cold_fields(L=16)
    for _ in range(5):
        rho, momx, momy, E, load, _ = step_melt(
            rho, momx, momy, E, load, p, rng, melt_frac=1.0)
    assert E.min() >= 0.0, "E must never go negative"


# --- run_melt_arc: arc-level contracts ---

def test_run_melt_arc_cold_control_never_tips():
    """hpc=0 with melt_frac=1: no heat from breakdown -> never tips."""
    p = BreakdownParams(hpc=0.0, drive_amount=1.0, drive_sites=4)
    r = run_melt_arc(L=20, steps=300, params=p, seed=0, melt_frac=1.0)
    assert r["tip_step"] is None
    assert r["max_temp"] < 2.5


def test_run_melt_arc_tips_after_quiet():
    """Fast drive with melt_frac=1.0: substrate tips after a quiet period, not at step 0."""
    p = BreakdownParams(hpc=0.2, drive_amount=2.0, drive_sites=8)
    r = run_melt_arc(L=20, steps=800, params=p, seed=0, melt_frac=1.0)
    assert r["tip_step"] is not None, "should tip under fast drive"
    assert r["tip_step"] > 3, "must be after a quiet period, not instant"
    assert r["max_temp"] >= 2.5


def test_run_melt_arc_conserves():
    """Energy conserved throughout a full melt arc."""
    r = run_melt_arc(L=30, steps=200, seed=0, melt_frac=1.0)
    assert abs(r["energy_residual"]) < 1e-6


def test_run_melt_arc_intensification_differs_from_no_melt():
    """The n_hot curve after tip should differ between melt_frac=0 and 1.0.

    We don't prescribe which is faster or slower — the model decides.
    We only assert the curves are measurably different (the latch has an effect).
    """
    p = BreakdownParams(hpc=0.15, drive_amount=1.5, drive_sites=6)
    r0 = run_melt_arc(L=40, steps=2500, params=p, seed=0, melt_frac=0.0)
    r1 = run_melt_arc(L=40, steps=2500, params=p, seed=0, melt_frac=1.0)
    # Both must tip
    assert r0["tip_step"] is not None and r1["tip_step"] is not None
    # The n_hot traces must differ at some point (melt latch has a measurable effect)
    nh0 = np.array(r0["n_hot_trace"])
    nh1 = np.array(r1["n_hot_trace"])
    # Use the shorter trace for comparison
    n = min(len(nh0), len(nh1))
    assert not np.allclose(nh0[:n], nh1[:n]), (
        "melt latch should produce a measurably different n_hot curve"
    )
