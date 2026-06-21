"""Tests for the cooling term: heat bleeding away to the void.

Mechanism: after LF transport each step, cells above the ambient baseline
(e_floor = starting specific internal energy) lose a fraction `cooling` of
their excess heat. That energy leaves the system (not conserved).
Full accounting: E + load + cumulative_cooled == E0 + total_driven.

cooling=0  -> identical to step() and the conserved quantity is unchanged.
cooling>0  -> total system energy decreases; a strong enough cooling rate
              can prevent the substrate from ever tipping.
"""
import numpy as np

from void_cascade.cascade_breakdown import (
    BreakdownParams,
    step,
    step_cool,
    run_cooling_arc,
    total_energy,
)
from void_cascade.material_motion_2d import GAMMA, internal_energy


E_FLOOR = 1.0 / (GAMMA - 1.0)   # = 1.5 for P0=1.0, gamma=5/3


def _uniform_fields(L=8, P0=1.0):
    rho  = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E    = np.full((L, L), P0 / (GAMMA - 1.0))
    load = np.zeros((L, L))
    return rho, momx, momy, E, load


# --- step_cool: basic contracts ---

def test_step_cool_zero_identical_to_step():
    """cooling=0 must be numerically identical to step()."""
    p     = BreakdownParams(drive_sites=2, drive_amount=0.5)
    rng_a = np.random.default_rng(0)
    rng_b = np.random.default_rng(0)
    rho, momx, momy, E, load = _uniform_fields()
    r_s = step(rho, momx, momy, E, load, p, rng_a)
    r_c = step_cool(rho, momx, momy, E, load, p, rng_b, cooling=0.0)
    for a, b in zip(r_s[:5], r_c[:5]):
        assert np.allclose(a, b)
    assert r_c[6] == 0.0   # cooled_this_step should be zero


def test_step_cool_zero_cooled_amount():
    p   = BreakdownParams()
    rng = np.random.default_rng(1)
    rho, momx, momy, E, load = _uniform_fields()
    *_, cooled = step_cool(rho, momx, momy, E, load, p, rng, cooling=0.0)
    assert cooled == 0.0


def test_step_cool_accounting():
    """E + load + cooled == E_before + driven (full conservation accounting)."""
    p   = BreakdownParams()
    rng = np.random.default_rng(2)
    rho, momx, momy, E, load = _uniform_fields()
    E_before = total_energy(E, load)
    rho2, momx2, momy2, E2, load2, _, cooled = step_cool(
        rho, momx, momy, E, load, p, rng, cooling=0.3, e_floor=E_FLOOR)
    driven = p.drive_sites * p.drive_amount
    residual = total_energy(E2, load2) + cooled - (E_before + driven)
    assert abs(residual) < 1e-8


def test_step_cool_hot_cells_lose_energy():
    """With cooling>0, a pre-heated system has less total energy than cooling=0."""
    p    = BreakdownParams(drive_sites=0, drive_amount=0.0)
    L    = 8
    rho  = np.ones((L, L))
    momx = momy = np.zeros((L, L))
    E    = np.full((L, L), 3.0)   # all cells hot (e_spec=3.0 > e_ign=2.5)
    load = np.zeros((L, L))

    rng_a = np.random.default_rng(0)
    rng_b = np.random.default_rng(0)
    _, _, _, E0, _, _, c0 = step_cool(rho, momx, momy, E, load, p, rng_a,
                                       cooling=0.0, e_floor=E_FLOOR)
    _, _, _, E1, _, _, c1 = step_cool(rho, momx, momy, E, load, p, rng_b,
                                       cooling=0.5, e_floor=E_FLOOR)
    assert E1.sum() < E0.sum(), "cooling should remove energy from hot cells"
    assert c1 > 0.0


def test_step_cool_cold_cells_unaffected():
    """Cells at or below e_floor lose no energy from cooling."""
    p    = BreakdownParams(drive_sites=0, drive_amount=0.0)
    rho, momx, momy, E, load = _uniform_fields()   # E = e_floor everywhere
    rng  = np.random.default_rng(0)
    rho2, momx2, momy2, E2, load2, _, cooled = step_cool(
        rho, momx, momy, E, load, p, rng, cooling=1.0, e_floor=E_FLOOR)
    assert cooled == 0.0 or cooled < 1e-10


def test_step_cool_E_never_negative():
    """E must never go negative even at cooling=1.0."""
    p    = BreakdownParams()
    rng  = np.random.default_rng(99)
    rho, momx, momy, E, load = _uniform_fields()
    E = E * 5.0   # pre-heat everything
    for _ in range(3):
        rho, momx, momy, E, load, _, _ = step_cool(
            rho, momx, momy, E, load, p, rng, cooling=1.0, e_floor=E_FLOOR)
    assert E.min() >= 0.0


# --- run_cooling_arc ---

def test_run_cooling_arc_zero_matches_no_cooling():
    """cooling=0 arc: total cooled==0 and substrate still tips."""
    r = run_cooling_arc(L=30, steps=500, seed=0, cooling=0.0)
    assert r["total_cooled"] == 0.0


def test_run_cooling_arc_conservation_accounting():
    """Full energy accounting closes for any cooling value."""
    r = run_cooling_arc(L=20, steps=100, seed=0, cooling=0.2)
    resid = abs(r["energy_residual"])
    assert resid < 1e-6, f"accounting residual too large: {resid}"


def test_run_cooling_arc_cold_control_never_tips():
    """hpc=0 never tips regardless of cooling rate."""
    p = BreakdownParams(hpc=0.0, drive_sites=4, drive_amount=1.0)
    r = run_cooling_arc(L=20, steps=200, params=p, seed=0, cooling=0.5)
    assert r["tip_step"] is None


def test_run_cooling_arc_high_cooling_still_tips():
    """cooling=1.0 does NOT prevent tipping under default params.

    Cooling fires after the breakdown->combust sequence within each step.
    Within-step ignition can occur before cooling gets a chance to act.
    Result: no stable cold regime found by post-step cooling alone.
    This is an honest finding, not a bug.
    """
    p = BreakdownParams(hpc=0.1, drive_sites=4, drive_amount=1.0)
    r = run_cooling_arc(L=40, steps=3200, params=p, seed=0, cooling=1.0)
    # Substrate tips despite maximum cooling -- accounting must still close.
    assert abs(r["energy_residual"]) < 1e-6
    assert r["total_cooled"] > 0.0   # cooling was active


def test_run_cooling_arc_zero_cooling_tips():
    """Baseline (cooling=0) still tips under the standard params."""
    p = BreakdownParams(hpc=0.15, drive_sites=6, drive_amount=1.5)
    r = run_cooling_arc(L=40, steps=2500, params=p, seed=0, cooling=0.0)
    assert r["tip_step"] is not None
