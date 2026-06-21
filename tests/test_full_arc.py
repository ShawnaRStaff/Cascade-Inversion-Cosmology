"""Tests for the full arc (continuous run) and the pressure-rigidity latch.

Full arc: cold substrate -> slow buildup -> tip -> eruption front, all in ONE
run with no seam between phases. The engine (cascade_breakdown) just keeps
running; we track what happens after the tip.

Rigidity latch: cold loaded cells express their stored load as elastic pressure
during LF transport (rigidity > 0). This is the pressure-resisting kind, NOT
velocity damping. When a cold cell crosses the ignition threshold during that
transport, it keeps the borrowed energy as thermal (crystal melts). With
rigidity=0 the result is numerically identical to plain step().
"""
import numpy as np
import pytest

from void_cascade.cascade_breakdown import (
    BreakdownParams,
    run_buildup_tip,
    run_full_arc,
    step,
    step_rigid,
    total_energy,
)

# Aggressive params: tips within ~100 steps on L=20; ~600 steps on L=40.
_FAST = BreakdownParams(hpc=0.2, drive_amount=2.0, drive_sites=8)


# ---------------------------------------------------------------------------
# Continuous arc
# ---------------------------------------------------------------------------


def test_full_arc_conservation():
    """Energy (E + load) is conserved across the full arc."""
    r = run_full_arc(L=20, steps_buildup=200, steps_after=80, params=_FAST, seed=0)
    assert abs(r["energy_residual"]) < 1e-4


def test_full_arc_tips_and_front_grows():
    """After the tip, the number of plasma cells grows -- the front is propagating."""
    r = run_full_arc(L=40, steps_buildup=700, steps_after=400, params=_FAST, seed=0)
    assert r["tip_step"] is not None, "FAST params must tip within 700 steps on L=40"
    post = [s for s in r["snapshots"] if s["t"] >= r["tip_step"]]
    assert len(post) >= 2
    assert post[-1]["n_hot"] >= post[0]["n_hot"]


def test_full_arc_cold_control_never_erupts():
    """hpc=0: no stress->heat, no tip, no plasma cells at any snapshot."""
    p = BreakdownParams(hpc=0.0, drive_amount=2.0, drive_sites=8)
    r = run_full_arc(L=20, steps_buildup=150, steps_after=50, params=p, seed=0)
    assert r["tip_step"] is None
    assert all(s["n_hot"] == 0 for s in r["snapshots"])


def test_full_arc_snapshots_have_required_keys():
    r = run_full_arc(L=20, steps_buildup=150, steps_after=50, params=_FAST, seed=0)
    for s in r["snapshots"]:
        for key in ("t", "n_hot", "max_temp", "max_speed", "mean_load", "phase"):
            assert key in s


def test_full_arc_buildup_tip_residuals_agree():
    """run_full_arc with steps_after=0 must agree with run_buildup_tip."""
    r1 = run_full_arc(L=20, steps_buildup=150, steps_after=0, params=_FAST, seed=7)
    r2 = run_buildup_tip(L=20, steps=150, params=_FAST, seed=7)
    assert abs(r1["energy_residual"] - r2["energy_residual"]) < 1e-6


# ---------------------------------------------------------------------------
# Rigidity latch
# ---------------------------------------------------------------------------


def test_rigid_zero_matches_plain_step():
    """step_rigid(rigidity=0) is numerically identical to step()."""
    L = 8
    rho = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)
    load = np.zeros((L, L))
    load[4, 4] = 3.0
    p = BreakdownParams()
    a = step(
        rho.copy(), momx.copy(), momy.copy(), E.copy(), load.copy(),
        p, np.random.default_rng(7),
    )
    b = step_rigid(
        rho.copy(), momx.copy(), momy.copy(), E.copy(), load.copy(),
        p, np.random.default_rng(7), rigidity=0.0,
    )
    for arr_a, arr_b in zip(a[:5], b[:5]):
        assert np.allclose(arr_a, arr_b), "rigidity=0 must match plain step exactly"


def test_rigid_step_conserves():
    """step_rigid conserves E + load for any rigidity in [0, 1]."""
    L = 16
    rho = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)
    rng0 = np.random.default_rng(0)
    load = rng0.uniform(0.0, 1.5, (L, L))
    p = BreakdownParams()
    driven = p.drive_sites * p.drive_amount
    before = total_energy(E, load)
    _, _, _, E2, load2, _ = step_rigid(
        rho, momx, momy, E, load, p, np.random.default_rng(1), rigidity=0.5,
    )
    assert abs(total_energy(E2, load2) - (before + driven)) < 1e-6


def test_full_arc_with_rigidity_conserves():
    """run_full_arc with rigidity > 0 still conserves E + load."""
    r = run_full_arc(
        L=20, steps_buildup=200, steps_after=80, params=_FAST, seed=0, rigidity=0.5,
    )
    assert abs(r["energy_residual"]) < 1e-3


def test_rigid_e_stays_non_negative():
    """After step_rigid, internal energy is non-negative in every cell."""
    L = 16
    rho = np.ones((L, L))
    momx = np.zeros((L, L))
    momy = np.zeros((L, L))
    E = np.full((L, L), 1.5)
    load = np.random.default_rng(0).uniform(0.0, 3.0, (L, L))
    p = BreakdownParams()
    _, _, _, E2, _, _ = step_rigid(
        rho, momx, momy, E, load, p, np.random.default_rng(0), rigidity=1.0,
    )
    from void_cascade.material_motion_2d import internal_energy
    e_int = internal_energy(rho, momx, momy, E2)
    assert float(e_int.min()) >= -1e-10
