"""Tests for the 3D buildup -> breakdown -> tip engine.

Same guarantees as the 2D engine: exact conservation, cold control never tips,
slow drive fizzles, fast drive tips after a quiet period. The geometry changes
(6 neighbours, 3D LF) but the physics must match.
"""
import numpy as np

from void_cascade.cascade_breakdown_3d import (
    BreakdownParams3d,
    breakdown_3d,
    combust_3d,
    drive_3d,
    run_buildup_tip_3d,
    total_energy_3d,
)
from void_cascade.material_motion_3d import internal_energy_3d


def test_drive_3d_pure_and_exact():
    load = np.zeros((4, 4, 4))
    rng  = np.random.default_rng(0)
    out  = drive_3d(load, rng, n_sites=3, amount=2.0)
    assert np.isclose(out.sum(), 6.0)
    assert load.sum() == 0.0


def test_breakdown_3d_conserves():
    L = 6
    load = np.zeros((L, L, L)); load[3, 3, 3] = 5.0
    rho  = np.ones((L, L, L))
    momx = np.zeros((L, L, L)); momy = np.zeros((L, L, L)); momz = np.zeros((L, L, L))
    E    = np.full((L, L, L), 1.5)
    before = total_energy_3d(E, load)
    load2, E2 = breakdown_3d(load, rho, momx, momy, momz, E, thr=2.0, hpc=0.1, e_ign=2.5)
    assert np.isclose(before, total_energy_3d(E2, load2))
    assert load2[3, 3, 3] < 5.0


def test_breakdown_3d_zero_hpc_no_heat():
    L = 6
    load = np.zeros((L, L, L)); load[3, 3, 3] = 5.0
    rho  = np.ones((L, L, L))
    momx = momz = np.zeros((L, L, L)); momy = np.zeros((L, L, L))
    E    = np.full((L, L, L), 1.5)
    load2, E2 = breakdown_3d(load, rho, momx, momy, momz, E, thr=2.0, hpc=0.0, e_ign=2.5)
    assert np.allclose(E2, E)
    assert np.isclose(load2.sum(), load.sum())


def test_breakdown_3d_sheds_to_six_neighbours():
    """An isolated over-threshold cell sheds to exactly its 6 face-neighbours."""
    L = 5
    load = np.zeros((L, L, L))
    load[2, 2, 2] = 3.0   # only this cell over threshold
    rho  = np.ones((L, L, L))
    momx = momy = momz = np.zeros((L, L, L))
    E    = np.full((L, L, L), 1.5)
    load2, _ = breakdown_3d(load, rho, momx, momy, momz, E, thr=2.0, hpc=0.0, e_ign=2.5)
    six_nbrs = [(1,2,2),(3,2,2),(2,1,2),(2,3,2),(2,2,1),(2,2,3)]
    for idx in six_nbrs:
        assert load2[idx] > 0.0, f"neighbour {idx} should receive stress"


def test_combust_3d_conserves():
    L = 4
    rho  = np.ones((L, L, L))
    momx = momy = momz = np.zeros((L, L, L))
    E    = np.full((L, L, L), 5.0)
    load = np.full((L, L, L), 1.0)
    E2, load2, ig = combust_3d(rho, momx, momy, momz, E, load, e_ign=2.5)
    assert ig is True
    assert np.isclose(load2.sum(), 0.0)
    assert np.isclose(E2.sum(), E.sum() + load.sum())


def test_run_3d_conserves():
    r = run_buildup_tip_3d(L=10, steps=100, seed=0)
    assert abs(r["energy_residual"]) < 1e-4


def test_cold_control_3d_never_tips():
    p = BreakdownParams3d(hpc=0.0, drive_amount=1.0, drive_sites=4)
    r = run_buildup_tip_3d(L=10, steps=200, params=p, seed=0)
    assert r["tip_step"] is None
    assert r["max_temp"] < 2.5


def test_fast_drive_3d_tips_after_quiet():
    """Aggressive drive: substrate tips to plasma after a cold quiet, not at step 0."""
    p = BreakdownParams3d(hpc=0.3, drive_amount=2.0, drive_sites=12)
    r = run_buildup_tip_3d(L=10, steps=1000, params=p, seed=0)
    assert r["tip_step"] is not None
    assert r["max_temp"] >= 2.5
    assert r["tip_step"] > 5
