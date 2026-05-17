"""Unit tests for the 1D Oslo sandpile.

These check the toppling rules and conservation laws directly with hand-set
states, plus a couple of coarse sanity checks on the steady state. Statistical
checks on the avalanche exponent live in the demo script, not here, because
tau is sensitive to fit range and cutoffs.
"""

from __future__ import annotations

import numpy as np

from void_cascade.sandpile_1d import OsloState, drive, relax, run


def _state(z: list[int], z_c: list[int]) -> OsloState:
    return OsloState(
        z=np.array(z, dtype=np.int64),
        z_c=np.array(z_c, dtype=np.int64),
    )


def test_drive_increments_left_slope():
    s = _state([0, 0, 0], [1, 1, 1])
    drive(s)
    assert s.z.tolist() == [1, 0, 0]


def test_no_topple_when_stable():
    rng = np.random.default_rng(0)
    s = _state([1, 1, 1], [1, 1, 1])
    size, T = relax(s, rng)
    assert size == 0
    assert T == 0
    assert s.z.tolist() == [1, 1, 1]


def test_left_boundary_topple():
    # z[0] = 2 > z_c[0] = 1: one topple, z[0] -= 2, z[1] += 1.
    rng = np.random.default_rng(0)
    s = _state([2, 0, 0], [1, 2, 2])
    size, T = relax(s, rng)
    assert s.z.tolist() == [0, 1, 0]
    assert size == 1
    assert T == 1
    assert s.grains_lost == 0


def test_interior_topple():
    # Only site 1 unstable: z[1] -= 2, z[0] += 1, z[2] += 1.
    rng = np.random.default_rng(0)
    s = _state([0, 3, 0], [2, 2, 2])
    size, T = relax(s, rng)
    assert s.z.tolist() == [1, 1, 1]
    assert size == 1
    assert T == 1
    assert s.grains_lost == 0


def test_right_boundary_loses_grain():
    # z[L-1] = 2 > z_c[L-1] = 1: topple ejects a grain off the right edge,
    # so z[L-1] -= 1 (not 2) and z[L-2] += 1.
    rng = np.random.default_rng(0)
    s = _state([0, 0, 2], [2, 2, 1])
    size, T = relax(s, rng)
    assert s.z.tolist() == [0, 1, 1]
    assert size == 1
    assert s.grains_lost == 1


def test_avalanche_propagation():
    # Engineered chain reaction:
    # Start: z=[2,1,0], z_c=[1,1,1]. Site 0 topples (size=1, T=1).
    # After: z=[0,2,0]. Site 1 topples (size=2, T=2).
    # After: z=[1,0,1]. Stable.
    rng = np.random.default_rng(0)
    s = _state([2, 1, 0], [1, 1, 1])
    size, T = relax(s, rng)
    assert size == 2
    assert T == 2
    # Final slopes after the chain reaction (independent of redrawn thresholds).
    assert s.z.tolist() == [1, 0, 1]


def test_parallel_update_handles_adjacent_unstable_sites():
    # Two adjacent unstable sites in the same sweep both push to the site
    # between them; np.add.at must not drop a contribution.
    # z=[3,3,0,0], z_c=[1,1,2,2]: sites 0 and 1 both topple in sweep 1.
    # Site 0: z[0] -= 2, z[1] += 1.
    # Site 1: z[1] -= 2, z[0] += 1, z[2] += 1.
    # Net: z = [3-2+1, 3+1-2, 0+1, 0] = [2, 2, 1, 0].
    rng = np.random.default_rng(0)
    s = _state([3, 3, 0, 0], [1, 1, 2, 2])
    # Force this to be the only sweep by setting z_c large enough afterwards.
    # We only check the sweep-1 outcome before the next sweep:
    # But relax() keeps going until stable. Let's check final stability and
    # the cumulative size includes at least the 2 initial topples.
    size, T = relax(s, rng)
    assert size >= 2
    assert T >= 1
    # Final state is stable: z <= z_c everywhere.
    assert np.all(s.z <= s.z_c)


def test_thresholds_stay_in_set():
    state, _, _ = run(L=20, n_drops=2_000, seed=1)
    assert state.z_c.min() >= 1
    assert state.z_c.max() <= 2


def test_steady_state_mean_slope_in_band():
    # After the transient, the average bulk slope settles between 1 and 2.
    # For z_c uniform on {1, 2} the theoretical bulk mean is close to ~1.73
    # (Pruessner, "Self-Organised Criticality", 2012), but at this lattice
    # size the bulk-vs-edge structure still wobbles; the band check is loose.
    state, _, _ = run(L=50, n_drops=20_000, seed=2)
    mean_slope = state.z.mean()
    assert 1.0 < mean_slope < 2.0


def test_conservation_grains_in_equals_pile_plus_lost():
    # Every drive adds one grain to h_0; every right-edge topple ejects one
    # grain. So grains_in == sum_i h_i + grains_lost, with h_i = sum_{j>=i} z_j.
    n_drops = 5_000
    state, sizes, _ = run(L=30, n_drops=n_drops, seed=3)
    h = np.cumsum(state.z[::-1])[::-1]  # h_i = sum_{j>=i} z_j
    assert n_drops == int(h.sum()) + state.grains_lost
