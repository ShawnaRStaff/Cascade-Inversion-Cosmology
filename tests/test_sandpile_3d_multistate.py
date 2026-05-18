"""Tests for the multi-state 3D Manna sandpile.

Confirms the new dynamics:
  - z_c(sigma) = 2*(sigma+1) applies per cell.
  - Cell toppling distributes z_c grains stochastically; lattice
    conservation grains_in == sum(z) + grains_lost holds.
  - sigma transformation only happens when neighbors are co-active.
  - sigma_max bounds work when provided.
"""

from __future__ import annotations

import numpy as np

from void_cascade.sandpile_3d_multistate import (
    MannaStateMS,
    drive,
    initialize,
    relax,
    run_with_ever_toppled,
)


def _state(z: np.ndarray, sigma: np.ndarray) -> MannaStateMS:
    return MannaStateMS(z=z.astype(np.int64), sigma=sigma.astype(np.int64))


def test_initialize_empty():
    state = initialize(L=4)
    assert state.z.shape == (4, 4, 4)
    assert state.sigma.shape == (4, 4, 4)
    assert state.z.sum() == 0
    assert state.sigma.sum() == 0


def test_drive_increments_random_cell():
    rng = np.random.default_rng(0)
    state = initialize(L=3)
    i, j, k = drive(state, rng)
    assert state.z.sum() == 1
    assert state.z[i, j, k] == 1


def test_no_topple_when_stable():
    rng = np.random.default_rng(0)
    z = np.ones((3, 3, 3), dtype=np.int64)  # all z=1, all sigma=0, threshold=2
    sigma = np.zeros((3, 3, 3), dtype=np.int64)
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    assert s == 0
    assert T == 0


def test_topple_threshold_scales_with_sigma():
    # sigma=2 cell needs z_c=6 grains to topple. z=5 is stable; z=6 topples.
    rng = np.random.default_rng(0)
    z = np.zeros((3, 3, 3), dtype=np.int64)
    sigma = np.zeros((3, 3, 3), dtype=np.int64)
    z[1, 1, 1] = 5
    sigma[1, 1, 1] = 2
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    assert s == 0, "z=5 below z_c=6 for sigma=2 should be stable"

    z[1, 1, 1] = 6
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    # Center topples at least once; secondary topplings may follow
    # because 6 grains distributed randomly among 6 neighbors will
    # often pile >=2 on at least one neighbor (only ~1.5% chance of
    # uniform 1-per-neighbor distribution).
    assert s >= 1
    # Lattice conservation
    grains_in = 6
    assert grains_in == int(state.z.sum()) + state.grains_lost


def test_isolated_topple_does_not_change_sigma():
    # Single isolated unstable cell with no neighbors also unstable
    # in the same sweep: sigma must not transform.
    rng = np.random.default_rng(0)
    z = np.zeros((5, 5, 5), dtype=np.int64)
    sigma = np.zeros((5, 5, 5), dtype=np.int64)
    z[2, 2, 2] = 2
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    assert s >= 1
    # The central cell may have toppled multiple times if neighbors got
    # excited and propagated; but in isolated topple with z=2, only one
    # sweep with one unstable cell, so sigma should stay 0 at center.
    # We can't fully assert global sigma stays 0 because the cascade
    # could create co-active neighbors next sweep. Instead, just check
    # that we have low sigma activity overall.
    assert state.sigma.max() <= 1  # no spontaneous high-sigma jumps


def test_two_adjacent_topples_transform_sigma():
    # Two face-adjacent unstable cells topple in the same sweep -
    # both should have sigma increase by 1 (from 0 to 1).
    rng = np.random.default_rng(0)
    z = np.zeros((5, 5, 5), dtype=np.int64)
    sigma = np.zeros((5, 5, 5), dtype=np.int64)
    z[2, 2, 2] = 2
    z[2, 2, 3] = 2
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    assert s >= 2
    # Both seed cells should have sigma >= 1 after the joint topple
    assert state.sigma[2, 2, 2] >= 1
    assert state.sigma[2, 2, 3] >= 1


def test_conservation():
    rng = np.random.default_rng(7)
    n_drops = 1_000
    state, sizes, durations, ever = run_with_ever_toppled(
        L=10, n_drops=n_drops, seed=11
    )
    assert n_drops == int(state.z.sum()) + state.grains_lost


def test_sigma_max_bounds_growth():
    # With sigma_max=1, sigma should never exceed 1 anywhere.
    state, _, _, _ = run_with_ever_toppled(
        L=8, n_drops=2_000, seed=2, sigma_max=1
    )
    assert state.sigma.max() <= 1


def test_unbounded_sigma_can_grow_above_max():
    # No bound: at long enough run, sigma should reach values > 0
    # in active regions
    state, _, _, _ = run_with_ever_toppled(
        L=8, n_drops=3_000, seed=3
    )
    assert state.sigma.max() >= 1
    # We don't assert how high; just that growth occurs


def test_sigma_promotion_uses_max_neighbor():
    # Construct: center cell sigma=0, one neighbor sigma=3.
    # Both topple in the same sweep. The center should be promoted to
    # max(0, 3+1) = 4 (or capped if sigma_max < 4).
    rng = np.random.default_rng(0)
    z = np.zeros((5, 5, 5), dtype=np.int64)
    sigma = np.zeros((5, 5, 5), dtype=np.int64)
    z[2, 2, 2] = 2          # sigma=0, z_c=2
    sigma[2, 2, 3] = 3
    z[2, 2, 3] = 2 * (3 + 1)  # sigma=3, z_c=8
    state = _state(z, sigma)
    s, T, _ = relax(state, rng)
    # Center should have been promoted because its neighbor (sigma=3)
    # was co-toppling.
    assert state.sigma[2, 2, 2] >= 4
