"""Unit tests for the 3D Manna sandpile.

Toppling and conservation are checked deterministically with hand-built
states. The cosmologically-relevant percolation check lives in
test_percolation.py; here we only confirm the sandpile dynamics are
correct.
"""

from __future__ import annotations

import numpy as np

from void_cascade.sandpile_3d import (
    MannaState3D,
    drive,
    initialize,
    relax,
    run,
    run_with_ever_toppled,
)


def _state(z: np.ndarray) -> MannaState3D:
    return MannaState3D(z=z.astype(np.int64))


def test_initialize_empty():
    state = initialize(L=4)
    assert state.z.shape == (4, 4, 4)
    assert state.z.sum() == 0
    assert state.grains_lost == 0


def test_drive_adds_one_grain():
    rng = np.random.default_rng(0)
    state = initialize(L=3)
    i, j, k = drive(state, rng)
    assert state.z.sum() == 1
    assert state.z[i, j, k] == 1


def test_no_topple_when_stable():
    rng = np.random.default_rng(0)
    z = np.ones((3, 3, 3), dtype=np.int64)
    state = _state(z)
    s, T, _ = relax(state, rng)
    assert s == 0
    assert T == 0
    assert state.z.sum() == 27


def test_interior_topple_conserves_six_local_neighborhood():
    # Single interior unstable site at center (1,1,1) loses 2 grains.
    # No diagonal contributions; all 2 grains land on face neighbors.
    rng = np.random.default_rng(0)
    z = np.zeros((3, 3, 3), dtype=np.int64)
    z[1, 1, 1] = 2
    state = _state(z)
    s, T, _ = relax(state, rng)
    assert s == 1
    assert T == 1
    assert state.z[1, 1, 1] == 0
    # The 6 face neighbors collectively gained 2 grains.
    nn = (
        state.z[0, 1, 1] + state.z[2, 1, 1]
        + state.z[1, 0, 1] + state.z[1, 2, 1]
        + state.z[1, 1, 0] + state.z[1, 1, 2]
    )
    assert nn == 2
    assert state.grains_lost == 0
    # Diagonal cells untouched (sandpile is face-connected).
    assert state.z[0, 0, 0] == 0
    assert state.z[2, 2, 2] == 0


def test_corner_topple_can_lose_grains():
    # Corner has only 3 face neighbors on-lattice. Three of the 6 directions
    # are off-lattice, so on average half the grains are lost. Conservation
    # must hold exactly.
    rng = np.random.default_rng(0)
    z = np.zeros((3, 3, 3), dtype=np.int64)
    z[0, 0, 0] = 2
    state = _state(z)
    grains_in = int(state.z.sum())
    s, _T, _ = relax(state, rng)
    assert s == 1
    assert int(state.z.sum()) + state.grains_lost == grains_in


def test_run_conservation():
    n_drops = 1_500
    state, _, _ = run(L=10, n_drops=n_drops, seed=7)
    assert n_drops == int(state.z.sum()) + state.grains_lost


def test_run_with_ever_toppled_consistency():
    state, sizes, durations, ever = run_with_ever_toppled(
        L=10, n_drops=500, seed=11
    )
    # Conservation still holds.
    assert 500 == int(state.z.sum()) + state.grains_lost
    # ever_toppled is the union of avalanche supports: it must contain
    # at least one True per non-zero-size avalanche.
    if sizes.max() > 0:
        assert ever.any()
    # Shape sanity.
    assert ever.shape == (10, 10, 10)
    assert ever.dtype == bool


def test_run_with_ever_toppled_callback_can_halt():
    halted_at: list[int] = []

    def stop_after_10(t, ever, sizes, durations):
        if t >= 10:
            halted_at.append(t)
            return True
        return False

    _state, sizes, durations, ever = run_with_ever_toppled(
        L=8, n_drops=10_000, seed=3, check_every=1,
        percolation_callback=stop_after_10,
    )
    assert halted_at
    # Sizes should be trimmed to the actually-executed range.
    assert sizes.size <= 12


def test_steady_state_density_in_band():
    # Quick sanity that the 3D Manna density settles somewhere below z_c=2.
    # Run is short to keep the test fast; the bound is loose.
    state, _, _ = run(L=12, n_drops=4_000, seed=13)
    rho = float(state.z.mean())
    assert 0.2 < rho < 1.5


def test_relax_handles_fully_primed_lattice():
    # A 4^3 lattice initialized at z=2 everywhere should release a large
    # avalanche and stabilize, with conservation intact.
    rng = np.random.default_rng(5)
    z = np.full((4, 4, 4), 2, dtype=np.int64)
    state = _state(z)
    grains_in = int(state.z.sum())
    s, T, _ = relax(state, rng)
    assert s > 0
    assert T > 0
    assert np.all(state.z < 2)
    assert int(state.z.sum()) + state.grains_lost == grains_in
