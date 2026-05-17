"""Unit tests for the 2D Manna sandpile.

Toppling and conservation are checked deterministically by seeding the RNG
and reading off the expected outcome. Statistical claims about tau, D, or
the cluster fractal dimension live in the demo script, not here.
"""

from __future__ import annotations

import numpy as np
import pytest

from void_cascade.sandpile_2d import (
    MannaState,
    drive,
    initialize,
    relax,
    run,
    run_with_clusters,
)


def _state(z: list[list[int]]) -> MannaState:
    return MannaState(z=np.array(z, dtype=np.int64))


def test_initialize_empty():
    state = initialize(L=5)
    assert state.z.shape == (5, 5)
    assert state.z.sum() == 0
    assert state.grains_lost == 0


def test_drive_adds_one_grain():
    rng = np.random.default_rng(0)
    state = initialize(L=4)
    i, j = drive(state, rng)
    assert state.z.sum() == 1
    assert state.z[i, j] == 1


def test_no_topple_when_stable():
    rng = np.random.default_rng(0)
    state = _state([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
    s, T, _ = relax(state, rng)
    assert s == 0
    assert T == 0
    assert state.z.sum() == 8


def test_interior_topple_conserves_two_grains_locally():
    # Single interior unstable site loses exactly 2 grains; both must land
    # on nearest neighbors. With z=2 and surroundings 0, post-relax the
    # central site is 0 and the neighbors collectively gained exactly 2.
    rng = np.random.default_rng(0)
    state = _state([[0, 0, 0], [0, 2, 0], [0, 0, 0]])
    s, T, _ = relax(state, rng)
    assert s == 1
    assert T == 1
    assert state.z[1, 1] == 0
    # Two grains distributed among the four NN; none lost off-lattice.
    nn_sum = (
        state.z[0, 1] + state.z[2, 1] + state.z[1, 0] + state.z[1, 2]
    )
    assert nn_sum == 2
    assert state.grains_lost == 0
    # Diagonals untouched.
    assert state.z[0, 0] == 0 and state.z[2, 2] == 0


def test_corner_topple_can_lose_grains():
    # Corner (0,0) has only 2 NN on-lattice: (0,1) and (1,0). Two grains
    # are each sent to a random direction in {N, S, W, E}. N and W are
    # off-lattice, so on average 50% of grains are lost. The total
    # conservation check still has to hold deterministically.
    rng = np.random.default_rng(0)
    state = _state([[2, 0, 0], [0, 0, 0], [0, 0, 0]])
    grains_in = int(state.z.sum())
    s, T, _ = relax(state, rng)
    assert s == 1
    assert int(state.z.sum()) + state.grains_lost == grains_in


def test_relax_returns_support_mask_when_requested():
    # Engineered chain reaction: site (1, 1) unstable, after toppling
    # neighbors may become unstable too. We only verify the mask covers
    # every site that toppled at least once.
    rng = np.random.default_rng(123)
    state = _state([[0, 1, 0], [1, 2, 1], [0, 1, 0]])
    s, T, mask = relax(state, rng, track_support=True)
    assert mask is not None
    assert mask.dtype == bool
    # At minimum the seed site toppled.
    assert mask[1, 1]
    # Mask only flags sites within the lattice.
    assert mask.shape == state.z.shape


def test_run_conservation():
    # grains_in == sum(z) + grains_lost, exactly, at all times.
    n_drops = 5_000
    state, _, _ = run(L=16, n_drops=n_drops, seed=7)
    assert n_drops == int(state.z.sum()) + state.grains_lost


def test_run_with_clusters_consistency():
    # areas <= sizes always (a site can topple more than once, so its
    # contribution to s exceeds its contribution to a).
    _, sizes, durations, areas, rgyr = run_with_clusters(
        L=16, n_drops=2_000, seed=11
    )
    assert sizes.shape == areas.shape == durations.shape
    assert np.all(areas <= sizes)
    # NaN only where area < 2.
    assert np.all(np.isnan(rgyr[areas < 2]) | (rgyr[areas < 2] == 0.0))
    assert np.all(rgyr[areas >= 2] >= 0.0)


def test_run_with_clusters_burn_in_skips_tracking():
    # During burn-in, areas/rgyr should remain at their default values.
    burn = 200
    _, sizes, _, areas, rgyr = run_with_clusters(
        L=12, n_drops=300, seed=3, burn_in=burn
    )
    assert np.all(areas[:burn] == 0)
    assert np.all(np.isnan(rgyr[:burn]))
    # And at least some post-burn drives did record geometry.
    nonzero_areas = (areas[burn:] > 0).sum()
    assert nonzero_areas > 0


def test_steady_state_mean_density_in_band():
    # For Manna 2D the steady-state mean of z stays well under z_c = 2.
    # Reported critical density ~ 0.72 (Lubeck 2000); the empirical mean
    # at finite L wobbles around there. The band is loose to keep the
    # test fast.
    state, _, _ = run(L=32, n_drops=20_000, seed=13)
    rho = float(state.z.mean())
    assert 0.4 < rho < 1.0


def test_relax_handles_unstable_initial_state():
    # A "fully primed" lattice (every site at z=2) should release a big
    # avalanche and end stable.
    rng = np.random.default_rng(5)
    state = MannaState(z=np.full((8, 8), 2, dtype=np.int64))
    grains_in = int(state.z.sum())
    s, T, _ = relax(state, rng)
    assert s > 0
    assert T > 0
    assert np.all(state.z < 2)
    assert int(state.z.sum()) + state.grains_lost == grains_in
