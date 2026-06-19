"""Tests for the irreversible-fracture sandpile.

The defining rule: a cell topples at most once, then is permanently
fractured. We test that rule directly, plus grain conservation in both
modes and the resulting hard ceiling on total topplings.
"""

import numpy as np
import pytest

from void_cascade.irreversible import (
    IrreversibleState,
    _apply_sweep,
    initialize,
    relax,
    run,
)


def _state(z: np.ndarray, fractured: np.ndarray | None = None) -> IrreversibleState:
    z = z.astype(np.int64)
    if fractured is None:
        fractured = np.zeros_like(z, dtype=bool)
    return IrreversibleState(z=z, fractured=fractured)


def test_fractured_cell_does_not_topple():
    # A cell that is already fractured stays put even when over threshold.
    z = np.zeros((3, 3, 3), dtype=np.int64)
    z[1, 1, 1] = 5
    fractured = np.zeros((3, 3, 3), dtype=bool)
    fractured[1, 1, 1] = True
    state = _state(z, fractured)
    n = _apply_sweep(state, np.random.default_rng(0), mode="sink")
    assert n == 0
    assert state.z[1, 1, 1] == 5


def test_unfractured_cell_topples_once_then_is_fractured():
    z = np.zeros((3, 3, 3), dtype=np.int64)
    z[1, 1, 1] = 2
    state = _state(z)
    n = _apply_sweep(state, np.random.default_rng(0), mode="sink")
    assert n == 1
    assert state.fractured[1, 1, 1]
    assert state.z[1, 1, 1] == 0  # lost its 2 grains


def test_conservation_sink():
    res = run(L=8, n_drops=20_000, mode="sink", seed=1)
    assert res["grains_accounted"] == res["grains_in"]


def test_conservation_hole():
    res = run(L=8, n_drops=20_000, mode="hole", seed=1)
    assert res["grains_accounted"] == res["grains_in"]


def test_total_topplings_never_exceeds_volume():
    # The whole point of irreversibility: each cell fractures at most once.
    res = run(L=8, n_drops=30_000, mode="sink", seed=2)
    assert res["total_topplings"] <= res["volume"]


def test_goes_silent_once_fully_fractured():
    # Small lattice driven hard: once everything is fractured, the
    # confirm-quiet tail must contain no topplings at all.
    res = run(L=6, n_drops=200_000, mode="sink", seed=3)
    assert res["fully_fractured_at_drop"] is not None
    assert res["events_after_full_fracture_max"] == 0


def test_unknown_mode_raises():
    z = np.zeros((3, 3, 3), dtype=np.int64)
    z[1, 1, 1] = 2
    with pytest.raises(ValueError):
        _apply_sweep(_state(z), np.random.default_rng(0), mode="nonsense")
