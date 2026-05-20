"""Tests for the checkpoint module.

Core property: a run that is interrupted, checkpointed, and resumed
must produce bit-identical output to an uninterrupted run with the
same seed and parameters.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from void_cascade.checkpoint import (
    CheckpointPayload,
    CheckpointPolicy,
    load_checkpoint,
    restore_rng,
    save_checkpoint,
)
from void_cascade.sandpile_3d import MannaState3D, drive, initialize, relax


def _run_segment(
    rng: np.random.Generator,
    state: MannaState3D,
    ever_toppled: np.ndarray,
    sizes: np.ndarray,
    durations: np.ndarray,
    start_drop: int,
    end_drop: int,
) -> None:
    """Drive (start_drop, end_drop] applying drive+relax each tick."""
    for t in range(start_drop, end_drop):
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True)
        sizes[t] = s
        durations[t] = T
        if mask is not None:
            ever_toppled |= mask


def _run_uninterrupted(L: int, seed: int, n_drops: int):
    rng = np.random.default_rng(seed)
    state = initialize(L)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    _run_segment(rng, state, ever_toppled, sizes, durations, 0, n_drops)
    return state, sizes, durations, ever_toppled


def _run_with_interruption(
    L: int, seed: int, n_drops: int, interrupt_at: int, ckpt_path: Path
):
    """Run to interrupt_at, save checkpoint, restart from checkpoint, finish."""
    # First segment.
    rng = np.random.default_rng(seed)
    state = initialize(L)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    _run_segment(rng, state, ever_toppled, sizes, durations, 0, interrupt_at)

    save_checkpoint(
        ckpt_path,
        CheckpointPayload(
            L=L,
            seed=seed,
            n_drops_max=n_drops,
            drop=interrupt_at,
            state=state,
            rng_state=rng.bit_generator.state,
            sizes=sizes,
            durations=durations,
            ever_toppled=ever_toppled,
            snapshots=[],
        ),
    )

    # Discard everything in memory to simulate process death.
    del rng, state, ever_toppled, sizes, durations

    payload = load_checkpoint(ckpt_path)
    assert payload is not None
    rng2 = restore_rng(payload.rng_state)
    state2 = payload.state
    ever_toppled2 = payload.ever_toppled
    sizes2 = payload.sizes
    durations2 = payload.durations

    _run_segment(
        rng2, state2, ever_toppled2, sizes2, durations2, payload.drop, n_drops
    )
    return state2, sizes2, durations2, ever_toppled2


def test_resume_is_bit_identical(tmp_path: Path) -> None:
    """Interrupted + resumed run must match uninterrupted run exactly."""
    L = 16
    seed = 12345
    n_drops = 2000
    interrupt_at = 700

    state_a, sizes_a, durations_a, ever_a = _run_uninterrupted(L, seed, n_drops)
    state_b, sizes_b, durations_b, ever_b = _run_with_interruption(
        L, seed, n_drops, interrupt_at, tmp_path / "ckpt.npz"
    )

    assert np.array_equal(state_a.z, state_b.z), "final z field differs"
    assert state_a.grains_lost == state_b.grains_lost, "grains_lost differs"
    assert np.array_equal(sizes_a, sizes_b), "sizes history differs"
    assert np.array_equal(durations_a, durations_b), "durations history differs"
    assert np.array_equal(ever_a, ever_b), "ever_toppled set differs"


def test_resume_multiple_interrupts(tmp_path: Path) -> None:
    """Multiple checkpoints in a single run still produce identical output."""
    L = 12
    seed = 999
    n_drops = 1500
    interrupts = [300, 600, 900, 1200]

    state_a, sizes_a, durations_a, ever_a = _run_uninterrupted(L, seed, n_drops)

    rng = np.random.default_rng(seed)
    state = initialize(L)
    ever_toppled = np.zeros((L, L, L), dtype=bool)
    sizes = np.zeros(n_drops, dtype=np.int64)
    durations = np.zeros(n_drops, dtype=np.int64)
    current = 0

    for interrupt_at in interrupts:
        _run_segment(rng, state, ever_toppled, sizes, durations, current, interrupt_at)
        ckpt = tmp_path / f"ckpt_{interrupt_at}.npz"
        save_checkpoint(
            ckpt,
            CheckpointPayload(
                L=L, seed=seed, n_drops_max=n_drops, drop=interrupt_at,
                state=state, rng_state=rng.bit_generator.state,
                sizes=sizes, durations=durations, ever_toppled=ever_toppled,
                snapshots=[],
            ),
        )
        # Discard and restore from the checkpoint we just wrote.
        del rng, state, ever_toppled, sizes, durations
        payload = load_checkpoint(ckpt)
        assert payload is not None
        rng = restore_rng(payload.rng_state)
        state = payload.state
        ever_toppled = payload.ever_toppled
        sizes = payload.sizes
        durations = payload.durations
        current = payload.drop

    _run_segment(rng, state, ever_toppled, sizes, durations, current, n_drops)

    assert np.array_equal(state_a.z, state.z)
    assert state_a.grains_lost == state.grains_lost
    assert np.array_equal(sizes_a, sizes)
    assert np.array_equal(durations_a, durations)
    assert np.array_equal(ever_a, ever_toppled)


def test_load_missing_returns_none(tmp_path: Path) -> None:
    assert load_checkpoint(tmp_path / "does_not_exist.npz") is None


def test_atomic_write_replaces_existing(tmp_path: Path) -> None:
    """Saving twice to the same path should replace the file atomically."""
    L = 8
    seed = 1
    ckpt_path = tmp_path / "ckpt.npz"

    rng = np.random.default_rng(seed)
    state = initialize(L)
    payload_a = CheckpointPayload(
        L=L, seed=seed, n_drops_max=100, drop=10,
        state=state, rng_state=rng.bit_generator.state,
        sizes=np.zeros(100, dtype=np.int64),
        durations=np.zeros(100, dtype=np.int64),
        ever_toppled=np.zeros((L, L, L), dtype=bool),
        snapshots=[{"drop": 10}],
    )
    save_checkpoint(ckpt_path, payload_a)
    first_size = ckpt_path.stat().st_size

    payload_b = CheckpointPayload(
        L=L, seed=seed, n_drops_max=100, drop=20,
        state=state, rng_state=rng.bit_generator.state,
        sizes=np.zeros(100, dtype=np.int64),
        durations=np.zeros(100, dtype=np.int64),
        ever_toppled=np.zeros((L, L, L), dtype=bool),
        snapshots=[{"drop": 10}, {"drop": 20}],
    )
    save_checkpoint(ckpt_path, payload_b)

    loaded = load_checkpoint(ckpt_path)
    assert loaded is not None
    assert loaded.drop == 20
    assert len(loaded.snapshots) == 2


def test_policy_fires_on_drop_interval() -> None:
    pol = CheckpointPolicy(every_drops=100, every_seconds=999999)
    assert not pol.should_save(50)
    assert pol.should_save(100)
    pol.mark_saved(100)
    assert not pol.should_save(150)
    assert pol.should_save(200)


def test_policy_fires_on_time_interval(monkeypatch: pytest.MonkeyPatch) -> None:
    fake_time = [1000.0]

    def mock_time():
        return fake_time[0]

    monkeypatch.setattr("void_cascade.checkpoint.time.time", mock_time)
    pol = CheckpointPolicy(every_drops=1_000_000, every_seconds=10.0)
    fake_time[0] = 1005.0
    assert not pol.should_save(100)
    fake_time[0] = 1011.0
    assert pol.should_save(100)


def test_snapshots_round_trip(tmp_path: Path) -> None:
    """Snapshots (list of dicts) survive save/load."""
    L = 8
    seed = 1
    rng = np.random.default_rng(seed)
    state = initialize(L)
    snapshots = [
        {"drop": 100, "p": 0.5, "max_s": 42},
        {"drop": 200, "p": 0.8, "max_s": 100},
    ]
    save_checkpoint(
        tmp_path / "ckpt.npz",
        CheckpointPayload(
            L=L, seed=seed, n_drops_max=300, drop=200,
            state=state, rng_state=rng.bit_generator.state,
            sizes=np.zeros(300, dtype=np.int64),
            durations=np.zeros(300, dtype=np.int64),
            ever_toppled=np.zeros((L, L, L), dtype=bool),
            snapshots=snapshots,
        ),
    )
    loaded = load_checkpoint(tmp_path / "ckpt.npz")
    assert loaded is not None
    assert loaded.snapshots == snapshots
