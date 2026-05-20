"""Smoke test for the FSS sweep driver.

Goal: confirm run_one completes end-to-end on a tiny configuration
and produces a valid final .npz output. Also confirms that an
interrupted job (via simulating a midpoint checkpoint) resumes
correctly through the sweep driver's wrapper, not just the
checkpoint module in isolation.

These tests are fast (a few seconds each). They are NOT the
scientific verification - those live elsewhere - just functional
plumbing checks.
"""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import numpy as np
import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT = REPO_ROOT / "scripts" / "run_milestone6_fss_sweep.py"
SRC = REPO_ROOT / "src"

if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


def _import_sweep_module():
    """Load the sweep script as a module (it lives outside src/).

    The dataclass decorator looks up cls.__module__ in sys.modules, so
    we have to register the module before exec_module runs.
    """
    if "fss_sweep" in sys.modules:
        return sys.modules["fss_sweep"]
    spec = importlib.util.spec_from_file_location("fss_sweep", SCRIPT)
    mod = importlib.util.module_from_spec(spec)
    sys.modules["fss_sweep"] = mod
    spec.loader.exec_module(mod)
    return mod


def test_run_one_completes_small_L(tmp_path: Path) -> None:
    """run_one finishes a small L=8, 2000-drop job and writes a valid npz."""
    mod = _import_sweep_module()
    ckpt = tmp_path / "ckpt.npz"
    out = tmp_path / "final.npz"
    log = tmp_path / "run.log"

    result = mod.run_one(
        L=8,
        seed=42,
        n_drops_max=2000,
        ckpt_path=ckpt,
        out_path=out,
        log_path=log,
    )

    assert out.exists(), "final output not written"
    assert log.exists(), "log file not written"
    assert not ckpt.exists(), "checkpoint should be cleaned up on success"
    assert result.L == 8
    assert result.seed == 42
    assert result.drops_executed == 2000
    assert 0.0 <= result.final_p <= 1.0
    assert result.peak_size >= 1

    # Verify the output file is readable and has the expected fields.
    data = np.load(out, allow_pickle=True)
    assert int(data["L"]) == 8
    assert int(data["seed"]) == 42
    assert data["sizes"].shape == (2000,)
    assert data["ever_toppled"].shape == (8, 8, 8)


def test_resume_through_sweep_driver(tmp_path: Path) -> None:
    """An interrupted run resumed via run_one matches an uninterrupted run."""
    mod = _import_sweep_module()
    # Reference: uninterrupted run.
    ref_out = tmp_path / "ref_final.npz"
    mod.run_one(
        L=8,
        seed=99,
        n_drops_max=3000,
        ckpt_path=tmp_path / "ref_ckpt.npz",
        out_path=ref_out,
        log_path=tmp_path / "ref.log",
    )
    ref = np.load(ref_out, allow_pickle=True)

    # Resumed run: write a hand-rolled checkpoint at mid-run, then
    # let run_one finish from it.
    ckpt_path = tmp_path / "resume_ckpt.npz"
    out_path = tmp_path / "resume_final.npz"
    log_path = tmp_path / "resume.log"

    # First half: replay 1500 drops manually using the same primitives
    # as the sweep, then save a checkpoint mimicking what the policy
    # would do.
    from void_cascade.checkpoint import CheckpointPayload, save_checkpoint
    from void_cascade.sandpile_3d import drive, initialize, relax

    rng = np.random.default_rng(99)
    state = initialize(8)
    ever_toppled = np.zeros((8, 8, 8), dtype=bool)
    sizes = np.zeros(3000, dtype=np.int64)
    durations = np.zeros(3000, dtype=np.int64)
    for t in range(1500):
        drive(state, rng)
        s, T, mask = relax(state, rng, track_support=True)
        sizes[t] = s
        durations[t] = T
        if mask is not None:
            ever_toppled |= mask
    save_checkpoint(
        ckpt_path,
        CheckpointPayload(
            L=8, seed=99, n_drops_max=3000, drop=1500,
            state=state, rng_state=rng.bit_generator.state,
            sizes=sizes, durations=durations,
            ever_toppled=ever_toppled, snapshots=[],
        ),
    )

    # Now call run_one which should detect the checkpoint and resume.
    mod.run_one(
        L=8,
        seed=99,
        n_drops_max=3000,
        ckpt_path=ckpt_path,
        out_path=out_path,
        log_path=log_path,
    )
    resumed = np.load(out_path, allow_pickle=True)

    # Bit-identical output expected.
    assert np.array_equal(ref["sizes"], resumed["sizes"])
    assert np.array_equal(ref["durations"], resumed["durations"])
    assert np.array_equal(ref["ever_toppled"], resumed["ever_toppled"])
    assert np.array_equal(ref["final_z"], resumed["final_z"])
    assert int(ref["grains_lost"]) == int(resumed["grains_lost"])


def test_seed_for_is_deterministic() -> None:
    mod = _import_sweep_module()
    assert mod.seed_for(64, 0) == 16_400
    assert mod.seed_for(64, 1) == 16_401
    assert mod.seed_for(128, 0) == 22_800
    # Distinct across L values.
    assert mod.seed_for(32, 0) != mod.seed_for(48, 0)


def test_build_job_list_matches_defaults(tmp_path: Path) -> None:
    mod = _import_sweep_module()
    jobs = mod.build_job_list(
        [32, 48], {32: 2, 48: 3}, {32: 100, 48: 200}, tmp_path
    )
    assert len(jobs) == 5  # 2 + 3
    assert jobs[0][0] == 32  # L
    assert jobs[2][0] == 48
    # Paths point into the provided sweep_dir.
    for job in jobs:
        assert tmp_path in job[3].parents
