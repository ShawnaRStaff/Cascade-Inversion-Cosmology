"""Checkpoint save/load for long-running 3D Manna simulations.

Designed for spot-instance compute (AWS spot, preemptible cloud VMs)
where the host can be reclaimed mid-run with only a few minutes of
notice. The contract:

  - Save the *minimum sufficient state* to resume bit-identically.
  - Atomic writes (tmp file + rename) so an interrupt during save
    never corrupts an existing checkpoint.
  - Save every N drops AND every M seconds, whichever fires first.
  - Resume: load the file, hand its state back to the caller, return
    the drop counter so the loop continues from there.

The state required for deterministic resume of the 3D Manna model is:
  - state.z          (the lattice grain field)
  - state.grains_lost (cumulative boundary loss, for conservation)
  - rng.bit_generator.state  (PCG64 internal state)
  - drop counter (which iteration to continue from)
  - sizes[:drop], durations[:drop]  (avalanche history so far)
  - ever_toppled (cumulative percolation set; needed because it's not
    derivable from z alone)
  - snapshots (the per-snapshot summary dicts the caller accumulates)
  - meta: L, n_drops_max, seed (for sanity check on resume)

Anything else (analysis, plots, etc.) is downstream of these and can be
recomputed.
"""

from __future__ import annotations

import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from .sandpile_3d import MannaState3D


@dataclass
class CheckpointPayload:
    """Everything a sim needs to resume bit-identically."""

    L: int
    seed: int
    n_drops_max: int
    drop: int  # next drop index to execute
    state: MannaState3D
    rng_state: dict  # rng.bit_generator.state at time of save
    sizes: np.ndarray  # per-event toppling count; full-length, zeros after drop
    durations: np.ndarray  # full-length, zeros after drop
    ever_toppled: np.ndarray
    snapshots: list[dict]
    # Per-event UNIQUE cells toppled at least once during that avalanche.
    # Distinct from sizes (which counts every topple, including re-topples).
    # unique_sizes[t] <= L^3 always; unique_sizes[t] / L^3 is the
    # "fraction of lattice involved in event t" metric.
    unique_sizes: np.ndarray = None  # type: ignore


def save_checkpoint(
    path: str | Path,
    payload: CheckpointPayload,
) -> None:
    """Atomically write the checkpoint to disk.

    Uses a tmp file + os.replace so a crash during save never
    corrupts an existing checkpoint file.
    """
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    # snapshots is a list[dict], save as object array. rng_state is a
    # nested dict; same treatment.
    # unique_sizes may be None for backward-compatible old checkpoints.
    # Save a sentinel zero-length array in that case.
    unique_sizes = (
        payload.unique_sizes
        if payload.unique_sizes is not None
        else np.zeros(0, dtype=np.int64)
    )
    np.savez_compressed(
        tmp,
        L=np.int64(payload.L),
        seed=np.int64(payload.seed),
        n_drops_max=np.int64(payload.n_drops_max),
        drop=np.int64(payload.drop),
        z=payload.state.z,
        grains_lost=np.int64(payload.state.grains_lost),
        rng_state=np.array([payload.rng_state], dtype=object),
        sizes=payload.sizes,
        durations=payload.durations,
        ever_toppled=payload.ever_toppled,
        snapshots=np.array(payload.snapshots, dtype=object),
        unique_sizes=unique_sizes,
    )
    # np.savez adds .npz to tmp; account for that.
    tmp_with_ext = tmp.with_suffix(tmp.suffix + ".npz") if tmp.suffix != ".npz" else tmp
    if tmp_with_ext.exists() and tmp_with_ext != tmp:
        os.replace(tmp_with_ext, path)
    else:
        os.replace(tmp, path)


def load_checkpoint(path: str | Path) -> CheckpointPayload | None:
    """Load a checkpoint. Returns None if path does not exist.

    Raises ValueError if the file is malformed or version-incompatible.
    Allow_pickle is required for the rng_state and snapshots fields.
    """
    path = Path(path)
    if not path.exists():
        return None
    data = np.load(path, allow_pickle=True)
    try:
        rng_state = data["rng_state"][0]  # unbox the 1-element object array
        snapshots = list(data["snapshots"])
        state = MannaState3D(
            z=data["z"].astype(np.int64),
            grains_lost=int(data["grains_lost"]),
        )
        # unique_sizes optional for backward compat with old checkpoints.
        unique_sizes = (
            data["unique_sizes"] if "unique_sizes" in data.files else None
        )
        if unique_sizes is not None and unique_sizes.size == 0:
            unique_sizes = None
        return CheckpointPayload(
            L=int(data["L"]),
            seed=int(data["seed"]),
            n_drops_max=int(data["n_drops_max"]),
            drop=int(data["drop"]),
            state=state,
            rng_state=rng_state,
            sizes=data["sizes"],
            durations=data["durations"],
            ever_toppled=data["ever_toppled"],
            snapshots=snapshots,
            unique_sizes=unique_sizes,
        )
    except KeyError as e:
        raise ValueError(f"Malformed checkpoint at {path}: missing {e}") from e


def restore_rng(rng_state: dict) -> np.random.Generator:
    """Reconstruct a numpy Generator from a saved bit_generator state.

    The state dict has 'bit_generator': 'PCG64' (or similar) and the
    internal state needed to continue the stream.
    """
    bg_name = rng_state["bit_generator"]
    # The bit-generator classes live on np.random; getattr by name.
    bg_cls = getattr(np.random, bg_name)
    bg = bg_cls()
    bg.state = rng_state
    return np.random.Generator(bg)


class CheckpointPolicy:
    """When to checkpoint.

    Fires when EITHER the drop interval OR the wall-clock interval
    is exceeded since the last checkpoint, whichever comes first.
    """

    def __init__(
        self,
        every_drops: int = 50_000,
        every_seconds: float = 600.0,
    ):
        self.every_drops = int(every_drops)
        self.every_seconds = float(every_seconds)
        self._last_drop = 0
        self._last_time = time.time()

    def should_save(self, current_drop: int) -> bool:
        now = time.time()
        if current_drop - self._last_drop >= self.every_drops:
            return True
        if now - self._last_time >= self.every_seconds:
            return True
        return False

    def mark_saved(self, current_drop: int) -> None:
        self._last_drop = current_drop
        self._last_time = time.time()
