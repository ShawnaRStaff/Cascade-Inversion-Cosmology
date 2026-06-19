"""3D Manna sandpile with IRREVERSIBLE fracture (scrutiny experiment).

The standard model self-heals: a cell empties when it topples, then
refills from neighbours and new input and topples again — forever. A real
fracture would not heal. Here a cell topples AT MOST ONCE. The first time
it becomes unstable it "fractures" and is then permanently changed:

  mode="sink": the fractured cell still receives grains (they accumulate)
               but it never topples again.
  mode="hole": the fractured cell is a hole; grains routed into it are
               lost from the system, like an internal open boundary.

Question this addresses: Milestone 5 found a permanent, never-ending
regime of huge events. Does that survive when fracture cannot heal, or
was the permanence an artifact of the recharge? We do not presuppose the
answer — we measure it.

Conservation (both modes): grains_in == sum(z) + grains_lost.
Irreversibility invariant: total topplings over all time <= L**3 (each
cell fractures at most once).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MODES = ("sink", "hole")


@dataclass
class IrreversibleState:
    """State of a 3D Manna sandpile in which fracture is permanent.

    Attributes
    ----------
    z : ndarray[int], shape (L, L, L)
        Local grain count per site.
    fractured : ndarray[bool], shape (L, L, L)
        True once a cell has toppled. A fractured cell never topples again.
    grains_lost : int
        Grains lost off the boundary plus (in hole mode) grains routed
        into fractured cells.
    """

    z: np.ndarray
    fractured: np.ndarray
    grains_lost: int = 0

    @property
    def L(self) -> int:
        return int(self.z.shape[0])


def initialize(L: int) -> IrreversibleState:
    """Empty L x L x L lattice, nothing fractured yet."""
    return IrreversibleState(
        z=np.zeros((L, L, L), dtype=np.int64),
        fractured=np.zeros((L, L, L), dtype=bool),
    )


def drive(state: IrreversibleState, rng: np.random.Generator) -> tuple[int, int, int]:
    """Add one grain at a uniformly random site."""
    L = state.L
    i = int(rng.integers(0, L))
    j = int(rng.integers(0, L))
    k = int(rng.integers(0, L))
    state.z[i, j, k] += 1
    return i, j, k


def _apply_sweep(state: IrreversibleState, rng: np.random.Generator, mode: str) -> int:
    """One parallel sweep. Only NOT-yet-fractured unstable cells topple.

    Each cell that topples is marked fractured (it has spent its one
    fracture) and so can never topple again. Returns the topple count.
    """
    zf = state.z.ravel()
    frf = state.fractured.ravel()
    L = state.z.shape[0]
    N = zf.size
    LL = L * L

    unstable_flat = np.flatnonzero((zf >= 2) & (~frf))
    n = unstable_flat.size
    if n == 0:
        return 0

    # These cells spend their one and only fracture now.
    frf[unstable_flat] = True

    i_idx = unstable_flat // LL
    rem = unstable_flat - i_idx * LL
    j_idx = rem // L
    k_idx = rem - j_idx * L

    zf[unstable_flat] -= 2

    dirs = rng.integers(0, 6, size=2 * n)
    ti = np.tile(i_idx, 2)
    tj = np.tile(j_idx, 2)
    tk = np.tile(k_idx, 2)
    ti[dirs == 0] += 1
    ti[dirs == 1] -= 1
    tj[dirs == 2] += 1
    tj[dirs == 3] -= 1
    tk[dirs == 4] += 1
    tk[dirs == 5] -= 1

    on_lattice = (
        (ti >= 0) & (ti < L) & (tj >= 0) & (tj < L) & (tk >= 0) & (tk < L)
    )
    state.grains_lost += int(2 * n - on_lattice.sum())

    valid_flat = ti[on_lattice] * LL + tj[on_lattice] * L + tk[on_lattice]

    if mode == "hole":
        # A grain routed into an already-fractured cell falls into the hole.
        into_hole = frf[valid_flat]
        state.grains_lost += int(into_hole.sum())
        valid_flat = valid_flat[~into_hole]
    elif mode != "sink":
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")

    if valid_flat.size:
        dz = np.bincount(valid_flat, minlength=N).astype(np.int64, copy=False)
        zf += dz

    return n


def relax(state: IrreversibleState, rng: np.random.Generator, mode: str) -> tuple[int, int]:
    """Topple until no not-yet-fractured cell is unstable. Returns (size, duration).

    Guaranteed to terminate: every sweep fractures at least one new cell,
    and there are finitely many cells.
    """
    s = 0
    T = 0
    while True:
        n = _apply_sweep(state, rng, mode)
        if n == 0:
            break
        s += n
        T += 1
    return s, T


def run(
    L: int,
    n_drops: int,
    mode: str,
    seed: int | None = None,
    confirm_quiet_drops: int = 10_000,
) -> dict:
    """Drive an empty lattice until n_drops, or until everything is fractured.

    Once every cell has fractured, no cell can ever topple again, so the
    system is provably silent. We stop early at that point but keep driving
    `confirm_quiet_drops` more so the silence is measured, not assumed.

    Returns a dict of traces and summary numbers.
    """
    rng = np.random.default_rng(seed)
    state = initialize(L)
    volume = L ** 3
    sizes = np.zeros(n_drops, dtype=np.int64)
    sample_every = max(1, n_drops // 400)
    frac_drops, frac_values = [], []
    fully_fractured_at = None
    drops_done = 0

    for t in range(n_drops):
        drive(state, rng)
        s, _ = relax(state, rng, mode)
        sizes[t] = s
        drops_done = t + 1
        if t % sample_every == 0:
            frac_drops.append(t)
            frac_values.append(float(state.fractured.mean()))
        if fully_fractured_at is None and bool(state.fractured.all()):
            fully_fractured_at = t
            # Keep driving a little to confirm everything stays silent.
            tail_end = min(n_drops, t + 1 + confirm_quiet_drops)
            for tt in range(t + 1, tail_end):
                drive(state, rng)
                s2, _ = relax(state, rng, mode)
                sizes[tt] = s2
                drops_done = tt + 1
            break

    sizes = sizes[:drops_done]
    tail = sizes[fully_fractured_at + 1:] if fully_fractured_at is not None else np.array([], dtype=np.int64)
    last_tenth = sizes[max(0, drops_done - drops_done // 10):]

    return {
        "mode": mode,
        "L": L,
        "volume": volume,
        "drops_done": int(drops_done),
        "total_topplings": int(sizes.sum()),
        "max_event_size": int(sizes.max()) if sizes.size else 0,
        "max_event_pct_of_volume": (float(sizes.max()) / volume * 100.0) if sizes.size else 0.0,
        "final_fractured_fraction": float(state.fractured.mean()),
        "fully_fractured_at_drop": fully_fractured_at,
        "events_after_full_fracture_max": int(tail.max()) if tail.size else 0,
        "mean_event_last_tenth": float(last_tenth.mean()) if last_tenth.size else 0.0,
        "grains_in": int(drops_done),
        "grains_accounted": int(state.z.sum() + state.grains_lost),
        "frac_trace_drops": frac_drops,
        "frac_trace_values": frac_values,
        "sizes": sizes.tolist(),
    }
