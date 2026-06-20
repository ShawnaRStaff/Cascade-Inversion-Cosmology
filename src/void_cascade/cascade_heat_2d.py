"""2D damage-coupled heat model (validated-SOC version).

Same idea as cascade_heat (1D) but on a 2D square lattice with the Manna
rule this project validated in M2 (threshold 2, shed 2 grains to random
neighbours of 4; tau~1.27). 1D Manna was degenerate (tau~0.78), so the 2D
version is what we trust for the heat numbers -- and it lets us ask the new
question: when the substrate tips, does the catastrophe spread as a FRONT
across the lattice, or pop everywhere at once?

Reuses the cell state (`HeatState`), params (`CascadeParams`) and the melt
`release` rule unchanged; only the lattice dimension of the avalanche and
the heat diffusion are new.
"""

from __future__ import annotations

import numpy as np

from void_cascade.cascade_heat import CascadeParams
from void_cascade.heat_gated import HeatState, release


def initialize_2d(L: int) -> HeatState:
    return HeatState(
        density=np.zeros((L, L), dtype=float),
        cracks=np.zeros((L, L), dtype=np.int64),
        heat=np.zeros((L, L), dtype=float),
        released=np.zeros((L, L), dtype=bool),
    )


def drive_2d(state: HeatState, p: CascadeParams, rng: np.random.Generator) -> None:
    L = state.density.shape[0]
    flat = state.density.ravel()
    sites = rng.integers(0, L * L, size=p.n_drive_sites)
    np.add.at(flat, sites, p.drive_amount)


def avalanche_2d(state: HeatState, p: CascadeParams, rng: np.random.Generator,
                 max_sweeps: int) -> int:
    """2D Manna relax: cells over threshold shed 2 grains to random neighbours
    of 4. Returns avalanche size (total fractures)."""
    L = state.density.shape[0]
    N = L * L
    dens = state.density.ravel()
    cr = state.cracks.ravel()
    ht = state.heat.ravel()
    rel = state.released.ravel()
    fd = p.fracture_density
    total = 0
    sweeps = 0
    while True:
        over = (dens >= fd) & (~rel)
        idx = np.flatnonzero(over)
        n = idx.size
        if n == 0:
            break
        sweeps += 1
        cr[idx] += 1
        ht[idx] += p.heat_per_crack
        dens[idx] -= 2.0
        # Two grains per fracture, each to a random one of the 4 neighbours.
        src = np.concatenate([idx, idx])
        i = src // L
        j = src % L
        dirs = rng.integers(0, 4, size=src.size)
        ni = i.copy()
        nj = j.copy()
        ni[dirs == 0] += 1
        ni[dirs == 1] -= 1
        nj[dirs == 2] += 1
        nj[dirs == 3] -= 1
        valid = (ni >= 0) & (ni < L) & (nj >= 0) & (nj < L)
        if valid.any():
            tgt = ni[valid] * L + nj[valid]
            dens += np.bincount(tgt, minlength=N)
        total += n
        if sweeps > max_sweeps:
            break
    return total


def diffuse_and_cool_2d(heat: np.ndarray, p: CascadeParams) -> np.ndarray:
    """5-point Laplacian with open (cold=0) boundaries, fractional cooling,
    floored at absolute zero. diffuse should be <= 0.25 for stability."""
    lap = -4.0 * heat
    lap[1:, :] += heat[:-1, :]
    lap[:-1, :] += heat[1:, :]
    lap[:, 1:] += heat[:, :-1]
    lap[:, :-1] += heat[:, 1:]
    out = heat + p.diffuse * lap
    out *= (1.0 - p.cooling)
    np.maximum(out, 0.0, out=out)
    return out


def step_2d(state: HeatState, p: CascadeParams, rng: np.random.Generator,
            max_sweeps: int) -> tuple[int, int]:
    drive_2d(state, p, rng)
    size = avalanche_2d(state, p, rng, max_sweeps)
    state.heat = diffuse_and_cool_2d(state.heat, p)
    n_rel, _ = release(state, p)
    return size, n_rel


def run_2d(L: int, n_steps: int, p: CascadeParams, seed: int | None = None) -> dict:
    """Drive the 2D lattice. Tracks the released region over time so we can
    see whether the catastrophe spreads as a front or pops everywhere."""
    rng = np.random.default_rng(seed)
    state = initialize_2d(L)
    max_sweeps = 50 * L * L
    cascade_sizes = np.zeros(n_steps, dtype=np.int64)
    sample_every = max(1, n_steps // 300)
    steps_axis, peak_heat_trace, n_released_trace = [], [], []
    peak_heat_overall = 0.0
    min_heat_overall = 0.0
    first_release_step: int | None = None
    snapshots: list[dict] = []
    snap_offsets = {0, 1, 2, 4, 8, 16}  # steps after first release to snapshot

    for t in range(n_steps):
        size, n_rel = step_2d(state, p, rng, max_sweeps)
        cascade_sizes[t] = size
        if first_release_step is None and n_rel > 0:
            first_release_step = t
        if first_release_step is not None and (t - first_release_step) in snap_offsets:
            snapshots.append({"step": t,
                              "after_tip": t - first_release_step,
                              "released": state.released.copy()})
        hmax = float(state.heat.max())
        hmin = float(state.heat.min())
        peak_heat_overall = max(peak_heat_overall, hmax)
        min_heat_overall = min(min_heat_overall, hmin)
        if t % sample_every == 0:
            steps_axis.append(t)
            peak_heat_trace.append(hmax)
            n_released_trace.append(int(state.released.sum()))

    return {
        "L": L,
        "n_steps": n_steps,
        "cooling": p.cooling,
        "ran_away": bool(state.released.any()),
        "first_release_step": first_release_step,
        "fraction_released_final": float(state.released.mean()),
        "peak_heat_overall": peak_heat_overall,
        "min_heat_overall": min_heat_overall,
        "max_cracks_on_a_cell": int(state.cracks.max()),
        "overall_max_cascade": int(cascade_sizes.max()),
        "cascade_sizes": cascade_sizes.tolist(),
        "steps_axis": steps_axis,
        "peak_heat_trace": peak_heat_trace,
        "n_released_trace": n_released_trace,
        "snapshots": snapshots,
    }
