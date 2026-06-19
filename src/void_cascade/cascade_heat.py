"""Damage-coupled heat model: fracturing that goes scattered -> cascading.

This fixes the flaw the earlier heat_gated model had: there, cells cracked
independently from a uniform drive, so heat was dribbled out thin and only a
zero-cooling case could ever tip. Here fracturing is **sandpile-coupled** --
a fracture pushes its load onto a neighbour, which can be pushed over
threshold too, producing an **avalanche**. Early on (sparse, low density)
avalanches are tiny and scattered; as the substrate fills toward critical,
avalanches grow and cluster. Clustered avalanches dump heat **concentrated**
in space and time, which is what can outrun cooling and tip the substrate.

This reunites the validated SOC avalanche dynamics (M1-M3) with the heat
gating (see docs/notes/heat_gated_model_design.md). The heat machinery
(`diffuse_and_cool`, `release`) and the cell state are reused from
`heat_gated` unchanged; only the fracturing is new.

Honest question, unchanged: stay cold forever, or run away? Built so strong
cooling stays cold -- not rigged to tip.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from void_cascade.heat_gated import HeatState, diffuse_and_cool, release


@dataclass(frozen=True)
class CascadeParams:
    fracture_density: float   # threshold to fracture (Manna z_c, use 2.0)
    heat_per_crack: float     # heat from one fracture
    diffuse: float            # neighbour heat sharing
    cooling: float            # fraction of heat lost per step (guardrail)
    melt_heat: float          # heat above which frozen cracks release
    release_factor: float     # heat burst per released crack
    drive_amount: float       # darkness added per driven site per step
    n_drive_sites: int        # how many sites get driven each step


def initialize(L: int) -> HeatState:
    return HeatState(
        density=np.zeros(L, dtype=float),
        cracks=np.zeros(L, dtype=np.int64),
        heat=np.zeros(L, dtype=float),
        released=np.zeros(L, dtype=bool),
    )


def drive(state: HeatState, p: CascadeParams, rng: np.random.Generator) -> None:
    """Absorb darkness at a few random sites (canonical sandpile drive)."""
    L = state.density.size
    sites = rng.integers(0, L, size=p.n_drive_sites)
    np.add.at(state.density, sites, p.drive_amount)


def avalanche(
    state: HeatState,
    p: CascadeParams,
    rng: np.random.Generator,
    max_sweeps: int,
) -> int:
    """Relax until stable. Cells over threshold fracture and push their load
    to a random neighbour (off-lattice = lost), which can cascade. Returns the
    avalanche size (total fractures this event).

    Stochastic single-neighbour targeting plus boundary loss makes avalanches
    finite (the Manna property); `max_sweeps` is a backstop only.
    """
    L = state.density.size
    fd = p.fracture_density
    total = 0
    sweeps = 0
    while True:
        over = (state.density >= fd) & (~state.released)
        idx = np.flatnonzero(over)
        n = idx.size
        if n == 0:
            break
        sweeps += 1
        state.cracks[idx] += 1
        state.heat[idx] += p.heat_per_crack
        # Manna rule: shed 2 units, one grain to each of two random neighbours.
        # A cell must accumulate toward the threshold, so a single received
        # grain only fractures a cell already near threshold -> bounded,
        # self-organising avalanches (not a lattice-spanning random walk).
        state.density[idx] -= 2.0
        directions = rng.integers(0, 2, size=2 * n) * 2 - 1  # -1 or +1
        targets = np.concatenate([idx, idx]) + directions
        valid = (targets >= 0) & (targets < L)
        if valid.any():
            incoming = np.bincount(targets[valid], minlength=L)  # 1 unit per grain
            state.density += incoming
        total += n
        if sweeps > max_sweeps:
            break
    return total


def step(
    state: HeatState,
    p: CascadeParams,
    rng: np.random.Generator,
    max_sweeps: int,
) -> tuple[int, int]:
    """One step: absorb darkness -> avalanche -> spread/cool heat -> maybe release."""
    drive(state, p, rng)
    size = avalanche(state, p, rng, max_sweeps)
    state.heat = diffuse_and_cool(state.heat, p)
    n_rel, _ = release(state, p)
    return size, n_rel


def run(L: int, n_steps: int, p: CascadeParams, seed: int | None = None) -> dict:
    """Drive the line for n_steps. Returns avalanche/heat traces + summary."""
    rng = np.random.default_rng(seed)
    state = initialize(L)
    max_sweeps = 50 * L
    cascade_sizes = np.zeros(n_steps, dtype=np.int64)
    sample_every = max(1, n_steps // 400)
    steps_axis, peak_heat_trace = [], []
    peak_heat_overall = 0.0
    min_heat_overall = 0.0
    first_release_step: int | None = None

    for t in range(n_steps):
        size, n_rel = step(state, p, rng, max_sweeps)
        cascade_sizes[t] = size
        if first_release_step is None and n_rel > 0:
            first_release_step = t
        hmax = float(state.heat.max())
        if hmax > peak_heat_overall:
            peak_heat_overall = hmax
        hmin = float(state.heat.min())
        if hmin < min_heat_overall:
            min_heat_overall = hmin
        if t % sample_every == 0:
            steps_axis.append(t)
            peak_heat_trace.append(hmax)

    q = max(1, n_steps // 4)
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
        "mean_cascade_first_quarter": float(cascade_sizes[:q].mean()),
        "mean_cascade_last_quarter": float(cascade_sizes[-q:].mean()),
        "cascade_sizes": cascade_sizes.tolist(),
        "steps_axis": steps_axis,
        "peak_heat_trace": peak_heat_trace,
    }
