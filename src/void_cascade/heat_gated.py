"""Heat-gated substrate: cold-frozen buildup vs heat runaway.

See docs/notes/heat_gated_model_design.md for the full reasoning. Short
version: cells start at absolute zero (heat = 0). They slowly absorb
"darkness" (density rises) and crack when dense enough. While cold, a crack
is *frozen in place* -- it counts as hidden damage (cracks accumulate) but
nothing else happens. Each crack gives off a little heat; heat spreads to
neighbours and is cooled away (lost to the surrounding cold). If breaking
ever makes heat faster than cooling removes it, a cell crosses the melt
point and its accumulated cracks *release at once* -- a heat burst that can
push neighbours over too. That is the thermal runaway.

The honest question: does it stay cold forever (cooling wins) or run away
(breaking wins)? Built so BOTH are possible -- strong cooling must keep it
frozen, or the model is rigged.

No bonds: "frozen" is not re-gluing, it's that below the melt point breaks
cannot express as separation. Absolute zero is a hard floor (heat >= 0).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class HeatParams:
    fracture_density: float  # density at which a cell cracks
    fracture_relief: float   # density removed by one crack
    drive_rate: float        # darkness absorbed per step (slow)
    heat_per_crack: float    # heat released by one crack
    diffuse: float           # neighbour heat sharing (0..0.5 for stability)
    cooling: float           # fraction of heat lost per step (the guardrail)
    melt_heat: float         # heat above which frozen cracks release
    release_factor: float    # heat burst per released crack


@dataclass
class HeatState:
    density: np.ndarray   # float per cell
    cracks: np.ndarray    # int, accumulated fractures (hidden damage)
    heat: np.ndarray      # float >= 0 (absolute zero floor)
    released: np.ndarray  # bool, has this cell let its cracks go


def initialize(L: int, seed: int | None = None, density_noise: float = 0.5) -> HeatState:
    """All cold (heat 0). Small random starting density so cells crack at
    staggered times rather than all at once."""
    rng = np.random.default_rng(seed)
    return HeatState(
        density=rng.uniform(0.0, density_noise, size=L),
        cracks=np.zeros(L, dtype=np.int64),
        heat=np.zeros(L, dtype=float),
        released=np.zeros(L, dtype=bool),
    )


def drive(state: HeatState, p: HeatParams) -> None:
    """Absorb a little darkness everywhere."""
    state.density += p.drive_rate


def fracture(state: HeatState, p: HeatParams) -> int:
    """Cells over the density threshold crack: +1 crack, density relieved, a
    little heat released. The crack is frozen in place (handled by release)."""
    over = state.density >= p.fracture_density
    n = int(over.sum())
    if n:
        state.cracks[over] += 1
        state.density[over] -= p.fracture_relief
        state.heat[over] += p.heat_per_crack
    return n


def diffuse_and_cool(heat: np.ndarray, p: HeatParams) -> np.ndarray:
    """Spread heat to neighbours, lose some to the surrounding cold, floor at 0.

    1D Laplacian with open (cold = 0) boundaries, then a fractional cooling
    loss everywhere. Cooling is the guardrail: enough of it and heat can never
    build to the melt point.
    """
    lap = np.empty_like(heat)
    lap[1:-1] = heat[:-2] + heat[2:] - 2.0 * heat[1:-1]
    lap[0] = heat[1] - 2.0 * heat[0]    # left side open to cold (0)
    lap[-1] = heat[-2] - 2.0 * heat[-1]  # right side open to cold (0)
    out = heat + p.diffuse * lap
    out *= (1.0 - p.cooling)
    np.maximum(out, 0.0, out=out)
    return out


def release(state: HeatState, p: HeatParams) -> tuple[int, int]:
    """Cells at/above the melt point let their accumulated cracks go, once.

    Returns (n_cells_released, total_cracks_freed). Each releasing cell dumps
    a heat burst proportional to the cracks it had been hiding.
    """
    hot = (state.heat >= p.melt_heat) & (~state.released) & (state.cracks > 0)
    n = int(hot.sum())
    freed = int(state.cracks[hot].sum())
    if n:
        state.heat[hot] += p.release_factor * state.cracks[hot]
        state.released[hot] = True
        state.cracks[hot] = 0
    return n, freed


def step(state: HeatState, p: HeatParams) -> tuple[int, int]:
    """One timestep: absorb -> crack -> spread/cool heat -> maybe release."""
    drive(state, p)
    n_frac = fracture(state, p)
    state.heat = diffuse_and_cool(state.heat, p)
    n_rel, _ = release(state, p)
    return n_frac, n_rel


def run(L: int, n_steps: int, p: HeatParams, seed: int | None = None) -> dict:
    """Drive the line for n_steps. Returns traces + summary.

    `ran_away` is True iff any cell ever crossed the melt point and released.
    """
    state = initialize(L, seed=seed, density_noise=p.fracture_density * 0.5)
    sample_every = max(1, n_steps // 400)
    steps_axis, peak_heat, n_released, total_cracks = [], [], [], []
    peak_heat_overall = 0.0
    min_heat_overall = 0.0
    first_release_step: int | None = None
    hidden_cracks_at_first_release = 0

    for t in range(n_steps):
        _n_frac, n_rel = step(state, p)
        if first_release_step is None and n_rel > 0:
            first_release_step = t
            # cracks already zeroed on released cells; the freed count is what
            # let go this step -- recompute from this step's release is lost, so
            # report the total cracks still frozen elsewhere as the hidden load.
            hidden_cracks_at_first_release = int(state.cracks.sum())
        hmax = float(state.heat.max())
        if hmax > peak_heat_overall:
            peak_heat_overall = hmax
        hmin = float(state.heat.min())
        if hmin < min_heat_overall:
            min_heat_overall = hmin
        if t % sample_every == 0:
            steps_axis.append(t)
            peak_heat.append(hmax)
            n_released.append(int(state.released.sum()))
            total_cracks.append(int(state.cracks.sum()))

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
        "hidden_cracks_at_first_release": hidden_cracks_at_first_release,
        "steps_axis": steps_axis,
        "peak_heat_trace": peak_heat,
        "n_released_trace": n_released,
        "total_cracks_trace": total_cracks,
    }
