"""Edgeless, GROWING substrate: the front never meets an edge.

We do not impose an edge (an edge would reflect energy back and smuggle in
"we are in a maker's box" -- and an edge must be proven, not assumed). So the
domain GROWS: whenever the catastrophe front nears the current border, we pad
fresh, pre-loaded substrate ahead of it. Energy never reflects, never wraps,
never returns -- it always advances into fresh substrate, matching the
premise of an endless pre-existing substrate.

Honest constraint still in force: this imposes only LOCAL connectivity (a
local cell-to-cell grid). It imposes NO edge and NO global geometry. Truly
imposing nothing (emergent geometry) is a separate, harder problem. Distances
are still local-grid (flat-local) -- flagged, not hidden.

The experiment: load + accumulate hidden damage (the eons), ignite one spot,
then propagate by heat-release cascade only (no driving) while growing the
domain. Question: does the front SUSTAIN and advance into fresh substrate, or
fizzle -- now with no edge to confuse the answer?

Reuses the 2D Manna avalanche, the heat diffusion, and the melt release rule.
"""

from __future__ import annotations

import numpy as np

from void_cascade.cascade_heat import CascadeParams
from void_cascade.cascade_heat_2d import avalanche_2d, diffuse_and_cool_2d, initialize_2d
from void_cascade.heat_gated import release


def near_border(mask: np.ndarray, margin: int) -> bool:
    return bool(mask[:margin].any() or mask[-margin:].any()
               or mask[:, :margin].any() or mask[:, -margin:].any())


def pad_domain(state, chunk: int, rng: np.random.Generator) -> None:
    """Grow symmetrically by `chunk` on all sides with FRESH pre-loaded
    substrate. Fresh cells are sampled from the current COLD (un-released)
    cells, so the new substrate is statistically identical to substrate that
    went through the same eons but the catastrophe has not yet reached. Heat
    starts at absolute zero there; nothing released yet. Old state is preserved
    in the centre."""
    d, cr, ht, rl = state.density, state.cracks, state.heat, state.released
    R, C = d.shape
    nR, nC = R + 2 * chunk, C + 2 * chunk
    cold = ~rl
    src_d = d[cold] if cold.any() else d.ravel()
    src_cr = cr[cold] if cold.any() else cr.ravel()
    idx = rng.integers(0, src_d.size, size=nR * nC)
    new_d = src_d[idx].reshape(nR, nC).astype(float)
    new_cr = src_cr[idx].reshape(nR, nC).astype(np.int64)
    new_ht = np.zeros((nR, nC), dtype=float)
    new_rl = np.zeros((nR, nC), dtype=bool)
    sl = (slice(chunk, chunk + R), slice(chunk, chunk + C))
    new_d[sl] = d
    new_cr[sl] = cr
    new_ht[sl] = ht
    new_rl[sl] = rl
    state.density, state.cracks, state.heat, state.released = new_d, new_cr, new_ht, new_rl


def ignite(state, p: CascadeParams, radius_frac: float = 0.05) -> None:
    """Push a central patch over the melt point to start the catastrophe."""
    R, C = state.heat.shape
    cy, cx = R // 2, C // 2
    r = max(1, int(radius_frac * min(R, C)))
    state.heat[cy - r:cy + r + 1, cx - r:cx + r + 1] = p.melt_heat * 2.0


def _fair_sites(n_cells: int, baseline: int = 48 * 48) -> int:
    return max(1, round(n_cells / baseline))


def front_radius(released: np.ndarray) -> float:
    ys, xs = np.nonzero(released)
    if xs.size == 0:
        return 0.0
    cy, cx = released.shape[0] / 2.0, released.shape[1] / 2.0
    return float(np.sqrt((xs - cx) ** 2 + (ys - cy) ** 2).max())


def run_grow(L0: int, accumulate_steps: int, propagate_steps: int, p: CascadeParams,
             seed: int, margin: int = 4, chunk: int = 20, max_size: int = 300) -> dict:
    """Phase 1: load + accumulate hidden cracks (heat off). Phase 2: ignite
    centre. Phase 3: propagate by heat-release cascade (no driving) while
    growing the domain so the front never meets an edge."""
    rng = np.random.default_rng(seed)
    st = initialize_2d(L0)
    st.density[:] = rng.uniform(0.0, 2.0, size=(L0, L0))
    avalanche_2d(st, p, rng, 50 * L0 * L0)

    # Phase 1: accumulate hidden damage with heat OFF.
    p_off = CascadeParams(**{**p.__dict__, "heat_per_crack": 0.0})
    for _ in range(accumulate_steps):
        n = _fair_sites(st.density.size)
        sites = rng.integers(0, st.density.size, size=n)
        np.add.at(st.density.ravel(), sites, p.drive_amount)
        avalanche_2d(st, p_off, rng, 50 * st.density.size)
    cracks_at_ignition = int(st.cracks.mean().round())

    # Phase 2: ignite.
    ignite(st, p)

    # Phase 3: propagate (heat-release cascade only) + grow.
    trace = []
    grow_events = 0
    capped = False
    for t in range(propagate_steps):
        st.heat = diffuse_and_cool_2d(st.heat, p)
        release(st, p)
        fr = front_radius(st.released)
        trace.append((t, int(st.released.sum()), fr, st.released.shape[0]))
        if near_border(st.released, margin):
            if st.released.shape[0] + 2 * chunk <= max_size:
                pad_domain(st, chunk, rng)
                grow_events += 1
            else:
                capped = True
                break

    fronts = [fr for _, _, fr, _ in trace]
    sizes_released = [n for _, n, _, _ in trace]
    sustained = grow_events > 0 and fronts[-1] > fronts[min(len(fronts) - 1, 10)]
    return {
        "L0": L0, "mean_cracks_at_ignition": cracks_at_ignition,
        "propagate_steps_run": len(trace),
        "grow_events": grow_events, "capped_at_max_size": capped,
        "final_size": int(st.released.shape[0]),
        "max_front_radius": float(max(fronts)) if fronts else 0.0,
        "final_released_cells": int(sizes_released[-1]) if sizes_released else 0,
        "front_sustained_and_grew": bool(sustained),
        "ever_touched_edge": False,  # by construction we grow before it can
        "trace": trace,
    }
