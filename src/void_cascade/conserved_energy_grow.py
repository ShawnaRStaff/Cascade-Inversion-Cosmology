"""The full picture: conserved energy on an edgeless growing substrate.

Combines the two principles:
  - energy is CONSERVED, only changing form (potential <-> kinetic), and
  - the substrate is EDGELESS (grows fresh, pre-loaded substrate ahead of the
    front), so energy never reaches an edge -> boundary loss -> 0.

Together: energy is *truly* conserved (no cooling sink, no boundary sink), and
the front -- driven by kinetic energy priming the loaded substrate ahead --
runs free, into ever-fresh substrate, with nothing to stop it. That is the
researcher's picture: the still-rolling event we are far behind and mistake
for an echo.

We grow when KINETIC energy (which runs ahead of the flip front, priming it)
nears the border, so energy is padded with fresh substrate before it could
leak. Fresh substrate carries its own stored potential (loaded during its own
eons), so padding ADDS energy -- accounted into energy_in to keep the
conservation law exact.

Cornerstone: residual ~ 0 AND boundary_lost ~ 0 (true conservation).
Honest caveat unchanged: flat-local connectivity; we ignite the front;
combustion is the activity wavefront, not material motion.
"""

from __future__ import annotations

import numpy as np

from void_cascade.conserved_energy import (
    EnergyParams, EnergyState, combust, conservation_residual,
    diffuse_kinetic, ignite,
)


def initialize_2d(L: int) -> EnergyState:
    return EnergyState(
        potential=np.zeros((L, L), dtype=float),
        kinetic=np.zeros((L, L), dtype=float),
        flipped=np.zeros((L, L), dtype=bool),
    )


def kinetic_near_border(kinetic: np.ndarray, margin: int, eps: float = 1e-9) -> bool:
    return bool((kinetic[:margin] > eps).any() or (kinetic[-margin:] > eps).any()
               or (kinetic[:, :margin] > eps).any() or (kinetic[:, -margin:] > eps).any())


def pad_energy_domain(state: EnergyState, chunk: int, rng: np.random.Generator) -> None:
    """Grow symmetrically with fresh pre-loaded substrate (potential sampled
    from current cold cells; kinetic 0; nothing flipped). The fresh potential
    is new energy -> added to energy_in so conservation stays exact."""
    pot, kin, fl = state.potential, state.kinetic, state.flipped
    R, C = pot.shape
    nR, nC = R + 2 * chunk, C + 2 * chunk
    cold = ~fl
    src = pot[cold] if cold.any() else pot.ravel()
    idx = rng.integers(0, src.size, size=nR * nC)
    new_pot = src[idx].reshape(nR, nC).astype(float)
    sl = (slice(chunk, chunk + R), slice(chunk, chunk + C))
    added = float(new_pot.sum() - new_pot[sl].sum())  # border (fresh) potential
    new_kin = np.zeros((nR, nC), dtype=float)
    new_fl = np.zeros((nR, nC), dtype=bool)
    new_pot[sl] = pot
    new_kin[sl] = kin
    new_fl[sl] = fl
    state.potential, state.kinetic, state.flipped = new_pot, new_kin, new_fl
    state.energy_in += added


def front_radius(flipped: np.ndarray) -> float:
    ys, xs = np.nonzero(flipped)
    if xs.size == 0:
        return 0.0
    cy, cx = flipped.shape[0] / 2.0, flipped.shape[1] / 2.0
    return float(np.sqrt((xs - cx) ** 2 + (ys - cy) ** 2).max())


def run_grow_energy(L0: int, propagate_steps: int, p: EnergyParams, seed: int,
                    margin: int = 5, chunk: int = 20, max_size: int = 320,
                    load_lo: float = 0.6, load_hi: float = 0.98) -> dict:
    rng = np.random.default_rng(seed)
    st = initialize_2d(L0)
    st.potential[:] = rng.uniform(load_lo * p.flip_threshold,
                                  load_hi * p.flip_threshold, size=(L0, L0))
    st.energy_in = float(st.potential.sum())
    ignite(st, p)

    trace = []
    grow_events = 0
    capped = False
    for t in range(propagate_steps):
        combust(st, p)
        diffuse_kinetic(st, p)
        trace.append((t, int(st.flipped.sum()), front_radius(st.flipped),
                      st.flipped.shape[0]))
        if kinetic_near_border(st.kinetic, margin):
            if st.flipped.shape[0] + 2 * chunk <= max_size:
                pad_energy_domain(st, chunk, rng)
                grow_events += 1
            else:
                capped = True
                break

    fronts = [fr for _, _, fr, _ in trace]
    sustained = grow_events > 0 and fronts[-1] > fronts[min(len(fronts) - 1, 10)]
    return {
        "L0": L0, "steps_run": len(trace), "grow_events": grow_events,
        "capped": capped, "final_size": int(st.flipped.shape[0]),
        "max_front_radius": float(max(fronts)) if fronts else 0.0,
        "final_flipped_fraction": float(st.flipped.mean()),
        "front_sustained_and_grew": bool(sustained),
        "conservation_residual": float(conservation_residual(st)),
        "boundary_lost": float(st.boundary_lost),
        "energy_in": float(st.energy_in),
        "trace": trace,
    }
