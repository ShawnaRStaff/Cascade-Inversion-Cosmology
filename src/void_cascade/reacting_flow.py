"""Reacting flow (1D detonation): the unification of motion + self-feeding front.

One energy-conserving model with BOTH:
  - a compressible fluid (rho, mom, E) -- material motion (the implosion/
    rebound physics), and
  - a fixed-in-space `fuel` field -- the pre-loaded substrate's stored energy.

Where the fluid gets hot enough (specific internal energy >= ignition), the
local fuel COMBUSTS: its stored energy is released into the fluid (fuel -> E),
raising pressure, driving the shock, which heats the fluid ahead and ignites
the fuel ahead -- a self-sustaining DETONATION front. This couples the
inversion (collapse/rebound, material motion) with the expansion (self-feeding
front) in a single model, instead of bolting two representations together.

Energy is conserved: E_total(fluid) + fuel(remaining) = const (periodic BC).
Fuel is fixed in space (the substrate is the medium the front runs through).
Honest scope: 1D reacting Euler (Lax-Friedrichs); abstract units; a detonation
is the bubble-wall / phase-transition-front family of the project's references.
"""

from __future__ import annotations

import numpy as np

from void_cascade.material_motion import (
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed,
)


def combust(rho, mom, E, fuel, e_ign: float):
    """Cells whose specific internal energy exceeds ignition burn their fuel,
    releasing it into E (energy conserved: fuel -> internal energy)."""
    e_specific = internal_energy(rho, mom, E) / rho
    fire = (e_specific >= e_ign) & (fuel > 0.0)
    released = np.where(fire, fuel, 0.0)
    return E + released, fuel - released, int(fire.sum())


def run_detonation(N=400, fuel0=3.0, e_ign=2.5, rho0=1.0, P0=1.0, ignite_E=3.0,
                   steps=500, cfl=0.4, dx=1.0) -> dict:
    """Loaded substrate (fuel everywhere), ignite the centre, watch a
    detonation front propagate outward through the fuel."""
    rho = np.full(N, rho0, dtype=float)
    mom = np.zeros(N, dtype=float)
    E = np.full(N, P0 / (GAMMA - 1.0), dtype=float)
    fuel = np.full(N, fuel0, dtype=float)
    c = N // 2
    E[c - 3:c + 4] += ignite_E  # hot spot to ignite the centre
    total0 = float(E.sum() + fuel.sum())

    t_trace, front, burned_frac = [], [], []
    min_P = np.inf
    for t in range(steps):
        E, fuel, _ = combust(rho, mom, E, fuel, e_ign)
        s = max_wave_speed(rho, mom, E)
        dt = cfl * dx / max(s, 1e-9)
        rho, mom, E = lax_friedrichs_step(rho, mom, E, dx, dt)
        _u, P = _primitives_P(rho, mom, E)
        min_P = min(min_P, float(P.min()))
        burned = fuel < fuel0 * 0.5
        idx = np.nonzero(burned)[0]
        fr = float(np.abs(idx - c).max()) if idx.size else 0.0
        t_trace.append(t)
        front.append(fr)
        burned_frac.append(float(burned.mean()))

    total1 = float(E.sum() + fuel.sum())
    front = np.asarray(front)
    # detonation speed: slope of front position vs step, over the growth region
    grew = front < (N / 2 - 3)
    if grew.sum() > 5:
        tt = np.asarray(t_trace)[grew]
        speed = float(np.polyfit(tt, front[grew], 1)[0])
    else:
        speed = 0.0
    return {
        "N": N, "fuel0": fuel0, "e_ign": e_ign, "steps": steps,
        "energy_residual": total1 - total0,
        "min_pressure": float(min_P),
        "final_burned_fraction": float(burned_frac[-1]),
        "max_front": float(front.max()),
        "front_speed": speed,
        "propagated": bool(burned_frac[-1] > 0.5),
        "t_trace": t_trace, "front": front.tolist(), "burned_fraction": burned_frac,
    }


def _primitives_P(rho, mom, E):
    u = mom / rho
    P = (GAMMA - 1.0) * (E - 0.5 * mom * u)
    return u, P
