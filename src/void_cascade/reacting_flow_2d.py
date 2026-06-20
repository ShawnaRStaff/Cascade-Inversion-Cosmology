"""End-to-end (stages 2->3): 2D implosion that ignites a detonation.

One continuous run, one model: a ring of substrate implodes, collapses to the
centre, heats to a plasma spike -- and that heat IGNITES the surrounding fuel
(the loaded substrate), launching a detonation front that expands outward. So
the inversion (collapse/rebound, material motion) and the expansion (self-
feeding front) happen in a single 2D reacting compressible fluid, with the
implosion *causing* the ignition.

= material_motion_2d (compressible Euler, Lax-Friedrichs) + a fuel field that
combusts when the fluid gets hot (fuel -> internal energy). Energy conserved
(E_fluid + fuel = const, periodic). Honest scope: 2D reacting Euler; abstract
units; LF diffusive; this is stages 2+3 -- the buildup (stage 1) is folded in
next, on the same grid.
"""

from __future__ import annotations

import numpy as np

from void_cascade.material_motion_2d import (
    GAMMA, internal_energy, lax_friedrichs_step, max_wave_speed, primitives,
)


def combust(rho, momx, momy, E, fuel, e_ign):
    """Cells hot enough (specific internal energy >= e_ign) burn their fuel,
    releasing it into E (energy conserved: fuel -> internal energy)."""
    e_specific = internal_energy(rho, momx, momy, E) / rho
    fire = (e_specific >= e_ign) & (fuel > 0.0)
    released = np.where(fire, fuel, 0.0)
    return E + released, fuel - released, int(fire.sum())


def implosion_plus_fuel_ic(L, u0, R0, fuel0, rho0=1.0, P0=1.0):
    """Converging ring (implosion) on top of a uniformly loaded fuel field."""
    cy = cx = L / 2.0
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    rx, ry = xx - cx, yy - cy
    r = np.sqrt(rx * rx + ry * ry)
    speed = np.where(r < R0, -u0 * np.sin(np.pi * np.clip(r / R0, 0, 1)), 0.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        ux = np.where(r > 0, speed * rx / r, 0.0)
        uy = np.where(r > 0, speed * ry / r, 0.0)
    rho = np.full((L, L), rho0)
    momx, momy = rho * ux, rho * uy
    E = P0 / (GAMMA - 1.0) + 0.5 * rho * (ux * ux + uy * uy)
    fuel = np.full((L, L), fuel0)
    return rho, momx, momy, E, fuel


def run_end_to_end(L=160, u0=2.5, R0=50, fuel0=3.0, e_ign=2.5, steps=500, cfl=0.4,
                   dx=1.0, band=4, snap_steps=(0,)) -> dict:
    rho, momx, momy, E, fuel = implosion_plus_fuel_ic(L, u0, R0, fuel0)
    c = L // 2
    sl = (slice(c - band, c + band), slice(c - band, c + band))
    total0 = float(E.sum() + fuel.sum())
    fuel_init = float(fuel.sum())

    t_trace, dens_c, heat_c, burned_frac, burn_radius = [], [], [], [], []
    ignite_step = None
    snaps = {}
    yy, xx = np.mgrid[0:L, 0:L].astype(float)
    rr = np.sqrt((xx - c) ** 2 + (yy - c) ** 2)
    for t in range(steps):
        if t in snap_steps:
            snaps[t] = rho.copy()
        E, fuel, nfire = combust(rho, momx, momy, E, fuel, e_ign)
        if ignite_step is None and nfire > 0:
            ignite_step = t
        s = max_wave_speed(rho, momx, momy, E)
        dt = cfl * dx / max(s, 1e-9)
        rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, dx, dt)
        burned = fuel < fuel0 * 0.5
        t_trace.append(t)
        dens_c.append(float(rho[sl].mean()))
        heat_c.append(float(internal_energy(rho, momx, momy, E)[sl].mean()))
        burned_frac.append(float(burned.mean()))
        burn_radius.append(float(rr[burned].max()) if burned.any() else 0.0)
    snaps[steps - 1] = rho.copy()

    dens_c = np.asarray(dens_c)
    return {
        "L": L, "u0": u0, "R0": R0, "fuel0": fuel0, "e_ign": e_ign, "steps": steps,
        "energy_residual": float(E.sum() + fuel.sum() - total0),
        "ignite_step": ignite_step,
        "peak_central_density": float(dens_c.max()),
        "peak_density_step": int(dens_c.argmax()),
        "final_burned_fraction": float(burned_frac[-1]),
        "fuel_burned_fraction": float(1.0 - fuel.sum() / fuel_init) if fuel_init else 0.0,
        "max_burn_radius": float(max(burn_radius)),
        "t_trace": t_trace, "central_density": dens_c.tolist(), "central_heat": heat_c,
        "burned_fraction": burned_frac, "burn_radius": burn_radius, "snapshots": snaps,
    }
