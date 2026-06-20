"""Material motion in 2D: the radial implosion -> plasma -> rebound shell.

The synthesis run. A ring of substrate rushes inward (an implosion), collapses
toward the centre, compresses and heats to a plasma-like spike, and the
pressure rebounds it outward as an expanding shell -- and that outward shell
IS the expansion front. One model, the whole inversion->expansion arc, with
real material motion (momentum) and conserved energy.

2D compressible Euler, Lax-Friedrichs (robust), periodic boundaries so mass /
momentum / energy are conserved to machine precision. Conserved variables per
cell: rho, momx, momy, E. Pressure P=(gamma-1)(E - 0.5 rho (u^2+v^2)),
gamma=5/3. Internal energy (heat) = E - 0.5 rho (u^2+v^2).

Honest scope: generic 2D compressible fluid (the continuum law for moving
material with pressure), connected to the theory as the substrate imploding
and rebounding. Flat-local grid; abstract units (no calibration yet).
"""

from __future__ import annotations

import numpy as np

GAMMA = 5.0 / 3.0


def primitives(rho, momx, momy, E):
    u = momx / rho
    v = momy / rho
    P = (GAMMA - 1.0) * (E - 0.5 * (momx * u + momy * v))
    return u, v, P


def internal_energy(rho, momx, momy, E):
    return E - 0.5 * (momx * momx + momy * momy) / rho


def max_wave_speed(rho, momx, momy, E) -> float:
    u, v, P = primitives(rho, momx, momy, E)
    c = np.sqrt(np.maximum(GAMMA * P / rho, 0.0))
    return float(np.max(np.sqrt(u * u + v * v) + c))


def _flux_x(rho, momx, momy, E):
    u, v, P = primitives(rho, momx, momy, E)
    return momx, momx * u + P, momx * v, (E + P) * u


def _flux_y(rho, momx, momy, E):
    u, v, P = primitives(rho, momx, momy, E)
    return momy, momy * u, momy * v + P, (E + P) * v


def lax_friedrichs_step(rho, momx, momy, E, dx, dt):
    Fx = _flux_x(rho, momx, momy, E)
    Fy = _flux_y(rho, momx, momy, E)
    U = (rho, momx, momy, E)

    def upd(Ui, Fxi, Fyi):
        avg = 0.25 * (np.roll(Ui, -1, 0) + np.roll(Ui, 1, 0)
                     + np.roll(Ui, -1, 1) + np.roll(Ui, 1, 1))
        dFx = np.roll(Fxi, -1, 0) - np.roll(Fxi, 1, 0)
        dFy = np.roll(Fyi, -1, 1) - np.roll(Fyi, 1, 1)
        return avg - (dt / (2.0 * dx)) * (dFx + dFy)

    return tuple(upd(U[i], Fx[i], Fy[i]) for i in range(4))


def radial_converging_ic(L, u0, R0, rho0=1.0, P0=1.0):
    """A ring of material (peak speed at R0/2) rushing inward toward the centre."""
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
    return rho, momx, momy, E


def run_radial_implosion(L=160, u0=1.5, R0=50, steps=500, cfl=0.4, dx=1.0,
                         band=5, snap_steps=(0,)) -> dict:
    rho, momx, momy, E = radial_converging_ic(L, u0, R0)
    c = L // 2
    sl = (slice(c - band, c + band), slice(c - band, c + band))
    mass0, mx0, my0, en0 = rho.sum(), momx.sum(), momy.sum(), E.sum()

    t_trace, dens_c, heat_c = [], [], []
    min_P, min_rho = np.inf, np.inf
    snaps = {}
    for t in range(steps):
        if t in snap_steps:
            snaps[t] = rho.copy()
        s = max_wave_speed(rho, momx, momy, E)
        dt = cfl * dx / max(s, 1e-9)
        rho, momx, momy, E = lax_friedrichs_step(rho, momx, momy, E, dx, dt)
        _u, _v, P = primitives(rho, momx, momy, E)
        min_P = min(min_P, float(P.min()))
        min_rho = min(min_rho, float(rho.min()))
        t_trace.append(t)
        dens_c.append(float(rho[sl].mean()))
        heat_c.append(float(internal_energy(rho, momx, momy, E)[sl].mean()))
    snaps[steps - 1] = rho.copy()

    dens_c = np.asarray(dens_c)
    peak_i = int(dens_c.argmax())
    return {
        "L": L, "u0": u0, "R0": R0, "steps": steps,
        "mass_residual": float(rho.sum() - mass0),
        "energy_residual": float(E.sum() - en0),
        "momentum_residual": float(abs(momx.sum() - mx0) + abs(momy.sum() - my0)),
        "min_pressure": float(min_P), "min_density": float(min_rho),
        "peak_central_density": float(dens_c.max()),
        "central_density_initial": float(dens_c[0]),
        "peak_density_step": peak_i,
        "central_heat_initial": float(heat_c[0]),
        "central_heat_at_peak": float(heat_c[peak_i]),
        "rebounded": bool(peak_i < steps - 1 and dens_c[-1] < dens_c[peak_i]),
        "t_trace": t_trace, "central_density": dens_c.tolist(), "central_heat": heat_c,
        "snapshots": snaps,
    }
